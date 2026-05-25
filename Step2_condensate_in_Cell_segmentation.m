% Batch Condensate Analysis with Watershed and Shape Filtering
% Author: debalina-datta
% Date: 2025-09
% Processes multiple TIFF images and corresponding TIFF masks.
% Uses adaptive+global thresholding, watershed splitting, and shape filtering.
%
% Pipeline overview:
%   1. For each image frame, load the raw fluorescence image and its
%      Cellpose-generated cell mask.
%   2. Iterate over every segmented cell in the mask.
%   3. Run three complementary condensate detection passes:
%        a) Standard condensates  — adaptive + Otsu thresholding + watershed
%        b) Super-saturated rings — bright ring structures near pixel saturation
%        c) Super-saturated spots — fully saturated point-like condensates
%   4. Apply multi-parameter shape filtering to reject noise and artefacts.
%   5. Save per-frame CSV results and annotated overlay images.

%% Enable parallel processing
% Start a local parallel pool if one is not already running.
% parfor below distributes frames across available CPU workers.
if isempty(gcp('nocreate'))
    parpool('local');
end

%% Input Parameters
clear all;
% close all;

originalImageFolder = 'cyto_K907A_slices'; % folder containing raw fluorescence TIFF slices
cellPoseMaskFolder  = 'cyto_K907A_masks';  % folder containing Cellpose 16-bit cell masks

outputFolderPrefix        = 'cyto_K907';   % prefix used for all output sub-folders
pixelToMicron             = 221.87 / 1024; % physical pixel size in µm — update from image metadata

% Derived output folder names
csvOutputFolder        = [outputFolderPrefix, '_csv_results'];
overlayOutputFolder    = [outputFolderPrefix, '_overlays'];
overlayNoIDOutputFolder = [outputFolderPrefix, '_overlays_no_id'];
mkdir(csvOutputFolder);
mkdir(overlayOutputFolder);
mkdir(overlayNoIDOutputFolder);

% Collect file listings — expects filenames matching 'frame*.tif'
originalImages = dir(fullfile(originalImageFolder, 'frame*.tif'));
cellPoseMasks  = dir(fullfile(cellPoseMaskFolder,  'frame*_mask.tif'));

if isempty(originalImages) || isempty(cellPoseMasks)
    error('No TIFF files found in the specified folders.');
end

%% Analysis — one frame per parallel worker
parfor fileIdx = 1:length(originalImages)

    % Derive the expected mask filename from the image filename
    [~, tiffBaseName, ~] = fileparts(originalImages(fileIdx).name);
    maskFileName = [tiffBaseName, '_mask'];
    maskFilePath = fullfile(cellPoseMaskFolder, [maskFileName, '.tif']);

    % Skip frames whose mask is missing
    if ~isfile(maskFilePath)
        disp(['Warning: Mask not found for ', originalImages(fileIdx).name, '. Skipping...']);
        continue;
    end

    % Load the raw image and the corresponding Cellpose cell mask
    originalImage = imread(fullfile(originalImages(fileIdx).folder, originalImages(fileIdx).name));
    cellPoseMask  = imread(maskFilePath);

    % Skip frames where Cellpose found no cells (mask is all background)
    if all(cellPoseMask(:) == 0)
        disp(['Warning: Mask file ', maskFilePath, ' is empty (all zeros). Skipping...']);
        continue;
    end

    % Convert RGB masks to grayscale if Cellpose saved a colour PNG
    if size(cellPoseMask, 3) > 1
        cellPoseMask = rgb2gray(cellPoseMask);
    end

    % Cast to double so integer cell IDs can be compared precisely
    cellPoseMask  = double(cellPoseMask);
    uniqueCellIDs = unique(cellPoseMask);
    uniqueCellIDs(uniqueCellIDs == 0) = []; % remove background label (0)

    % Preallocate per-frame accumulators
    results          = [];                          % will hold one row per condensate
    condensateOverlay = zeros(size(originalImage)); % accumulates all accepted condensate masks
    condensateID     = 0;                           % running condensate counter across all cells

    % =========================================================
    % PASS 1 — Standard Condensate Detection
    %   Strategy: combine adaptive (local) and Otsu (global) thresholds
    %   so that both dim diffuse condensates and bright focal ones are captured.
    %   Watershed then separates touching objects before shape filtering.
    % =========================================================
    for i = 1:length(uniqueCellIDs)
        cellID   = uniqueCellIDs(i);
        cellMask = (cellPoseMask == cellID);

        % Normalise intensity to [0,1] within the cell region only;
        % pixels outside the cell are zeroed so they don't influence thresholds.
        cellImage = mat2gray(originalImage) .* cellMask;

        % Whole-cell summary measurements (used in the CSV output)
        cellAreaPixels    = sum(cellMask(:));
        cellNumPixels     = cellAreaPixels;
        cellTotalIntensity = sum(originalImage(cellMask));
        cellMeanIntensity  = cellTotalIntensity / cellNumPixels;
        cellAreaMicron     = cellAreaPixels * pixelToMicron^2;

        % --- 1a. Dual thresholding ---
        % Adaptive threshold: captures dim condensates that fall below the global Otsu level.
        % Sensitivity 0.55 and 21x21 neighbourhood — tune based on typical condensate size.
        localThreshold = adaptthresh(cellImage, 0.55, 'NeighborhoodSize', [21,21]);
        binaryLocal    = imbinarize(cellImage, localThreshold);

        % Global Otsu threshold: reliably detects bright, high-contrast condensates.
        globalThresh  = graythresh(cellImage);
        binaryGlobal  = imbinarize(cellImage, globalThresh);

        % Union of both masks — keeps any object flagged by either method.
        % bwareaopen removes isolated 1-pixel noise specks (< 2 px).
        binaryCondensates = binaryLocal | binaryGlobal;
        binaryCondensates = bwareaopen(binaryCondensates, 2);

        % --- 1b. Watershed splitting ---
        % Distance transform peaks correspond to condensate centres.
        % imextendedmax suppresses shallow local maxima (h=2) to avoid over-splitting.
        % Watershed then cuts along the ridge lines between adjacent peaks.
        distanceTransform = bwdist(~binaryCondensates);
        localMax          = imextendedmax(distanceTransform, 2);
        imposedMin        = imimposemin(-distanceTransform, localMax);
        watershedLabels   = watershed(imposedMin);

        % Watershed lines (label == 0) lie on object boundaries — remove them
        % so touching condensates become separate labelled regions.
        binaryCondensates(watershedLabels == 0) = 0;

        % --- 1c. Label regions and measure shape/intensity properties ---
        [condensateLabels, numCondensates] = bwlabel(binaryCondensates);

        % Shape filtering thresholds — adjust to match your condensate biology
        minSize        = 4;                        % minimum area in pixels
        maxSize        = 500;                      % maximum area in pixels
        minIntensity   = 0.1 * max(cellImage(:)); % must be at least 10% of cell peak
        maxEccentricity = 0.95;                   % rejects highly elongated objects
        maxAspectRatio  = 5;                       % major/minor axis ratio cap
        minSolidity     = 0.8;                     % rejects irregular / fragmented shapes

        stats = regionprops(condensateLabels, cellImage, ...
            'Area', 'Perimeter', 'MeanIntensity', 'PixelList', 'Centroid', 'PixelValues', ...
            'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Solidity');

        foundCondensate = false; % tracks whether any condensate passed filters in this cell

        for regionIdx = 1:numCondensates
            area         = stats(regionIdx).Area;
            perimeter    = stats(regionIdx).Perimeter;
            meanIntensity = stats(regionIdx).MeanIntensity;
            pixelValues  = stats(regionIdx).PixelValues;

            maxPixelValue  = max(pixelValues);
            maxToMeanRatio = maxPixelValue / mean(pixelValues); % >1 indicates a bright focal peak
            stdIntensity   = std(double(pixelValues));          % spread of intensity within the region
            circularity    = (4 * pi * area) / (perimeter^2);  % 1 = perfect circle; <1 = elongated/irregular
            eccentricity   = stats(regionIdx).Eccentricity;
            aspectRatio    = stats(regionIdx).MajorAxisLength / max(1, stats(regionIdx).MinorAxisLength);
            solidity       = stats(regionIdx).Solidity;

            % Reject region if any filter criterion fails;
            % maxToMeanRatio >= 1.6 ensures the region has a genuine bright centre
            % rather than being a diffuse haze of roughly uniform intensity.
            if area < minSize || area > maxSize         || ...
               meanIntensity   < minIntensity           || ...
               maxToMeanRatio  < 1.6                    || ...
               eccentricity    > maxEccentricity        || ...
               aspectRatio     > maxAspectRatio         || ...
               solidity        < minSolidity
                % Erase rejected region from the binary mask
                binaryCondensates(condensateLabels == regionIdx) = 0;
            else
                % Region passed all filters — record it
                foundCondensate = true;
                condensateID    = condensateID + 1;
                % Feret diameter = longest distance between any two boundary pixels
                feretDiameter = max(pdist(stats(regionIdx).PixelList));
                results = [results; cellID, condensateID, cellAreaMicron, ...
                           area * pixelToMicron^2, ...
                           circularity, meanIntensity, ...
                           feretDiameter * pixelToMicron, ...
                           maxToMeanRatio, stdIntensity, ...
                           cellNumPixels, cellTotalIntensity, cellMeanIntensity];
            end
        end

        % Cells with no valid condensate still get a row so they appear in the CSV
        % (all condensate columns set to NaN to distinguish from detected values).
        if ~foundCondensate
            results = [results; cellID, NaN, cellAreaMicron, NaN, NaN, NaN, NaN, NaN, NaN, ...
                       cellNumPixels, cellTotalIntensity, cellMeanIntensity];
        end

        % Accumulate accepted condensate pixels for the overlay image
        condensateOverlay = condensateOverlay + binaryCondensates;
    end


    % =========================================================
    % PASS 2 — Super-Saturated Ring Condensate Detection
    %   Targets ring-shaped condensates whose centre pixel intensity
    %   approaches the detector saturation limit.  Standard thresholding
    %   can miss these because the ring interior may be dim relative to
    %   the bright rim.  Strategy: threshold at 95% of per-cell maximum,
    %   then dilate to bridge any gap across the ring centre.
    % =========================================================
    for i = 1:length(uniqueCellIDs)
        cellID   = uniqueCellIDs(i);
        cellMask = (cellPoseMask == cellID);

        % Work in raw (non-normalised) intensity space to preserve absolute brightness
        maxIntensityValue = double(max(originalImage(cellMask)));
        highThresh        = 0.95 * maxIntensityValue; % adjust threshold fraction if rings are missed

        % Select pixels brighter than the threshold and inside the cell
        brightMask = (double(originalImage) >= highThresh) & cellMask;

        % Dilate by a small disk to bridge the dim gap in the ring interior.
        % Increase disk radius if rings have a wide dim centre.
        se          = strel('disk', 2);
        dilatedMask = imdilate(brightMask, se);
        dilatedMask = bwareaopen(dilatedMask, minSize); % remove noise specks

        [ringLabels, numRings] = bwlabel(dilatedMask);
        ringStats = regionprops(ringLabels, originalImage, ...
            'Area', 'Perimeter', 'PixelList', 'Centroid', 'PixelValues', ...
            'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Solidity');

        for rrIdx = 1:numRings
            area        = ringStats(rrIdx).Area;
            perimeter   = ringStats(rrIdx).Perimeter;
            circularity = (4 * pi * area) / (perimeter^2);
            eccentricity = ringStats(rrIdx).Eccentricity;
            aspectRatio = ringStats(rrIdx).MajorAxisLength / max(1, ringStats(rrIdx).MinorAxisLength);
            solidity    = ringStats(rrIdx).Solidity;

            % Ring condensates are accepted with a relaxed circularity floor (0.1)
            % since dilation can distort the shape slightly.
            % The maxToMeanRatio and stdIntensity filters are not applied here
            % because intensity is nearly uniform near saturation.
            if area >= minSize && area <= maxSize       && ...
               circularity  >= 0.1                     && ...
               eccentricity <= maxEccentricity          && ...
               aspectRatio  <= maxAspectRatio           && ...
               solidity     >= minSolidity
                condensateID  = condensateID + 1;
                feretDiameter = max(pdist(ringStats(rrIdx).PixelList));
                % MaxToMeanRatio set to 1, StdIntensity to 0 — not meaningful at saturation
                results = [results; cellID, condensateID, ...
                           cellAreaPixels * pixelToMicron^2, ...
                           area * pixelToMicron^2, ...
                           circularity, maxIntensityValue, ...
                           feretDiameter * pixelToMicron, ...
                           1, 0, cellNumPixels, cellTotalIntensity, cellMeanIntensity];
                condensateOverlay = condensateOverlay + (ringLabels == rrIdx);
            end
        end
    end


    % =========================================================
    % PASS 3 — Super-Saturated Spot Detection
    %   Catches small, fully saturated point condensates that reach the
    %   global image maximum (i.e., pixel value == max across the entire
    %   frame).  These may be under-detected by adaptive thresholding
    %   if their neighbourhood is uniformly bright.
    % =========================================================
    for i = 1:length(uniqueCellIDs)
        cellID   = uniqueCellIDs(i);
        cellMask = (cellPoseMask == cellID);
        cellImage = mat2gray(originalImage) .* cellMask;

        % Global image maximum — pixels at this value are considered fully saturated
        maxIntensityValue = double(max(originalImage(:)));
        superSatMask      = (cellImage == maxIntensityValue) & cellMask;
        superSatMask      = bwareaopen(superSatMask, minSize);

        [superSatLabels, numSuperSat] = bwlabel(superSatMask);
        superSatStats = regionprops(superSatLabels, cellImage, ...
            'Area', 'Perimeter', 'PixelList', 'Centroid', 'PixelValues', ...
            'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Solidity');

        for ssIdx = 1:numSuperSat
            area        = superSatStats(ssIdx).Area;
            perimeter   = superSatStats(ssIdx).Perimeter;
            circularity = (4 * pi * area) / (perimeter^2);
            eccentricity = superSatStats(ssIdx).Eccentricity;
            aspectRatio = superSatStats(ssIdx).MajorAxisLength / max(1, superSatStats(ssIdx).MinorAxisLength);
            solidity    = superSatStats(ssIdx).Solidity;

            % Same relaxed acceptance criteria as ring condensates (Pass 2)
            if area >= minSize && area <= maxSize       && ...
               circularity  >= 0.1                     && ...
               eccentricity <= maxEccentricity          && ...
               aspectRatio  <= maxAspectRatio           && ...
               solidity     >= minSolidity
                condensateID  = condensateID + 1;
                feretDiameter = max(pdist(superSatStats(ssIdx).PixelList));
                results = [results; cellID, condensateID, cellAreaMicron, ...
                           area * pixelToMicron^2, ...
                           circularity, maxIntensityValue, ...
                           feretDiameter * pixelToMicron, ...
                           1, 0, cellNumPixels, cellTotalIntensity, cellMeanIntensity];
                condensateOverlay = condensateOverlay + (superSatLabels == ssIdx);
            end
        end
    end

    % =========================================================
    % Save Results
    % =========================================================

    % --- CSV: one row per condensate (NaN row where no condensate was found) ---
    resultsTable = array2table(results, ...
        'VariableNames', {'CellID', 'CondensateID', 'CellArea', 'CondensateSize', ...
                          'Circularity', 'MeanIntensity', 'FeretDiameter', ...
                          'MaxToMeanRatio', 'StdIntensity', ...
                          'CellNumPixels', 'CellTotalIntensity', 'CellMeanIntensity'});

    csvFileName = fullfile(csvOutputFolder, [tiffBaseName, '_results.csv']);
    writetable(resultsTable, csvFileName);

    % --- Overlay with condensate boundaries and numeric IDs ---
    % Green outlines + yellow ID labels make it easy to cross-reference the CSV.
    figure('Visible','off');
    imshow(originalImage, []);
    hold on;
    boundaries   = bwboundaries(condensateOverlay > 0);
    statsOverlay = regionprops(condensateOverlay > 0, 'Centroid');
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1.5);
        if k <= length(statsOverlay)
            centroid = statsOverlay(k).Centroid;
            text(centroid(1), centroid(2), num2str(k), 'Color', 'y', ...
                'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
    end
    hold off;
    overlayFileName = fullfile(overlayOutputFolder, [tiffBaseName, '_overlay.png']);
    saveas(gcf, overlayFileName);
    close;

    % --- Overlay without IDs — cleaner version for figures / publications ---
    figure('Visible','off');
    imshow(originalImage, []);
    hold on;
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1.5);
    end
    hold off;
    overlayNoIDFileName = fullfile(overlayNoIDOutputFolder, [tiffBaseName, '_overlay_no_id.png']);
    saveas(gcf, overlayNoIDFileName);
    close;

    disp(['Processed ', tiffBaseName, ': CSV saved to ', csvFileName, ...
          ', overlays saved to ', overlayFileName, ' and ', overlayNoIDFileName]);
end

disp('Batch analysis complete! All files processed.');
