# Cell Cluster Segmentation

A two-step pipeline for automated segmentation and quantitative analysis of protein condensates (clusters) within individual cells from fluorescence microscopy images.

**Step 1 — Cell Segmentation (Python):** Batch-processes TIFF/PNG/JPG images using a custom-trained [Cellpose 2.0](https://github.com/MouseLand/cellpose) model to generate per-cell binary masks.

**Step 2 — Condensate Analysis (MATLAB):** Detects, segments, and quantifies fluorescent condensates inside each segmented cell using adaptive thresholding, watershed splitting, and multi-parameter shape filtering. Exports per-condensate measurements to CSV and produces annotated overlay images.

---

## Pipeline Overview

```
Raw fluorescence images (TIFF/PNG/JPG)
            │
            ▼
 ┌─────────────────────────────┐
 │  Step 1 · Python (Cellpose) │
 │  step1_cellpose_segmentation│
 │  _batch.py                  │
 └─────────────────────────────┘
            │  16-bit TIFF cell masks
            ▼
 ┌─────────────────────────────┐
 │  Step 2 · MATLAB            │
 │  Step2_condensate_in_Cell_  │
 │  segmentation.m             │
 └─────────────────────────────┘
            │
     ┌──────┼──────────────┐
     ▼      ▼              ▼
   CSV    Overlay PNGs   Overlay PNGs
 results  (with IDs)    (without IDs)
```

---

## Prerequisites

### Step 1 — Python environment

| Dependency | Tested version |
|---|---|
| Python | ≥ 3.8 |
| [Cellpose](https://github.com/MouseLand/cellpose) | 2.x |
| Pillow | any recent |
| tifffile | any recent |
| NumPy | any recent |
| PyTorch | ≥ 1.10 (with optional CUDA for GPU) |

Install all Python dependencies with:

```bash
pip install cellpose pillow tifffile numpy torch
```

> **GPU support:** If a CUDA-capable GPU is available, the script detects it automatically via `torch.cuda.is_available()` and runs Cellpose on GPU. No manual flag is required.

### Step 2 — MATLAB environment

| Requirement | Notes |
|---|---|
| MATLAB | R2020b or later recommended |
| Image Processing Toolbox | Required (`adaptthresh`, `watershed`, `regionprops`, etc.) |
| Parallel Computing Toolbox | Optional but recommended for `parfor` speedup |

---

## Repository Structure

```
Cell-Cluster-Segmentation/
├── step1_cellpose_segmentation_batch.py   # Python: cell segmentation via Cellpose
├── Step2_condensate_in_Cell_segmentation.m # MATLAB: condensate analysis & quantification
├── README.md
└── LICENSE
```

---

## Usage

### Step 1 · Cell Segmentation (Python)

**1. Prepare your images**

Place all raw fluorescence image files (`.tif`, `.tiff`, `.png`, `.jpg`, or `.jpeg`) in a single input folder.

**2. Configure parameters** inside `step1_cellpose_segmentation_batch.py`:

```python
input_folder  = './cyto_K907A_slices'   # folder containing raw images
output_folder = './cyto_K907A_masks'    # folder where masks will be saved
model_type    = 'CP_optoirek907a'       # custom Cellpose 2.0 model name
diameter      = 100                     # estimated cell diameter (pixels)
flow_threshold      = 0.4              # Cellpose flow threshold
cellprob_threshold  = 0.0              # Cellpose cell probability threshold
```

**3. Run the script**

```bash
python step1_cellpose_segmentation_batch.py
```

**Output:** One 16-bit TIFF mask file per input image saved to `output_folder`, named `<original_name>_mask.tif`. Each pixel in the mask encodes a unique integer cell ID (0 = background).

The script processes all images in parallel using `ProcessPoolExecutor` to maximize throughput.

---

### Step 2 · Condensate Analysis (MATLAB)

**1. Configure parameters** at the top of `Step2_condensate_in_Cell_segmentation.m`:

```matlab
originalImageFolder = 'cyto_K907A_slices';  % raw TIFF image folder
cellPoseMaskFolder  = 'cyto_K907A_masks';   % mask folder (Step 1 output)
outputFolderPrefix  = 'cyto_K907';          % prefix for all output folders
pixelToMicron       = 221.87 / 1024;        % microns per pixel (update from image metadata)
```

**2. Run the script** in MATLAB:

```matlab
run('Step2_condensate_in_Cell_segmentation.m')
```

or open the file in the MATLAB Editor and press **Run**.

**Output folders created automatically:**

| Folder | Contents |
|---|---|
| `<prefix>_csv_results/` | One CSV per image with per-condensate measurements |
| `<prefix>_overlays/` | Overlay PNGs with condensate boundaries and numeric IDs |
| `<prefix>_overlays_no_id/` | Overlay PNGs with condensate boundaries only (no labels) |

---

## Output CSV Columns

Each row in the CSV corresponds to one detected condensate (or a NaN row if no condensate passed filters in a given cell).

| Column | Units | Description |
|---|---|---|
| `CellID` | — | Integer ID of the parent cell (from Cellpose mask) |
| `CondensateID` | — | Unique condensate counter within the image |
| `CellArea` | µm² | Area of the parent cell |
| `CondensateSize` | µm² | Area of the condensate |
| `Circularity` | — | 4π·Area / Perimeter² (1 = perfect circle) |
| `MeanIntensity` | a.u. | Mean pixel intensity within the condensate |
| `FeretDiameter` | µm | Maximum caliper (Feret) diameter |
| `MaxToMeanRatio` | — | Ratio of maximum to mean pixel intensity |
| `StdIntensity` | a.u. | Standard deviation of pixel intensities |
| `CellNumPixels` | pixels | Total pixel count of the parent cell |
| `CellTotalIntensity` | a.u. | Sum of all pixel intensities in the parent cell |
| `CellMeanIntensity` | a.u. | Mean pixel intensity of the parent cell |

---

## Condensate Detection Logic (Step 2)

The MATLAB script runs three complementary detection passes for each cell:

1. **Standard condensates** — Combined adaptive (neighborhood 21×21, sensitivity 0.55) and Otsu global thresholding, followed by watershed splitting of touching objects.
2. **Super-saturated ring condensates** — Detects bright ring-like structures by thresholding at 95 % of the per-cell maximum intensity and dilating to capture dim ring centers.
3. **Super-saturated spot condensates** — Finds pixels at the global image maximum, indicating fully saturated bright spots.

All detected regions are passed through a multi-parameter shape filter before being accepted:

| Filter | Threshold |
|---|---|
| Area | 4 – 500 px |
| Mean intensity | ≥ 10 % of cell maximum |
| Eccentricity | ≤ 0.95 |
| Aspect ratio | ≤ 5 |
| Solidity | ≥ 0.80 |
| Max-to-mean ratio | ≥ 1.6 (standard pass only) |

---

## Key Parameters to Tune

| Parameter | Location | Effect |
|---|---|---|
| `diameter` | Step 1 Python | Expected cell diameter in pixels — the most impactful Cellpose parameter |
| `flow_threshold` | Step 1 Python | Lower values produce more conservative (fewer) cell boundaries |
| `cellprob_threshold` | Step 1 Python | Adjusts sensitivity to dim cells |
| `pixelToMicron` | Step 2 MATLAB | Must match your microscope's pixel size |
| Adaptive threshold sensitivity (`0.55`) | Step 2 MATLAB | Controls condensate segmentation sensitivity; tune based on cluster size |
| Adaptive neighborhood size (`[21,21]`) | Step 2 MATLAB | Should scale with average condensate diameter |
| `minSize` / `maxSize` | Step 2 MATLAB | Pixel area range for valid condensates |

---


---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Author

**Debalina Datta**
