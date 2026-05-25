"""
utils.py
Shared helper functions for the IRE1-mNG / ER-Tracker co-localization pipeline.

Channel convention (0-indexed):
    CH0 = ER-Tracker Red
    CH1 = FL-IRE1-mNG
    CH2 = NucBlue
"""

import numpy as np
import tifffile
from scipy import ndimage
from skimage.filters import gaussian, threshold_otsu
from skimage.morphology import remove_small_objects
from skimage.segmentation import watershed


# ---------------------------------------------------------------------------
# TIFF loading
# ---------------------------------------------------------------------------

def load_tiff(path):
    """
    Load a multi-channel TIFF and return a (C, Y, X) uint16 array.

    Handles both 3-D (C, Y, X) and 4-D (T, C, Y, X) inputs.
    For 4-D files, only the first time frame (T=0) is returned.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    np.ndarray, shape (C, Y, X), dtype uint16
    """
    arr = tifffile.imread(str(path))
    if arr.ndim == 4:
        arr = arr[0]          # take first time frame
    if arr.ndim != 3:
        raise ValueError(f"Unexpected array shape {arr.shape} for {path}")
    return arr.astype(np.uint16)


# ---------------------------------------------------------------------------
# Nucleus segmentation
# ---------------------------------------------------------------------------

def nucleus_segment(ch_nuc, min_size_px=500, sigma=1.0):
    """
    Segment nuclei from a single-channel nuclear stain image (NucBlue / DAPI).

    Steps:
        1. Gaussian smoothing (sigma=1 px)
        2. Otsu global threshold
        3. Distance-transform watershed to separate touching nuclei
        4. Remove objects smaller than min_size_px

    Parameters
    ----------
    ch_nuc : np.ndarray, shape (Y, X)
        Nuclear stain channel (uint16).
    min_size_px : int
        Minimum nucleus area in pixels. Objects smaller than this are discarded.
    sigma : float
        Gaussian smoothing sigma in pixels.

    Returns
    -------
    labels : np.ndarray, shape (Y, X), dtype int32
        Integer label array; 0 = background, 1..N = individual nuclei.
    """
    smoothed = gaussian(ch_nuc.astype(np.float32), sigma=sigma)
    thresh   = threshold_otsu(smoothed)
    binary   = smoothed > thresh

    # Distance-transform watershed
    distance = ndimage.distance_transform_edt(binary)
    local_max = distance > (0.3 * distance.max())
    markers, _ = ndimage.label(local_max)
    labels = watershed(-distance, markers, mask=binary)

    # Remove small objects
    for lbl in np.unique(labels):
        if lbl == 0:
            continue
        if np.sum(labels == lbl) < min_size_px:
            labels[labels == lbl] = 0

    # Re-label contiguously
    labels, _ = ndimage.label(labels > 0)
    return labels.astype(np.int32)


# ---------------------------------------------------------------------------
# Background subtraction
# ---------------------------------------------------------------------------

def bg_subtract(channel, cyto_mask_all):
    """
    Subtract a per-image background estimate from a single channel.

    Background is defined as the median intensity of all cytoplasm pixels
    that do not belong to any segmented cell (i.e., pixels in the image
    field that are outside all cell masks but within the cytoplasm region).
    If no background pixels are available, the global image median is used.

    Parameters
    ----------
    channel : np.ndarray, shape (Y, X)
        Single fluorescence channel (uint16 or float).
    cyto_mask_all : np.ndarray, shape (Y, X), bool
        Boolean mask marking all segmented cytoplasm pixels across all cells.

    Returns
    -------
    np.ndarray, shape (Y, X), dtype float32
        Background-subtracted channel, clipped to >= 0.
    """
    bg_pixels = channel[~cyto_mask_all]
    if bg_pixels.size > 0:
        bg_val = np.median(bg_pixels.astype(np.float32))
    else:
        bg_val = np.median(channel.astype(np.float32))

    result = channel.astype(np.float32) - bg_val
    return np.clip(result, 0, None)


# ---------------------------------------------------------------------------
# MOC computation
# ---------------------------------------------------------------------------

def compute_moc(ch1_pixels, ch2_pixels):
    """
    Compute the Manders' Overlap Coefficient (MOC) for two co-registered
    fluorescence channels within a region of interest.

    MOC = sum(CH1 * CH2) / sqrt(sum(CH1^2) * sum(CH2^2))

    MOC ranges from 0 (no overlap) to 1 (complete overlap).

    Parameters
    ----------
    ch1_pixels : np.ndarray, 1-D
        Pixel intensities for channel 1 (ER-Tracker Red) within the ROI.
    ch2_pixels : np.ndarray, 1-D
        Pixel intensities for channel 2 (FL-IRE1-mNG) within the ROI.

    Returns
    -------
    float or np.nan
        MOC value, or np.nan if the denominator is zero.
    """
    num   = np.sum(ch1_pixels * ch2_pixels)
    denom = np.sqrt(np.sum(ch1_pixels ** 2) * np.sum(ch2_pixels ** 2))
    if denom == 0:
        return np.nan
    return float(num / denom)
