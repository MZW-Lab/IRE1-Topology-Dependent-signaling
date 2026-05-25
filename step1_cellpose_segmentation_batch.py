#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  3 14:41:58 2025

@author: debalina datta
"""

import os
from cellpose import models, io
from PIL import Image
import tifffile
import numpy as np
import torch
from concurrent.futures import ProcessPoolExecutor

# --- Parameters ---
input_folder = './cyto_K907A_slices'  # Path to the folder containing the images
output_folder = './cyto_K907A_masks'  # Output folder for segmented images
model_type = 'CP_optoirek907a'  # Cellpose 2.0 model
diameter = 100 #diameter as estimated by cellpose training
flow_threshold = 0.4
cellprob_threshold = 0.0

# --- Ensure output directory exists ---
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# --- GPU Check ---
use_gpu = torch.cuda.is_available()  # Check if GPU is available

if use_gpu:
    print("GPU is available. Running on GPU.")
else:
    print("GPU is not available. Running on CPU.")

# --- Initialize Cellpose Model ---
model = models.CellposeModel(gpu=use_gpu, model_type=model_type)

# --- Function to Process a Single Image ---
def process_image(filename):
    input_path = os.path.join(input_folder, filename)
    
    # Skip non-image files
    if not filename.lower().endswith(('.png', '.jpg', '.jpeg', '.tif', '.tiff')):
        print(f"Skipping non-image file: {filename}")
        return

    print(f"Processing {filename}...")

    # Load the image
    try:
        img = io.imread(input_path)
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return

    # Segment the image
    try:
        masks, flows, _ = model.eval(
            img,
            diameter=diameter,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
            channels=[0, 0]  # Assume grayscale images
        )
    except Exception as e:
        print(f"Error during segmentation of {filename}: {e}")
        return

    # Save the segmented mask as a TIFF file
    output_path = os.path.join(output_folder, f"{os.path.splitext(filename)[0]}_mask.tif")
    try:
        tifffile.imwrite(output_path, masks.astype(np.uint16))  # Save as 16-bit TIFF
        print(f"Saved segmented mask to {output_path}")
    except Exception as e:
        print(f"Error saving mask for {filename}: {e}")

# --- Process All Images in Parallel ---
if __name__ == "__main__":
    with ProcessPoolExecutor() as executor:
        executor.map(process_image, os.listdir(input_folder))

print("Segmentation completed for all images.")