#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import imageio
import pandas as pd

# Exercise 1: Loading Image Data
# Define constants
GENES = ["APEX1", "PIM2", "POLR2B", "SRSF1"]
FIELDS = ["field0", "field1"]
CHANNELS = ["DAPI", "PCNA", "nascentRNA"]
IMAGE_PATH = "/Users/cmdb/Desktop/week10/image/"  # Ensure path ends with "/"
OUTPUT_FILE = "Nuclei_data.csv"
DEFAULT_IMAGE_SHAPE = (1024, 1024)  # Default shape if no images are found

# Dictionary to store images
images = {}

# Load images into a dictionary
for gene in GENES:
    for field in FIELDS:
        for channel in CHANNELS:
            file_path = os.path.expanduser(f"{IMAGE_PATH}{gene}_{field}_{channel}.tif")
            if os.path.exists(file_path):
                try:
                    img = imageio.v3.imread(file_path).astype(np.uint16)
                    images[f"{gene}_{field}_{channel}"] = img
                except Exception as e:
                    print(f"Error loading image {file_path}: {e}")
            else:
                print(f"File not found: {file_path}")

# Determine image shape
IMAGE_SHAPE = next(iter(images.values()), np.zeros(DEFAULT_IMAGE_SHAPE)).shape

# Combine channels into multi-channel arrays
image_arrays = []
for gene in GENES:
    for field in FIELDS:
        combined_img = np.zeros((*IMAGE_SHAPE, len(CHANNELS)), dtype=np.uint16)
        for idx, channel in enumerate(CHANNELS):
            key = f"{gene}_{field}_{channel}"
            if key in images:
                combined_img[:, :, idx] = images[key]
        image_arrays.append(combined_img)

# Exercise 2: Cell Identification and Filtering

# Step 2.1: Create binary masks from the DAPI channel
masks = [(img[:, :, 0] >= np.mean(img[:, :, 0])) for img in image_arrays]

# Step 2.2: Label the masks
def find_labels(mask):
    from scipy.ndimage import label
    labeled, _ = label(mask)
    return labeled

label_arrays = [find_labels(mask) for mask in masks]

# Step 2.3: Filter labels by size
def filter_by_size(labels, min_size=100, max_size=None):
    sizes = np.bincount(labels.ravel())
    max_size = max_size or sizes.max()
    for label_id, size in enumerate(sizes):
        if size < min_size or size > max_size:
            labels[labels == label_id] = 0
    return labels

label_arrays = [filter_by_size(labels) for labels in label_arrays]

# Further filter labels using mean ± SD
def filter_by_size_msd(labels):
    sizes = np.bincount(labels.ravel())[1:]  # Exclude background (label 0)
    if len(sizes) == 0:
        return labels  # If no nuclei are labeled, return unchanged labels
    mean_size = np.mean(sizes)
    std_size = np.std(sizes)
    lower, upper = mean_size - std_size, mean_size + std_size
    return filter_by_size(labels, min_size=lower, max_size=upper)

label_arrays = [filter_by_size_msd(labels) for labels in label_arrays]

# Exercise 3: Signal Measurement
# Collect nucleus signal measurements
data = {
    "gene": [],
    "field": [],
    "nucleus": [],
    "PCNA": [],
    "nascentRNA": [],
    "log2_ratio": []
}

for idx, labels in enumerate(label_arrays):
    gene = GENES[idx // len(FIELDS)]
    field = FIELDS[idx % len(FIELDS)]
    for nucleus_id in range(1, labels.max() + 1):  # Exclude background (label 0)
        mask = (labels == nucleus_id)
        pcna_signal = np.mean(image_arrays[idx][:, :, 1][mask]) if np.any(mask) else 0
        nascent_rna_signal = np.mean(image_arrays[idx][:, :, 2][mask]) if np.any(mask) else 0
        log2_ratio = (
            np.log2(nascent_rna_signal / pcna_signal)
            if pcna_signal > 0 else 0
        )
        # Append data
        data["gene"].append(gene)
        data["field"].append(field)
        data["nucleus"].append(nucleus_id)
        data["PCNA"].append(pcna_signal)
        data["nascentRNA"].append(nascent_rna_signal)
        data["log2_ratio"].append(log2_ratio)

# Convert to DataFrame and save
nuclei_data = pd.DataFrame(data)
nuclei_data.to_csv(OUTPUT_FILE, index=False)
print(f"Data saved to {OUTPUT_FILE}")
