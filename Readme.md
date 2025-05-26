"# Spatial Data Splitting Tool" 
"R script for balanced spatial data splitting using cluster-based methods." 
# Spatial Data Splitting Tool

This R script provides a balanced cluster-based splitting method for spatial data, ensuring proper spatial distribution between training and testing datasets.

## Features

- **Cluster-based splitting**: Uses k-means clustering to maintain spatial representativeness
- **Balanced distribution**: Ensures even distribution of test points across all spatial clusters
- **Spatial metrics**: Calculates intra-group and inter-group distances
- **Visualization**: Generates plots showing the spatial distribution of train/test splits
- **Multiple outputs**: Saves results in TXT, PNG, TIFF, and PDF formats

## Requirements

- R (version 4.0 or higher)
- Required packages:
  - FNN
  - ggplot2
  - svglite
  - dplyr
  - terra
  - sf

## Usage

1. Prepare your data frame with coordinates in columns named 'x' and 'y'
2. Adjust parameters in the `params` list as needed
3. Run the script

## Outputs

- `balanced_cluster_split_plot.png`/.tiff/.pdf: Visualization of the split
- `trainAria3.txt` and `testAria3.txt`: Split datasets
- `splitting_report.txt`: Detailed metrics about the split

## Parameters

Key adjustable parameters:
- `test_size`: Fraction of data for testing (default 0.2)
- `n_clusters`: Number of spatial clusters (default 12)
- `max_points_per_cluster`: Maximum test points per cluster (default 2)
- Various weights for optimization criteria

## Example

```r
# Load your spatial data
pts <- read.csv("your_data.csv")

# Ensure it has x and y columns
names(pts)[names(pts) %in% c("longitude", "lat")] <- c("x", "y")

# Run the splitting function
split_result <- balanced_cluster_split(pts, params)