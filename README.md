# LidR_EcoTraits: Comprehensive Crown Metrics Analysis for Tree Canopies Using Lidar Data

## Overview
This R script aims to facilitate an extensive analysis of crown metrics using LiDAR data, focusing on tree canopies. The functionalities encompass package installations, LiDAR data reading, preprocessing, and the computation of various metrics such as canopy ratios, Effective Number of Layers (ENL), box dimensions, alpha shapes, among others.


## Context and Background
Originally developed as part of a PhD thesis, this script has been documented and published in a specific academic publication. Its primary purpose is to serve as an initial framework for subsequent projects, offering insights and guidance on deriving ecologically significant traits for individual trees.

## Installation
To ensure all required packages are installed, run the following R code:
```{r}
# List of packages
packages <- c("lidR", "lidRmetrics", "Rdimtools", "Lmoments", "fitdistrplus", 
              "alphashape3d", "tidyr", "purrr", "concaveman", "sf", "viewshed3d", 
              "pracma", "ggplot2", "fsbrain", "rayshader", "rgl", "terra")

# Function to check and install packages
check_and_install <- function(pkg){
  if(!require(pkg, character.only = TRUE)){
    install.packages(pkg, dependencies = TRUE)
  }
}

# Apply the function to the list of packages
sapply(packages, check_and_install)
```
Note: Some specific packages like Rdimtools and viewshed3d may necessitate the developer version available on GitHub or additional dependencies. Refer to respective package documentation for detailed instructions.

## How to Use
1. Data Preparation: Ensure your lidar data is in .laz or .las format. The script reads a sample file named MixedConifer.laz from the lidR package.

2. Run Script: Execute the entire R script. This will:

- Load necessary libraries.
- Read the LiDAR data.
- Perform data cleaning and filtering.
- Calculate various crown metrics.
- Save the computed metrics to an output CSV file (`output.csv``).

3. Visualisation: The script includes visualisation functions to plot various metrics such as tree visibility, leaf area index, and more.

## Functions overview
- 'canopy_ratio': Calculates the canopy ratio based on height percentiles.
- 'ENL': Computes the Effective Number of Layers (ENL) for 2D data.
- 'bd_fun': Estimates the Box Dimension of the point cloud.
- 'ashape': Computes the Alpha shape metrics.
- 'canopyShape': Calculates 2D shape metrics based on Owen et al. (2021).

## Output
The script produces an output CSV file (output.csv) containing the computed metrics for each tree in the LiDAR dataset.


## Authors and acknowledgment
The script credits its authors and acknowledges package contributors. A special acknowledgment is extended to Francois Du Toit for their contribution towards the alpha shape function.

## MIT License

Copyright (c) 2024 Leonard Hambrecht

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


## Project status
This project is archival in nature and is not actively maintained. It is intended to act as a foundational resource for subsequent projects, providing guidance on deriving ecologically meaningful traits for individual trees. Users are cautioned that the functions may become obsolete or break due to updates in associated packages. Last successfully executed in September 2023.
