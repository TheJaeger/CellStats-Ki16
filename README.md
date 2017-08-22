# CellStats-Ki16
MATLAB scripts to generate MST-based statistics on cells in control vs. diseased samples

## Files
The main folder contains the following files:
* **__main.m__**: Main script to generate data and create global comparisons
* **__getData.m__**: Functin script which extracts cells and extracts parameters from each image
* **__ROIdetection.m__**: Extracts ROI from image to create a mask
* **__distancePoints.m__**: Part of MatGeom toolbox to compute pair-wise Euclidean distances
* **__drawGraph.m__**: Plots MST and nodes

## Running Instructions
Running __main.m__ in MATLAB will prompt the the user to select folders containing images of control and diseases samples
Upon confirming folders, the script analyses each images and output all data into the **Output** folder.



