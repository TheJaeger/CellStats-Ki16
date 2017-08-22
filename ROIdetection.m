function ROI = ROIdetection(imageArray)
%   ROIdetection.m computes circular ROIs and stores them in the 3D
%   matrix 'ROI' defined by MxNxP. MxN represents the image size and P
%   represents the number of ROIs discovered.
%
%   ------------------------------------------------------------------
%   Author: Siddhartha Dhiman
%   e-mail: sdhiman@buffalo.edu

%%  Extract blue channel
IB = imbinarize(imageArray(:,:,3));
IB_fill = imfill(IB, 'holes');

%% Locate ROIs
regsIB = regionprops(IB_fill, 'Centroid', 'Area');

%% Find Radius of Each ROI
for r = 1:size(regsIB, 1)
    R(r) = sqrt(regsIB(r).Area/pi);
end

%% Reform ROIs
imageSizeX = size(IB,2);
imageSizeY = size(IB,1);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
for s = 1:size(regsIB, 1)
    ROI(:,:,s) = (rowsInImage - regsIB(s).Centroid(2)).^2 ... 
        + (columnsInImage - regsIB(s).Centroid(1)).^2 <= R(s).^2;
end

clearvars variables -except ROI





