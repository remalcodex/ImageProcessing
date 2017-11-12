clc;
clear all;
close all;
imtool close all;
addpath(genpath('Functions/'));
Im = imread('Inputs/LakePowell1G.jpg');

%Step 1: Image Colorization.
%Step 1.1-Histogram Equalization.
histCDF = histogramEqualize(Im);
fig1 = figure(1);
bar(histCDF);
xlabel('Intensity')
ylabel('cdf')
saveas(fig1, 'Outputs/hitogramCDF', 'png')

%Step 1.2- Intensity transoformation.
equalizedIm = intensityTransformation(Im, histCDF);
equalizedImS = uint8(equalizedIm);
imwrite(equalizedImS, 'Outputs/EqualizedeImage.png');
histEqualized = histogram(equalizedIm, 256, 0, 255);
fig2 = figure(2);
bar(histEqualized)
title('Equalized Histogram')
xlabel('Intensity value')
ylabel('Correspoding Intensity')
saveas(fig2, 'Outputs/equalizedHitogram1', 'png');

%Step 1.3 -Inverse Histogram Equalization.
mapping = inverseHistogramEq(Im);
fig3 = figure(3);
bar(mapping);
title('Mapping from equalized to grayscale image.');
xlabel('Equalized intensity');
ylabel('Grascale intensity');
saveas(fig3, 'Outputs/mapping', 'png');

%Step 1.4-Inverse Histogram equalizaiton and coloring the image.
Im = imread('Inputs/LakePowell1G.jpg');
ImRGB = imread('Inputs/LakePowell2.jpg');
newImRGB = colorizeImage(Im, ImRGB);
imwrite(newImRGB, 'Outputs/ColoredImage.png');

%Part 2: Counting cells
lakeFileName = 'Inputs/LakePowell1G.jpg';
%Step 2.1-Otsu Threshold
%Passing the filename because iminfo needs filename to read the file info. 
%iminfo is then used to get the bit depth of the image.
otsuIm = uint8(otsuThreshold(lakeFileName));
imwrite(otsuIm,'Outputs/otsuThreshold.png');

%Step 2.2-Adaptive Thesholding.
adaptiveIm = adaptiveThresholding(lakeFileName,[200,200],1000);%334/2,770/2,1000
adaptiveIm = uint8(adaptiveIm);
imwrite(adaptiveIm,'Outputs/adaptiveThreshold.png');

%Step 2.3-Counting cells.
cellFileName = 'Inputs/CellImage.tif';
countCells(cellFileName);

%Remove close all if you wish to see intermediate outputs.
close all;