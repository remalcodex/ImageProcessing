clc;
clear all;
close all;
imtool close all;
addpath(genpath('Functions/'));

Im = imread('Inputs/shapes_noise.tif');
%Im = imread('Inputs/sunset.png');
%Im = imread('Inputs/turkeys.tif');

%Step 1 -- Histogram
hist = histogram(Im, 256 , 0, 255);
f1 = figure (1);
bar(1:256, hist);
xlabel('Intensity');
ylabel('No of pixels');
saveas(f1, 'Outputs/histogram1' ,'png');

%Step 2.1 -- Threshold
thresholdIm = threshold(Im, 120, 256);
thresholdImS = mat2gray(thresholdIm);
%imtool(thresholdImS);
imwrite(thresholdImS, 'Outputs/Threshold1.png');

%Step 2.2 -- Flood Fill
seed = [262,103]; % Seed for shapes_noise and turkeys image.
%seed = [50,50]; % Seed for sunset image.
labelImage = zeros(size(thresholdIm));
[labelImage, area] = floodFill(thresholdIm, seed, 200, labelImage);
labelImageS = mat2gray(labelImage);
%imtool(labelImageS);
imwrite(labelImageS, 'Outputs/floodFillImage.png');

%Step 2.3 -- Connected Component.
ccLabelImage = zeros(size(thresholdIm));
ccLabelImage = connectedComponent(thresholdIm, 1, ccLabelImage);
ccLabelImageS = label2rgb(ccLabelImage);
%imtool(ccLabelImageS);
imwrite(ccLabelImageS, 'Outputs/connectedComponent.png');

%Step 3 -- Denoising
ImDenoise = imread('Inputs/brain.tif');
ImD1 = threshold(ImDenoise, 0, 45);
ImD2 = threshold(ImDenoise, 46, 65);
ImD3 = threshold(ImDenoise, 66,255);

% ImDenoise = imread('Inputs/sunset.png');
% ImD1 = threshold(ImDenoise, 0, 30);
% ImD2 = threshold(ImDenoise, 31, 110);
% ImD3 = threshold(ImDenoise, 111,255);

ccImD1 = zeros(size(ImD1));
ccImD2 = zeros(size(ImD2));
ccImD3 = zeros(size(ImD3));
[ccImD1, ccImD1Info] = connectedComponent(ImD1, 1, ccImD1);
[ccImD2, ccImD2Info] = connectedComponent(ImD2, size(ccImD1Info,1)+1, ccImD2);
[ccImD3, ccImD3Info] = connectedComponent(ImD3, size(ccImD1Info,1)+size(ccImD2Info,1)+1, ccImD3);

f2 = figure (2)
subplot(1,3,1);
imshow(ccImD1);
title('Threshold 0-45');
subplot(1,3,2);
imshow(ccImD2);
title('Threshold 46-65');
subplot(1,3,3);
imshow(ccImD3);
title('Threshold 66-255');
saveas(f2, 'Outputs/DenoiseThresholds', 'png');

newIm = ccImD1 + ccImD2 + ccImD3;
newImInfo = [ccImD1Info; ccImD2Info; ccImD3Info];
newImRGB=label2rgb(newIm);
%imtool(newImRGB);
imwrite(newImRGB, 'Outputs/BeforeDenoising.png');
    
finalIm = denoising(newIm, newImInfo, 200);
finalImRGB = label2rgb(finalIm);
%imtool(finalImRGB);
imwrite(finalImRGB, 'Outputs/AfterDenoising.png');

%Step 4 -- Motion Detection.
ImMotionStatic = imread('Inputs/houndog2S.png');
ImMotionDynamic = imread('Inputs/houndog1S.png');

ImMotionSub = subtractImages(ImMotionStatic, ImMotionDynamic);
grayIm = mat2gray(ImMotionSub);
%imtool(grayIm);
imwrite(grayIm, 'Outputs/MotionDetectionGray.png');

ImT1 = threshold(ImMotionSub, 25, 255);
%imtool(ImT1);
imwrite(ImT1, 'Outputs/MotionDetectionThreshold.png');

ccIm1 = zeros(size(ImT1));
[ccIm1, ccIm1Info] = connectedComponent(ImT1, 1, ccIm1);
ccIm1RGB = label2rgb(ccIm1);
%imtool(ccIm1RGB);
imwrite(ccIm1RGB, 'Outputs/MotionDetectionCCL.png');

finalIm = denoising(ccIm1, ccIm1Info, 10);
finalImRGB = label2rgb(finalIm);
%imtool(finalImRGB);
imwrite(finalImRGB, 'Outputs/MotionDetectionDenoised.png');

close all; %This is to close the two figures that pop up for histogram and threshold.