function countCells(cellFileName)
    Im = imread(cellFileName);
    %Applying otsu thresholding.
    otsuIm = otsuThreshold(cellFileName);
    otsuImS = uint8(otsuIm);
    imwrite(otsuImS, 'Outputs/cellsOtsuThreshold.png');
    
    %Applying connected components to otsu image.
    labeledIm = zeros(size(otsuIm));
    [labeledIm, imInfo] = connectedComponent(otsuIm, 1, labeledIm);
    labeledImS = label2rgb(labeledIm);
    imwrite(labeledImS, 'Outputs/cellsOtsuCCLImage.png');
     
    %Denoising the otsu image.    
    [denoisedIm, denoisedImInfo] = denoising(labeledIm, imInfo, 10);
    finalImS = label2rgb(denoisedIm);
    imwrite(finalImS, 'Outputs/cellsOtsuDenoisedImage.png');
    disp(strcat('Otsu cell count:',  num2str(size(denoisedImInfo,1))));
    
    %---------------------------------------------------------------------%
    %Applying adaptive thresholding.
    adaptiveIm = adaptiveThresholding(cellFileName, [60,60], 200);
    adaptiveImS = uint8(adaptiveIm);
    imwrite(adaptiveImS, 'Outputs/cellsAdaptiveThreshold.png');
    
    %Applying connected components to adaptive image.
    labeledIm = zeros(size(adaptiveIm));
    [labeledIm, imInfo] = connectedComponent(adaptiveIm, 1, labeledIm);
    labeledImS = label2rgb(labeledIm);
    imwrite(labeledImS, 'Outputs/cellsAdaptiveCCLImage.png');
    
    %Denoising the adaptive thresholding image.
    [denoisedIm, denoisedImInfo] = denoising(labeledIm, imInfo, 10);
    finalImS = label2rgb(denoisedIm);
    imwrite(finalImS, 'Outputs/cellsAdaptiveDenoisedImage.png');
    disp(strcat('Adaptive Threshold cell count:',  num2str(size(denoisedImInfo,1))));
end