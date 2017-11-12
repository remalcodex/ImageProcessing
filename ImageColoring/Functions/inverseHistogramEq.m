function A = inverseHistogramEq(inImGray)
    maxIntensity = 256;
    sGray = histogramEqualize(inImGray);   
   
    %Equalizing image.
    grayImEqualized = intensityTransformation(inImGray, sGray); 
    sGrayEq = histogramEqualize(grayImEqualized);
    
    %Mapping cdf of equlized image to cdf of original image.
    mapping = zeros(1, maxIntensity);
    for i=1:size(sGrayEq,2)
        val = sGrayEq(1,i);
        closestJ = 1;
        closestVal = sGray(1,closestJ);
        for j=1:size(sGray,2)
            if abs(closestVal - val) > abs(sGray(1,j) - val)
                closestJ = j;
                closestVal = sGray(1,j);
            end
        end
        mapping(1,i) = closestJ;
    end
    
    A = mapping;
    
    %Reproducing the original grayscale image.
    ImCheck = zeros(size(grayImEqualized), 'uint8');
    for i=1:size(grayImEqualized, 1)
        for j=1:size(grayImEqualized,2)
            ImCheck(i,j) = mapping((grayImEqualized(i,j)+1));
        end
    end
    imwrite(ImCheck, 'Outputs/reproducedImage.png');
end