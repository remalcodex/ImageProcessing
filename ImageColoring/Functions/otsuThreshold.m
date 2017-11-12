function A = otsuThreshold(cellFileName)
    Im = imread(cellFileName);
    info = imfinfo(cellFileName);    
    maxIntensity = 2^(info.BitDepth);
    
    %Creating pdf.    
    hist = histogram(Im, maxIntensity, 0, maxIntensity-1);
    imSize = double(numel(Im));
    pdf = zeros(size(hist), 'double');
    for i=1:size(pdf,2)
        pdf(1,i) = double(hist(1,i) / imSize);
    end
    
    %Calculating otsu threshold.
    maxVariance = -1;
    maxThreshold = -1;
    for threshold=1:size(hist,2)
        [pdfWeight1, pdfWeight2] = calculateWeight(pdf, threshold);
        [pdfMean1, pdfMean2] = calculateMeans(pdf, threshold, pdfWeight1, pdfWeight2);
        variance = pdfWeight1*pdfWeight2*((pdfMean1-pdfMean2)^2);
        if variance > maxVariance
            maxVariance = variance;
            maxThreshold = threshold;
        end
    end    
    
    newIm = zeros(size(Im));
    for i=1:size(Im,1)
        for j=1:size(Im,2)
            if Im(i,j) >= maxThreshold
                newIm(i,j) = 255;            
            end
        end
    end
    A = newIm;
end

function [A,B] = calculateWeight(inPDF, inThresh)
    weight1 = 0;
    weight2 = 0;
    for i=1:size(inPDF,2)
        if i < inThresh
            weight1 = weight1 + inPDF(1,i);
        else
            weight2 = weight2 + inPDF(1,i);
        end
    end    
    A = weight1;
    B = weight2;
end

function [A,B] = calculateMeans(inPDF, inThresh, inWeight1, inWeight2)
    mean1 = 0;
    mean2 = 0;
    for i=1:size(inPDF,2)
        if i < inThresh
            mean1 = mean1 + (i*inPDF(1,i));
        else
            mean2 = mean2 + (i*inPDF(1,i));
        end
    end
    mean1 = mean1/inWeight1;
    mean2 = mean2/inWeight2;
    A = mean1;
    B = mean2;
end