function A = histogramEqualize(inIm)
    
    hist = histogram(inIm, 256, 0, 255);
    %bar(1:256, hist);
    maxIntensity = 256;
    imSize = numel(inIm);
    
    %Creating PDF.
    inPDF = createPdf(hist, imSize);

    %calculate cdf.
    cdf = double(zeros(1, maxIntensity));
    previousCDF = double(0.0);
    for i=1:size(cdf, 2)
        cdf(1,i) = previousCDF + inPDF(1,i);
        previousCDF = cdf(1,i);
    end    
    
    s = cdf * double(maxIntensity-1);    
    %Rounding the cdf, s.
    for i=1:size(s, 2)
        s(1,i) = round(s(1,i));
    end
    
    A = s;
end

%Remove unnecessary doubles later.
function B = createPdf(inHist, inSize)
    outHist = double(zeros(size(inHist)));
    for i=1:size(inHist,2)
        outHist(1,i) = double(inHist(1,i) / double(inSize));
    end
    B = outHist;
end