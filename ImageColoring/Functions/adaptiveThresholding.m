function A = adaptiveThresholding(cellFileName, bSize, minVariance)
    Im = imread(cellFileName);
    info = imfinfo(cellFileName);    
    maxIntensity = 2^(info.BitDepth);
    
    nLoops = ceil(size(Im)./bSize);
    newIm = zeros(size(Im));
    for i=0:(nLoops(1,1)-1)
        for j=0:(nLoops(1,2)-1)
            %Calculating indices.
            startI = i*bSize(1,1)+1;
            startJ = j*bSize(1,2)+1;
            endI = startI + bSize(1,1)-1;
            endJ = startJ + bSize(1,2)-1;
            if endI > size(Im,1)
                endI = size(Im,1);
            end
            if endJ > size(Im,2)
                endJ = size(Im,2);
            end
            
            blockIm = Im(startI:endI, startJ:endJ);
            newIm(startI:endI, startJ:endJ) = thresholdBlock(blockIm, minVariance, maxIntensity);
        end
    end
    A = newIm;
end

function A = thresholdBlock(Im, minVariance, maxIntensity)
    hist = histogram(Im, maxIntensity, 0, maxIntensity-1);
    imSize = double(numel(Im));
    pdf = zeros(size(hist), 'double');
    for i=1:size(pdf,2)
        pdf(1,i) = double(hist(1,i) / imSize);
    end
    
    meanB = double(0.0);
    for i=1:size(pdf,2)
        meanB = meanB + (i*pdf(1,i));
    end
    
    varianceB = double(0.0);
    for i=1:size(pdf,2)
        varianceB = varianceB + ((i-meanB)^2)*pdf(1,i);
    end
    
    newIm = zeros(size(Im));
    if varianceB > minVariance
        for i=1:size(Im,1)
            for j=1:size(Im,2)            
                if Im(i,j)>(meanB+sqrt(varianceB))
                    newIm(i,j) = 255;
                end
            end
        end
    end    
    A = newIm;
end