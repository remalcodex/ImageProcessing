function A = histogram(Im, nBins, minI, maxI)
    Im = double(Im);
    noBinValues = ceil((maxI-minI+1)/nBins);
    A = zeros(1, nBins);
    for i = 1:size(Im,1)
        for j = 1:size(Im,2)
            if ((Im(i,j) >= minI) && (Im(i,j) <= maxI))
                pixVal = Im(i,j) + 1 - minI; %(Normalizing): Since if min is 200 then 200->x will be bin 1.
                binNo = ceil(pixVal/noBinValues);
                A(binNo) = A(binNo)+1;
            end
        end
    end    
end