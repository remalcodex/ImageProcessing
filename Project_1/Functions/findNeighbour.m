function A = findNeighbour(Im, ImInfo)

    seed = [ImInfo(1,1), ImInfo(1,2)];
    if Im(seed(1,1), seed(1,2)) == 0
        disp('Couldnot find neighbour');
        A = -1;
        return;
    end
    
    neighbourPixels = zeros(1,256);
    
    stack = [seed(1,1), seed(1,2)];
    while ~(isempty(stack))
        i = stack(1,1);
        j = stack(1,2);        
        stack(1, :) = [];       
        
        if (boundCheck(Im, [i, j]) == 0)
            A = labelOut;
            B = area;
            continue;
        end
        
        if labelOut(i, j) ~= 0
            A = labelOut;
            B = area;
            continue;
        end
        
        seed1 = [i, j] + [+1, 0];
        seed2 = [i, j] + [-1, 0];
        seed3 = [i, j] + [0, +1];
        seed4 = [i, j] + [0, -1];
        
        labelOut(i, j) = labelVal;
        area = area + 1;
        
        if (boundCheck(Im,seed1)~=0) 
            pixVal = Im(seed1(1,1), seed1(1,2));
            if (Im(i, j) == pixVal)
                stack = [[seed1(1,1), seed1(1,2)]; stack];
            else
                neighbourPixels(1, pixVal) = neighbourPixels(1, pixVal) + 1;
            end
        end
        if (boundCheck(Im,seed2)~=0)
            pixVal = Im(seed2(1,1), seed2(1,2));
            if (Im(i, j) == pixVal)
                stack = [[seed2(1,1), seed2(1,2)]; stack];
            else
                neighbourPixels(1, pixVal) = neighbourPixels(1, pixVal) + 1;
            end
        end
        if (boundCheck(Im,seed3)~=0) 
            pixVal = Im(seed3(1,1), seed3(1,2));
            if (Im(i, j) == pixVal)
                stack = [[seed3(1,1), seed3(1,2)]; stack];
            else
                neighbourPixels(1, pixVal) = neighbourPixels(1, pixVal) + 1;
            end
        end
        if (boundCheck(Im,seed4)~=0)
            pixVal = Im(seed4(1,1), seed4(1,2));
            if (Im(i, j) == pixVal)
                stack = [[seed4(1,1), seed4(1,2)]; stack];
            else
                neighbourPixels(1, pixVal) = neighbourPixels(1, pixVal) + 1;
            end
        end
    end
 
    A = labelOut; 
    B = area;
    C = max(neighbourPixels);
end