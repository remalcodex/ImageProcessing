function [A, B, C, D] = floodFill(Im, seed, labelVal, labelOut)    

    area = 0;
    %Return if seed value is 0;
    if Im(seed(1,1), seed(1,2)) == 0
        disp('Couldnot run flood fill');
        A = labelOut;
        B = area;
        return;
    end
    
    neighbourList = [];
    visitedNodes = zeros(size(Im));
    
    stack = [seed(1,1), seed(1,2)];
    while ~(isempty(stack))
        i = stack(1,1);
        j = stack(1,2);        
        stack(1, :) = [];       
        
        if (boundCheck(Im, [i, j]) == 0)
            A = labelOut;
            B = area;
            D = visitedNodes;
            continue;
        end
        
        if visitedNodes(i,j) ~= 0
            A = labelOut;
            B = area;
            D = visitedNodes;
            continue;
        end
        
        seed1 = [i, j] + [+1, 0];
        seed2 = [i, j] + [-1, 0];
        seed3 = [i, j] + [0, +1];
        seed4 = [i, j] + [0, -1];
        
        visitedNodes(i,j) = 'V'; %Marking the current pixel as visited.
        labelOut(i, j) = labelVal;
        area = area + 1;
        
        if (boundCheck(Im,seed1)~=0) 
            pixVal = Im(seed1(1,1), seed1(1,2));
            if (Im(i, j) == pixVal)
                stack = [[seed1(1,1), seed1(1,2)]; stack];
            else
                neighbourList = addNeighbour(neighbourList, pixVal);
            end
        end
        if (boundCheck(Im,seed2)~=0)
            pixVal = Im(seed2(1,1), seed2(1,2));
            if (Im(i, j) == pixVal)
                stack = [[seed2(1,1), seed2(1,2)]; stack];
            else
                neighbourList = addNeighbour(neighbourList, pixVal);
            end
        end
        if (boundCheck(Im,seed3)~=0) 
            pixVal = Im(seed3(1,1), seed3(1,2));
            if (Im(i, j) == pixVal)
                stack = [[seed3(1,1), seed3(1,2)]; stack];
            else
                neighbourList = addNeighbour(neighbourList, pixVal);
            end
        end
        if (boundCheck(Im,seed4)~=0)
            pixVal = Im(seed4(1,1), seed4(1,2));
            if (Im(i, j) == pixVal)
                stack = [[seed4(1,1), seed4(1,2)]; stack];
            else
                neighbourList = addNeighbour(neighbourList, pixVal);
            end
        end
    end
 
    A = labelOut; 
    B = area;
    
    %Taking the maximum of the neighbours.
    maxVal = max(neighbourList(:,2));
    maxIndex = find(neighbourList(:,2) == maxVal);
    maxIndex = maxIndex(1,1);
    C = neighbourList(maxIndex, 1);
    
    D = visitedNodes;
end

%This checks if the image is in the bounds of the image.
function Check = boundCheck(Im, seed)
    if ((seed(1,1) <= 0) || (seed(1,1) > size(Im, 1)))
        Check = 0;
        return;
    end
    
    if ((seed(1,2) <= 0) || (seed(1,2) > size(Im, 2)))
        Check = 0;
        return;
    end
    
    Check = 1;
end

% This function adds the val to the neighbour list. It either increments
% the count of the label already present in the list or creates a new entry
% with count as 1.
function A = addNeighbour(neighbourList, val)
    if isempty(neighbourList)
        neighbourList = [neighbourList; val 1];
        A = neighbourList;
        return;
    end

    values = (neighbourList(:, 1) == val);
    if max(values) == 0
        neighbourList = [neighbourList; val, 1];
    else
        index = find(values == max(values));
        neighbourList(index, 2) = neighbourList(index, 2) + 1;
    end
    A = neighbourList;
end