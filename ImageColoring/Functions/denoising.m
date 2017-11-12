function [finalIm, finalImInfo] = denoising(Im, ImInfo, threshold)
    
    finalIm = Im;
    newImInfo = [];
    for k = 1:size(ImInfo, 1)
        if ImInfo(k,3) <= threshold
            %Finding the maximum neighbour.
            dummyOut = zeros(size(Im));
            maxNeighbour = -1;
            [dummyOut, area, maxNeighbour] = floodFill(Im, ImInfo(k, 1:2), ImInfo(k, 4), dummyOut);
            ImInfo(k, 5) = maxNeighbour;
            newImInfo = [newImInfo; ImInfo(k, :)];
        end
    end   
    
    %Since we know the maximum neighbour, we run floodfill again to fill in
    %the pixels with new label.    
    for k = 1:size(newImInfo)
        dummyOut = zeros(size(Im));
        [dummyOut, area, maxNeighbour, visitedNodes] = floodFill(Im, newImInfo(k, 1:2), newImInfo(k, 5), dummyOut);
        
        %Placing the new label to only those pixels that were visited.
        for i = 1:size(finalIm,1)
            for j = 1:size(finalIm, 2)
                if visitedNodes(i,j) ~= 0
                    finalIm(i,j) = newImInfo(k, 5);
                end
            end
        end
        
        %Removing the info from main info.
        ImInfo(newImInfo(k,4)-(k-1),:) = [];
    end   
    finalImInfo = ImInfo;
end