function [A, B] = connectedComponent(Im, startLabel, labeledOut)
    labelNo = startLabel - 1;
    labelInfo = [];
    
    %visitedNodesG keeps track of all the pixels that have been visited.
    visitedNodesG = zeros(size(Im));
    for i = 1:size(Im, 1)
        for j = 1:size(Im, 2)
            pixVal = Im(i,j);
            if (pixVal ~= 0) && (visitedNodesG(i,j) == 0)
                labelNo = labelNo + 1;
                seed = [i, j];
                [labeledOut, area, maxNeighbor, visitedNodes] = floodFill(Im, seed, labelNo, labeledOut);             
                visitedNodesG = visitedNodesG + visitedNodes;
                %labelInfo keeps track of connected component in the image.
                %With properties likee start position, area, labelNo,
                %maximum neighbour.
                labelInfo = [labelInfo ; i, j, area, labelNo, maxNeighbor];               
            end
        end
    end
    A = labeledOut;
    B = labelInfo;
end