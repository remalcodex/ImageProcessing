function A = intensityTransformation(inIm, inEqualizedHist)
    newIm = zeros(size(inIm));
    for i=1:size(inIm,1)
        for j=1:size(inIm,2)
            newIm(i,j) = inEqualizedHist(1,(inIm(i,j)+1));
        end
    end
    A= newIm;
end