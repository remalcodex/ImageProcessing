function A = subtractImages(Im1, Im2)
    outIm = zeros(size(Im1));
    for i=1:size(Im1,1)
        for j=1:size(Im2,2)
            outIm(i,j) = abs(Im1(i,j) - Im2(i,j));
        end
    end
    A = outIm;
end