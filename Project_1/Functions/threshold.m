function A = threshold(Im, minT, maxT)
    A = zeros(size(Im));
    for i = 1:size(Im, 1)
        for j = 1:size(Im, 2)
            if (Im(i,j) >= minT && Im(i,j) <= maxT)
                A(i,j) = 256;
        end
    end
end