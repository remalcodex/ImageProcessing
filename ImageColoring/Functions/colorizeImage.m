function A = colorizeImage(ImG, ImRGB)
    histEqualized = histogramEqualize(ImG);
    invHistEqualizationR = histogramEqualize(ImRGB(:,:,1));%Red channel
    invHistEqualizationG = histogramEqualize(ImRGB(:,:,2));%Green channel
    invHistEqualizationB = histogramEqualize(ImRGB(:,:,3));%Blue channel
    %invHistEqualization = cat(3, invHistEqualizationR, invHistEqualizationG, invHistEqualizationB);
    newIm = zeros([size(ImG),3], 'uint8');

    for ch=1:3
        %Picking the right histogram.
        invHistEqualizationA = [];
        if ch==1
            invHistEqualizationA = invHistEqualizationR;
        elseif ch==2
            invHistEqualizationA = invHistEqualizationG;
        elseif ch==3
            invHistEqualizationA = invHistEqualizationB;
        end
        
        for i=1:size(ImG, 1)
            for j=1:size(ImG,2)        
                sk = histEqualized(1,(ImG(i,j)+1));
                smallestZq = 1;
                smallestGzq = invHistEqualizationA(1,smallestZq);
                %Calculating the closest value to sk.
                for q=1:size(invHistEqualizationA,2)
                    if sk - smallestGzq > abs(sk - invHistEqualizationA(1,q))
                        smallestGzq = invHistEqualizationA(1,q);
                        smallestZq = q;                        
                    end
                end
                newIm(i,j,ch) = uint8(smallestZq); %Adding for RGB 8 bits each.
            end
        end
    end
    
    A = newIm;
end