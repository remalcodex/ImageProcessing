function phaseCorrelation(inputFolder)
    fileNames = dir(strcat(inputFolder, '*.png'));
    
    %Do preprecessing to make image sizes same.
    [images, imagesReal] = padImages(inputFolder, fileNames);
    imagesGray = [];
    
    %Creating low pass filtering.
    lpIm = images{1};
    lpSize = size(lpIm);
    %For 0 that is low pass filtering i found sigma to 70 the best for toy
    %input.
    lowFilter = createLowPassFilterG(lpSize(1,1:2), 1, 50, 1);
    %figure;imshow(lowFilter, [])
    
    %Making fourier images.
    fourierIm = [];
    for k = 1:numel(images)
        %Fourier transform
        Im = images{k};
        if numel(size(Im)) == 3
            Im = rgb2gray(Im);
        end
        Im = histeq(Im, 256);
        %figure;imshow(Im)
        Im = double(Im); %Dont forget to double the image before fourier transform.
        imagesGray{k} = Im;
        fIm = fft2(Im);
        FIm = fftshift(fIm);
        fourierIm{k} = FIm;
        %F = abs(F);
        %F = log(F+1);        
        %Applying low pass filtering.
        %G = FIm .* lowFilter; % Do this in next step when you are done with convolution.              
        %G = ifftshift(G); %Inverse Fourier transform.
        %G = ifft2(G);
        %figure;imshow(uint8(abs(G)))        
    end
    
    %Making the canvas ready with 1 image.
    imagesOnCanvas = [];
    firstIm = imagesGray{1};
    imageSize = size(firstIm);
    canvas = 0;
    if numel(size(images{1})) == 3
        canvas = zeros([imageSize(1,1:2)*numel(imagesGray)*2 - imageSize(1,1:2), 3], 'uint8');
    else
        canvas = zeros(imageSize*numel(imagesGray)*2 - imageSize, 'uint8');
    end
    cs = size(canvas);
    cX = cs(2) / 2;
    cY = cs(1) / 2;    
    y1 = cY - (imageSize(1)/2);
    y2 = cY + (imageSize(1)/2)-1;
    x1 = cX - (imageSize(2)/2);
    x2 = cX + (imageSize(2)/2)-1;
    canvas(y1:y2, x1:x2, :) = uint8(images{1});
    imagesOnCanvas = [1, y1, x1];
    
    %Checking the best mosaic with peeks now.
    for i = 1:(numel(fourierIm)-1)
        for j = i+1:numel(fourierIm)
            Ga = fourierIm{i};
            Gb = fourierIm{j};
            Gb = conj(Gb);
            G = Ga .* Gb;
            G = G ./ abs(G);
            G = G .* lowFilter;
            G = ifftshift(G);
            G = ifft2(G);
            realG = real(G);
            
            %Finding the peak.
            [maxIndex, maxVal] = findPeak(realG);            
            %mesh(realG)            
            
            %Checking where this image fits. i.e in which quadrant.
            if maxVal > 0.002
                imHeight = size(imagesGray{1}, 1);
                imWidth = size(imagesGray{1}, 2);
                xI = maxIndex(2);
                yI = maxIndex(1);
                
                Im1 = imagesGray{i};
                Im2 = imagesGray{j};
                maxArray = [0,0,0,0];
                
                q11 = Im1(1:yI, 1:xI);
                q12 = Im2((imHeight - yI + 1):end, (imWidth - xI + 1):end);
                %lpFilter = createLowPassFilter(size(q11));
                %maxVal1 = calculateCorrelation(q11, q12, lpFilter);
                maxVal11 = calculateDIC(q11, q12);
                maxArray(1) = maxVal11;
                
                q21 = Im1(1:yI, xI:end);
                q22 = Im2((imHeight - yI + 1):end, 1:(imWidth - xI + 1));
                %lpFilter = createLowPassFilter(size(q21));
                %maxVal2 = calculateCorrelation(q21, q22, lpFilter);
                maxVal21 = calculateDIC(q21, q22);
                maxArray(2) = maxVal21;
                
                q31 = Im1(yI+1:end, 1:xI);
                q32 = Im2(1:(imHeight - yI), (imWidth - xI + 1):end);
                %lpFilter = createLowPassFilter(size(q31));
                %maxVal3 = calculateCorrelation(q31, q32, lpFilter);
                maxVal31 = calculateDIC(q31, q32);
                maxArray(3) = maxVal31;
                
                q41 = Im1(yI+1:end, xI+1:end);
                q42 = Im2(1:(imHeight - yI), 1:(imWidth - xI));
                %lpFilter = createLowPassFilter(size(q41));
                %maxVal4 = calculateCorrelation(q41, q42, lpFilter);
                maxVal41 = calculateDIC(q41, q42);
                maxArray(4) = maxVal41;
                
                %Dont forget to add to the canvas list.
                %Try doing the conjigate reverse.
                [mVal, mIndex] = max(maxArray);
                iIndex = isPresentOnCanvas(imagesOnCanvas , i);
                jIndex = isPresentOnCanvas(imagesOnCanvas, j);
                %Do mosaicing.
                if mIndex == 1
                    if iIndex ~= 0 && jIndex == 0
                        imY = imagesOnCanvas(iIndex, 2);
                        imX = imagesOnCanvas(iIndex, 3);
                        newY = imY - (imageSize(1)-yI);
                        newX = imX - (imageSize(2)-xI);                        
                        canvas(newY:newY + size(imagesReal{j},1)-1, newX:newX + size(imagesReal{j},2)-1) = imagesReal{j};
                        imagesOnCanvas = [imagesOnCanvas; [j, newY, newX]];
                    end
                    if iIndex == 0 && jIndex ~= 0
                        imY = imagesOnCanvas(jIndex, 2);
                        imX = imagesOnCanvas(jIndex, 3);
                        newY = imY - (imageSize(1)-yI);
                        newX = imX - (imageSize(2)-xI);
                        canvas(newY:newY + size(imagesReal{i},1)-1, newX:newX + size(imagesReal{i},2)-1) = imagesReal{i};
                        imagesOnCanvas = [imagesOnCanvas; [i, newY, newX]];
                    end                    
                end
                if mIndex == 2
                    if iIndex ~= 0 && jIndex == 0
                        imY = imagesOnCanvas(iIndex, 2);
                        imX = imagesOnCanvas(iIndex, 3);
                        newY = imY - (imageSize(1)-yI);
                        newX = imX + xI;
                        %Add averaging later on. Maybe???
                        canvas(newY:newY + size(imagesReal{j},1)-1, newX:newX + size(imagesReal{j},2)-1) = imagesReal{j};
                        imagesOnCanvas = [imagesOnCanvas; [j, newY, newX]];
                    end
                    if iIndex == 0 && jIndex ~= 0
                        imY = imagesOnCanvas(jIndex, 2);
                        imX = imagesOnCanvas(jIndex, 3);
                        newY = imY - (imageSize(1)-yI);
                        newX = imX + xI;
                        %Add averaging later on. Maybe???
                        canvas(newY:newY + imageSizesize(imagesReal{i}, 1)-1, newX:newX + size(imagesReal{i},2)-1) = imagesReal{i};
                        imagesOnCanvas = [imagesOnCanvas; [i, newY, newX]];
                    end                    
                end
                if mIndex == 3
                    if iIndex ~= 0 && jIndex == 0
                        imY = imagesOnCanvas(iIndex, 2);
                        imX = imagesOnCanvas(iIndex, 3);
                        newY = imY + yI;
                        newX = imX - (imageSize(2)-xI);
                        %Add averaging later on. Maybe???
                        canvas(newY:newY + size(imagesReal{j}, 1)-1, newX:newX + size(imagesReal{j}, 2)-1) = imagesReal{j};
                        imagesOnCanvas = [imagesOnCanvas; [j, newY, newX]];
                    end
                    if iIndex == 0 && jIndex ~= 0
                        imY = imagesOnCanvas(jIndex, 2);
                        imX = imagesOnCanvas(jIndex, 3);
                        newY = imY + yI;
                        newX = imX - (imageSize(2)-xI);
                        %Add averaging later on. Maybe???
                        canvas(newY:newY + size(imagesReal{i}, 1)-1, newX:newX + size(imagesReal{i}, 2)-1) = imagesReal{i};
                        imagesOnCanvas = [imagesOnCanvas; [i, newY, newX]];
                    end                    
                end
                if mIndex == 4                   
                    if iIndex ~= 0 && jIndex == 0
                        imY = imagesOnCanvas(iIndex, 2);
                        imX = imagesOnCanvas(iIndex, 3);
                        newY = imY + yI - 1;
                        newX = imX + xI - 1;
                        %Add averaging later on. Maybe???
                        canvas(newY:newY + size(imagesReal{j}, 1)-1, newX:newX + size(imagesReal{j}, 2)-1, :) = imagesReal{j};
                        imagesOnCanvas = [imagesOnCanvas; [j, newY, newX]];
                    end
                    if iIndex == 0 && jIndex ~= 0
                        imY = imagesOnCanvas(jIndex, 2);
                        imX = imagesOnCanvas(jIndex, 3);
                        newY = imY + yI - 1;
                        newX = imX + xI - 1;
                        %Add averaging later on. Maybe???
                        canvas(newY:newY + size(imagesReal{i}, 1)-1, newX:newX + size(imagesReal{i}, 2)-1) = imagesReal{i};
                        imagesOnCanvas = [imagesOnCanvas; [i, newY, newX]];
                    end
                end
            end            
        end
    end
    %imshow(canvas);
    canvasResized = reduceCanvasSize(canvas, imagesOnCanvas, imagesReal);
    imwrite(canvasResized, 'Outputs/mosaic.png');
end

%Another way to implement DIC.
function [A,B] = calculateCorrelation(Im1, Im2, lowPassFilter)
    %Expecting the images in double and histogram equalized.
    fIm1 = fft2(Im1);
    Ga = fftshift(fIm1);    
    Ga = conj(Ga);
    
    fIm2 = fft2(Im2);
    Gb = fftshift(fIm2);
    
    G = Ga .* Gb;
    G = G ./ abs(G);
    G = G .* lowPassFilter;
    G = ifftshift(G);
    G = ifft2(G);
    realG = real(G);
    [maxIndex, maxVal] = findPeak(realG);
    
    A = maxVal;
    B = maxIndex;
end

function [A] = calculateDIC(Im1, Im2)
    mean1 = mean(mean(Im1));
    mean2 = mean(mean(Im2));
    Im1 = Im1 - mean1;
    Im2 = Im2 - mean2;
    Im1S = Im1 .^ 2;
    Im2S = Im2 .^ 2;
    A = sum(sum(Im1 .* Im2))/sqrt(sum(sum(Im1S)) * sum(sum(Im2S)));    
end

function [A] = reduceCanvasSize(canvas, imagesOnCanvas, imagesReal)    
    xMin = min(imagesOnCanvas(:,3));
    [xMax, xInd] = max(imagesOnCanvas(:,3));
    yMin = min(imagesOnCanvas(:,2));
    [yMax, yInd] = max(imagesOnCanvas(:,2));
    xMax = xMax + size(imagesReal{xInd}, 2);
    yMax = yMax + size(imagesReal{yInd}, 1);
    A = canvas(yMin:yMax, xMin:xMax, :);
end

function [A] = createLowPassFilterG(imSize, kernelNo, sigma, power)
    lowFilter = zeros(imSize);
    centreX = imSize(1,2)/2;
    centreY = imSize(1,1)/2;
    for i = 1:imSize(1,2)
        for j = 1:imSize(1,1)
            distance = (i-centreX)^2 + (j-centreY)^2;
            if kernelNo == 1
                lowFilter(j,i) = exp(-(distance / (2 * sigma^2)));
            elseif kernelNo == 2
                distance = sqrt(distance);
                lowFilter(j,i) = 1/(1+(distance/sigma)^(2*power));
            elseif kernelNo == 0
                distance = sqrt(distance);
                if distance <= sigma
                    lowFilter(j,i) = 1;
                else
                    lowFilter(j,i) = 0;
                end
            end
        end
    end
    A = lowFilter;
end

function [A] = isPresentOnCanvas(list, imageNo)    
    [val, loc] = ismember(imageNo, list(1:end,1));
    A = loc;
end

function [A, B] = findPeak(im)
    maxVal = im(1,1);
    maxIndex = [1,1];
    for i = 1:size(im,1)
        for j = 1:size(im,2)
            if maxVal < im(i,j)
                maxVal = im(i,j);
                maxIndex = [i, j];
            end
        end
    end
    
    %Finding a connected region around the peak.    
    maxVal1 = maxVal - (20*maxVal/100);
    binIm = im > maxVal1;    
    [labelout, count, neighbours, avgV] = floodFill(binIm, maxIndex, 100, zeros(size(binIm)));
    avgV = avgV / count;
    avgV = round(avgV);
    
    A = avgV;
    B = maxVal;
end

function [A, B] = padImages(inputFolder, fileNames)
    images = [];
    imagesReal = [];
    maxWidth = 0;
    maxHeight = 0;
    fileNames = {fileNames.name};
    for k = 1:numel(fileNames)
        file = fileNames{k};        
        Im = imread(strcat(inputFolder, file));        
        images{k} = Im;
        imagesReal{k} = Im;
        if maxHeight < size(Im, 1)
            maxHeight = size(Im, 1);
        end
        if maxWidth < size(Im, 2)
            maxWidth = size(Im, 2);
        end
    end
    
    %Padding the images.
    for k = 1:numel(images)
        Im = images{k};
        if size(Im, 1) < maxHeight
            Im = padarray(Im, [(maxHeight - size(Im, 1)), 0], 'replicate', 'post');
            %Im = padarray(Im, [(maxHeight - size(Im, 1)), 0], 0, 'post');
        end
        if size(Im, 2) < maxWidth
            Im = padarray(Im, [0, (maxWidth - size(Im, 2))], 'replicate', 'post');
            %Im = padarray(Im, [0, (maxWidth - size(Im, 2))], 0, 'post');
        end
        images{k} = Im;
    end
    A = images;
    B = imagesReal;
end