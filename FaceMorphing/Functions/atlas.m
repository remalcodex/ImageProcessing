function morph(inParamsFileName)    
    [controlPoints, controlPointsMean, stdParam, kernelNo, inputFiles, outputFileName] = readParams(inParamsFileName);
    
    %Calculating the canvas size and creating the canvas image.
    xRange = [[0,0];[0,0];[0,0];[0,0]];
    yRange = [[0,0];[0,0];[0,0];[0,0]];
    for t = 1:size(inputFiles,2)
        inputFile = imread(strcat('Inputs/', inputFiles{t}));        
        bVector = [0;0;0;controlPointsMean(:,1);0;0;0;controlPointsMean(:,2)];
        
        %Calculating A matrix.
        cp = controlPoints(:,(t*2-1):(t*2));
        aSubMatrix = createSubMatrix(cp, stdParam, kernelNo);
        aMatrix = zeros(size(aSubMatrix)*2);
        aSubMatrixSize = size(aSubMatrix);
        aMatrix(1:aSubMatrixSize(1), 1:aSubMatrixSize(2)) = aSubMatrix;
        aMatrix(aSubMatrixSize(1)+1:end, aSubMatrixSize(2)+1:end) = aSubMatrix;
        
        xVector = pinv(aMatrix) * bVector;
        [xRangeT, yRangeT] = doForwardWarping(cp, inputFile, xVector, stdParam, kernelNo);
        xRange(t, :) = xRangeT;
        yRange(t, :) = yRangeT;
    end
    xMin = min(xRange(:,1));
    xMax = max(xRange(:,2));
    yMin = min(yRange(:,1));
    yMax = max(yRange(:,2));
    xRange = [xMin, xMax];
    yRange = [yMin, yMax];
    xRange = ceilFix(xRange);
    yRange = ceilFix(yRange);
    imageSize = [yRange(1,2) - yRange(1,1) + 1, xRange(1,2) - xRange(1,1) + 1];    
    imageOffset = [1-xRange(1,1), 1-yRange(1,1)]; %Adding 1 since index in matlab starts with 1.
    canvasImage = zeros(imageSize); %This is the intermediate image.
    
    %Offseting the control point.
    controlPointsMean = bsxfun(@plus, controlPointsMean, imageOffset);
        
    %Creating the a inverse matrix.
    aInvSubMatrix = createSubMatrix(controlPointsMean, stdParam, kernelNo);
    aInvMatrix = zeros(size(aInvSubMatrix)*2);
    aInvSubMatrixSize = size(aInvSubMatrix);
    aInvMatrix(1:aInvSubMatrixSize(1), 1:aInvSubMatrixSize(2)) = aInvSubMatrix;
    aInvMatrix(aInvSubMatrixSize(1)+1:end, aInvSubMatrixSize(2)+1:end) = aInvSubMatrix;
    
    for t = 1:size(inputFiles,2)
        %Debugging purpose only.
        %canvasImage = zeros(imageSize);
        
        inputFile = imread(strcat('Inputs/', inputFiles{t}));
        %Doing the contrast adjustment with adaptive histogram
        %equalization.
        inputFile = adapthisteq(inputFile);
        %Inverse weight.
        cp = controlPoints(:,(t*2-1):(t*2));
        xInvVector = pinv(aInvMatrix) * [0;0;0;cp(:,1);0;0;0;cp(:,2)];        
        
        for i=1:size(canvasImage, 1) % i is y.
            for j=1:size(canvasImage, 2) % j is x.                
                pixelCoords = [i,j];
                
                %Calculating A matrix.
                phiMatrixP = calculateBasis(controlPointsMean, pixelCoords, stdParam, kernelNo);
                subMatrixP = [phiMatrixP, i, j, 1];
                aMatrixP = zeros(size(subMatrixP)*2);
                subMatrixSizeP = size(subMatrixP);
                aMatrixP(1:subMatrixSizeP(1), 1:subMatrixSizeP(2)) = subMatrixP;
                aMatrixP(subMatrixSizeP(1)+1:end, subMatrixSizeP(2)+1:end) = subMatrixP;
                
                finalPos = aMatrixP * xInvVector;
                val1 = 0;
                if finalPos(1,1) >= 1 && finalPos(1,1) <= size(inputFile,2)
                    if finalPos(2,1) >= 1 && finalPos(2,1) <= size(inputFile,1)                        
                        val1 = bilinearInterpolation(inputFile, finalPos);                      
                    end
                end
                canvasImage(i,j) = canvasImage(i,j) + val1;
                
                %Averaging in the end.
                if t == size(inputFiles,2)
                    canvasImage(i,j) = canvasImage(i,j) / t;
                end
            end
        end
    end
    canvasImage = uint8(canvasImage);
    imwrite(canvasImage, strcat('Outputs/', outputFileName));    
end

%Interpolation for gray image.
function A = bilinearInterpolation(inputImage, inputPoint)
    inX = inputPoint(1,1);
    inY = inputPoint(2,1);
    
    x1 = floor(inX); x2 = ceil(inX);
    y1 = floor(inY); y2 = ceil(inY);
    
    f11 = double(inputImage(y1, x1)); %TODO: Fix this in morphing as well.
    f12 = double(inputImage(y2, x1));
    f21 = double(inputImage(y1, x2));
    f22 = double(inputImage(y2, x2));
    
    A = (f11 * (x2 - inX) * (y2 - inY)) + (f21 * (inX - x1) * (y2 - inY)) + (f12 * (x2 - inX) * (inY - y1)) + (f22 * (inX - x1) * (inY - y1));
    A = A/((x2-x1) * (y2-y1));
end

function A = ceilFix(x)
    A = ceil(abs(x)).*sign(x);
end

%This function returns the range of x and y in the forward transformation.
function [xRange, yRange] = doForwardWarping(controlPoints, inputImage, xVector, stdParam, kernelNo)
    xArray = [0,0];
    yArray = [0,0];
    firstPass = 0;
    for i=1:size(inputImage,1) % i is y.
        for j=1:size(inputImage,2) % j is x.
            pixelCoords = [i,j];
            phiMatrix1 = calculateBasis(controlPoints, pixelCoords, stdParam, kernelNo);
            subMatrix = [phiMatrix1, i, j, 1];
            aMatrix = zeros(size(subMatrix)*2);
            subMatrixSize = size(subMatrix);
            aMatrix(1:subMatrixSize(1), 1:subMatrixSize(2)) = subMatrix;
            aMatrix(subMatrixSize(1)+1:end, subMatrixSize(2)+1:end) = subMatrix;
            finalPos = aMatrix * xVector;
            if firstPass == 0
                firstPass = 1;
                xArray(1,1) = finalPos(1,1);
                xArray(1,2) = finalPos(1,1);
                yArray(1,1) = finalPos(2,1);
                yArray(1,2) = finalPos(2,1);
            end
        
            %Calculating new dimension of the array
            if finalPos(1,1) < xArray(1,1)
                xArray(1,1) = finalPos(1,1);
            end
            if finalPos(1,1) > xArray(1,2)
                xArray(1,2) = finalPos(1,1);
            end
            if finalPos(2,1) < yArray(1,1)
                yArray(1,1) = finalPos(2,1);
            end
            if finalPos(2,1) > yArray(1,2)
                yArray(1,2) = finalPos(2,1);
            end
        end
    end
    
    xRange = xArray;
    yRange = yArray;
    return;
end

%Calculates the phi matrix.
function phiMatrix = calculateBasis(controlPoints, newPoint, stdParam, kernelNo)
    phiMatrix = zeros(1, size(controlPoints,1));    
    tempMatrix = bsxfun(@minus, controlPoints, newPoint); %Has both x and y points in each row.
    
    % Gaussian Method.
    if (kernelNo == 1)
        tempMatrix = tempMatrix .^ 2;
        tempMatrix = tempMatrix(:,1) + tempMatrix(:,2);
        varianceT = stdParam^2;
        varianceT = (-2) * varianceT;
        tempMatrix = tempMatrix / varianceT;
        tempMatrix = exp(tempMatrix);
        tempMatrix = sqrt(tempMatrix);
    end
    
    %Spline Method.
    if (kernelNo == 2)
        tempMatrix1 = tempMatrix .^ 2;
        tempMatrix1 = tempMatrix1(:,1) + tempMatrix1(:,2);
        tempMatrix2 = sqrt(tempMatrix1);        
        for i=1:size(tempMatrix2,1)
            if tempMatrix2(i,1) ~= 0
                tempMatrix2(i,1) = log(tempMatrix2(i,1));
            end
        end
        tempMatrix = tempMatrix1 .* tempMatrix2;
    end
    
    %Inverse Quadratic Method.
    if (kernelNo == 3)
        tempMatrix = tempMatrix .^ 2;
        tempMatrix = tempMatrix(:,1) + tempMatrix(:,2); %This is already squared.
        tempMatrix = (stdParam.^2) * tempMatrix;        
        tempMatrix = tempMatrix + 1;
        tempMatrix = 1./tempMatrix;
    end
    
    phiMatrix(1, :) = tempMatrix;
end

function A = createSubMatrix(controlPoints, stdParam, kernelNo)
    noOfPoints = size(controlPoints, 1);
    phiMatrix = zeros(noOfPoints, noOfPoints);
    aSubMatrix = zeros(noOfPoints + 3, noOfPoints + 3);
    for i = 1:noOfPoints     
        aSubMatrix(1,i) = controlPoints(i,1);
        aSubMatrix(2,i) = controlPoints(i,2);
        aSubMatrix(3+i,end-1) = controlPoints(i,1);
        aSubMatrix(3+i,end-2) = controlPoints(i,2);
        
        %Calculating phi Matrix. Replace this with calclate basis method.        
        tempMatrix = calculateBasis(controlPoints, controlPoints(i,:), stdParam, kernelNo);
        phiMatrix(:, i) = tempMatrix;
    end
    
    aSubMatrix(3,1:noOfPoints) = 1;
    aSubMatrix(4:end,end) = 1;
    aSubMatrix(4:end,1:size(phiMatrix,2)) = phiMatrix;
    A = aSubMatrix;    
end

function [C, CM, stdParam, kernelNo, inputFiles, outputFileName] = readParams(inParamsFileName)
    fileId = fopen(inParamsFileName);
    if fileId == -1
        disp('cant open the params file.')
        return
    end
    
    stdParam = 0;
    counter = 0;
    noOfPoints = 0;
    controlPoints = 0; %This will be the full array of all control points.  
    noOfImages = 0;
    inputFiles = [];
    outputFileName = '';    
    
    line = fgetl(fileId);
    while ischar(line)
        if numel(line) ~= 0 && ~strcmpi(line(1), '/') && ~strcmpi(line(2), '/')
            counter = counter + 1;
            if counter == 1
                strList = strsplit(line);                
                noOfPoints = str2double(strList(1))/2;
                noOfImages = str2double(strList(2));
                controlPoints = zeros(noOfPoints, noOfImages*2);
                line = fgetl(fileId);
                continue;
            end
            if (counter-1) <= (noOfPoints*2)
                strList = strsplit(line);                
                if(mod(counter, 2) == 0) %Processing x coordinate
                    index = floor(counter/2);
                    for i = 1:noOfImages
                       controlPoints(index,(i*2-1)) = str2double(strList(i));
                    end                    
                else %Processing y coordinate
                    index = floor(counter/2);
                    for i = 1:noOfImages
                        controlPoints(index,(i*2)) = str2double(strList(i));
                    end
                end
                line = fgetl(fileId);
                continue;
            end
            if counter <= 1 + (noOfPoints*2) + noOfImages
                inputFiles{end+1} = line;
                line = fgetl(fileId);
                continue;
            end
            if counter <= 1 + (noOfPoints*2) + noOfImages + 1
                outputFileName = line;
                line = fgetl(fileId);
                continue;
            end
            if counter <= 1 + (noOfPoints*2) + noOfImages + 1 + 1
                stdParam = str2double(line);
                line = fgetl(fileId);
                continue;
            end
            if counter <= 1 + (noOfPoints*2) + noOfImages + 1 + 1 + 1
                kernelNo = str2double(line);
                line=fgetl(fileId);
                continue;
            end
            %disp(line)
        end
        line = fgetl(fileId);
    end
    
    %Take mean of all controlPoints.
    controlPointsMean = zeros(noOfPoints,2);
    for i = 1:noOfImages
        controlPointsMean(:, 1) = controlPointsMean(:, 1) + controlPoints(:,(2*i-1));
        controlPointsMean(:,2) = controlPointsMean(:,2) + controlPoints(:,(2*i));        
    end
    controlPointsMean = controlPointsMean / noOfImages;
    
    C = controlPoints;
    CM = controlPointsMean;
    fclose(fileId);
end

%Code for testing.
function test(controlPoints, inputImage, xVector, stdParam, kernelNo)
    intermediateImage = zeros(size(inputImage));
    for i=1:size(inputImage,1) % i is y.
        for j=1:size(inputImage,2) % j is x.
            pixelCoords = [i,j];
            phiMatrix1 = calculateBasis(controlPoints, pixelCoords, stdParam, kernelNo);
            subMatrix = [phiMatrix1, i, j, 1];
            aMatrix = zeros(size(subMatrix)*2);
            subMatrixSize = size(subMatrix);
            aMatrix(1:subMatrixSize(1), 1:subMatrixSize(2)) = subMatrix;
            aMatrix(subMatrixSize(1)+1:end, subMatrixSize(2)+1:end) = subMatrix;
            finalPos = aMatrix * xVector;
            finalPos = round(finalPos);
            finalPos = uint16(finalPos);
            if(finalPos(1,1) ~= j || finalPos(2,1) ~= i)
             disp('damn!') 
            end
            
            intermediateImage(finalPos(2,1), finalPos(1,1)) = double(inputImage(i,j));
        end
    end
    intermediateImage = uint8(intermediateImage);
    imshow(intermediateImage);
end