function morph(inParamsFileName)    
    [aMatrix1, aMatrix2, controlPoints1, controlPoints2, step, stdParam, kernelNo, inputFiles, outputFileName] = readParams(inParamsFileName);
    
    inputFile1 = imread(strcat('Inputs/', inputFiles{1}));
    inputFile2 = imread(strcat('Inputs/', inputFiles{2}));       
    
    controlPointsB = zeros(size(controlPoints1));    
    for t = 0.0:step:1.0
        %t = 0.5;
        controlPointsB = (1-t)*controlPoints1 + t*controlPoints2;
        bVector = [0;0;0;controlPointsB(:,1);0;0;0;controlPointsB(:,2)];
        %Calculating x Vector weights.
        xVector1 = pinv(aMatrix1) * bVector;
        xVector2 = pinv(aMatrix2) * bVector;
        
        %Code for testing the warp. Ignore this.
        %test(controlPoints1, inputFile1, xVector1, stdParam);
        
        %Transforming the image. Calculating x and y min/max.
        [x1Range, y1Range] = doForwardWarping(controlPoints1, inputFile1, xVector1, stdParam, kernelNo);
        [x2Range, y2Range] = doForwardWarping(controlPoints2, inputFile2, xVector2, stdParam, kernelNo);
    
        %Calculating new min max.
        xRange = [0,0];
        if x1Range(1,1) < x2Range(1,1)
            xRange(1,1) = x1Range(1,1);
        else
            xRange(1,1) = x2Range(1,1);
        end
        if x1Range(1,2) > x2Range(1,2)
            xRange(1,2) = x1Range(1,2);
        else
            xRange(1,2) = x2Range(1,2);
        end
        
        yRange = [0,0];
        if y1Range(1,1) < y2Range(1,1)
            yRange(1,1) = y1Range(1,1);
        else
            yRange(1,1) = y2Range(1,1);
        end
        if y1Range(1,2) > y2Range(1,2)
            yRange(1,2) = y1Range(1,2);
        else
            yRange(1,2) = y2Range(1,2);
        end
        
        %Creating canvas for inverse mapping now.        
        xRange = ceilFix(xRange);
        yRange = ceilFix(yRange);
        imageSize = [yRange(1,2) - yRange(1,1) + 1, xRange(1,2) - xRange(1,1) + 1];
        if numel(size(inputFile1)) == 3
            imageSize = [imageSize, 3];
        end
        imageOffset = [1-xRange(1,1), 1-yRange(1,1)]; %Adding 1 since index in matlab starts with 1.
        intermediateImage = zeros(imageSize); %This is the intermediate image.
        
        %Offseting the control point.
        controlPointsB = bsxfun(@plus, controlPointsB, imageOffset);
        
        %Creating the a inverse matrix 1.
        aInvSubMatrix = createSubMatrix(controlPointsB, stdParam, kernelNo);
        aInvMatrix = zeros(size(aInvSubMatrix)*2);
        aInvSubMatrixSize = size(aInvSubMatrix);
        aInvMatrix(1:aInvSubMatrixSize(1), 1:aInvSubMatrixSize(2)) = aInvSubMatrix;
        aInvMatrix(aInvSubMatrixSize(1)+1:end, aInvSubMatrixSize(2)+1:end) = aInvSubMatrix;
        
        %Inverse weights.
        xInvVector1 = pinv(aInvMatrix) * [0;0;0;controlPoints1(:,1);0;0;0;controlPoints1(:,2)];
        xInvVector2 = pinv(aInvMatrix) * [0;0;0;controlPoints2(:,1);0;0;0;controlPoints2(:,2)];
                
        for i=1:size(intermediateImage, 1) % i is y.
            for j=1:size(intermediateImage, 2) % j is x.                
                pixelCoords = [i,j];
                phiMatrixP = calculateBasis(controlPointsB, pixelCoords, stdParam, kernelNo);
                subMatrixP = [phiMatrixP, i, j, 1]; %Maybe this is the problem for guassian.
                aMatrixP = zeros(size(subMatrixP)*2);
                subMatrixSizeP = size(subMatrixP);
                aMatrixP(1:subMatrixSizeP(1), 1:subMatrixSizeP(2)) = subMatrixP;
                aMatrixP(subMatrixSizeP(1)+1:end, subMatrixSizeP(2)+1:end) = subMatrixP;
                
                finalPos1 = aMatrixP * xInvVector1;
                finalPos2 = aMatrixP * xInvVector2;                
            
                val1 = 0;
                val2 = 0;
                if numel(size(inputFile1)) == 3
                    val1 = [0,0,0];
                    val2 = [0,0,0];
                end
                
                %Checking if the transformed point is inside the range of
                %the image.
                if finalPos1(1,1) >= 1 && finalPos1(1,1) <= size(inputFile1,2)
                    if finalPos1(2,1) >= 1 && finalPos1(2,1) <= size(inputFile1,1)
                        if numel(size(inputFile1)) == 3
                            val1 = bilinearInterpolation3(inputFile1, finalPos1);
                        else
                            val1 = bilinearInterpolation(inputFile1, finalPos1);
                        end
                    end
                end
                if finalPos2(1,1) >= 1 && finalPos2(1,1) <= size(inputFile2,2)
                    if finalPos2(2,1) >= 1 && finalPos2(2,1) <= size(inputFile2,1)
                        if numel(size(inputFile1)) == 3
                            val2 = bilinearInterpolation3(inputFile2, finalPos2);
                        else
                            val2 = bilinearInterpolation(inputFile2, finalPos2);
                        end
                    end
                end
                
                if numel(size(inputFile1)) == 3
                    intermediateImage(i,j,1) = (1-t)*val1(1) + t*val2(1);
                    intermediateImage(i,j,2) = (1-t)*val1(2) + t*val2(2);
                    intermediateImage(i,j,3) = (1-t)*val1(3) + t*val2(3);
                else
                    intermediateImage(i,j) = (1-t)*val1 + t*val2;
                end
            end
        end
        %Uncomment these if you wish to see the intermediate images.
        %figure;
        %imshow(intermediateImage,[]);
        intermediateImage = uint8(intermediateImage);
        imwrite(intermediateImage,strcat('Outputs/', outputFileName,'_' , num2str(uint16(t * 1.0/step)), '.png'));
    end    
end

%Interpolation for gray image.
function A = bilinearInterpolation(inputImage, inputPoint)
    inX = inputPoint(1,1);
    inY = inputPoint(2,1);
    
    x1 = floor(inX); x2 = ceil(inX);
    y1 = floor(inY); y2 = ceil(inY);
    
    f11 = inputImage(y1, x1);
    f12 = inputImage(y2, x1);
    f21 = inputImage(y1, x2);
    f22 = inputImage(y2, x2);
    
    A = (f11 * (x2 - inX) * (y2 - inY)) + (f21 * (inX - x1) * (y2 - inY)) + (f12 * (x2 - inX) * (inY - y1)) + (f22 * (inX - x1) * (inY - y1));
    A = A/((x2-x1) * (y2-y1));
    
    if (x2-x1) == 0 || (y2-y1)==0
        A = inputImage(uint16(round(inY)),uint16(round(inX)));
    end
end

%Interpolation for color image.
function A = bilinearInterpolation3(inputImage, inputPoint)
    A = [0,0,0];
    for i = 1:3
        inX = inputPoint(1,1);
        inY = inputPoint(2,1);
    
        x1 = floor(inX); x2 = ceil(inX);
        y1 = floor(inY); y2 = ceil(inY);
    
        f11 = inputImage(y1, x1, i);
        f12 = inputImage(y2, x1, i);
        f21 = inputImage(y1, x2, i);
        f22 = inputImage(y2, x2, i);
    
        A(i) = (f11 * (x2 - inX) * (y2 - inY)) + (f21 * (inX - x1) * (y2 - inY)) + (f12 * (x2 - inX) * (inY - y1)) + (f22 * (inX - x1) * (inY - y1));
        A(i) = A(i)/((x2-x1) * (y2-y1));
        
        if (x2-x1) == 0 || (y2-y1)==0
            A(i) = inputImage(uint16(round(inY)),uint16(round(inX)));
        end
    end
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

function [A1, A2, C1, C2, step, stdParam, kernelNo, inputFiles, outputFileName] = readParams(inParamsFileName)
    fileId = fopen(inParamsFileName);
    if fileId == -1
        disp('cant open the params file.')
        return
    end
    
    stdParam = 0;
    counter = 0;
    noOfPoints = 0;    
    controlPoints1 = 0;
    controlPoints2 = 0;
    noOfImages = 0;
    inputFiles = [];
    outputFileName = '';
    tValue = 0;
    
    line = fgetl(fileId);
    while ischar(line)
        if numel(line) ~= 0 && ~strcmpi(line(1), '/') && ~strcmpi(line(2), '/')
            counter = counter + 1;
            if counter == 1
                strList = strsplit(line);                
                noOfPoints = str2double(strList(1))/2;
                noOfImages = str2double(strList(2));                
                controlPoints1 = zeros(noOfPoints, 2);
                controlPoints2 = zeros(noOfPoints, 2);
                line = fgetl(fileId);
                continue;
            end
            if (counter-1) <= (noOfPoints*2)
                strList = strsplit(line);                
                if(mod(counter, 2) == 0) %Processing x coordinate
                    index = floor(counter/2);                                    
                    controlPoints1(index,1) = str2double(strList(1));
                    controlPoints2(index,1) = str2double(strList(2));
                else %Processing y coordinate
                    index = floor(counter/2);                                  
                    controlPoints1(index,2) = str2double(strList(1));
                    controlPoints2(index,2) = str2double(strList(2));
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
                tValue = str2double(line);
                line = fgetl(fileId);
                continue;
            end
            if counter <= 1 + (noOfPoints*2) + noOfImages + 1 + 1 + 1
                stdParam = str2double(line);
                line = fgetl(fileId);
                continue;
            end
            if counter <= 1 + (noOfPoints*2) + noOfImages + 1 + 1 + 1 + 1
                kernelNo = str2double(line);
                line=fgetl(fileId);
                continue;
            end
            %disp(line)
        end
        line = fgetl(fileId);
    end
    
    %Creating transform matrix A from c1.
    aSubMatrix1 = createSubMatrix(controlPoints1, stdParam, kernelNo);
    aMatrix1 = zeros(size(aSubMatrix1)*2);
    aSubMatrixSize1 = size(aSubMatrix1);
    aMatrix1(1:aSubMatrixSize1(1), 1:aSubMatrixSize1(2)) = aSubMatrix1;
    aMatrix1(aSubMatrixSize1(1)+1:end, aSubMatrixSize1(2)+1:end) = aSubMatrix1;
    
    %Creating transform matrix A from c2.
    aSubMatrix2 = createSubMatrix(controlPoints2, stdParam, kernelNo);
    aMatrix2 = zeros(size(aSubMatrix2)*2);
    aSubMatrixSize2 = size(aSubMatrix2);
    aMatrix2(1:aSubMatrixSize2(1), 1:aSubMatrixSize2(2)) = aSubMatrix2;
    aMatrix2(aSubMatrixSize2(1)+1:end, aSubMatrixSize2(2)+1:end) = aSubMatrix2;
    
    A1 = aMatrix1;
    A2 = aMatrix2;
    C1 = controlPoints1;
    C2 = controlPoints2;
    step = 1.0/tValue;    
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