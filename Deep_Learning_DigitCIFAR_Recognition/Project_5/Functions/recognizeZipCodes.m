function recognizeZipCodes(net)
    for i=1:2
        fname = strcat('Inputs/ZipCodes-Test/test', num2str(i), '.jpg');
        im_test = imread(fname);
        sub_images = preprocessImage(im_test, i);
        
        recognizedNumber = '';
        for j=1:numel(sub_images)
            sub_image = sub_images{j};
            im_test = im2single(sub_image);
            %Why is this code working without the mean subtraction.!!!!
            %data = bsxfun(@minus, im_test, dataMean);
    
            res = vl_simplenn(net, im_test) ;
            %     for i=1:size(res(end).x,2)
            %         [score(i),pred(i)] = max(squeeze(res(end).x(1,i,:))) ;
            %     end
            [score, pred] = max(squeeze(res(end).x(1,1,:)));
            %score
            %pred = pred-1
            pred = pred-1;
            recognizedNumber = strcat(recognizedNumber, int2str(pred));
        end        
        fprintf('Recognized nnumber: %s\n', recognizedNumber)
    end
end

function [A] = preprocessImage(im, im_number)
    im = rgb2gray(im);
    imThresh = graythresh(im);
    im = im < imThresh*256;
    imLabeled = bwlabel(im, 4);
    siz=size(im);
    stats = regionprops(im, 'BoundingBox');
    ObjCell={};
    for i = 1:numel(stats)
        bbox = stats(i).BoundingBox;
        if bbox(3)*bbox(4) < 50
            continue
        end
        idx_x=[bbox(1)-2 bbox(1)+bbox(3)+2];
        idx_y=[bbox(2)-2 bbox(2)+bbox(4)+2];
        if idx_x(1)<1, idx_x(1)=1; end
        if idx_y(1)<1, idx_y(1)=1; end
        if idx_x(2)>siz(2), idx_x(2)=siz(2); end
        if idx_y(2)>siz(1), idx_y(2)=siz(1); end
        % Crop the object and write to ObjCell
        %im=L==i;
        subImage = uint8(im(idx_y(1):idx_y(2),idx_x(1):idx_x(2)) * 255);
        subImage = imresize(subImage, [28 28]);
        ObjCell{end+1} = subImage;
    end
    
    f1 = figure;
    for i=1:numel(ObjCell)        
        subplot(1,numel(ObjCell),i);
        imshow(ObjCell{i});     
    end
    saveas(f1, strcat('Outputs/Digits', int2str(im_number), '.png'));

    figure
    imshow(im);
    im = im *255;
    A = ObjCell;
end