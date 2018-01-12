function imdb = getMNISTData(doOddEven)
% Preparing the imdb structure, returns image data with mean image subtracted
f=fopen('Inputs/train-images.idx3-ubyte','r') ;
x1=fread(f,inf,'uint8');
fclose(f) ;
x1=permute(reshape(x1(17:end),28,28,60e3),[2 1 3]) ;

f=fopen('Inputs/t10k-images.idx3-ubyte','r') ;
x2=fread(f,inf,'uint8');
fclose(f) ;
x2=permute(reshape(x2(17:end),28,28,10e3),[2 1 3]) ;

f=fopen('Inputs/train-labels.idx1-ubyte','r') ;
y1=fread(f,inf,'uint8');
fclose(f) ;
y1=double(y1(9:end)')+1 ;

f=fopen('Inputs/t10k-labels.idx1-ubyte','r') ;
y2=fread(f,inf,'uint8');
fclose(f) ;
y2=double(y2(9:end)')+1 ;

if doOddEven == 1
    for y1_i = 1:numel(y1)
        if mod(y1(y1_i), 2) == 0
            y1(y1_i) = 2;
        else
            y1(y1_i) = 1;
        end
    end
    
    for y2_i = 1:numel(y2)
        if mod(y2(y2_i), 2) == 0
            y2(y2_i) = 2;
        else
            y2(y2_i) = 1;
        end
    end
end

set = [ones(1,numel(y1)) 3*ones(1,numel(y2))];
data = single(reshape(cat(3, x1, x2),28,28,1,[]));
dataMean = mean(data(:,:,:,set == 1), 4);
data = bsxfun(@minus, data, dataMean) ;

imdb.images.data = data ;
imdb.images.data_mean = dataMean;
imdb.images.labels = cat(2, y1, y2) ;
imdb.images.set = set ;
imdb.meta.sets = {'train', 'val', 'test'} ;
if doOddEven == 1
    imdb.meta.classes = arrayfun(@(x)sprintf('%d',x),0:1,'uniformoutput',false) ;
else
    imdb.meta.classes = arrayfun(@(x)sprintf('%d',x),0:9,'uniformoutput',false) ;
end