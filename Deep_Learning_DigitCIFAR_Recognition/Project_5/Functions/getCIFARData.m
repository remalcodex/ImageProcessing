function imdb = getCIFARData()
files = {};
for i=1:5
    files{i} = strcat('Inputs/cifar-10-batches-mat/data_batch_', int2str(i), '.mat');
end
files{6} = strcat('Inputs/cifar-10-batches-mat/test_batch.mat');
file_set = uint8([ones(1, 5), 3]);

data = cell(1, numel(files));
labels = cell(1, numel(files));
sets = cell(1, numel(files));
for fi = 1:numel(files)
  fd = load(files{fi}) ;
  data{fi} = permute(reshape(fd.data',32,32,3,[]),[2 1 3 4]) ;
  labels{fi} = fd.labels' + 1; % Index from 1
  sets{fi} = repmat(file_set(fi), size(labels{fi}));
end

set = cat(2, sets{:});
data = single(cat(4, data{:}));

% Removing mean
dataMean = mean(data(:,:,:,set == 1), 4);
data = bsxfun(@minus, data, dataMean);

% Read it in Andrew Ng's  website. Try settting them false.
contrastNormalization = true;
if contrastNormalization == true
  z = reshape(data,[],60000) ;
  z = bsxfun(@minus, z, mean(z,1)) ;
  n = std(z,0,1) ;
  z = bsxfun(@times, z, mean(n) ./ max(n, 40)) ;
  data = reshape(z, 32, 32, 3, []) ;
end

% Try with setting them false.
whitenData = true;
if whitenData == true
  z = reshape(data,[],60000) ;
  W = z(:,set == 1)*z(:,set == 1)'/60000 ;
  [V,D] = eig(W) ;
  % the scale is selected to approximately preserve the norm of W
  d2 = diag(D) ;
  en = sqrt(mean(d2)) ;
  z = V*diag(en./max(sqrt(d2), 10))*V'*z ;
  data = reshape(z, 32, 32, 3, []) ;
end

clNames = load(fullfile('Inputs/cifar-10-batches-mat', 'batches.meta.mat'));

imdb.images.data = data ;
imdb.images.labels = single(cat(2, labels{:})) ;
imdb.images.set = set;
imdb.meta.sets = {'train', 'val', 'test'} ;
imdb.meta.classes = clNames.label_names;