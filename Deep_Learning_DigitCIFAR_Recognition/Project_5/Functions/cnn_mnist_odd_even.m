function [net, info, dataMean] = cnn_mnist_odd_even()
opts.network = [] ;
opts.expDir = fullfile(pwd, 'Outputs', 'oddEvenModel');
opts.train = struct();
opts.train.gpus = [];

% ----Preparing the data.
net = cnn_mnist_init_odd_even() ;
imdb = getMNISTData(1);
net.meta.classes.name = arrayfun(@(x)sprintf('%d',x),1:2,'UniformOutput',false) ;

% ------Training the network
 trainfn = @cnn_train ;
[net, info] = trainfn(net, imdb, getBatch(), 'expDir', opts.expDir, ...
  net.meta.trainOpts, opts.train, 'val', find(imdb.images.set == 3)) ;
dataMean = imdb.images.data_mean;

function fn = getBatch()
fn = @(x,y) getSimpleNNBatch(x,y) ;

function [images, labels] = getSimpleNNBatch(imdb, batch)
images = imdb.images.data(:,:,:,batch) ;
labels = imdb.images.labels(1,batch) ;
