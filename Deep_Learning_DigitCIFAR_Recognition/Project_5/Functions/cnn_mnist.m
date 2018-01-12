function [net, info, dataMean] = cnn_mnist()
    opts.networkType = 'simplenn' ;
    opts.expDir = fullfile(pwd, 'Outputs', 'firstModel');

    opts.train = struct();
    opts.train.gpus = [];

    % Preparing model and data.
    net = cnn_mnist_init() ;
    imdb = getMNISTData(0);
    net.meta.classes.name = arrayfun(@(x)sprintf('%d',x),1:10,'UniformOutput',false) ;

    % Training
    trainfn = @cnn_train;
    [net, info] = trainfn(net, imdb, getBatch(), 'expDir', opts.expDir, ...
      net.meta.trainOpts, opts.train, 'val', find(imdb.images.set == 3)) ;
    dataMean = imdb.images.data_mean;


function fn = getBatch()
    fn = @(x,y) getSimpleNNBatch(x,y) ;

function [images, labels] = getSimpleNNBatch(imdb, batch)
    images = imdb.images.data(:,:,:,batch) ;
    labels = imdb.images.labels(1,batch) ;
