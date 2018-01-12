function [net, info] = cnn_cifar()    
    opts.expDir = fullfile(pwd, 'Outputs', 'cifarModel');   
    opts.train = struct() ;
    opts.train.gpus = [];

    % -------- Preparing model and data.
    net = cnn_cifar_init();
    imdb = getCIFARData();
    net.meta.classes.name = imdb.meta.classes(:)' ;

    % ----------------------------- Training.
    [net, info] = cnn_train(net, imdb, getBatch(), 'expDir', opts.expDir, ...
      net.meta.trainOpts, opts.train, 'val', find(imdb.images.set == 3)) ;

function fn = getBatch()
    fn = @(x,y) getNNBatch(x,y) ;

function [images, labels] = getNNBatch(imdb, batch)
    images = imdb.images.data(:,:,:,batch) ;
    labels = imdb.images.labels(1,batch) ;
    if rand > 0.5, images=fliplr(images) ; end
