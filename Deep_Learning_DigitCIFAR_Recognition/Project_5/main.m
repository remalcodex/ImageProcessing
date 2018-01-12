clc;
close all;
clear;
imtool close all;
addpath(genpath('Functions/'));
addpath(genpath('Functions/mnistHelper'));

%run the vl setup path.
run('matconvnet-1.0-beta25/matlab/vl_setupnn.m')

%Step 1.1
[net, info, dataMean] = cnn_mnist();
save('Outputs/model.mat', '-struct', 'net') ;

%Step 1.2
[net, info, dataMean] = cnn_mnist_odd_even();

%Step 2.1
net = load('Outputs/model.mat') ;
net.layers{end}.type = 'softmax';
recognizeZipCodes(net)

%Step 3.1
cnn_cifar();

close all;
