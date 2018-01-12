function identifyEvenOdd(net)
    %------------------ Odd-Even recognition.
    %Modifying the netowrk to add another FCN layer with one node for odd-even
    %numbers.
    %maybe add bias.
    %net.layers{end+1} = struct('type', 'fcn', ...
     %                          'weights', {{f*randn(10,1,1,1, 'single'),zeros(1,50,'single')}}, ...
      %                         'stride', 1, ...
       %                        'pad', 0) ;
    % layer = fullyConnectedLayer(1);
    % layer.Weights = [0,1,0,1,0,1,0,1,0,1]
    % net.layers{end+1} = layer;
    dump = randn(1,1,10,1, 'single');
    dump(:,:,1,:) = 0;
    dump(:,:,2,:) = 0;
    dump(:,:,3,:) = 1;
    dump(:,:,4,:) = 0;
    dump(:,:,5,:) = 1;
    dump(:,:,6,:) = 0;
    dump(:,:,7,:) = 1;
    dump(:,:,8,:) = 0;
    dump(:,:,9,:) = 1;
    dump(:,:,10,:) = 0;

    net.layers{end+1} = struct('type', 'conv', ...
                               'weights', {{dump, zeros(1,1,'single')}}, ...
                               'stride', 1, ...
                               'pad', 0) ;
    net = vl_simplenn_tidy(net);
    vl_simplenn_display(net) ;
    for i = 1:9
        fname = strcat('Inputs/Odd-Even-Test/im', num2str(i), '.png');
        im_test = imread(fname);    
        im_test = im2single(im_test);
        %Why is this code working without the mean subtraction.!!!!
        %data = bsxfun(@minus, im_test, dataMean);

         res = vl_simplenn(net, im_test) ;
    %     for i=1:size(res(end).x,2)
    %         [score(i),pred(i)] = max(squeeze(res(end).x(1,i,:))) ;
    %     end
        [score, pred] = max(squeeze(res(end).x(1,1,:)));
        score
        pred = pred-1

        % accumulate errors
    %     error = sum([error, [...
    %       sum(double(gather(res(end).x))) ;
    %       reshape(params.errorFunction(params, labels, res),[],1) ; ]],2) ;
    end
end