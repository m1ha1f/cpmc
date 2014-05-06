function testCMF()
    addpath('../external_code')
    addpath('../external_code/CMF v1.0');
    
	%img_name = '2010_000238'; % airplane and people   
    img_name = '2007_009084'; % dogs, motorbike, chairs, people    
    %img_name = '2010_002868'; % buses   
    %img_name = '2010_003781'; % cat, bottle, potted plants

    I = imread(['../data/JPEGImages/' img_name '.jpg']);
    I = single(I)/255;

    rows = size(I, 1);
    cols = size(I, 2);

    width = 20;
    height = 20;

    foregroundSeeds = extractSeeds(I, 150, 165, width, height);

    backgroundSeeds = extractSeeds(I, 1, 1, cols, 1);
    backgroundSeeds = [backgroundSeeds; extractSeeds(I, 1, 1, 1, rows)];
    backgroundSeeds = [backgroundSeeds; extractSeeds(I, 1, cols, 1, rows)];
    backgroundSeeds = [backgroundSeeds; extractSeeds(I, rows, 1, cols, 1)];

%     backgroundSeeds = extractSeeds(I, 204, 420, width, height);
%     backgroundSeeds = [backgroundSeeds; extractSeeds(I, 335, 186, width, height)];
%     backgroundSeeds = [backgroundSeeds; extractSeeds(I, 300, 20, width, height)];
%     
    
%     backgroundSeeds = backgroundSeeds(randsample(size(backgroundSeeds, 1), 25), :);

    Cs = computeCapacity(I, foregroundSeeds);
    Ct = computeCapacity(I, backgroundSeeds);
%     f = 0.12;
%     b = 0.75;
%     
%     Cs = abs(I - f);
%     Ct = abs(I - b);

%     penalty = 0.5*ones(rows,cols);

    CG = colgrad(I);
    penalty = 0.5*(CG + 0.5).^-1;



    varParas = [rows; cols; 300; 1e-4; 0.3; 0.16];
%     para 0,1 - rows, cols of the given image
%     para 2 - the maximum number of iterations
%     para 3 - the error bound for convergence
%     para 4 - cc for the step-size of augmented Lagrangian method
%     para 5 - the step-size for the graident-projection of p

    [uu, erriter,num,tt] = CMF_GPU(single(penalty), single(Cs), single(Ct), single(varParas));
    
%     uu = im2bw(uu, level);
    
    figure, imagesc(uu);

end

