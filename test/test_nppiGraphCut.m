function test_nppiGraphCut()
    addpath('../external_code');
    addpath('../external_code/npp');
   
    addpath('../external_code/vlfeats/toolbox/kmeans/');
    addpath('../external_code/vlfeats/toolbox/mex/mexa64/');
    addpath('../external_code/vlfeats/toolbox/mex/mexglx/');
    
    
	%img_name = '2010_000238'; % airplane and people   
    img_name = '2007_009084'; % dogs, motorbike, chairs, people    
    %img_name = '2010_002868'; % buses   
    %img_name = '2010_003781'; % cat, bottle, potted plants

    I = imread(['../data/JPEGImages/' img_name '.jpg']);
    I = single(I)/255;

    rows = size(I, 1);
    cols = size(I, 2);

    width = 15;
    height = 15;
    
    prow = 170;
    pcol = 150;

    foregroundSeeds = extractSeeds(I, prow, pcol, width, height);
    
    [ix iy] = meshgrid(1:height, 1:width);
    seedids = [ix(:)+prow-1 iy(:)+pcol-1];
    seedids = sub2ind([rows cols], seedids(:, 1), seedids(:, 2));

    backgroundSeeds = [];
%     backgroundSeeds = [backgroundSeeds; extractSeeds(I, 1, 1, cols, 1)];
    backgroundSeeds = [backgroundSeeds; extractSeeds(I, 1, 1, 1, rows)];
    backgroundSeeds = [backgroundSeeds; extractSeeds(I, 1, cols, 1, rows)];
%     backgroundSeeds = [backgroundSeeds; extractSeeds(I, rows, 1, cols, 1)];

    backgroundids = sub2ind([rows cols], 1:rows, ones(1, rows));
    backgroundids = [backgroundids sub2ind([rows cols], 1:rows, ones(1, rows)*cols)];

%     backgroundSeeds = extractSeeds(I, 204, 420, width, height);
%     backgroundSeeds = [backgroundSeeds; extractSeeds(I, 335, 186, width, height)];
%     backgroundSeeds = [backgroundSeeds; extractSeeds(I, 300, 20, width, height)];
%     
    
%     backgroundSeeds = backgroundSeeds(randsample(size(backgroundSeeds, 1), 25), :);

    Cs = computeCapacity(I, foregroundSeeds);
    Ct = computeCapacity(I, backgroundSeeds);

%       rows = 10;
%       cols = 15;
% 
%       Cs = zeros(rows, cols);
%       Cs(3:7, 3:7) = 1;
%       
%       Ct = zeros(rows, cols);
%       Ct(:, 1) = 1;
%       Ct(:, cols) = 1;
      

    CG = colgrad(I);  

    CONTRAST_SENSITIVE_WEIGHT = 1;
    POTTS_WEIGHT = 1/1000;
    SIGMA = 0.001;

    leftTranspose = CG(1: rows, 2:cols) - CG(1:rows, 1:cols-1);
    leftTranspose = abs(leftTranspose);
    leftTranspose = normVal(leftTranspose, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
    leftTranspose = [zeros(rows, 1) leftTranspose];
    leftTranspose = leftTranspose';

    rightTranspose = CG(1:rows, 1:cols-1) - CG(1:rows, 2:cols);
    rightTranspose = abs(rightTranspose);
    rightTranspose = normVal(rightTranspose, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
    rightTranspose = [rightTranspose zeros(rows, 1)];
    rightTranspose = rightTranspose';

    top = CG(2:rows, :) - CG(1:rows-1, :);
    top = abs(top);
    top = normVal(top, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
    top = [zeros(1, cols); top];

    bottom = CG(1:rows-1, :) - CG(2:rows, :);
    bottom = abs(bottom);
    bottom = normVal(bottom, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
    bottom = [bottom; zeros(1, cols)];
    
    CTerminals = single(Cs-Ct);
    

    for i = -2:2
        ct = CTerminals;
        lambda = i/10;
        ct(seedids) = ct(seedids) + lambda;
        ct(backgroundids) = ct(backgroundids) - lambda;
        labels = nppiGraphcut_32f8u_mex(cols, rows, ct, leftTranspose, rightTranspose, top, bottom); 
        figure, imagesc(labels);
    end

    unique(labels)

%     uu = im2bw(uu, level);
    
    

end

