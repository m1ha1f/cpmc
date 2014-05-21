function test_nppiGraphCut()
    addpath('../external_code');
    addpath('../external_code/npp');
    addpath('../external_code/CMF v1.0');
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

    width = 5;
    height = 5;

    foregroundSeeds = extractSeeds(I, 170, 150, width, height);

    backgroundSeeds = [];
%     backgroundSeeds = [backgroundSeeds; extractSeeds(I, 1, 1, cols, 1)];
    backgroundSeeds = [backgroundSeeds; extractSeeds(I, 1, 1, 1, rows)];
    backgroundSeeds = [backgroundSeeds; extractSeeds(I, 1, cols, 1, rows)];
%     backgroundSeeds = [backgroundSeeds; extractSeeds(I, rows, 1, cols, 1)];

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

    leftTranspose = CG(1: nrows, 2:ncols) - CG(1:nrows, 1:ncols-1);
    leftTranspose = abs(leftTranspose);
    leftTranspose = obj.normVal(leftTranspose, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
    leftTranspose = [zeros(nrows, 1) leftTranspose];
    leftTranspose = leftTranspose';

    rightTranspose = CG(1:nrows, 1:ncols-1) - CG(1:nrows, 2:ncols);
    rightTranspose = abs(rightTranspose);
    rightTranspose = obj.normVal(rightTranspose, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
    rightTranspose = [rightTranspose zeros(nrows, 1)];
    rightTranspose = rightTranspose';

    top = CG(2:nrows, :) - CG(1:nrows-1, :);
    top = abs(top);
    top = obj.normVal(top, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
    top = [zeros(1, ncols); top];

    bottom = CG(1:nrows-1, :) - CG(2:nrows, :);
    bottom = abs(bottom);
    bottom = obj.normVal(bottom, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
    bottom = [bottom; zeros(1, ncols)];
    CCD = single(Cs-Ct);
    
    labels = nppiGraphcut_32f8u_mex(cols, rows, CCD, CLT, CRT, CT, CB); 
    
    unique(labels)

%     uu = im2bw(uu, level);
    
    figure, imagesc(labels);

end

