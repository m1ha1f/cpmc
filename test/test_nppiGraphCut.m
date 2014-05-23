function test_nppiGraphCut()
    addpath('../code');
    addpath('../external_code');
    addpath('../external_code/npp');
    addpath('../external_code/paraFmex');
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
    foregroundids = [ix(:)+prow-1 iy(:)+pcol-1];
    foregroundids = sub2ind([rows cols], foregroundids(:, 1), foregroundids(:, 2));

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
    

    for i = 0
        ct = CTerminals;
        lambda = i/10;
        ct(foregroundids) = ct(foregroundids) + lambda;
        ct(backgroundids) = ct(backgroundids) - lambda;
        labels = nppiGraphcut_32f8u_mex(cols, rows, ct, leftTranspose, rightTranspose, top, bottom); 
        figure, imagesc(labels);
    end

    N = sparse(rows*cols+2, rows*cols+2);
    s = rows*cols+1;
    t = rows*cols+2;

    Srows = [];
    Scols = [];
    Svals = [];
    for i = 1 : rows
        for j = 1: cols
            ci = sub2ind([rows cols], i, j);
            if i > 1
                ct = sub2ind([rows cols], i-1, j);
%                 N(ci, ct) = top(i, j);
                Srows = [Srows ci];
                Scols = [Scols ct];
                Svals = [Svals top(i, j)];
            end
            if i < rows
                cb = sub2ind([rows cols], i+1, j);
%                 N(ci, cb) = bottom(i, j);
                Srows = [Srows ci];
                Scols = [Scols cb];
                Svals = [Svals bottom(i, j)];
            end
            if j > 1
                cl = sub2ind([rows cols], i, j-1);
%                 N(ci, cl) = leftTranspose(j, i);
                Srows = [Srows ci];
                Scols = [Scols cl];
                Svals = [Svals leftTranspose(j, i)];
            end
            if j < cols
                cr = sub2ind([rows cols], i, j+1);
%                 N(ci, cr) = rightTranspose(j, i);
                Srows = [Srows ci];
                Scols = [Scols cr];
                Svals = [Svals rightTranspose(j, i)];
            end
            
            Srows = [Srows s];
            Scols = [Scols ci];
            Svals = [Svals Cs(i, j)];
            
            Srows = [Srows ci];
            Scols = [Scols t];
            Svals = [Svals Ct(i, j)];
        end
    end 
    
    N = sparse(Srows, Scols, double(Svals), rows*cols+2, rows*cols+2);

    lambda_edges = [ones(width*height, 1)*s foregroundids];
    lambda_edges = [lambda_edges; backgroundids' ones(size(backgroundids, 2), 1)*t];

    lambda_weights = ones(size(lambda_edges, 1), 1);
    lambda_offsets = Cs(foregroundids);
    lambda_offsets = double([lambda_offsets; Ct(backgroundids)']);

    [cuts, lambdas] = hoch_pseudo_par_mex(N, lambda_edges, lambda_weights, lambda_offsets, s, t, -1, 150, 20);
    
    for i = 1:size(cuts, 2)
        img = reshape(cuts(1:rows*cols, i), rows, cols);
        figure, imagesc(img);
    end
    

end

