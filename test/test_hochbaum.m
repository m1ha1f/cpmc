function [cuts, lambdas] = test_hochbaum(I, K, foregroundids, foregroundSeeds, backgroundids, backgroundSeeds, l, u)
    rows = size(I, 1);
    cols = size(I, 2);

    Cs = computeCapacity(I, foregroundSeeds);
    Ct = computeCapacity(I, backgroundSeeds);

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
    
    foregroundids = int32(foregroundids);
    backgroundids = int32(backgroundids);

    s = rows*cols+1;
    t = rows*cols+2;

    n_edges = 2*(3*rows*cols-rows-cols);
       
    [tx1, ty1] = meshgrid(2:rows, 1:cols);
    [tx2, ty2] = meshgrid(1:rows-1, 1:cols);
    t1 = sub2ind([rows cols], tx1(:), ty1(:));
    t2 = sub2ind([rows cols], tx2(:), ty2(:));
    tv = top(2:rows, :);
    tv = tv(:);
    
    [bx1, by1] = meshgrid(1:rows-1, 1:cols);
    [bx2, by2] = meshgrid(2:rows, 1:cols);
    b1 = sub2ind([rows cols], bx1(:), by1(:));
    b2 = sub2ind([rows cols], bx2(:), by2(:));
    bv = bottom(1:rows-1, :);
    bv = bv(:);
    
    [lx1, ly1] = meshgrid(1:rows, 2:cols);
    [lx2, ly2] = meshgrid(1:rows, 1:cols-1);
    l1 = sub2ind([rows cols], lx1(:), ly1(:));
    l2 = sub2ind([rows cols], lx2(:), ly2(:));
    left = leftTranspose';
    lv = left(:, 2:cols);
    lv = lv(:);
    
    [rx1, ry1] = meshgrid(1:rows, 1:cols-1);
    [rx2, ry2] = meshgrid(1:rows, 2:cols);
    r1 = sub2ind([rows cols], rx1(:), ry1(:));
    r2 = sub2ind([rows cols], rx2(:), ry2(:));
    right = rightTranspose';
    rv = right(:, 1:cols-1);
    rv = rv(:);
    Srows = [t1; b1; l1; r1; s*ones(rows*cols, 1);    (1:rows*cols)'];
    Scols = [t2; b2; l2; r2; (1:rows*cols)';          t*ones(rows*cols, 1)];
    Svals = [tv; bv; lv; rv; Cs(:);                   Ct(:) ];
    
    
    N = sparse(Srows, Scols, double(Svals), rows*cols+2, rows*cols+2);

    lambda_edges = [ones(numel(foregroundids), 1)*s foregroundids];
    lambda_edges = [lambda_edges; backgroundids' ones(size(backgroundids, 2), 1)*t];
    lambda_edges = double(lambda_edges);

    lambda_weights = ones(size(lambda_edges, 1), 1);
    lambda_offsets = Cs(foregroundids);
    lambda_offsets = double([lambda_offsets; Ct(backgroundids)']);

%     t0 = tic();
    [cuts, lambdas] = my_hoch_pseudo_par_mex(N, lambda_edges, lambda_weights, lambda_offsets, s, t, l, u, K);
%     t = toc(t0);


    cuts = cuts(1:rows*cols, :);
    cuts = ~cuts;
    % fprintf('Hochbaum: %f\n', t);
    
    % for i = 1:size(cuts, 2)
    %     img = reshape(cuts(1:rows*cols, i), rows, cols);
    %     figure, imagesc(~img);
    % end
end