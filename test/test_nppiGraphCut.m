function [cuts, lambdas] = test_nppiGraphCut(I, K, foregroundids, foregroundSeeds, backgroundids, backgroundSeeds, l, u)
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
    
    CTerminals = single(Cs-Ct);
    
    lambda_range =  [logspace(log10(l), log10(u), K)];
%     lambda_range(lambda_range > (u+eps)) = [];
%     lambda_range(l > lambda_range) = [];
%     lambda_range = [0 lambda_range];
    % lambda_range = single(lambda_range)/1000;
    lambda_range = int32(round(lambda_range*1000));
    
    foregroundids = int32(foregroundids);
    backgroundids = int32(backgroundids);

    CTerminals = int32(round(CTerminals*1000));
    leftTranspose = int32(round(leftTranspose*1000));
    rightTranspose = int32(round(rightTranspose*1000));
    top = int32(round(top*1000));
    bottom = int32(round(bottom*1000));
    
    t0 = tic();
%     [cuts, lambdas] = nppiGraphcut_32s8u_multi_mex(cols, rows, CTerminals, leftTranspose, rightTranspose, top, bottom, ...
%         K, lambda_range, numel(foregroundids), foregroundids', numel(backgroundids), backgroundids');
    [cuts, lambdas] = nppiGraphcut_32s8u_multi_mex(CTerminals, leftTranspose, rightTranspose, top, bottom, ...
                        lambda_range, foregroundids, backgroundids);
    t = toc(t0);
    % fprintf('nppiGraphCut: %f\n', t);
end

