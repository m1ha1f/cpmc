% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function pred = predict_regressor(Feats, w)
    assert(~any(isnan(w)));
    
    pred = w' * [ones(1,size(Feats,2)); Feats];
    % inverse logit transformation
    %pred = 1 - 1./(exp(pred) + 1);
end
