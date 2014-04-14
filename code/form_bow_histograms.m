% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function H = form_bow_histograms(I, N_WORDS)
%  n_feats_so_far = 0;
%   H = sparse(N_WORDS, length(n_features_per_image));
%   
%   for j=1:length(n_features_per_image)
%     H(:,j) = vl_ikmeanshist(N_WORDS,I(n_feats_so_far+1:n_feats_so_far+n_features_per_image(j)));
%     n_feats_so_far = n_feats_so_far + n_features_per_image(j);
%   end
    H = cell(1,length(I));
    %parfor(k=1:length(I), 8)
    for k=1:length(I)
        H{k} = int32(vl_ikmeanshist(N_WORDS,I{k}));       
    end
end
