% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

%function SvmSegm_show_best_segments(I, Q, masks)
% I is the image
% Q has the qualities for each ground truth segment
% masks are the computed segments
function [best, max_score] = SvmSegm_show_best_segments(I, Q, masks)
    if(iscell(Q))
        Q = Q{1};
    end
    
    if(iscell(masks))
        masks = masks{1};
    end
    
    for i=1:numel(Q)
        [max_score(i), best_seg_id(i)] = max(Q(i).q);
        tit{i} = sprintf('%f', Q(i).q(best_seg_id(i)));
    end
    best = subplot_auto_transparent(masks(:,:,best_seg_id), I, tit);
    %figure;
    %sc(seg, 'rand')
end