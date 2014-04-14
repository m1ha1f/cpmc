% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function Q = SvmSegm_segment_quality(img_name, exp_dir, masks, segm_quality_type)             
    if iscell(masks)
        masks = cell2mat(masks);
    end

    name = [exp_dir 'SegmentationObject/' img_name '.png'];

    ground_truth_obj_segs = imread(name);

    un = unique(ground_truth_obj_segs);
    un(un==0) = [];
    un(un==255) = [];
    
    care = (ground_truth_obj_segs~=255);    
    parfor k=1:numel(un)
        ground_truth_k = zeros(size(ground_truth_obj_segs));
        ground_truth_k(ground_truth_obj_segs == un(k)) = 1;

        [duh1, duh2, this_Q] = myCalcCandScoreFigureGroundAll(masks,ground_truth_k, segm_quality_type, care);
        Q(k).q = this_Q;
    end
end
