% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function [masks] = cpmc_masks(exp_dir, img_name, segm_pars)    
    DefaultVal('*segm_pars', '[]');
    
    img_folder = [exp_dir 'JPEGImages/'];   
    
    % extract initial pool
    dir_name = [exp_dir 'MySegmentsMat/' segm_pars.name];
    if(~exist([exp_dir 'MySegmentsMat/' segm_pars.name], 'dir'))
        mkdir(dir_name);
    end
    
    filename = [dir_name '/' img_name '.mat'];
    if(~exist(filename, 'file'))
        [masks] = SvmSegm_extract_segments_nosave(img_folder, {img_name}, segm_pars);
        masks = masks{1};
        save(filename, 'masks');
    else
        var = load(filename);
        masks = var.masks;
    end
end