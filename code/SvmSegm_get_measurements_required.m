% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function [meas_req, meas_need_shape, figure_ground, the_pars] = SvmSegm_get_measurements_required(class_features)
   the_pars = cell(numel(class_features),1);
   for i=1:numel(class_features)
       switch(class_features{i})
           case 'signature'
               meas_req{i} = 'signature';
               meas_need_shape(i) = true;               
               figure_ground{i} = 'figure';               
           case 'tps_hog'
               meas_req{i} = 'tps_hog';
               meas_need_shape(i) = true;               
               figure_ground{i} = 'figure';
               the_pars{i}.n_ori = 10;
               the_pars{i}.n_levels = 3;
               the_pars{i}.withpb = true;
               the_pars{i}.isbow = false;
           case 'mask_phog_scale_inv_10_orientations_3_levels'
               meas_req{i} = 'mask_phog_scale_inv';
               meas_need_shape(i) = true;               
               figure_ground{i} = 'figure';
               the_pars{i}.n_ori = 10;
               the_pars{i}.n_levels = 3;
               the_pars{i}.withpb = true;
               the_pars{i}.isbow = false;
           case 'mask_phog_scale_inv_20_orientations_2_levels'
               meas_req{i} = 'mask_phog_scale_inv';
               meas_need_shape(i) = true;               
               figure_ground{i} = 'figure';
               the_pars{i}.n_ori = 20;
               the_pars{i}.n_levels = 2;
               the_pars{i}.withpb = true;          
               the_pars{i}.isbow = false;
           case 'back_mask_phog_nopb_20_orientations_3_levels'
               meas_req{i} = 'back_mask_phog';
               meas_need_shape(i) = true;
               figure_ground{i} = 'figure';
               the_pars{i}.n_ori = 20;
               the_pars{i}.n_levels = 3;    
               the_pars{i}.withpb = false;
               the_pars{i}.isbow = false;
           case 'back_mask_phog_20_orientations_3_levels'
               meas_req{i} = 'back_mask_phog';
               meas_need_shape(i) = true;
               figure_ground{i} = 'figure';
               the_pars{i}.n_ori = 20;
               the_pars{i}.n_levels = 3;    
               the_pars{i}.withpb = true;
               the_pars{i}.isbow = false;               
           case 'bbox_phog_scale_inv_nopb_10_orientations_3_levels'
               meas_req{i} = 'bbox_phog_scale_inv';
               meas_need_shape(i) = true;
               figure_ground{i} = 'figure';
               the_pars{i}.n_ori = 10;
               the_pars{i}.n_levels = 3;
               the_pars{i}.withpb = false;
               the_pars{i}.isbow = false;                 
           case 'bbox_phog_nopb_10_orientations_3_levels'
               meas_req{i} = 'bbox_phog';
               meas_need_shape(i) = true;
               figure_ground{i} = 'figure';
               the_pars{i}.n_ori = 10;
               the_pars{i}.n_levels = 3;
               the_pars{i}.withpb = false;
               the_pars{i}.isbow = false;
           case 'bow_dense_sift_4_scales_figure_300'
               meas_req{i} = 'dense_sift_4_scales';
               meas_need_shape(i) = false;
               figure_ground{i} = 'figure';
               the_pars{i}.isbow = true;
           case 'bow_dense_sift_4_scales_ground_300'
               meas_req{i} = 'dense_sift_4_scales';
               meas_need_shape(i) = false;
               figure_ground{i} = 'ground';
               the_pars{i}.isbow = true;
           case 'bow_dense_color_sift_3_scales_figure_300'
               meas_req{i} = 'dense_color_sift_3_scales';
               meas_need_shape(i) = false;
               figure_ground{i} = 'figure';
               the_pars{i}.isbow = true;               
           case 'bow_dense_color_sift_3_scales_ground_300'
               meas_req{i} = 'dense_color_sift_3_scales';
               meas_need_shape(i) = false;
               figure_ground{i} = 'ground';
               the_pars{i}.isbow = true;
           case 'mask_local_shape_contexts'
               meas_req{i} = 'local_shape_contexts_boundary';
               meas_need_shape(i) = true;
               figure_ground{i} = 'boundary';
               the_pars{i}.isbow = false;
               the_pars{i}.codebook = 'kmeans_SegmenterLongRangeSegmenter_mask_local_shape_contexts_300_words';
           case 'idsc_boundary'
               meas_req{i} = 'idsc_boundary';
               meas_need_shape(i) = true;
               figure_ground{i} = 'boundary';      
               the_pars{i}.isbow = false;
       end
   end
end