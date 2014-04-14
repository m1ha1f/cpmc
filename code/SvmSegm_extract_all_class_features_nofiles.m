% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function D = SvmSegm_extract_all_class_features_nofiles(exp_dir, img_name, mask_type, class_feats, delete_local_feats)
   DefaultVal('delete_local_feats', 'true');
   
   img_name = {img_name};
   
   [meas_req, needs_shape, figure_ground, the_pars] = SvmSegm_get_measurements_required(class_feats);   
   
   %
   % Histograms of gradients
   %
   
   parfor i=1:numel(meas_req)
       filename = [exp_dir 'MyMeasurements/' mask_type '_' class_feats{i} '/' img_name{1} '.mat'];
       if(~needs_shape(i))
           filename = [exp_dir 'MyMeasurements/' meas_req{i} '/' img_name{1} '.mat'];
           if(~exist(filename, 'file'))
               SvmSegm_extract_measurements_new(meas_req{i}, img_name, exp_dir, inf);
           end
       else
           if(~exist(filename, 'file'))
               SvmSegm_extract_measurements_on_masks(mask_type, meas_req{i}, the_pars{i}, img_name, exp_dir);
           end
       end
   end
   
   
   %
   % Bag of words
   %
   parfor i=1:numel(class_feats)
       if(the_pars{i}.isbow)
           SvmSegm_get_bow_like_histograms_multiple_segments(mask_type, meas_req{i}, {figure_ground{i}}, img_name, exp_dir);
       end
   end
   
   if(delete_local_feats)
       for i=1:numel(class_feats)
           if(the_pars{i}.isbow)
               filename = [exp_dir 'MyMeasurements/' meas_req{i} '/' img_name{1} '.mat'];
               if(exist(filename, 'file'))
                   delete([exp_dir 'MyMeasurements/' meas_req{i} '/' img_name{1} '.mat']);
               end
           end
       end
   end
   
   for i=1:numel(class_feats)
       var = load([exp_dir 'MyMeasurements/' mask_type '_' class_feats{i} '/' img_name{1} '.mat']);
       D{i} = var.D;
   end
end
