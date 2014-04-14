% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

%function SvmSegm_get_bow_like_histograms_multiple_segments(mask_type, measurement_type, figure_or_ground, img_names, exp_dir)
%
% Goes through each segment, finds local features inside and projects them
% into a codebook. Then forms an histogram. 
function SvmSegm_get_bow_like_histograms_multiple_segments(mask_type, measurement_type, figure_or_ground, img_names, exp_dir, codebook_selection)
    codebook_dir = [exp_dir 'MyCodebooks/'];
    representation = 'bow';
    
    default_codebook_selection = ['kmeans_' measurement_type '_300_words.mat'];
    DefaultVal('codebook_selection', 'default_codebook_selection');    
    codebook_file = [codebook_dir codebook_selection];
   
    var = load(codebook_file, 'codebook');
    N_WORDS = size(var.codebook,2);
    codebook = var.codebook;
        
    for g=1:length(figure_or_ground)
        figure_ground = figure_or_ground{g};                                
        
        output_dir = [exp_dir 'MyMeasurements/' mask_type '_' representation '_' measurement_type '_' figure_ground '_' int2str(N_WORDS) '/'];
        if(~exist(output_dir, 'dir'))
            mkdir(output_dir);
        end
        
        %parfor (i=1:length(img_names),8)
        for i=1:length(img_names)
            i            
            
            measurement_file = [output_dir img_names{i} '.mat'];

            %if(exist(measurement_file, 'file'))
            %    continue;
            %end
            
            % 1. compute projections for all descriptors in every image
            var2 = myload([exp_dir 'MyMeasurements/' measurement_type '/' img_names{i} '.mat']);
            
            I = vl_ikmeanspush(uint8(var2.D), int32(codebook));
            
            % 2. count occurrences in masks
            masks_dir = [exp_dir 'MySegmentsMat/' mask_type '/'];
            
            var3 = load([masks_dir img_names{i} '.mat']); % masks
            if(iscell(var3.masks))
                masks = var3.masks{1};
            else
                masks = var3.masks;
            end
            
            if(strcmp(figure_ground, 'ground'))
                masks = ~masks;
            end
            
            filtI = cell(size(masks,3),1);
            %parfor (k=1:length(filtI), 8)
            for k=1:length(filtI)
                [filtI{k}] = filter_feats_outside_roi(I,var2.F,masks(:,:,k),'points');
            end
            Isegms = filtI;
            
            % 3. form histograms
            if(strcmp(representation, 'bow')) % bag of words
                D = cell2mat(form_bow_histograms(Isegms, N_WORDS));
            elseif(strcmp(representation, 'spatial_pyramid')) % spatial pyramid
                error('not ready');
                %Meas{i} = form_spatial_pyramid_histograms(F, I, n_features_per_image_mask, N_WORDS, img_dims{i});
            end
            
            mysave(measurement_file, D);
        end
    end
end

function v = myload(thepath)
    v = load(thepath, 'D', 'F');
end

function mysave(measurement_file,D)
    save([measurement_file], 'D');
end
