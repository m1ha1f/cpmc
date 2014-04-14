% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function SvmSegm_extract_measurements_on_masks(mask_type, type, pars, img_names, exp_dir, to_save, MAX_DIM, overwrite)    
    DefaultVal('*overwrite', 'false', '*MAX_DIM', 'inf', 'to_save', 'true');  
    %disp('starting extract features');
    
    codebook_dir = [exp_dir 'MyCodebooks/'];
    assert(~iscell(type));
    
    if(strcmp(type, 'mask_phog') || strcmp(type, 'back_mask_phog')) % aspect ratio invariant
       type_func = @run_exp_database_do_phog;        
        if(~pars.withpb)
            nopb_str = 'nopb_';
        else 
            nopb_str = '';
        end
        dir_name = [mask_type '_' type '_'  nopb_str int2str(pars.n_ori) '_orientations_' int2str(pars.n_levels) '_levels/'];
    elseif(strcmp(type, 'back_mask_phog_scale_inv') || strcmp(type, 'mask_phog_scale_inv')) % aspect ratio variant
        type_func = @run_exp_database_do_scale_inv_phog;
        dir_name = [mask_type '_' type '_'  int2str(pars.n_ori) '_orientations_' int2str(pars.n_levels) '_levels/'];
    elseif(strcmp(type, 'bbox_phog_scale_inv'))
        type_func = @run_exp_database_do_scale_inv_phog;
        dir_name = [mask_type '_' type '_' 'nopb_' int2str(pars.n_ori) '_orientations_' int2str(pars.n_levels) '_levels/'];
    elseif(strcmp(type, 'bbox_phog'))        
        type_func = @run_exp_database_do_phog;
        dir_name = [mask_type '_' type '_' 'nopb_' int2str(pars.n_ori) '_orientations_' int2str(pars.n_levels) '_levels/'];
    elseif(strcmp(type, 'mask_regionprops'))
        type_func = @run_exp_database_do_simple;
    elseif(strcmp(type, 'figure_color_histogram'))
        type_func = @run_exp_database_do_ch;
        pars.figure_ground = 'figure';
        dir_name = [mask_type '_' type '/'];
    elseif(strcmp(type, 'ground_color_histogram'))
        type_func = @run_exp_database_do_ch;
        pars.figure_ground = 'ground';
        dir_name = [mask_type '_' type '/'];
    elseif(strcmp(type, 'figure_texton_histogram'))
        type_func = @run_exp_database_do_texton_hist;
        pars.figure_ground = 'figure';
        dir_name =  [mask_type '_' type '/'];
    elseif(strcmp(type, 'simple_segment_feats'))
        type_func = @run_exp_database_do_simple_segment_feats;
        dir_name = [mask_type '_' type '/'];
        pars.pb_dir = [exp_dir 'PB/'];
    elseif(strcmp(type, 'extended_segment_feats'))
        type_func = @run_exp_database_do_extended_segment_feats;
        dir_name = [mask_type '_' type '/'];
        pars.pb_dir = [exp_dir 'PB/'];
    elseif(strcmp(type, 'composition_segment_feats'))
        type_func = @run_exp_database_do_composition_segment_feats;
        dir_name = [mask_type '_' type '/'];

        pars.local_feats_path = [exp_dir 'MyMeasurements'];
        pars.local_feats{1} = 'dense_color_sift_3_scales';
        pars.local_feats{2} = 'dense_sift_4_scales';

        pars.region_feats_path = [exp_dir 'MyFeatures/' mask_type '/'];
        pars.region_feats{1}.figure_file = 'bow_dense_sift_4_scales_figure_300';
        pars.region_feats{1}.ground_file = 'bow_dense_sift_4_scales_ground_300';
        pars.region_feats{2}.figure_file = 'bow_dense_color_sift_3_scales_figure_300';
        pars.region_feats{2}.ground_file = 'bow_dense_color_sift_3_scales_ground_300';
        pars.region_feats{3}.figure_file = 'figure_color_histogram';
        pars.region_feats{3}.ground_file = 'ground_color_histogram';

      error('fix this, don''t want to pass img_list as parameter');
      for i=1:length(pars.region_feats)
          load([pars.region_feats_path train_test_val '__' pars.region_feats{i}.figure_file], 'Meas');
          Meas = Meas(img_list);
          pars.region_feats{i}.figure = Meas;

          load([pars.region_feats_path train_test_val '__' pars.region_feats{i}.ground_file], 'Meas');
          Meas = Meas(img_list);
          pars.region_feats{i}.ground = Meas;
          clear Meas;
      end
    elseif(strcmp(type, 'back_mask_local_shape_contexts') || strcmp(type, 'local_shape_contexts_boundary'))
%        codebook_file = [codebook_dir 'kmeans_mask_local_shape_contexts_300_words.mat'];
%         if(exist(codebook_file, 'file'))
%             vars = load(codebook_file, 'codebook');
%             pars.codebook = vars.codebook;
%         else
%             pars.codebook = [];
%         end
        if(~isempty(pars.codebook))
            var = load([exp_dir 'MyCodebooks/' pars.codebook]);            
            pars.codebook = var.codebook;
            dir_name = [mask_type '_' 'bow_' type '/'];
        else
            dir_name = [mask_type '_' type '/'];
        end
        
        type_func = @run_exp_database_do_local_shape_contexts;
        

    else       
        error('no such type defined');
    end
    
    
    t = tic();    
    
    % parfor
    for i=1:numel(img_names)     
    %parfor i=1:numel(img_names)     
        %t_image = tic();
        img_name = img_names{i};
        
        if ~overwrite
            if(to_save && exist([exp_dir 'MyMeasurements/' dir_name img_name '.mat'], 'file'))
                continue;
            end
        end
        
        img_file = [exp_dir 'JPEGImages/' img_name];
        I = imread([img_file '.jpg']);
        
        max_dim = max(size(I));
        if(max_dim>MAX_DIM)
            MAX_DIM
            img_name
            I = imresize(I, MAX_DIM/max_dim);
        end
        
        masks_file = [exp_dir 'MySegmentsMat/'  mask_type '/' img_name '.mat'];
        pb_path = [exp_dir 'PB/'  img_name '_PB.mat'];
    
        var = load(masks_file, 'masks');
        masks = var.masks;
        if(iscell(masks)) % if there are alternatives get the first one
            masks = masks{1};
        end
        
        the_pars = add_j_to_pars(pars, 1, img_name); % is this ok ? need to check it out
        
        my_all(exp_dir, I, type, masks, pb_path, the_pars, dir_name, img_name, type_func, to_save);   
        %toc(t_image);
    end
    %feats_on_masks_time_all_imgs = toc(t);
end

function [F, D] = my_all(exp_dir, I, type, masks, pb_path, the_pars, dir_name, img_name, type_func, to_save)
    if(isempty(masks))
        F = [];
        D = [];
    else    
        [F,D] = type_func(I, type, masks, pb_path, the_pars);
    end
    
    %%%% if save is false it assumes you only passed one image!
    if(to_save)
        if(~exist([exp_dir 'MyMeasurements/' dir_name], 'dir'))
            mkdir([exp_dir 'MyMeasurements/' dir_name]);
        end
        save([exp_dir 'MyMeasurements/' dir_name img_name '.mat'], 'F', 'D');
    else
        error('not ready for this');
        assert(numel(img_names) == 1);
    end           
end

% to be able to use parfor
function pars = add_j_to_pars(pars,j, img_name)
    pars.j = j;
    pars.img_name = img_name;	
end


function [the_bboxes] = square_boxes_from_masks(masks, I)
    n_masks = size(masks,3);
    the_bboxes = zeros(4, n_masks);    
    
    for k=1:n_masks    
        % creates fixed aspect ratio bounding boxes (setting it as  a square here)
        % bbox -- (ytop,ybottom,xleft,xright)
        
        %center = [(bbox(1)+bbox(2)) (bbox(3)+bbox(4))]./2;
        %imshow(I); hold on;plot(center(:,1), center(:,2), 'o')
        
        props = regionprops(double(masks(:,:,k)), 'BoundingBox');     
        if(isempty(props))
            bbox = [1 2 3 4];
        else            
            bbox(1) = props.BoundingBox(2); %ymin
            bbox(2) = bbox(1) + props.BoundingBox(4); %ymax
            bbox(3) = props.BoundingBox(1); % xmin
            bbox(4) = bbox(3) + props.BoundingBox(3); % xmax
            bbox = round(bbox);
        end
        
        % adds some extra space
        MARGIN = [10 10 10 10];        
        bbox(1) = max(bbox(1) - MARGIN(1), 1);
        bbox(2) = min(bbox(2) + MARGIN(2), size(I,1));
        bbox(3) = max(bbox(3) - MARGIN(3), 1);
        bbox(4) = min(bbox(4) + MARGIN(4), size(I,2));
             
        yaxis = bbox(2) - bbox(1);
        xaxis = bbox(4) - bbox(3);
        
        to_add = ceil(abs(xaxis-yaxis)/2);
        if(xaxis>yaxis)
            bbox(2) = bbox(2)+to_add;
            bbox(1) = bbox(1)-to_add;
        else
            bbox(4) = bbox(4)+to_add;
            bbox(3) = bbox(3)-to_add;
        end
        
        the_bboxes(:,k) = bbox';
    end
end

function [F,D] = run_exp_database_do_scale_inv_phog(I, type, masks, pb_path, pars)
    n_bins = pars.n_ori;
    angle = 180; % no gradient direction
    n_pyramid_levels = pars.n_levels;

    if(size(I,3) ~= 1)
        I = rgb2gray(I);
    end
    I = (double(I)/255.0);  %I = (double(I)/255.0).*mask;

    pb_file = pb_path;
         
    if iscell(masks)
        masks = cell2mat(masks);
    end
    
    [the_bboxes] = square_boxes_from_masks(masks, I);
    %assert(all([the_bboxes(2,:) - the_bboxes(1,:)] == [the_bboxes(4,:) - the_bboxes(3,:)]));

    if(strcmp(type, 'back_mask_phog_scale_inv') || strcmp(type, 'bbox_phog_scale_inv'))
        WITH_BBOX = false;
        if(strcmp(type, 'bbox_phog_scale_inv'))            
            WITH_BBOX = true;            
        end
        [D] = phog_backmasked(I, n_bins, angle, n_pyramid_levels, the_bboxes, masks, pb_file, WITH_BBOX);
    elseif(strcmp(type,'mask_phog_scale_inv'))
        [D] = phog_masked(masks, n_bins, angle, n_pyramid_levels, the_bboxes, pb_file);        
    end
    D = single(D);
    F = the_bboxes;
    F = single(F);
end

function [F,D] = run_exp_database_do_phog(I, type, masks, pb_path, pars)          
    n_bins = pars.n_ori;
    angle = 180; % no gradient direction
    n_pyramid_levels = pars.n_levels;

    if(size(I,3) ~= 1)
        I = rgb2gray(I);
    end
    I = (double(I)/255.0);  %I = (double(I)/255.0).*mask;

    pb_file = pb_path;
         
    if(isempty(masks))
        n_masks = 0;
    else
        n_masks = size(masks,3);
    end
    
    the_bboxes = zeros(4, n_masks);    
       
    MARGIN = [10 10 10 10];
    for k=1:n_masks    
        props = regionprops(double(masks(:,:,k)), 'BoundingBox');
        if(isempty(props))
            bbox(1:4) = [1 2 3 4];
        else
            bbox(1) = props.BoundingBox(2); %ymin
            bbox(2) = bbox(1) + props.BoundingBox(4); %ymax
            bbox(3) = props.BoundingBox(1); % xmin
            bbox(4) = bbox(3) + props.BoundingBox(3); % xmax
        end
        bbox = round(bbox);
        bbox(1) = max(bbox(1) - MARGIN(1), 1);
        bbox(2) = min(bbox(2) + MARGIN(2), size(I,1));
        bbox(3) = max(bbox(3) - MARGIN(3), 1);
        bbox(4) = min(bbox(4) + MARGIN(4), size(I,2));
                
        the_bboxes(:,k) = bbox';
    end

    if(strcmp(type, 'mask_phog'))
        [D] = phog_masked(masks, n_bins, angle, n_pyramid_levels, the_bboxes);
    elseif(strcmp(type, 'back_mask_phog') || strcmp(type, 'back_mask_phog_nopb') || strcmp(type, 'bbox_phog'))
        bbox = false;
        if(strcmp(type, 'bbox_phog'))
            bbox = true;
        end
        
        if(pars.withpb)
            [D] = phog_backmasked(I, n_bins, angle, n_pyramid_levels, the_bboxes, masks, pb_file, bbox);
        else
            [D] = phog_backmasked(I, n_bins, angle, n_pyramid_levels, the_bboxes, masks, [], bbox);
        end
    else
        error('no such type');
    end
    
    D = single(D);
    F = the_bboxes;
    F =  single(F);
end

function [F,D] = run_exp_database_do_ch(I, type, masks, pb_path, pars)
    % color histogram
    N_BINS = 14;
    n_masks = size(masks,3);
    
    if(size(I,3) == 1)
        I(:,:,2) = I(:,:,1);
        I(:,:,3) = I(:,:,1);
    end
    Ir = double(I(:,:,1));
    Ig = double(I(:,:,2));
    Ib = double(I(:,:,3));
    
    %newI = rgb2hsv(I);
    %H = newI(:,:,3);
  
    if(strcmp(pars.figure_ground, 'ground'))
        masks = ~masks;
    end

    D = zeros(3*N_BINS, n_masks);
    for k=1:n_masks
        mask = masks(:,:,k);
        D(1:N_BINS,k) = hist(Ir(mask), N_BINS);
        D(N_BINS+1:2*N_BINS,k) = hist(Ig(mask), N_BINS);
        D((2*N_BINS) +1:3*N_BINS,k) = hist(Ib(mask), N_BINS);
    end
    F = [];
end

function [F,D] = run_exp_database_do_texton_hist(I, type, masks, pb_path, pars)
    % texton histogram
    N_BINS = 64;
    n_masks = size(masks,3);
    
    % load pb
    textons = myload(pb_path, 'textons');
    
    D = zeros( N_BINS, n_masks);
    
    if(strcmp(pars.figure_ground, 'ground'))
        masks = ~masks;
    end

    for k=1:n_masks
        mask = masks(:,:,k);
        D(:,k) = hist(double(textons(mask)), N_BINS); 
    end
    D = single(D);
    F = [];
end

function [F,D] = run_exp_database_do_local_shape_contexts(I, type, masks, pb_path, pars)
    if(size(I,3) >1)
        I = rgb2gray(I);
    end
   
    %load(pb_path); % not doing back_mask_local_shape_contexts yet

    nbins_theta = pars.theta; % 20
    nbins_r = pars.r;    % 8
    r_inner = pars.r_inner; %0.05, 0.01
    r_outer = pars.r_outer; %0.5 0.4
    
    %if(isempty(pars.codebook))
    %    error('not ready for that');
    %end
    
    if(~isempty(pars.codebook))        
        D = zeros(size(pars.codebook,2), size(masks,3));
    else
        D = [];
    end
    
    % load quality
    if isfield(pars,'quality_dir')
      load([pars.quality_dir pars.img_name '.mat']);
      [quals inds] = sort(Quality.q,'descend');
    else
      inds = 1:size(masks,3);
    end
    
%     pb = load(pb_path);
%     GrayI = pb.gPb_thin>10;
    for i=1:size(masks,3)
        %I_this = I.*uint8(masks(:,:,i));        
        
        if(strcmp(type, 'back_mask_local_shape_contexts'))
            % sample from outer contour, but consider inner edges
            NumPts = pars.NumPts; 
            
            
%             EdgeI = masks(:,:,inds(i)) .* GrayI; % apply mask to pb
            GrayI = I.*(uint8(masks(:,:,inds(i)))); % apply mask
            GrayI(GrayI == 0) = 255 ; % don't know why
            
            % canny
            EdgeI = edge(GrayI,'canny'); % using PB so no need for it
            [PRows PCols] = find(EdgeI == 1);
            NumPixels = size(PRows,1);            
            
            Interval = NumPixels/NumPts ;
            interV = 1 ;
            for iCount=1:NumPts            
                Bsamp(1,iCount) = PRows(floor(interV),1) ;
                Bsamp(2,iCount) = PCols(floor(interV),1) ;
                interV = interV + Interval ;            
            end
            

            % Now Compute the Shape Context Histogram
            out_vec = zeros(1, NumPts) ;
            Tsamp = zeros(1, NumPts) ;
            mean_dist_global = [] ;
            [desc, mean_dist_l] = sc_compute(Bsamp,Tsamp,mean_dist_global,nbins_theta,nbins_r, r_inner,r_outer,out_vec);

            if(~isempty(pars.codebook))
              proj = vl_ikmeanspush(uint8(desc)', int32(pars.codebook));
              D(:,i) = vl_ikmeanshist(length(pars.codebook), proj);
            else
              % if we don't have a codebook just retrieve some of the shape contexts 
                D = [D desc'];
                if(i==10) % this means only 10 masks will be used
                  break;
                end
            end
            F = [];
        elseif(strcmp(type, 'local_shape_contexts_boundary'))
          coords = bwboundaries(masks(:,:,inds(i)));
          coords = coords{1}';
          coords = coords(:,1:pars.SAMPLING_RATE:end);

          
          Bsamp = coords;
          Tsamp = zeros(1,size(Bsamp,2));    
          out_vec = Tsamp;
            % sample from outer contour, don't consider inner edges
            desc = sc_compute(Bsamp,Tsamp,[],nbins_theta,nbins_r,r_inner,r_outer, out_vec);
            if(~isempty(pars.codebook))
                proj = vl_ikmeanspush(uint8(desc)', int32(pars.codebook));
                D(:,i) = vl_ikmeanshist(length(pars.codebook), proj);
            else
                % if we don't have a codebook just retrieve some of the shape contexts 
                D = [D desc'];
                %if(i==10) % up to 10
                %    break;
                %end
            end
            F = [];
        else
            error('no such type');
        end
        %imshow(I_this);
        %pause;
    end
end

% features for learning to segment
function [F, D] = run_exp_database_do_simple_segment_feats(I, type, masks, pb_path, pars)
    % A. pairwise affinity features:
    % 1 - cut ratio ( sum of edges across cut, divided by their number )
    % 2 - cut
    % 3 - normalized cut
    % 4 - unbalanced normalized cut
    % 
    % B. Area features
    % from regionprops (17 feats)
           
    img_name = pars.img_name 
    %img_name = pb_path(strfind(pb_path, 'PB/')+3:end);
    %img_name = pb_path(1:end-7);
    
    
    
    
    %s = LongRangeSegmenter(I,img_name);    
    %s = s.set_pb_path(pars.pb_dir);
    %s = s.initialize();
 
    %img_dgraph = s.p_hgraph.pairwise_to_dgraph();
    resh_masks = reshape(masks, size(masks,1) * size(masks,2), size(masks,3));

    no_separate = true;
    S = SegmentProcessor([], resh_masks, I, pb_path,  25, no_separate);
    get_all = true;
    S = S.compute_energies(get_all);

    %
    %%%% Cut features %%%%
    %
    D_cut = zeros(8, size(masks,3));
    D_cut(1,:) = [S.energies(:).cut_ratio];
    D_cut(2,:) = [S.energies(:).cut];
    D_cut(3,:) = [S.energies(:).normalized_cut];
    D_cut(4,:) = [S.energies(:).unbalanced_normalized_cut];
    tmp_var = [S.energies(:).fraction_of_healthy_boundary];
    D_cut(5:end,:) = reshape(tmp_var, 4,size(tmp_var,2)/4);
    
    %
    %%%% Coarse Region Shape and Location features %%%%%%
    %
    D_crsl = zeros(19, size(masks,3)); 
    
    % absolute quantities
    for i=1:size(masks,3)
        props = regionprops(masks(:,:,i), 'Area', 'Centroid', 'BoundingBox', 'MajorAxisLength', ...
                                            'MinorAxisLength', 'Eccentricity', 'Orientation', 'ConvexArea', 'EulerNumber', 'EquivDiameter', 'Solidity', 'Extent', 'Perimeter');
        props = props(1);  % might be more than one, describe the first
        D_crsl(1:17,i) = struct2array(props(1)); 
                
        % convexity
        D_crsl(18, i) = props.Area / props.ConvexArea; % solidity is the same as this one
        % absolute distance to center of image
        D_crsl(19, i) = norm([size(I,1) size(I,2)]./2 - props.Centroid);
    end    
    
    
    D = [D_cut; D_crsl];
    F = size(I)';
end

function [F, D] = run_exp_database_do_extended_segment_feats(I, type, masks, pb_path, pars)
    % 
    im = double(I)/255.0;

    var = load(pb_path, 'textons');
    textons = double(var.textons);

%    [ ...
%    textons, ...
%    bg_r3, bg_r5,  bg_r10,  cga_r5, cga_r10, cga_r20, cgb_r5, cgb_r10, cgb_r20, tg_r5,  tg_r10,  tg_r20...
%    ] = mex_pb_parts_final_selected(im(:,:,1),im(:,:,2),im(:,:,3));

    %
    %%% histograms of textons
    %
    t = reshape(textons, numel(textons), 1);
    m = reshape(masks,numel(textons), size(masks,3));
    N = 65;
    hist_fg = zeros(N, size(m,2));
    hist_bg = zeros(N, size(m,2));
    t = t+1; % to avoid t==0, that has problem with int_hist
    for i=1:size(m,2)
        hist_fg(:,i) = int_hist(t(m(:,i)), N)';
        hist_bg(:,i) = int_hist(t(~m(:,i)), N)';
    end
    hist_fg = scale_data(hist_fg, 'norm_1');
    hist_bg = scale_data(hist_bg, 'norm_1');
    D = chi2_mex(single(hist_fg), single(hist_bg));
    
    %
    %%% histograms of brightness
    %
    N_BINS = 256;
    hist_fg_b = zeros(N_BINS, size(m,2));
    hist_bg_b = zeros(N_BINS, size(m,2));
    Igray = rgb2gray(I)+1;    
    Igray_r = double(reshape(Igray, numel(Igray), 1));
    for i=1:size(m,2)        
        hist_fg_b(:,i) = int_hist(Igray_r(m(:,i)), N_BINS)';
    end
    hist_fg_b = scale_data(double(hist_fg_b), 'norm_1');
    hist_bg_b = scale_data(double(hist_bg_b), 'norm_1');
    D_b = chi2_mex(single(hist_fg_b), single(hist_bg_b));          
    
    %%% contour energy
    load(pb_path); % gPb_thin
    gPb_thin_c = reshape(gPb_thin, numel(gPb_thin),1);
    
    all_bw = imdilate(masks, ones(3,3)) & ~masks;
    n = zeros(1,size(m,2));
    s = n;
    s_intra = n;
    for i=1:size(m,2)
%         t = tic();
%         bw{i} = bwboundaries(masks(:,:,i));
%         n(i) = size(bw{i}{1}, 1);
%         ids = bw{i}{1};
%         time_bwbound = toc(t)
        
        %t = tic();
        [ids1, ids2] = find(all_bw(:,:,i));
        ids = [ids1 ids2];
        n(i) = size(ids,1);
        %time_dilate = toc(t)
        
        s(i) = sum(gPb_thin( sub2ind(size(gPb_thin), ids(:,1), ids(:,2))));
        s_intra(i) = sum(gPb_thin_c(m(:,i)));         
    end
    %gPb = mean(gPb_orient,3);
    
    %%% curvilinear continuity 
    %%%( curvature using bruckstein discrete approximation )
    the_lines = cell(size(m,2),1);
    for i=1:size(m,2)
        %BREAK_FACTOR = 0.9;
        %thresh = BREAK_FACTOR*0.05*log(size(bw{i}{1},1))/log(1.1);
        thresh = 1;
        lines = lineseg({ids}, thresh);
        the_lines{i} = seglist2segarray(lines);        
        %jl_plot_lines(the_lines{i});
        %pause;
    end
    
    sum_curvature = zeros(size(m,2),1);
    for i=1:size(m,2) 
        if(size(the_lines{i},2) == 1)
            sum_curvature(i) = 0;
            continue;
        end
        
        curvatures = zeros(size(the_lines{i},2)-1,1);
        curvatures(1) = angle_between_linesegs(the_lines{i}(:,end), the_lines{i}(:,1));
        for j=1:size(the_lines{i},2)-1
            curvatures(j) = (angle_between_linesegs(the_lines{i}(:,j), the_lines{i}(:,j+1)))^2;      
            len_1 = norm(the_lines{i}(1:2,j) - the_lines{i}((3:4),j));
            len_2 = norm(the_lines{i}((1:2),j+1) - the_lines{i}(((3:4)),j+1));
            curvatures(j) = curvatures(j)/(min(len_1,len_2));
        end
        sum_curvature(i) = sum(curvatures); 
    end
        
    % 1. inter-region dissimilarity ( big is good )
    D_chi2_fg_bg = diag(D); % this seems a nice a feature, could grow it, by doing a new feature with the relative value to the best of the image
    
    % 2. intra-region similarity (simplicity, small is good)       
    D_fg_simplicity = sum(hist_fg>(1/300))'; 
    
    % 3. inter region brightness similarity
    D_chi2_fg_bg_b = diag(D_b);
    
    % 3. b) brigh
    % 4. intra-region brightness similarity (simplicity, small is good)
    D_fg_b_simplicity = sum(hist_fg_b>(1/N_BINS))';
    
    % 5. inter-region contour energy (similar to ratio cut)
    D_fg_bg_contour_energy = (s./n)';
    
    % 6. intra-region contour_energy
    D_fg_bg_intra_contour_energy = (s_intra ./sum(m))';
    
    % 7. curvilinear continuity
    D_cc = sum_curvature;
            
    F = [];
    D = [D_chi2_fg_bg D_fg_simplicity D_chi2_fg_bg_b D_fg_b_simplicity D_fg_bg_contour_energy D_fg_bg_intra_contour_energy D_cc]';
end
