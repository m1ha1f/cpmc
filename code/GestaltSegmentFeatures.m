% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef GestaltSegmentFeatures < SegmFeatures
    properties 
        pb_file
        gPb_thin
        resize_factor
    end
    
    methods
        function obj = GestaltSegmentFeatures(I, segments, pb_file)
            obj = obj@SegmFeatures(I,segments);
            obj.resize_factor = 0.5;
            obj.original_I = I;            
            
            obj.pb_file = pb_file;            
            var =load(obj.pb_file, 'gPb_thin');
            obj.gPb_thin = var.gPb_thin;            
        end                       
        
        function feats = get_feats(obj)            
            feat_names = fieldnames(obj.features);
            for i=1:numel(feat_names)
                n_times(i) = size(obj.features.(feat_names{i}),1);
            end
            n_feats = sum(n_times);
            n_segms = size(obj.segments,3);
            
            feats = zeros(n_feats, n_segms);
            
            counter = 1;
            for i=1:numel(feat_names)                
                feats(counter:(counter+n_times(i) -1), :) = full(obj.features.(feat_names{i}));
                counter = counter+ n_times(i);
            end
            
            feats = single(feats);
        end
        
        %%%%%%%%%%%%%%%%%%%%% feature computation %%%%%%%%%%%%%%%%%%%%

        function obj = compute_all_feats(obj)
            masks = obj.segments;
            n_masks = size(masks,3);

            % adapt neighborhood for the unary case!                                    
            obj.nb_ids = (1:n_masks)';
            obj.nb_ids = [obj.nb_ids obj.nb_ids+n_masks];                                                           
            
            %
            %%%% bla bla features %%%%
            %
                            
            m = reshape(masks,size(masks,1)*size(masks,2), size(masks,3));            

            if(isfield(obj.features, 'inter_texton_sim') || isfield(obj.features, 'intra_texton_sim'))                
                var = load(obj.pb_file, 'textons');
                textons = double(var.textons);
                %
                %%% histograms of textons
                %
                t = reshape(textons, numel(textons), 1);               
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
                D = chi2_mex(single(hist_fg), single(hist_bg), false);
                
                if(isfield(obj.features, 'inter_texton_dissim'))
                    % 1. inter-region dissimilarity ( big is good )
                    obj.features.inter_texton_dissim = diag(D)'; % this seems a nice a feature, could grow it, by doing a new feature with the relative value to the best of the image
                end
                
                if(isfield(obj.features, 'intra_texton_sim'))
                    % 2. intra-region similarity (simplicity, small is good)
                    obj.features.intra_texton_sim = sum(hist_fg>(1/300));
                end
            end
            
            if(isfield(obj.features, 'inter_brightness_sim') || isfield(obj.features, 'intra_brightness_sim'))
                %
                %%% histograms of brightness
                %
                N_BINS = 256;
                hist_fg_b = zeros(N_BINS, size(m,2));
                hist_bg_b = zeros(N_BINS, size(m,2));
                Igray = rgb2gray(obj.I)+1;    
                Igray_r = double(reshape(Igray, numel(Igray), 1));
                for i=1:size(m,2)        
                    hist_fg_b(:,i) = int_hist(Igray_r(m(:,i)), N_BINS)';
                    hist_bg_b(:,i) = int_hist(Igray_r(~m(:,i)), N_BINS)';
                end
                hist_fg_b = scale_data(double(hist_fg_b), 'norm_1');
                hist_bg_b = scale_data(double(hist_bg_b), 'norm_1');
                %D_b = pdist2(single(hist_fg_b)', single(hist_bg_b)', 'L1');          
                D_b = chi2_mex(single(hist_fg_b), single(hist_bg_b), false);          

                if(isfield(obj.features, 'inter_brightness_dissim'))
                    % 3. inter region brightness similarity
                    obj.features.inter_brightness_dissim = diag(D_b)';
                end
                
                if(isfield(obj.features, 'intra_brightness_sim'))
                    % 4. intra-region brightness similarity (simplicity, small is good)
                    obj.features.intra_brightness_sim = sum(hist_fg_b>(1/N_BINS))/N_BINS;
                end
            end                        
            
            if(isfield(obj.features, 'inter_contour_energy') || isfield(obj.features, 'intra_contour_energy'))
                %%% contour energy            
                load(obj.pb_file); % gPb_thin
                gPb_thin_c = reshape(gPb_thin, numel(gPb_thin),1);

                area = sum(reshape(masks,size(masks,1)*size(masks,2), size(masks,3)),1);
                
                % for robustness (?) grow the edges
                all_bw_len=false(size(masks));
                all_bw=false(size(masks));
                parfor i=1:size(masks,3)
                    all_bw_len(:,:,i) = imdilate(masks(:,:,i), ones(3,3)) & ~masks(:,:,i);
                    all_bw(:,:,i) = imdilate(masks(:,:,i), ones(5,5)) & imdilate(~masks(:,:,i), ones(5,5));
                end
                
                %all_bw = imdilate(masks, ones(5,5)) & imdilate(~masks, ones(5,5));
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
                    
                    a = find(all_bw_len(:,:,i));
                    n(i) = size(a,1);
                    %time_dilate = toc(t)
                    
                    %%% visualize ( ids are mislocalized by 1 or 2 pixels!
                    %sc(gPb_thin); hold on; plot(ids(:,2), ids(:,1), 'o')
                    s(i) = sum(gPb_thin( sub2ind(size(gPb_thin), ids(:,1), ids(:,2))));
                    s_intra(i) = sum(gPb_thin_c(m(:,i)))/area(i);         
                end
                %gPb = mean(gPb_orient,3);


                
                obj.features.inter_contour_energy = s;
                obj.features.inter_contour_energy_norm = (s./n);
                obj.features.intra_contour_energy = s_intra;
                obj.features.intra_contour_energy_norm = (s_intra ./sum(m));                
                %obj.features.ratio_intra_inter_contour_energy = s_intra./s;
            end
                                 
            if(isfield(obj.features, 'curv_continuity'))                     
                EXCLUDE_FRAME = true;
                obj.boundaries = obj.compute_boundaries(masks, EXCLUDE_FRAME);                
                SAMPLING_RATE = 15;
                INTERVAL = 10;
                mask_c = zeros(size(masks,3),1);
                parfor (i=1:size(masks,3),4)
                    n_conn = numel(obj.boundaries{i});
                    n_intervals = zeros(n_conn,1);
                    mask_c_cc = zeros(n_conn,1);
                    for j=1:n_conn          
                        this_bndr = obj.boundaries{i}{j};
                        bndr_length = size(this_bndr,1);

                        n_intervals(j) = numel(1:SAMPLING_RATE:bndr_length) - 1;
                        c = zeros(n_intervals(j), 1);
                        circle_center = zeros(2,n_intervals(j));
                        for k=1:n_intervals(j)           
                            range = round(max(1, (k*(SAMPLING_RATE) - (INTERVAL/2))):min(bndr_length,(k*(SAMPLING_RATE) + (INTERVAL/2))));
                            x = this_bndr(range,1);
                            y = this_bndr(range,2);
                            %[x,y] = ind2sub(sz, obj.adj_contour_conn_comp{i}.PixelIdxList{j}(range));

                            [out, curvature] = curvature_2d([x y]); 
                            circle_center(:,k) = out(:,[1 2]);
                            c(k) = curvature;                                                                       
                        end
                        mask_c_cc(j) = sum(c)/max(1,numel(c));
                    end     
                    mask_c(i) = sum(mask_c_cc)/max(1,numel(mask_c_cc));
                end
%                 %%% curvilinear continuity 
%                 %%%( curvature using bruckstein discrete approximation )
%                 the_lines = cell(size(m,2),1);
%                 for i=1:size(m,2)
%                     [ids1, ids2] = find(all_bw_len(:,:,i));
%                     ids = [ids1 ids2];
%                     
%                     %BREAK_FACTOR = 0.9;
%                     %thresh = BREAK_FACTOR*0.05*log(size(bw{i}{1},1))/log(1.1);
%                     thresh = 5;
%                     lines = lineseg({ids}, thresh);
%                     the_lines{i} = seglist2segarray(lines);        
%                     jl_plot_lines(the_lines{i});
%                     pause;
%                 end
% 
%                 sum_curvature = zeros(size(m,2),1);
%                 for i=1:size(m,2) 
%                     if(size(the_lines{i},2) == 1)
%                         sum_curvature(i) = 0;
%                         continue;
%                     end
% 
%                     curvatures = zeros(size(the_lines{i},2)-1,1);
%                     curvatures(1) = angle_between_linesegs(the_lines{i}(:,end), the_lines{i}(:,1));
%                     for j=1:size(the_lines{i},2)-1
%                         curvatures(j) = (angle_between_linesegs(the_lines{i}(:,j), the_lines{i}(:,j+1)))^2;
%                         %jl_plot_lines(the_lines{i}(:,j:j+1)); hold on;
%                         len_1 = norm(the_lines{i}(1:2,j) - the_lines{i}((3:4),j));
%                         len_2 = norm(the_lines{i}((1:2),j+1) - the_lines{i}(((3:4)),j+1));
%                         curvatures(j) = curvatures(j)/(min(len_1,len_2));
%                         %curvatures(j)
%                      end
%                     sum_curvature(i) = sum(curvatures)/numel(curvatures); 
%                 end                
                
                % 8. curvilinear continuity
                obj.features.curv_continuity = mask_c';
            end
            
            
            % visualization
            %
            % feats = fieldnames(obj.features);
            % for i=1:numel(feats)
            %   [val, id] = sort(obj.features.(feats{i})(1,:), 'descend');
            %   subplot_auto_transparent(obj.segments(:,:,id), obj.I, val);
            %   feats{i}             
            % end
        end                                   
        
        %%%%%%%%%%%%%%%% Features %%%%%%%%%%%%%%%%%%%%%%
    end
    methods (Static)
        function features = create_features_struct()
            features.inter_texton_dissim = []; % ~ok
            features.intra_texton_sim= []; % ~ok
            features.inter_brightness_dissim = []; % ~ok
            features.intra_brightness_sim = []; % ~ok
            features.inter_contour_energy = [];
            features.inter_contour_energy_norm = [];
            features.intra_contour_energy = []; 
            features.intra_contour_energy_norm = []; 
            features.curv_continuity = []; 
        end
        
        function print_feat_vals(vals)
            feats = GestaltSegmentFeatures.create_features_struct();
            feat_names = fieldnames(feats);
            assert(size(vals,1) == numel(feat_names));
            for i = 1:numel(feat_names)
                fprintf(['%s: %f\n'], feat_names{i}, vals(i));
            end
        end    
    end       
end
