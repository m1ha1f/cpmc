% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef Segmenter   
    properties
        I
        img_name
        params
        pb_path
        windows_path
        
        Pixels                
        p_hgraph    % pixel pairwise affinities
        
        L
        Superpixels
        sp_hgraph  % superpixel pairwise affinities
        Sp_centers % k-means color centers
        
        Background
        Foreground        
        
        P                  % GraphProb
                
        Segments      
        
        % early segment filtering
        MIN_NPIXELS
        MAX_ENERGY
        MIN_DISSIMILARITY
        DONT_BREAK_DISCONNECTED
        filter_segments
        randomize_N 
        
        % seed hypothesis stuff
        arrangement
        subframe_dims 
        frame_width
        RECT_DIMS
        grid_dims
        frame_grid_dims
        frame_types
        SEED_FRAME_WEIGHT
        img_frame
        windows
        window_gen_pars
        
        upper_bp
        
        resize_factor
        orig_I
        initialized
        
        gPb_thin
    end
    
    methods
        function obj = Segmenter(I, img_name)
            vl_twister('STATE',1313);
            
            obj.I = I;
            obj.img_name = img_name;
            
            obj.Foreground = 1;
            obj.Background = 2;            
            
            % segment filtering parameters
            obj.filter_segments = true;
            obj.MIN_NPIXELS = 100;
            %obj.MAX_ENERGY = 0.15;   % 200
            obj.MAX_ENERGY = 100; % 100
            %obj.MAX_ENERGY = inf;            
            obj.MIN_DISSIMILARITY = 0.05;
            obj.DONT_BREAK_DISCONNECTED = true;              
            obj.randomize_N = 1000;
            
            % seed hypothesis parameters
            obj.RECT_DIMS = [20 20]; % original [5 5]
            obj.frame_width = 1;
            obj.arrangement = 'grid';
            %obj.arrangement = 'grid';
            obj.grid_dims = 1;
            obj.frame_types = {'all'};  %'all_but_down', 'horiz', 'vert'};            
            obj.SEED_FRAME_WEIGHT = 400; % 400
            
            %obj.params.local_dist_type = 'pb';
            
            obj.upper_bp = 150; % 15 before
            
            obj.resize_factor = 1;
            obj.initialized = false;
            
            obj.params = JustSegmGen.create_pars();
            obj.params.CONTRAST_SENSITIVE_SIGMA = []; % others will check this
            
        end

        function obj = initialize(obj)  
            obj = obj.set_params();
            obj.orig_I = obj.I;
            if(obj.resize_factor ~= 1)
                obj.orig_I = obj.I;
                obj.I = imresize(obj.I, obj.resize_factor);
                obj.RECT_DIMS = ceil(obj.RECT_DIMS * obj.resize_factor);
            end            
            
            obj.P = GraphProb(obj.I);
            obj.P.upper_bp = obj.upper_bp;
            
            p_gen = JustSegmGen(obj.params);            

            if(strcmp(obj.params.local_dist_type, 'pb'))
                % compute pb on original img, get graph on resized version
                if(~exist([obj.pb_path obj.img_name  '_PB.mat'], 'file'))
                    t_pb = tic();
                    [gPb_orient, gPb_thin, textons] = globalPb_new_nofiles(obj.orig_I, [], 0.5);
                    save([obj.pb_path obj.img_name  '_PB.mat'], 'gPb_thin', 'gPb_orient', 'textons');
                    toc(t_pb)
                else
                    %load([obj.pb_path obj.img_name  '_PB.mat'], 'gPb_thin', 'gPb_orient', 'textons');
                    %obj.gPb_thin = imresize(gPb_thin, obj.resize_factor);                    
                end
                
                obj.p_hgraph = p_gen.apply(obj.I, [obj.pb_path obj.img_name  '_PB.mat']);
            else
                obj.p_hgraph = p_gen.apply(obj.I);
            end
            
%             obj.p_hgraph.values{1} = round(obj.p_hgraph.values{1}*obj.params.CONSTANT);
            obj.Pixels = obj.p_hgraph.nodes;
            
            obj.img_frame = obj.P.generate_img_frame('pixels', 'all', obj.frame_width);            
            
            % compute superpixels
%             [obj.L] = vgg_segment_gb(obj.I, 0.005, 50, 100, true);
%             %sc(L,'rand')
%             un = unique(obj.L);
%             obj.Superpixels = cell(numel(un), 1);            
%             pix = reshape(obj.I, size(obj.I,1)* size(obj.I,2), 3);    
%             K = 5;
%             for i=1:numel(un)
%                 obj.Superpixels{i} = find(obj.L == un(i));                
%                 obj.Sp_centers{i} = vl_ikmeans(pix(obj.Superpixels{i},:)', min(K, numel(obj.Superpixels{i})), 'Method', 'lloyd'); 
%             end
%             
%              extractor = SuperPixelHgraphExtractor(true, false);
%              img_centers = extractor.get_superpixel_centers(obj.L);
%              feats = extractor.get_superpixel_color_means(obj.I, obj.L);            
%              obj.sp_hgraph = exp(-0.0005*pdist2(double(feats)', double(feats)'));
             %sc(obj.sp_hgraph)
                             
                           
            %if(strcmp(class(obj), 'LongRangeSegmenter'))
            %    scales = [0.4 0.8] * obj.resize_factor;
            %    disp('computing huehistogram on patches!!!!');
            %    [obj.F,obj.D] = csift(obj.I, 3, scales, 'rgbhistogram');            
            %    obj.D = uint8(obj.D*255.0);
            %end
            obj.initialized = true;
        end
        
        function obj = compute_windows(obj)
            %%% compute putative windows (sort of regions of interest) %%%
            kind = 'sliding_window_detection'; % default type
            if(isfield(obj.window_gen_pars, 'kind'))
                kind = obj.window_gen_pars.kind;
            end
            
            if(~isempty(obj.windows_path))
                wingen = WindowGen(obj.orig_I, kind, size(obj.orig_I), [obj.windows_path '/' obj.img_name '_W.mat']);
            else
                wingen = WindowGen(obj.orig_I, kind, size(obj.orig_I), [obj.img_name '_W.mat']);
            end
            obj.windows = wingen.apply(obj.window_gen_pars.det_classes);            
            rsz_factor = size(obj.I, 1)/size(obj.orig_I,1);
            obj.windows = obj.windows*rsz_factor;
        end
        
        function obj = set_params(obj, par)
            if(nargin==1)                
                
                % params
                obj.params.DEBUG_MODE = false;
                obj.params.CONSTANT = 1000;
                
                if(isempty(obj.params.local_dist_type))
                    obj.params.local_dist_type = 'pb';
                end
                
                if(~isempty(obj.params.local_dist_type) && strcmp(obj.params.local_dist_type, 'color'))
                    obj.params.CONTRAST_SENSITIVE_SIGMA = 0.05; % 0.1
                else
                    obj.params.CONTRAST_SENSITIVE_SIGMA = 0.1; % 0.1
                end                
                
                obj.params.CONTRAST_SENSITIVE_WEIGHT =1;
                obj.params.PAIRWISE_OFFSET = 2/obj.params.CONSTANT;
                obj.params.n_neighbors =  4;
                
            else
                obj.params = par;
            end
        end        
        
         function obj = set_pb_path(obj, pb_path)
             obj.pb_path = pb_path;
         end         
        
        function obj = set_windows_path(obj, windows_path, window_gen_pars)
            obj.windows_path = windows_path;           
            if(exist('window_gen_pars', 'var'))
                obj.window_gen_pars = window_gen_pars;
            else
                obj.window_gen_pars = [];
            end
        end
        
        function obj = compute_segments(obj, max_n_segms)
            if(~obj.initialized)
                error('Run initialize() method before this');
            end
            
            if(~exist('max_n_segms', 'var'))
                max_n_segms = inf;
            end
            
            [obj.P] = obj.P.add_atomic_vars('pixels', obj.p_hgraph);          
            [obj.P] = obj.P.add_atomic_vars('classes', Hgraph(2));
            
            obj = obj.add_hyp_conns(); % this one gets specialized by subclasses
            
            %obj.display_subframes();
            if(~isempty(obj.P.hypConn))
                % obj.P = obj.P.solve('classes', obj.Foreground, obj.Background);
                t_runFlowTotal = tic();
                obj.P = obj.P.solveInstances('classes', obj.Foreground, obj.Background, obj.p_hgraph.leftTranspose, obj.p_hgraph.rightTranspose, obj.p_hgraph.top, obj.p_hgraph.bottom);
                t_runFlowTotal = toc(t_runFlowTotal)
            else
                return;
            end
            
            %obj.P.show_results();
            min_n_pixels = obj.resize_factor*obj.MIN_NPIXELS; % MIN_NPIXELS should be given in the size of the original image.
            S = SegmentProcessor(obj.P.prob_dgraph, obj.P.get_type_solution('pixels'), obj.P.I, [obj.pb_path obj.img_name '_PB'], min_n_pixels, obj.DONT_BREAK_DISCONNECTED);

            t_en = tic();
            S = S.compute_energies();
            t_energy = toc(t_en)
            
            min_n_segments = 5;
            
            %[val, id] = sort([S.energies.(S.deciding_energy_type)]);            
            %subplot_auto_transparent(S.segments(:,id(1:100)), obj.I, val(1:100));
            %subplot_auto_transparent(S.segments(:,id(end-100:end)), obj.I, val(end-100:end));
            %subplot_auto_transparent(S.segments(:,id(1:end)), obj.I, val(1:end))
            
            n_before_energy = size(S.segments,2);
            S = S.filter_segments('max_energy', obj.MAX_ENERGY, min_n_segments);
            %n_after_energy = size(S.segments,2);
            %fprintf('N segments after energy filtering: %d\n', n_after_energy);
            
            %subplot_auto_transparent(S.segments(:, [S.energies.cut_ratio] >= 100), obj.I)

            if(obj.filter_segments)   
                if(size(S.segments,2)> obj.randomize_N)
                    randn = randperm(size(S.segments,2));                    
                    S.segments = S.segments(:,randn(1:obj.randomize_N));
                    S.energies = S.energies(randn(1:obj.randomize_N));
                end
                    
                t = tic();
                %n_before = size(S.segments,2)
                %S = S.filter_segments('min_dissimilarity', 0.05);
                %n_after = size(S.segments,2)

                S = S.remove_repeated_segments();
                n_before_clustering = size(S.segments,2);
               
                if(n_before_clustering ~= 1)
                    S = S.filter_segments('min_dissimilarity', 0.05, min_n_segments, max_n_segms);
                end
                n_after_clustering = size(S.segments,2);
                %fprintf('N segments after clustering: %d\n', n_after_clustering);

                
                time_removing_similar = toc(t);
                obj.Segments = S.segments;
            else
                obj.Segments = S.segments; % obj.P.get_type_solution('pixels');
            end
            
            if(obj.resize_factor ~=1)
                orig_dims = [size(obj.orig_I,1) size(obj.orig_I,2)];
                n = size(obj.Segments,2);
                
                segm = reshape(obj.Segments, size(obj.I,1), size(obj.I,2), n);
                segm = imresize(segm, orig_dims, 'nearest'); % faster
                new_segms = reshape(segm, prod(orig_dims), n);
                
                obj.Segments = new_segms;
            end
        end        
        
        
        function obj = compute_segments_interactive(obj)            
            [obj.P] = obj.P.add_atomic_vars('pixels', obj.p_hgraph);          
            [obj.P] = obj.P.add_atomic_vars('classes', Hgraph(2));                                   
            [hyp_conns, types] = obj.generate_interactive_hyp(); % this one gets specialized by subclasses                     
            
            t0 = tic();
            obj.P = obj.P.add_hyp_connections(types, hyp_conns);
            obj.P = obj.P.solve('classes', obj.Foreground, obj.Background);
            time_solve = toc(t0)
            
            if(obj.filter_segments)
                min_n_segments = 5;
                %S = SegmentProcessor(obj.P.prob_dgraph, obj.P.get_type_solution('pixels'), obj.P.I, obj.MIN_NPIXELS, false);
                S = SegmentProcessor(obj.P.prob_dgraph, obj.P.get_type_solution('pixels'), obj.P.I, [obj.pb_path obj.img_name '_PB'], obj.MIN_NPIXELS, obj.DONT_BREAK_DISCONNECTED);
                S = S.compute_energies();

                S = S.filter_segments('max_energy', obj.MAX_ENERGY, min_n_segments); 
                S = S.filter_segments('min_dissimilarity', obj.MIN_DISSIMILARITY);

                obj.Segments = S.segments;            
            else
                obj.Segments = obj.P.get_type_solution('pixels');
            end
             
            if(obj.resize_factor ~=1)
                orig_dims = [size(obj.orig_I,1) size(obj.orig_I,2)];
                n = size(obj.Segments,2);
                
                segm = reshape(obj.Segments, size(obj.I,1), size(obj.I,2), n);
                segm = imresize(segm, orig_dims, 'nearest'); % faster
                new_segms = reshape(segm, prod(orig_dims), n);
                
                obj.Segments = new_segms;
            end
        end
        % unary potentials
        
        function val = compute_unary_values(obj, seed_region)             
            if(size(obj.I, 3) ~=3)
                obj.I(:,:,2) = obj.I(:,:,1);
                obj.I(:,:,3) = obj.I(:,:,1);
            end            

            %val = obj.compute_unary_values_patch_rgb(seed_region);
            val = obj.compute_unary_values_rgb(seed_region);          
        end
          
        function val = compute_unary_values_rgb(obj, seed_region)             
             K = 5; % 10
             
             %%%%%% RGB space %%%%%%%%%%%
             pix = reshape(obj.I, size(obj.I,1)* size(obj.I,2), 3);           
             
             %%%%%%% HSV space %%%%%%%%%%%
             %pix = uint8(255*reshape(rgb2hsv(obj.I), size(obj.I,1)* size(obj.I,2), 3));  
             
             %%%%%%% LAB space %%%%%%%%%%%%
%               pix = reshape(rgb2lab(obj.I), size(obj.I,1)* size(obj.I,2), 3);          
%               pix(:,1) = pix(:,1)*2.55;
%               pix(:,[2 3]) = 255*((pix(:,[2 3]) + 110) / 220);
%               pix = uint8(pix);
             %t = tic();
                                       
             if(K>size(seed_region,1))
                 seed_region = [seed_region; seed_region];
             end
             
             centers = vl_ikmeans(pix(seed_region,:)', K, 'Method', 'lloyd');  % centers in columns!             
             %toc(t);
             %obj.show_centers_colors(centers);
             
             %center = mean(double(pix_rgb(seed_region,:)));
             %dists = sum(abs(pix_rgb - center(ones(size(obj.I,1)* size(obj.I,2),1),:)),2);
             
             %ext_pix_rgb = pix(:,:,ones(K,1));
             
             dists = inf(size(pix,1),1);
             pix = single(pix);
             
             grower = ones(size(obj.I,1)* size(obj.I,2),1);
             for i=1:K
                 center = single(centers(:,i)');
                 %t = tic();
                 dists = min(dists, sum(abs(pix - center(grower,:)),2));
                 %toc(t)
             end
             %0.12 works best 600x400 images, 0.3 with 300x200, se
             val = double(0.3* obj.params.CONSTANT  * exp(-dists*0.07));             
             %val = double(120*exp(-dists*0.07));             
             %val = double(120*exp(-dists*0.07));             
             %sc(reshape(uint8(val),size(obj.I,1), size(obj.I,2)))             
             %pause;
             
             %val = obj.smooth_spherically(seed_region, max([size(obj.I,1) size(obj.I,2)])/10, val);
        end
        
        function val = compute_unary_values_patch_rgb(obj, seed_region)    
            error('not working');
            K = 5; % 10
            
            %%%%%% RGB space %%%%%%%%%%%
            pix = reshape(obj.I, size(obj.I,1)* size(obj.I,2), 3);
            
            f_lin = sub2ind(size(obj.I), round(obj.F(2,:)), round(obj.F(1,:)));
            [new_seed_region, seed_ids] = intersect(f_lin, seed_region);
            
            if(K>numel(new_seed_region))
                rp = randperm(numel(f_lin));
                seed_ids = rp(1:K); 
            end
                      
            [centers, ass] = vl_ikmeans(obj.D(:, seed_ids), min(K, numel(seed_ids)) , 'Method', 'lloyd');  % centers in columns!
            unique_center_ids = unique(ass);
            centers = centers(:, unique_center_ids);
            %toc(t);
            %show_centers_colors(centers);
            
            %center = mean(double(pix_rgb(seed_region,:)));
            %dists = sum(abs(pix_rgb - center(ones(size(obj.I,1)* size(obj.I,2),1),:)),2);
            
            %ext_pix_rgb = pix(:,:,ones(K,1));
            
            dists = inf(numel(f_lin),1);
            pix = single(obj.D)';
            
            grower = ones(numel(f_lin),1);
            for i=1:numel(unique_center_ids)
                center = single(centers(:,i)');
                %t = tic();
                dists = min(dists, sum(abs(pix - center(grower,:)),2));
                %toc(t)
            end
            %0.12 works best 600x400 images, 0.3 with 300x200, se
            val = double(obj.params.CONSTANT  * exp(-dists*0.009));
            
            % interpolate as stella yu did in imncut %            
            
%             new_val = zeros(size(obj.I,1), size(obj.I,2));
%             new_val(f_lin) = val;
%             sc(reshape(new_val,size(obj.I,1), size(obj.I,2)))
%              
            coords = round(obj.F(1:2,:));
            [x, y] = meshgrid(unique(coords(2,:)), unique(coords(1,:)));
            [xi, yi] = meshgrid(1:size(obj.I,1), 1:size(obj.I,2));

            if(size(x,1)*size(x,2) ~= numel(val))
                val = [val; 0];
            end
            
            ZI = INTERP2(x, y, reshape(val, size(x,1), size(x,2)), xi, yi)';
            subplot(1,2,1), sc(obj.I);
            subplot(1,2,2), sc(ZI);
            pause;
            val = reshape(ZI, size(obj.I,1) * size(obj.I,2), 1);
        end
            
        function [fg_val, bg_val] = compute_unary_values_iterative(obj, bbox)
            if(size(obj.I, 3) ~=3)
                obj.I(:,:,2) = obj.I(:,:,1);
                obj.I(:,:,3) = obj.I(:,:,1);
            end
            
            [fg_val, bg_val] = obj.compute_unary_values_rgb_bbox_iterative( bbox);            
        end
        
        function [fg_val, bg_val] = compute_unary_values_edges(obj, bbox) 
            error('just testing');
             Icrop = obj.P.crop_img_to_frame(obj.I, bbox);                
             
             ed = edge(histeq(rgb2gray(Icrop)), 0.5);
             fg_val = reshape(ed, size(ed,1)* size(ed,2), 1);                          
             fg_val = single(fg_val);
             
             b = bwdist(ed);
             bg_val = double(0.03* obj.params.CONSTANT  *(1-fg_val));
             fg_val =  double(1* obj.params.CONSTANT * fg_val);
        end
        
        function [fg_val, bg_val] = compute_unary_values_rgb_bbox_iterative(obj, bbox)
             K_bg = 20; % 10
             K_fg = 20;
             
             Icrop = obj.P.crop_img_to_frame(obj.I, bbox);
             pix = reshape(Icrop, size(Icrop,1)* size(Icrop,2), 3);                          
             spix = single(pix);
             
             zero_one_img = ones(size(Icrop,1), size(Icrop,2));             
             
             %if(numel(pix) == (size(obj.I,1)*size(obj.I,2)))
             subframe_perimeter = logical(bwperim(zero_one_img));
             %else             
             %   subframe_perimeter = setdiff();
             %end
             
             
             [centers_bg, ass] = vl_ikmeans(pix(subframe_perimeter,:)', K_bg, 'Method', 'lloyd');  % centers in columns!              
             %hist_bg = hist(double(ass), size(centers_bg,2));
             %hist_bg = (double(hist_bg)/ sum(hist_bg))';             
             [dists_bg, ass] = obj.dist_to_centers(spix, centers_bg);    
             %bg_val = hist_bg(ass);
             
             [dists_bg_sorted, ids_sort] = sort(dists_bg, 'descend');             
             min_n = K_fg;
             pass_test = (dists_bg_sorted > 100);
             fg_ids = ids_sort(1:max(sum(pass_test),min_n));              
             
             %%% version 1 %%%%
             %%% fg potentials %%%      
%              best_collision = inf;
%              for i=1:5               
%                  [centers_fg, ass] = vl_ikmeans(pix(fg_ids,:)', min(numel(fg_ids),K_fg), 'Method', 'lloyd');  % centers in columns!
%                  dists_fg = obj.dist_to_centers(spix, centers_fg);             
%                  d = sqrt(pdist2(double(centers_fg'), double(centers_bg')));
%                  d_min = min(d,[],2);
%                  [val, id] = sort(d_min, 'ascend');               
%                                                                                         
%                  collision = sum(val);
%                  if(collision < best_collision)
%                     best_collision = collision;
%                     the_fg_dist = dists_fg;
%                  end
%                  
%                  c = obj.assign_to_codebooks(spix, centers_bg, centers_fg);
%                  fg_ids = fg_ids(c);
%              end
%             the_bg_dist = dists_bg;

             %%% version 2 (the best) %%%%              
              [centers_fg, ass] = vl_ikmeans(pix(fg_ids,:)', min(numel(fg_ids),K_fg), 'Method', 'lloyd');  
              %hist_fg = hist(double(ass), size(centers_fg,2));
              %hist_fg = (double(hist_fg)/ sum(hist_fg))';
              [fg_dist, ass] = obj.dist_to_centers(spix, centers_fg); 
              %fg_val = hist_fg(ass);
              
              the_bg_dist = dists_bg;
              the_fg_dist = fg_dist;
         
%%%% version 3 -- with superpixels %%%%
%%% find superpixels fully contained inside subframe
%              inside_frame= [];
%              Lcrop = obj.P.crop_img_to_frame(obj.L, bbox);
%              unique_L = unique(obj.L);
%              id_translator = zeros(numel(unique_L),1);
%              un = unique(Lcrop);
%              new_ids = 1:numel(un);
%              id_translator(un) = new_ids;
%              
%              new_Sp = cell(numel(un),1);
%              subframe_perimeter_ids = find(subframe_perimeter);
%              for i=1:numel(un)
%                  new_Sp{i} = find(Lcrop == un(i));
%                  if(isempty(intersect(new_Sp{i}, subframe_perimeter_ids)))
%                      inside_frame = [inside_frame; un(i)];
%                  end
%              end
%              
%              outside_frame = setdiff(1:numel(obj.Superpixels), inside_frame);
%              
%              frame_to_others_sim = max(obj.sp_hgraph(outside_frame, inside_frame));
%              [val, sorted_ftos] = sort(frame_to_others_sim,'ascend');
%              ids = sorted_ftos(1:min(numel(sorted_ftos), 5));
%              
%              the_ones = inside_frame(ids);
%              
%              %id_translator
%              if(~isempty(the_ones))
%                  fg_ids = new_Sp{id_translator(the_ones)};
%              end
%              [centers_fg] = vl_ikmeans(pix(fg_ids,:)', min(numel(fg_ids),K_fg), 'Method', 'lloyd');
%              the_fg_dist = obj.dist_to_centers(spix, centers_fg);
%              the_bg_dist = dists_bg;
                                                       
                          
             % ending, shared by all version %             
             bg_val = double(0.2* obj.params.CONSTANT  *exp(-the_bg_dist*0.01));            
             fg_val =  double(0.2* obj.params.CONSTANT  *exp(-the_fg_dist*0.01));
             
             %bg_val = histeq(double(reshape(bg_val,size(Icrop,1), size(Icrop,2))));
             %fg_val = histeq(reshape(fg_val,size(Icrop,1), size(Icrop,2)));
             
             %bg_val = double(0.1* obj.params.CONSTANT  * bg_val);
             %fg_val =  double(0.1* obj.params.CONSTANT  * fg_val);     
             
                                       
             %[sum_sorted, ids_sum_sort] = sort(fg_val+bg_val, 'ascend');             
%              figure;
%              subplot(1,3,1), obj.show_centers_colors(centers_bg);
%              subplot(1,3,2), obj.show_centers_colors(centers_fg);
%              [c] = obj.assign_to_codebooks(spix, centers_bg, centers_fg);
%              subplot(1,3,3), imshow(reshape(c*0.5, size(Icrop,1), size(Icrop,2)))
%              
%              figure;             
%              subplot(1,3,1), sc(Icrop);
%              subplot(1,3,2), obj.show_centers_colors(centers_bg); title('boundary');
%              subplot(1,3,3), obj.show_centers_colors(centers_fg); title('rest');
             
             %figure;             
             % subplot(1,4,1), sc(Icrop);
             % subplot(1,4,2), sc(reshape(double(fg_val),size(Icrop,1), size(Icrop,2))); title('Foreground');
             % subplot(1,4,3), sc(reshape(double(bg_val),size(Icrop,1), size(Icrop,2))); title('Background');
             % subplot(1,4,4), sc(reshape(double(fg_val-bg_val),size(Icrop,1), size(Icrop,2))); title('Sum of likelihood');

%              figure;
%              to_display = zeros(size(Lcrop));
%              for i=1:numel(the_ones)
%                  to_display(Lcrop == the_ones(i)) = i;
%              end
%              subplot(1,2,1), sc(to_display);
%              subplot(1,2,2), sc(Lcrop, 'rand');
%              pause;
        end
        
        function [c, dists_1, dists_2] = assign_to_codebooks(obj, spix, c1, c2)
            dists_1 = obj.dist_to_centers(spix, c1);
            dists_2 = obj.dist_to_centers(spix, c2);         
            c = (dists_1 > dists_2) + 1;
        end
        
        function [dists, ass] = dist_to_centers(obj, pix, centers)
            grower = ones(size(pix,1),1);
            dists = inf(size(pix,1),1);
            ass = zeros(size(pix,1),1);
           
            for i=1:size(centers,2)
                center = single(centers(:,i)');
                [new_dists] = min(dists, sum(abs(pix - center(grower,:)),2));
                changed_dists = (new_dists ~=  dists);
                ass(changed_dists) = i;
                dists = new_dists;
            end
        end
        
        function show_centers_colors(obj, centers)
            for i=1:size(centers,2)
                theI(i,1,1:3) = uint8(centers(1:3,i));
            end
            sc(theI);
        end
        
        
        function val = compute_unary_values_gray(obj, seed_region)
             function show_centers_colors(centers)
                 for i=1:size(centers,2)
                     theI(i,1,1:3) = uint8(centers(1:3,i));                     
                 end          
                 imshow(theI);
             end
             
             K = 10;
             pix_gray = reshape(obj.I, size(obj.I,1)* size(obj.I,2), 1);
             %t = tic();
             centers = vl_ikmeans(pix_gray(seed_region)', K, 'Method', 'lloyd');  % centers in columns!             
             %toc(t);
             %show_centers_colors(centers);
             
             %center = mean(double(pix_rgb(seed_region,:)));
             %dists = sum(abs(pix_rgb - center(ones(size(obj.I,1)* size(obj.I,2),1),:)),2);
             
             %ext_pix_rgb = pix_rgb(:,:,ones(K,1));
             
             dists = inf(size(pix_gray,1),1);
             pix_gray = single(pix_gray);
             for i=1:K
                 center = single(centers(:,i)');
                 dists = min(dists, sum(abs(pix_gray - center(ones(size(obj.I,1)* size(obj.I,2),1),:)),2));
             end
                          
             dists = dists*3;
             val = double(80*exp(-dists*0.07));
             %hist(val)
             %imshow(reshape(uint8(val),size(obj.I,1), size(obj.I,2)))
        end        
        
         %%%%%%%%%%%%%
         % Visualization  %
         %%%%%%%%%%%%%
         
         function display_grids(obj)
             our_rect = obj.P.generate_rectangle_coords(obj.RECT_DIMS);
             grid = obj.P.generate_img_grid('pixels', our_rect, obj.external_grid_dims);

             theI = obj.I;
             if(size(theI,3) == 1)
                theI(:,:,2) = theI(:,:,1);
                theI(:,:,3) = theI(:,:,1);
             end
             
             newI = reshape(theI, size(theI,1)*size(theI,2), 3);
             for i=1:length(grid)
                 hyp = grid{i};
                 the_color = rand(1,3);
                 newI(hyp, :) = 255 * the_color(ones(length(hyp),1), :);
             end
             newI = reshape(newI, size(theI,1), size(theI,2), size(theI,3));
             imshow(newI);
             title('external_growth_grid');
         end
         
         function display_subframes(obj)
             [x,y] = cartprod_mex((-obj.subframe_dims(1)/2:obj.subframe_dims(1)/2)', (-obj.subframe_dims(2)/2:obj.subframe_dims(2)/2)');
            our_rect = [x y];
            fg_grid_hyp = obj.P.generate_img_grid('pixels', our_rect, obj.frame_grid_dims);
            
            imshow(obj.I);
            for i=1:numel(fg_grid_hyp)
                [x,y]  = ind2sub([size(obj.I,1) size(obj.I,2)], fg_grid_hyp{i});
                min_x = min(x);
                min_y = min(y);
                h = max(x) - min_x;
                w = max(y) - min_y;
               rectangle('Position', [min_y min_x w h]);
               pause;
            end          
         end
         
         function display_segments(obj, order)
             if(nargin == 2 && order)
                [val, id] = sort(sum(obj.Segments));
                obj.Segments = obj.Segments(:,id);
             end
             
             if(~isempty(obj.Segments))
                subplot_auto_transparent(obj.Segments, obj.I);
             end
         end
                  
        %%%%%%%%%%%%%%%%%%%%
        % hypothesis generation %
        %%%%%%%%%%%%%%%%%%%%                

        
        %%%%% for superpixels (very old code, no longer works) %%%%        
        function [hyp_conns, types] = generate_external_growth_hyps_sp(obj)
            %our_rect = obj.P.generate_rectangle_coords(obj.RECT_DIMS);
            %fg_grid_hyp = obj.P.generate_img_grid('superpixels', our_rect, obj.external_grid_dims);
            region = obj.P.get_vars_interactive('superpixels');
            fg_grid_hyp = {region};
            
            for i=1:length(fg_grid_hyp)
                hyp_conns{i,1} = {fg_grid_hyp{i}, obj.Foreground, inf, 0};
                other_sp = setdiff(obj.Superpixels, [fg_grid_hyp{i}]);
                hyp_conns{i,2} = {other_sp, obj.Background, 0, 1};
            end
            types{1,2} = {'superpixels', 'classes'};
            types{1,1} = {'superpixels', 'classes'};
        end
        
        function [hyp_conns, types] = generate_internal_growth_hyps_sp(obj)
            our_rect = obj.P.generate_rectangle_coords(obj.RECT_DIMS);
            region = obj.P.get_vars_interactive('superpixels');
            fg_grid_hyp = {region};
            bg_frame = obj.P.generate_img_frame('superpixels', 'all');
            
            for i=1:length(fg_grid_hyp)
                hyp_conns{i,1} = {bg_frame, obj.Background, obj.SEED_FRAME_WEIGHT, 0};
                hyp_conns{i,2} = {fg_grid_hyp{i}, obj.Foreground, inf, 0};
                hyp_conns{i,3} = {setdiff(obj.Superpixels, [fg_grid_hyp{i} bg_frame]), obj.Foreground, 0, 1};
            end
            types{1,1} = {'superpixels', 'classes'};
            types{1,2} = {'superpixels', 'classes'};
            types{1,3} = {'superpixels', 'classes'};
        end
        
        function [hyp_conns, types] = generate_internal_growth_hyps_sp_no_internal(obj)
            bg_frame = obj.P.generate_img_frame('superpixels', 'all');
            
            hyp_conns{1,1} = {bg_frame, obj.Background, obj.SEED_FRAME_WEIGHT, 0};
            hyp_conns{1,2} = {setdiff(obj.Superpixels, bg_frame), obj.Foreground, 0, 1};
            types{1,1} = {'superpixels', 'classes'};
            types{1,2} = {'superpixels', 'classes'};
        end
        
        function [hyp_conns, types] = generate_external_growth_hyps_mixed(obj)
            our_rect = obj.P.generate_rectangle_coords(obj.RECT_DIMS);
            fg_grid_hyp_sp = obj.P.generate_img_grid('superpixels', our_rect, obj.external_grid_dims);
            fg_grid_hyp_p = obj.P.generate_img_grid('pixels', our_rect, obj.external_grid_dims);
            
            mapAggs = obj.P.mapAggs('superpixels');
            for i=1:length(fg_grid_hyp_sp)
                hyp_conns{i,1} = {fg_grid_hyp_sp{i},obj.Foreground, inf, 0};
                pixels = cell2mat(mapAggs(fg_grid_hyp_sp{i}));
                hyp_conns{i,2} = {pixels, obj.Foreground, inf, 0};
                
                other_sp = setdiff(obj.Superpixels, [fg_grid_hyp_sp{i}]);
                other_pixels = cell2mat(mapAggs(other_sp));
                hyp_conns{i,3} = {other_pixels, obj.Background, 0, 1};
            end
            %types{1,1} = {'superpixels', 'classes'};
            types{1,1} = {'superpixels', 'classes'};
            types{1,2} = {'pixels', 'classes'};
            types{1,3} = {'pixels', 'classes'};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%% pixels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % uniform unary potentials
        
        function [hyp_conns, types] = generate_growth_hyps(obj, internal_external)
            types = {};
            
            frame_pixels = obj.img_frame;
            inside_frame = setdiff(obj.Pixels, frame_pixels);
            
            our_rect = obj.P.generate_rectangle_coords(obj.RECT_DIMS);
            [fg_grid_hyp] = obj.P.generate_seeds('pixels', our_rect, obj.arrangement, obj.grid_dims);
            
            hyp_conns = [];
            for i=1:numel(fg_grid_hyp)
                if(strcmp(internal_external, 'internal'))
                    [hyp_conn, types] = internal(obj, fg_grid_hyp{i});
                    hyp_conns = [hyp_conns; hyp_conn];
                else
                    [hyp_conn, types] = subframe_plus_internal(obj, fg_grid_hyp{i}, inside_frame);
                    hyp_conns = [hyp_conns; hyp_conn];
                end
            end
        end 
                
        % non-uniform unary potentials                     
        
        function [hyp_conns, types] = generate_growth_hyps_unary(obj, internal_external)            
            types = {};
            our_rect = obj.P.generate_rectangle_coords(obj.RECT_DIMS);
            [fg_grid_hyp] = obj.P.generate_seeds('pixels', our_rect, obj.arrangement, obj.grid_dims);            
           
            frame_pixels = obj.img_frame; %obj.P.generate_img_frame('pixels', 'all', obj.frame_width);
            inside_frame = setdiff(obj.Pixels, frame_pixels);
            hyp_conns = [];
            for i=1:length(fg_grid_hyp)        
                if(strcmp(internal_external, 'internal'))
                    [hyp_conn, types] = obj.internal_unary(fg_grid_hyp{i});
                else
                    [hyp_conn, types] = obj.subframe_plus_internal_unary(fg_grid_hyp{i}, inside_frame);
                end
                hyp_conns = [hyp_conns; hyp_conn];
            end            
        end        
        
        %
        % with subframes     
        %        
                
        % uniform, without fg seed
        function [hyp_conns, types] = generate_subframes_growth_hyps(obj)
            [x,y] = cartprod_mex((-obj.subframe_dims(1)/2:obj.subframe_dims(1)/2)', (-obj.subframe_dims(2)/2:obj.subframe_dims(2)/2)');
            our_rect = [x y];
            
            if(isempty(obj.frame_grid_dims)) 
                obj.frame_grid_dims = [1 1]; % means it's a frame seed, our favorite
                counter = 1;
                for frame_kind= {'all', 'all_but_down', 'horiz', 'vert', 'down', 'up'};                    
                    frame = obj.P.generate_img_frame( 'pixels', frame_kind, 1);
                    fg_grid_hyp{counter} = setdiff(1:size(obj.I,1)*size(obj.I,2),frame);
                    counter = counter + 1;
                end
            else            
                [fg_grid_hyp] = obj.P.generate_seeds('pixels', our_rect, obj.arrangement,obj.frame_grid_dims, [],  obj.windows);
            end
            
            hyp_conns = [];
            for i=1:numel(fg_grid_hyp)
                [hyp_conn, types] = obj.subframe(fg_grid_hyp{i});
                hyp_conns = [hyp_conns; hyp_conn];
            end
        end
        
        % uniform, with fg seed
        function [hyp_conns, types] = generate_subframe_plus_internal_growth_hyps(obj)
            hyp_conns = {};
            types = {};
                        
            [x,y] = cartprod_mex((-obj.subframe_dims(1)/2:obj.subframe_dims(1)/2)', (-obj.subframe_dims(2)/2:obj.subframe_dims(2)/2)');
            our_frame_rect = [x y];
            
            [x,y] = cartprod_mex((-obj.RECT_DIMS(1)/2:obj.RECT_DIMS(1)/2)', (-obj.RECT_DIMS(2)/2:obj.RECT_DIMS(2)/2)');
            our_rect = [x y];
            
            % frames
            [frame_grid_hyp] = obj.P.generate_seeds('pixels', our_frame_rect, obj.arrangement, obj.frame_grid_dims, [],  obj.windows);
           
            hyp_conns = [];
            types = [];
            %parfor i=1:numel(frame_grid_hyp)             
            for i=1:numel(frame_grid_hyp)             
                duh = obj.P.generate_seeds('pixels', our_rect, 'grid', obj.grid_dims, frame_grid_hyp{i});
                internal_grid_hyp = duh;                
                for j=1:numel(internal_grid_hyp)
                    internal_grid_hyp{j} = intersect(internal_grid_hyp{j}, frame_grid_hyp{i}); % without this can crash, when there's something in the border                    
                    [hyp_conn, type] = obj.subframe_plus_internal(internal_grid_hyp{j}, frame_grid_hyp{i});                    
                    hyp_conns = [hyp_conns; hyp_conn];    
                    types = [types; type];
                end
            end
            types = types(1,:);
        end
        
        % unary, without fg seed
        function [hyp_conns, types] = generate_subframes_growth_hyps_unary(obj, iterative)
            hyp_conns = {};
            types = {};
                    
            % unary potentials connecting only to the frame                    
            [x,y] = cartprod_mex((-obj.subframe_dims(1)/2:obj.subframe_dims(1)/2)', (-obj.subframe_dims(2)/2:obj.subframe_dims(2)/2)');
            our_frame_rect = [x y];
            
            % frames
            [frame_grid_hyp] = obj.P.generate_seeds('pixels', our_frame_rect, obj.arrangement, obj.frame_grid_dims, [],  obj.windows);
            
            hyp_conns = []; 
            types = [];
            
            %parfor
            for i=1:numel(frame_grid_hyp)             
                if(iterative)
                    [hyp_conn, type] = obj.subframe_unary_iterative(frame_grid_hyp{i});
                else
                    [hyp_conn, type] = obj.subframe_unary(frame_grid_hyp{i});
                end
                types = [types; type];
                hyp_conns = [hyp_conns; hyp_conn];
            end              
        end        
        
        % unary, with fg seed
        function [hyp_conns, types] = generate_subframe_plus_internal_growth_hyps_unary(obj)
            hyp_conns = {};
            types = {};
            
            % unary potentials connecting only to the frame
            [x,y] = cartprod_mex((-obj.subframe_dims(1)/2:obj.subframe_dims(1)/2)', (-obj.subframe_dims(2)/2:obj.subframe_dims(2)/2)');
            our_frame_rect = [x y];
            
            [x,y] = cartprod_mex((-obj.RECT_DIMS(1)/2:obj.RECT_DIMS(1)/2)', (-obj.RECT_DIMS(2)/2:obj.RECT_DIMS(2)/2)');
            our_rect = [x y];
            
            % frames
            [frame_grid_hyp] = obj.P.generate_seeds('pixels', our_frame_rect, obj.arrangement, obj.frame_grid_dims, [], obj.windows);
            
            hyp_conns = []; %[FIX THIS PART, add grid of squares inside these frames]
            types = [];
            %parfor i=1:numel(frame_grid_hyp)             
            for i=1:numel(frame_grid_hyp)                           
                duh = obj.P.generate_seeds('pixels', our_rect, 'grid', obj.grid_dims, frame_grid_hyp{i});
                internal_grid_hyp = duh;
                for j=1:numel(internal_grid_hyp)
                    internal_grid_hyp{j} = intersect(internal_grid_hyp{j}, frame_grid_hyp{i}); % without this can crash, when there's something in the border      
                    [hyp_conn, type] = obj.subframe_plus_internal_unary(internal_grid_hyp{j}, frame_grid_hyp{i});                    
                    hyp_conns = [hyp_conns; hyp_conn];    
                    types = [types; type];
                end
            end
           types = types(1,:);
        end
                
        %%%% not converted to new code %%%%
        function [hyp_conns, types] = generate_frame_from_bbox(obj)   
            bbox_dimx = obj.bbox(3)-obj.bbox(1)+1;
            bbox_dimy = obj.bbox(4)-obj.bbox(2)+1;
            [pixel_coords_x, pixel_coords_y] = cartprod_mex((3:bbox_dimx-2)', (3:bbox_dimy-2)');         
            inside_the_frame = sub2ind([bbox_dimy bbox_dimx], pixel_coords_y, pixel_coords_x);  
            outside_the_frame = setdiff(1:bbox_dimx*bbox_dimy, inside_the_frame);

            hyp_conns{1,1} = {outside_the_frame, obj.Background, obj.SEED_FRAME_WEIGHT, 0};
            hyp_conns{1,2} = {inside_the_frame, obj.Foreground, 0, 1};            
            
            types{1,1} = {'pixels', 'classes'};    
            types{1,2} = {'pixels', 'classes'};   
        end
        
        function [hyp_conns, types] = generate_frame_from_bbox_unary(obj)
            % don't add hyp_conns with pixels outside bbox!!
            % it has the frame of coordinates in the bbox
            CONSTANT = 0.1;
            bbox_dimx = obj.bbox(3)-obj.bbox(1)+1;
            bbox_dimy = obj.bbox(4)-obj.bbox(2)+1;
            [pixel_coords_x, pixel_coords_y] = cartprod_mex((3:bbox_dimx-2)', (3:bbox_dimy-2)');         
            inside_the_frame = sub2ind([bbox_dimy bbox_dimx], pixel_coords_y, pixel_coords_x);  
            on_the_frame = setdiff(1:bbox_dimx*bbox_dimy, inside_the_frame);
            
            [pixel_coords_x_est, pixel_coords_y_est] = cartprod_mex((obj.bbox(1):obj.bbox(3))', (obj.bbox(2):obj.bbox(4))');
            ids_for_unary_estimation = setdiff(1:size(obj.I,1)*size(obj.I,2), sub2ind([size(obj.I,1) size(obj.I,2)], pixel_coords_y_est, pixel_coords_x_est));
            region = zeros(size(obj.I,1), size(obj.I,2));
            region(ids_for_unary_estimation) = 1;
            ids_for_estimation = find(bwperim(~region));
            %ids_for_estimation = find(region);
            unary_values = CONSTANT*obj.compute_unary_values(ids_for_estimation); %#ok<FNDSB>
            
            our_rect = obj.P.generate_rectangle_coords(obj.RECT_DIMS);            
            [fg_grid_hyp] = obj.P.generate_seeds('pixels', our_rect, obj.arrangement, obj.grid_dims);
                
            for i=1:prod(obj.grid_dims)
                hyp_conns{i,1} = {on_the_frame, obj.Background, obj.SEED_FRAME_WEIGHT, 0};                                
                hyp_conns{i,2} = {inside_the_frame, obj.Background, unary_values(inside_the_frame), 0};                     
                hyp_conns{i,3} = {inside_the_frame, obj.Foreground, 0, 1};        
                hyp_conns{i,4} = {fg_grid_hyp{i}, obj.Foreground, inf, 0};        
                types{i,1} = {'pixels', 'classes'};    
                types{i,2} = {'pixels', 'classes'}; 
                types{i,3} = {'pixels', 'classes'}; 
                types{i,4} = {'pixels', 'classes'}; 
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Interactive
        
       function [hyp_conns, types] = generate_interactive_subframe_plus_internal(obj)
            disp('insert inner seed');
            the_rectangle = obj.P.get_vars_interactive('pixels');  
            disp('insert outer subframe');
            the_frame = obj.P.get_vars_interactive('pixels');
            
            [hyp_conns, types] = obj.subframe_plus_internal(the_rectangle, the_frame);
       end       
       
       function [hyp_conns, types] = generate_interactive_subframe_plus_internal_unary(obj)
            disp('insert inner seed');
            region = obj.P.get_vars_interactive('pixels');  
            disp('insert outer subframe');
            frame = obj.P.get_vars_interactive('pixels');                            
                       
            [hyp_conns, types] = obj.subframe_plus_internal_unary(region, frame);
       end           
         
        function [hyp_conns, types] = generate_interactive_growth_hyps(obj)
            region = obj.P.get_vars_interactive('pixels');            
            frame = obj.P.generate_img_frame('pixels', 'all', 1);            
            all_but_frame = setdiff(obj.Pixels, frame);
            
            [hyp_conns, types] = subframe_plus_internal(obj, region, all_but_frame);
         end
        
        function [hyp_conns, types] = generate_interactive_growth_hyps_unary(obj)
            region = obj.P.get_vars_interactive('pixels');            
            frame = obj.P.generate_img_frame('pixels', 'all',1);
            all_but_frame = setdiff(obj.Pixels, frame);
            
            [hyp_conns, types] = subframe_plus_internal_unary(obj, region, all_but_frame);
        end             
       
        %
        % basic functionalities, shared by interactive and non-interactive modes       
        %
        
        function [hyp_conns, types] = internal_unary(obj, the_rectangle)
            unary_values = obj.compute_unary_values(the_rectangle);
            bground_pixels = setdiff(obj.Pixels, the_rectangle);
            hyp_conns{1,1} = {the_rectangle, obj.Foreground, inf, 0};
            hyp_conns{1,2} = {bground_pixels, obj.Foreground, unary_values(bground_pixels), 0};            
            hyp_conns{1,3} = {bground_pixels, obj.Background, 0, 1};
        
            types{1,1} = {'pixels', 'classes'};
            types{1,2} = {'pixels', 'classes'};                   
            types{1,3} = {'pixels', 'classes'};        
        end
        
        function [hyp_conns, types] = internal(obj, the_rectangle)
            bground_pixels = setdiff(obj.Pixels, the_rectangle);
            hyp_conns{1,1} = {the_rectangle, obj.Foreground, inf, 0};
            hyp_conns{1,2} = {bground_pixels, obj.Background, 0, 1};
            
            types{1,1} = {'pixels', 'classes'};
            types{1,2} = {'pixels', 'classes'};
        end
        
        function [hyp_conns, types] = subframe(obj, the_frame)            
            bground_pixels = setdiff(setdiff(obj.Pixels, the_frame), obj.img_frame);
                
            hyp_conns{1,1} = {bground_pixels, obj.Background, inf, 0};
            hyp_conns{1,2} = {the_frame, obj.Foreground, 0, 1};
            hyp_conns{1,3} = {intersect(obj.img_frame, the_frame), obj.Background, obj.SEED_FRAME_WEIGHT, 0};

            types{1,1} = {'pixels', 'classes'};    
            types{1,2} = {'pixels', 'classes'};     
            types{1,3} = {'pixels', 'classes'};     
                        
            %if(isempty(bground_pixels))
            %    types(1) = [];
            %    hyp_conns(1) = [];
            %end
        end
        
        function [hyp_conns, types] = subframe_unary(obj, the_frame)
            CONSTANT = 1;
            CONSTANT_DIVISION = 10;
            bground_pixels = setdiff(setdiff(obj.Pixels, the_frame), obj.img_frame);
            
            region = zeros(size(obj.I,1), size(obj.I,2));
            region(the_frame) = 1;
            ids_for_estimation = find(bwperim(region));
            unary_values = CONSTANT*obj.compute_unary_values(ids_for_estimation); %#ok<FNDSB>
            hyp_conns{1,1} = {bground_pixels, obj.Background, inf, 0};            
            hyp_conns{1,2} = {obj.img_frame, obj.Background, obj.SEED_FRAME_WEIGHT, 0};
            hyp_conns{1,3} = {the_frame, obj.Background, unary_values(the_frame)/CONSTANT_DIVISION, 0};
            hyp_conns{1,4} = {the_frame, obj.Foreground, 0, 1};

            types{1,1} = {'pixels', 'classes'};            
            types{1,2} = {'pixels', 'classes'};
            types{1,3} = {'pixels', 'classes'};
            types{1,4} = {'pixels', 'classes'};
            
            %if(isempty(bground_pixels))
            %    types(1) = [];
            %    hyp_conns(1) = [];
            %end
        end
                
        function [hyp_conns, types] = subframe_unary_iterative(obj, frame)
            bground_pixels = setdiff(setdiff(obj.Pixels, frame), obj.img_frame); 
            
            % external seed - topological constraint            
            hyp_conns{1,1} = {bground_pixels, obj.Background, inf(numel(bground_pixels),1), 0};                              
            
            [unary_fg, unary_bg] = obj.compute_unary_values_iterative(frame); 
            hyp_conns{1,2} = {frame, obj.Background, unary_bg, 0};
            hyp_conns{1,3} = {frame, obj.Foreground, unary_fg, 0};

            % parametric potential
            hyp_conns{1,4} = {frame, obj.Foreground, 0, 1};           
            
            % always have some frame soft bias to belong in the background
            hyp_conns{1,5} = {obj.img_frame, obj.Background, obj.SEED_FRAME_WEIGHT, 0};
            
            % parametric potential
            hyp_conns{1,6} = {frame, obj.Background, 0, 1};           
            
            types{1,1} = {'pixels', 'classes'};
            types{1,2} = {'pixels', 'classes'};            
            types{1,3} = {'pixels', 'classes'};
            types{1,4} = {'pixels', 'classes'};            
            types{1,5} = {'pixels', 'classes'};
            types{1,6} = {'pixels', 'classes'};
        end
        
        function [hyp_conns, types] = subframe_plus_internal(obj, the_rectangle, the_frame)
            frame_but_rect = setdiff(the_frame, the_rectangle);
            bground_pixels = setdiff(setdiff(obj.Pixels, the_frame), obj.img_frame); 
            hyp_conns{1,1} = {the_rectangle, obj.Foreground, inf, 0};
            hyp_conns{1,2} = {bground_pixels, obj.Background,inf, 0};          
            hyp_conns{1,3} = {frame_but_rect, obj.Foreground, 0, 1};                        
            frame_cost = obj.frame_bg_bias(the_rectangle, size(obj.I));
            hyp_conns{1,4} = {obj.img_frame, obj.Background, frame_cost*obj.SEED_FRAME_WEIGHT, 0};
            %hyp_conns{1,5} = {frame_but_rect, obj.Background, 0, 1};      
            hyp_conns{1,5} = {frame_but_rect, obj.Background, 0, 1};      

            types{1,1} = {'pixels', 'classes'};
            types{1,2} = {'pixels', 'classes'};
            types{1,3} = {'pixels', 'classes'};            
            types{1,4} = {'pixels', 'classes'};
            %types{1,5} = {'pixels', 'classes'};
            types{1,5} = {'pixels', 'classes'};
            
%             if(isempty(bground_pixels))
%                 types(2) = [];
%                 hyp_conns(2) = [];
%             end
        end
        
        function [hyp_conns, types] = subframe_plus_internal_unary(obj, region, frame)
            [non_seed_pixels, ids] = setdiff(frame, region);
            bground_pixels = setdiff(setdiff(obj.Pixels, frame), obj.img_frame); 
            
            % internal seed - topological constraint
            hyp_conns{1,1} = {region, obj.Foreground, inf, 0}; 
            % external seed - topological constraint            
            hyp_conns{1,2} = {bground_pixels, obj.Background, inf , 0};            
            % unary potential reflecting similarity to bground pixels            
            zero_one_img = zeros(size(obj.I,1), size(obj.I,2));
            zero_one_img(frame) = 1;
            pixels_for_estimation = find(logical(bwperim(zero_one_img)));
            unary_values = obj.compute_unary_values(pixels_for_estimation); %#ok<FNDSB>
            hyp_conns{1,3} = {non_seed_pixels, obj.Background, unary_values(non_seed_pixels), 0};
            % unary potential reflecting similarity to region pixels          
            unary_values = obj.compute_unary_values(region);
            hyp_conns{1,4} = {non_seed_pixels, obj.Foreground, unary_values(non_seed_pixels), 0};
            % parametric potential
            hyp_conns{1,5} = {non_seed_pixels, obj.Foreground, 0, 1};                                    
            frame_cost = obj.frame_bg_bias(region, size(obj.I));
            %frame_cost = 1;
            % always have some frame soft bias to belong in the background
            hyp_conns{1,6} = {obj.img_frame, obj.Background, frame_cost*obj.SEED_FRAME_WEIGHT, 0};
            % parametric potential
            hyp_conns{1,7} = {non_seed_pixels, obj.Background, 0, 1};
            
            types{1,1} = {'pixels', 'classes'};
            types{1,2} = {'pixels', 'classes'};            
            types{1,3} = {'pixels', 'classes'};
            types{1,4} = {'pixels', 'classes'};            
            types{1,5} = {'pixels', 'classes'};            
            types{1,6} = {'pixels', 'classes'};
            types{1,7} = {'pixels', 'classes'};


%% grow from the outside
%             [non_seed_pixels, ids] = setdiff(frame, region);
%             bground_pixels = setdiff(setdiff(obj.Pixels, frame), obj.img_frame); 
%             hyp_conns{1,1} = {region, obj.Foreground, inf, 0}; 
%             hyp_conns{1,2} = {bground_pixels, obj.Background, inf , 0};    
%             unary_values = obj.compute_unary_values(region);
%             hyp_conns{1,3} = {non_seed_pixels, obj.Foreground, unary_values(non_seed_pixels), 0};
% %             hyp_conns{1,4} = {non_seed_pixels, obj.Foreground, 0, 1};
%             
%             types{1,1} = {'pixels', 'classes'};
%             types{1,2} = {'pixels', 'classes'};
%             types{1,3} = {'pixels', 'classes'};
%             types{1,4} = {'pixels', 'classes'};
        end        
        
        %
        % Code for smoothing unary terms spatially (e.g. away from seed )
        % 
        
        % for spatially weighting smooth negative seeds ( non-inf potentials for the frame )
        function [frame_cost] = frame_bg_bias(obj, internal_region, img_sz)
            [ir, jr] = ind2sub(img_sz, internal_region);
            [ifr, jfr] = ind2sub(img_sz, obj.img_frame);
            
            center = mean([ir jr]);            
            dist = [ifr jfr] - center(ones(numel(ifr),1),:);
            frame_cost = sqrt(sum(dist.*dist,2));
            frame_cost = frame_cost/max(frame_cost);            
        end
        
        function val = smooth_spherically(obj, seed_pixels, smoothing_radius, val)
            img_sz = [size(obj.I,1) size(obj.I,2)];
            [ir, jr] = ind2sub(img_sz, seed_pixels);
            [ifr, jfr] = ind2sub(img_sz, 1:prod(img_sz));            
            center = mean([ir jr]);         
            dists = [ifr; jfr]' - center(ones(numel(ifr),1),:);
            dist = sqrt(sum(dists.*dists,2));
            val = val.*(1 - scale_data(dist', 'zero_one')');
        end
    end
    
    methods (Abstract)
        
       [obj] = add_hyp_conns(obj)                    
                
       [hyp_conns, types] = generate_interactive_hyp(obj)       
       
    end   
end
