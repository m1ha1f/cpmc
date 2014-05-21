% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function [the_masks, the_additional] = SvmSegm_extract_segments_nosave(img_dir, the_img_names, options)
	the_additional = []; % we might not return anything here

    the_segm_pars = options;    
    window_gen_parms = options.window_gen_parms;
    windows_path = options.windows_folder;
    min_n_pixels = options.min_n_pixels;
    pb_path = options.pb_folder;
        
    if(isempty(pb_path))
        pb_path = img_dir;
    end
    
    if(isempty(windows_path))
        windows_path = img_dir;
    end
    
    for j=1:length(the_img_names)
        t = tic();
        masks = logical([]);
        
        for i=1:length(the_segm_pars.segm_methods)
            img_name = the_img_names{j};
            
            img_type = get_img_type(img_dir, the_img_names{j});
            
            I = imread([img_dir img_name img_type]);

            if(strcmp(the_segm_pars.segm_methods{i}, 'UniformSegmenter'))
                fprintf('Extracting segments using an uniform unary term...\n');
                s = UniformSegmenter(I,img_name); % without region model
            elseif(strcmp(the_segm_pars.segm_methods{i}, 'LongRangeSegmenter'))
                fprintf('Extracting segments using an color-based unary term...\n');
                s = LongRangeSegmenter(I,img_name); % with color-based region model
            elseif(strcmp(the_segm_pars.segm_methods{i}, 'GridOfFramesSegmenter'))
                fprintf('Extracting segments inside subframes...\n');
                s = GridOfFramesSegmenter(I,img_name);
            else
                error('no such method');
            end
            
            s = s.set_pb_path(pb_path);
            s = s.set_windows_path(windows_path, window_gen_parms);
            
            if(isfield(options, 'resize_factor') && ~isempty(options.resize_factor))
                s.resize_factor = options.resize_factor;
            end
            if(isfield(options, 'dont_break_disconnected') && ~isempty(options.dont_break_disconnected{i}))
                dont_break_disconnected = options.dont_break_disconnected{i};
                s.DONT_BREAK_DISCONNECTED = dont_break_disconnected;
            end
            if(isfield(options, 'sigma') && ~isempty(options.sigma{i}))
                sigma = options.sigma{i};
                s.params.CONTRAST_SENSITIVE_SIGMA = sigma;
            end
                        
            if(isfield(options, 'grid_dims') && ~isempty(options.grid_dims))
                s.grid_dims = options.grid_dims{i};
            end
            
            if(isfield(options, 'local_dist_type') && ~isempty(options.local_dist_type))
                s.params.local_dist_type =  options.local_dist_type;
            end
            
             if(isfield(options, 'filter_segments') && ~isempty(options.filter_segments{i}))
                 s.filter_segments = options.filter_segments{i};
             end

             if(isfield(options, 'max_energy') && ~isempty(options.filter_segments{i}))
                 s.MAX_ENERGY= options.max_energy{i};
             end

             if(isfield(options, 'randomize_N'))
                 s.randomize_N= options.randomize_N;
             end
             
             if(isfield(options, 'min_n_pixels'))
                 s.MIN_NPIXELS = min_n_pixels;
             end
             
            s = s.initialize();            

            if(isfield(options, 'max_n_segms') && ~isempty(options.filter_segments{i}))
                max_n_segms = options.max_n_segms(i);
            else
                max_n_segms = [];
            end
            s = s.compute_segments(max_n_segms);

            resh_masks = reshape(s.Segments, size(I,1), size(I,2), size(s.Segments,2));
            filtered_masks = resh_masks;
            if(isfield(options, 'morph_open') && options.morph_open)
                for k=1:size(resh_masks,3)
                    filtered_masks(:,:,k) = bwmorph(resh_masks(:,:,k), 'open');
                end
            end
            filtered_masks = reshape(filtered_masks, size(I,1)*size(I,2), size(s.Segments,2));
            
            %subplot_auto_transparent(s.Segments,I)
            masks = cat(2, masks, filtered_masks);           
        end
        
        if((numel(the_segm_pars.segm_methods) ~= 1) && (sum(options.max_n_segms) ~= inf))
            % filter repeated segments
            max_n_segms = sum(options.max_n_segms);
            pb_file = [pb_path '/' the_img_names{j} '_PB.mat'];
            resize = 1;
            MIN_PIXELS = min_n_pixels*options.resize_factor;
            MIN_DISSIMILARITY = 0.05;
            min_n_segments = 5;
            no_separate = false;
            S = SegmentProcessor([], masks, I, pb_file, MIN_PIXELS, no_separate, resize);
            t_en = tic();
            S = S.compute_energies();
            t_energy = toc(t_en)
            t_diss = tic();
            n_before = size(S.segments,2);
            
            
            S = S.remove_repeated_segments();
            S = S.filter_segments('min_dissimilarity', MIN_DISSIMILARITY, min_n_segments, max_n_segms);
            n_after = size(S.segments,2);
            time_removing_similar = toc(t_diss);
            obj.Segments = S.segments;
            masks = S.segments;
        end
        time_segmentation = toc(t);
        
        masks = logical(reshape(masks, size(I,1), size(I,2), size(masks,2)));     
        the_masks{j} = masks;		
    end
end

function img_type = get_img_type(img_dir, img_name)
    files = dir([img_dir img_name '*']);
    [name_len, id] = sort(arrayfun(@(a) numel(a.name), files), 'ascend');
    img_type = files(id(1)).name(end-3:end);
end