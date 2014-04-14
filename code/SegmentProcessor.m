% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef SegmentProcessor
    properties
        segments
        energies
        deciding_energy_type
        graph
        overlap_matrix
        I        
        gPb_fat
        gPb_thin
        resize_factor        
        remaining_orig_ids
        
        boundary_pixel_ids
    end
    methods
        function [obj] = SegmentProcessor(graph, segments, I, pb_file, MIN_NPIXELS, no_separate, resize_factor) 
            assert(ndims(segments) == 2); % each segment in a column
            
            obj.I = I;          
            if(isempty(graph))
                p_gen = JustSegmGen();
                p_gen.CONTRAST_SENSITIVE_SIGMA = 1.2;
                p_gen.local_dist_type = 'pb_fat'; % 'pb'
                g = p_gen.apply(I, pb_file); % it's not using pb by default!!
                obj.graph = g.pairwise_to_dgraph();
            else                
                obj.graph = graph;
            end
            obj.segments = segments;            
            if(~islogical(obj.segments))
                obj.segments = logical(obj.segments);
            end

            if(exist('pb_file', 'var') && (exist(pb_file, 'file') || exist([pb_file '.mat'], 'file')))
                load(pb_file);
                obj.gPb_fat = mean(gPb_orient,3);
                obj.gPb_thin = gPb_thin;
            else
                [gx, gy] = gradient(double(rgb2gray(I)));
                grad_map = (max(abs(gx),abs(gy)));
                obj.gPb_fat = grad_map/50.0;
                gPb_thin = obj.gPb_fat;
                gPb_thin(0.05>gPb_thin) = 0;                
                obj.gPb_thin = gPb_thin .* bwmorph(gPb_thin, 'skel', inf);
            end
                        
            if(exist('resize_factor', 'var'))
                obj.gPb_fat = imresize(obj.gPb_fat, resize_factor);                
                obj.gPb_thin = imresize(obj.gPb_thin, resize_factor);
                s = reshape(obj.segments, size(obj.I,1), size(obj.I,2), size(obj.segments,2));
                obj.segments = imresize(s, resize_factor, 'nearest'); % this one needs to be fast
                obj.segments = reshape(obj.segments, size(obj.gPb_fat,1)*size(obj.gPb_fat,2), size(segments,2));
                obj.I = imresize(obj.I, resize_factor);
                obj.resize_factor = resize_factor;
            else
                obj.resize_factor = 1;
            end
            
            %%% in some cases adding the reverse is good, but generally will add lots of
            %%% background segments, it might be that we can filter them
            %%% out with regression            
            %obj.segments = [obj.segments ~obj.segments];
              
            
            if(~(nargin > 3))
                MIN_NPIXELS = 25;
            end
            % separate into connected components
            
            if(~exist('no_separate', 'var'))
                no_separate = true;
            end
            
            if(no_separate) % don't separate means keep all non-connected components in a single segment                
                %segms = obj.segments;
                %[obj, nconn] = obj.separate_conn_comp(MIN_NPIXELS);               
                %obj.segments = [segms(:,nconn~=1) obj.segments];
            else % separate and only keep separated versions
                [obj, nconn] = obj.separate_conn_comp(MIN_NPIXELS);               
            end                
            
            obj.deciding_energy_type = 'cut_ratio';
            
            obj.remaining_orig_ids = 1:size(obj.segments,2);
            
            obj.boundary_pixel_ids = frame_pixel_ids(size(obj.I,1), size(obj.I,2), 1, 'all');
        end
        
        
        function [obj] = set_deciding_energy(obj, type)
            assert(strcmp(type, 'cut') || strcmp(type, 'cut_ratio') || strcmp(type, 'normalized_cut') || strcmp(type, 'unbalanced_normalized_cut') || strcmp(type, 'inter_region_contour_energy'));
            obj.deciding_energy_type = type;
        end
        
        function [obj] = compute_energies(obj, get_all)            
            if(exist('get_all', 'var')) && get_all
                feats_to_compute ={'cut', 'cut_ratio', 'normalized_cut', 'unbalanced_normalized_cut', 'fraction_of_healthy_boundary', 'inter_region_contour_energy'};
            else
                feats_to_compute = {obj.deciding_energy_type};
            end
            
            max_value = max(max(obj.graph.D));
                        
            obj.energies.cut = zeros(size(obj.segments,2),1);
            
            the_graph = obj.graph;
            the_segments = obj.segments;
            the_energies = obj.energies;
            
                        
            resh_segments = reshape(obj.segments, size(obj.I,1), size(obj.I,2), size(obj.segments,2));            
                           
            %parfor (i=1:size(obj.segments,2),8)
            for i=1:size(obj.segments,2)
                for j=1:numel(feats_to_compute)            
                %for i=1:size(obj.segments,2)
                    % 1. cut (from 0 to 1)        
                    links_across = the_graph.D(~the_segments(:,i), the_segments(:,i));                
                    cut = sum(sum(links_across));
                    if(strcmp(feats_to_compute{j}, 'cut'))
                        the_energies(i).cut = cut;
                    elseif(strcmp(feats_to_compute{j}, 'cut_ratio'))
                        % 2. cut ratio (minimizing it is called the sparsest cut problem)                         
                        n_edges_across = sum(sum(links_across>0));
                        the_energies(i).cut_ratio = cut./n_edges_across;
                        %cut
                        %n_edges_across
                        %the_energies(i).cut_ratio
                        %subplot_auto_transparent(resh_segments(:,:,i), obj.I);
                        %pause;
                    elseif(strcmp(feats_to_compute{j}, 'normalized_cut'))
                        % 3. normalized cut (cut/edges_0v) + (cut/edges_1v)
                        edges_00 = sum(sum(the_graph.D(~the_segments(:,i),~the_segments(:,i))));
                        edges_11 = sum(sum(the_graph.D(the_segments(:,i),the_segments(:,i))));
                        edges_0v = cut + edges_00;
                        edges_1v = cut + edges_11;
                        the_energies(i).normalized_cut = cut * ((1/edges_0v) + (1/edges_1v));
                    elseif(strcmp(feats_to_compute{j}, 'unbalanced_normalized_cut'))
                        % 4. unbalanced normalized cut (cut/edges_11)
                        the_energies(i).unbalanced_normalized_cut = min(10000, cut / (edges_11));
                    elseif(strcmp(feats_to_compute{j}, 'fraction_of_healthy_boundary'))
                        % 5. number of healthy edges, normalized by number of
                        % edges crossing it, for different thresholds             
                        vals = full(links_across(links_across>0));
                        thresholds = max_value.*[0.05 0.1 0.2 0.4];
                    
                        %size(vals)
                        %size(repmat(thresholds, size(vals,1), 1))
                        %size(repmat(vals,1, size(thresholds,2)))
                        %disp('over with this one')
                        if(size(vals,2) > size(vals,1))
                          vals = vals';
                        end
                        if(isempty(vals))
                          the_energies(i).fraction_of_healthy_boundary = [0 0 0 0];
                        else
                            the_energies(i).fraction_of_healthy_boundary = sum( repmat(thresholds, size(vals,1), 1) > repmat(vals,1, size(thresholds,2)))/n_edges_across;
                        end
                    elseif(strcmp(feats_to_compute{j}, 'inter_region_contour_energy'))
                        bw = bwboundaries(resh_segments(:,:,i));       
                        bw = cell2mat(bw);
                        
                        %imshow(resh_segments(:,:,i)); hold on;
                        %plot(bw(:,2),bw(:,1),'g','LineWidth',2);
                        
                        % get the reunion here                                                
                        pixel_ids =  sub2ind(size(obj.I), bw(:,1), bw(:,2));
                        pixel_ids = setdiff(pixel_ids, obj.boundary_pixel_ids);
                        n = numel(pixel_ids);
                        s = sum(obj.gPb_thin(pixel_ids));
                        the_energies(i).inter_region_contour_energy = -s/n;
                    else
                        error('no such type of energy');
                    end
                end
            end
            
            obj.energies = the_energies;
        end
        
        function [obj, n_conn_comp] = separate_conn_comp(obj, min_npixels)
            if(nargin == 1)
                min_npixels = 25;
            end
            
            segms = full(obj.segments);
            n_dims = size(segms,1);
            n_conn_comp = zeros(size(segms,2),1);
            
            for i=1:size(segms,2)
                segm = reshape(segms(:,i), size(obj.I,1), size(obj.I,2));
                bw(i) = bwconncomp(segm);

                % remove the very small ones immediately
                n_cc = numel(bw(i).PixelIdxList);
                to_remove = false(n_cc,1);
                for j=1:n_cc
                    if(min_npixels > numel(bw(i).PixelIdxList{j})) % smaller than min_npixels is rejected
                        to_remove(j) = true;
                    end
                end               
                
                bw(i).PixelIdxList(to_remove) = [];
                bw(i).NumObjects = bw(i).NumObjects - sum(to_remove);
                n_conn_comp(i) = bw(i).NumObjects;
            end                        
            
            clear segms;
            n_comp = sum([bw(:).NumObjects]);
            new_segms = false(n_dims, n_comp);
            remaining_n_comps = n_comp;
            
            while (remaining_n_comps > 0)
                if(bw(1).NumObjects>0)
                    for j=1:bw(1).NumObjects
                        new_segms(bw(1).PixelIdxList{j},n_comp-remaining_n_comps+1) = true;
                        remaining_n_comps = remaining_n_comps - 1;
                    end
                end
                bw(1) = [];
            end

            obj.segments = new_segms;
        end
        
        function [obj] = filter_diversified_energy(obj, n_to_keep)
            scores = [obj.energies(:).(obj.deciding_energy_type)];
            obj = obj.compute_overlap_matrix();
            new_scores = diversify_ranking(-scores, obj.overlap_matrix);            
            [sorted_new_scores, ids] = sort(new_scores, 'descend');
            winners = ids(1:min(numel(ids), n_to_keep));
            obj.segments = obj.segments(:,winners);
            obj.remaining_orig_ids = obj.remaining_orig_ids(winners);
        end
        
        function [obj] = filter_segments(obj, type, threshold, min_n_segms, max_n_segms) 
            if(isempty(obj.segments))
                return;
            end
            
            if(exist('min_n_segms', 'var'))             
                min_n_segms = min(min_n_segms, size(obj.segments,2));
            else
                min_n_segms = 0;
            end
                      
            if(~exist('max_n_segms', 'var'))
                max_n_segms= inf;
            end
            
            if(isempty(obj.segments))
                return;
            end
            
            if(nargin == 2)
                switch type
                    case 'min_npixels'
                        threshold = 150;
                    case 'max_energy'
                        threshold = 200;
                    case 'min_dissimilarity'
                        threshold = 0.05;
                    otherwise
                        error('no such threshold type');
                end
            end
            
            if(strcmp(type, 'min_pixels'))
                if(obj.resize_factor ~=1)
                    threshold = threshold*obj.resize_factor;
                end
            end
            
            switch type
                case 'min_npixels' 
                    % remove small connected components
                    sums = sum(obj.segments);
                    to_remove = (threshold > sums);
                case 'max_energy'
                    if(isempty(obj.energies))
                        error('you need to compute energies first');
                    end
                    to_remove = ([obj.energies(:).(obj.deciding_energy_type)] > threshold);
                    % visualization
                    %subplot_auto_transparent(reshape(obj.segments(:,to_remove), size(obj.I,1), size(obj.I,2), sum(to_remove)), obj.I);
                    [sorted_energies,indices] = sort([obj.energies(:).(obj.deciding_energy_type)], 'ascend');
                    if(min_n_segms > sum(~to_remove))                        
                        to_remove(indices(1:min_n_segms)) = 0;
                        to_remove(indices(min_n_segms+1:end)) = 1;
                    end
                                       
                    to_keep = ~to_remove;
                    if(numel(to_keep) > max_n_segms)
                        to_keep = indices(1:max_n_segms);
                        to_remove = ~to_keep;
                    end
               case 'min_dissimilarity' % requires energies to be computed
                   %subplot_auto_transparent(reshape(obj.segments,  size(obj.I,1), size(obj.I,2), size(obj.segments,2)), obj.I)
                   all_to_remove = [];
                   segm_backup = obj.segments;
                   segm = reshape(obj.segments, size(obj.I,1), size(obj.I,2), size(obj.segments,2));
                   segm = imresize(segm, 0.6, 'nearest');  
                   obj.segments = reshape(segm, size(segm,1)*size(segm,2), size(obj.segments,2));
                   
                   %%% splitting the set of segments into chunks to
                   %%% speed up clustering.  version without chunks is
                   %%% below, and probably works better
                   
                    %CHUNK_SIZE = 100;
                    CHUNK_SIZE = 500;
                    
                    if(~isempty(obj.overlap_matrix))
                        % it's already computed, no need for speedup                        
                        n_chunks = 1;
                        
                        chunked_segments = {obj.segments};
                        chunked_energies = {[obj.energies.(obj.deciding_energy_type)]};
                    else                                            
                        n_chunks = ceil(size(obj.segments,2)/CHUNK_SIZE);
                        chunked_segments = chunkify(obj.segments, n_chunks);
                        chunked_energies = chunkify([obj.energies.(obj.deciding_energy_type)], n_chunks);
                    end
                    assert(sum(cellfun(@(a) size(a,2), chunked_segments)) == size(obj.segments,2));
                    to_keep = [];
                    counter = 0;
                    for i=1:n_chunks
                        if(size(obj.overlap_matrix,1) == 1)
                            return;
                        end
                        
                        obj.segments = chunked_segments{i};
                        these_energies = chunked_energies{i};
                        if(isempty(obj.overlap_matrix) || (n_chunks ~= 1))
                            obj = obj.compute_overlap_matrix();
                        end
                        
                        % linkage wants a dissimilarity matrix
                        Z = linkage(squareform(-(obj.overlap_matrix-1),'tovector'));
                        % dendrogram(Z)
                        
                        ids = cluster(Z, 'cutoff', threshold, 'Criterion','distance');
                        
                        % for those in a same cluster, remove the ones with
                        % smaller energy
                        cluster_ids = unique(ids);
                            
                        for j=1:length(cluster_ids)                            
                            cluster_id = cluster_ids(j);
                            in_cluster = (ids==cluster_id);
                            
                            [min_val, min_id] = min(these_energies(in_cluster));
%                             subplot_auto_transparent( chunked_segments{i}(:,in_cluster), imresize(obj.I, 0.6, 'nearest'), these_energies(in_cluster));                                                        
%                            min_id
%                             pause;
                            
                            ids_in_cluster = find(in_cluster);
                            to_keep = [to_keep (ids_in_cluster(min_id) + counter)];
                        end                  
                        counter = counter + size(chunked_segments{i},2);
                    end
                    to_remove = setdiff(1:length(obj.energies), to_keep);
                    all_to_remove = [all_to_remove to_remove];
                    
                    % do the final round
                    obj.segments = cell2mat(chunked_segments);
                    obj.segments(:,to_remove) = [];                        
                    obj.energies(to_remove) = [];                  
                    
                    if(size(obj.segments,2) > 1)
                        obj = obj.compute_overlap_matrix();
                        to_remove_further = SegmentProcessor.cluster_distance(-(obj.overlap_matrix-1), threshold, max_n_segms, [obj.energies.(obj.deciding_energy_type)]);
                        all_to_remove = [all_to_remove to_keep(to_remove_further)];
                        
                        obj.overlap_matrix(to_remove_further, :) = [];
                        obj.overlap_matrix(:,to_remove_further) = [];

                        obj.segments(:, to_remove_further) = [];
                        obj.energies(to_remove_further) = [];
                    end
                    
                    segm_backup(:, all_to_remove) = [];
                    obj.segments = segm_backup;                    
                    obj.remaining_orig_ids(all_to_remove) = [];
                    
                    assert(size(obj.segments,2) == numel(obj.energies));                    
                    return;

%%% Version without chunks
%                 case 'min_dissimilarity'
%                    segm_backup = obj.segments;
%                    segm = reshape(obj.segments, size(obj.I,1), size(obj.I,2), size(obj.segments,2));
%                    segm = imresize(segm, 0.6, 'nearest');  
%                    obj.segments = reshape(segm,size(segm,1)*size(segm,2), size(obj.segments,2));
%                     if(isempty(obj.overlap_matrix))
%                         obj = obj.compute_overlap_matrix();
%                        %error('you need to compute overlap matrix first!');
%                     end
%                     
%                     % linkage wants a dissimilarity matrix
%                     if(size(obj.overlap_matrix,1) == 1)
%                         return;
%                     end
%                        
%                     %obj.overlap_matrix = obj.overlap_matrix - eye(size(obj.overlap_matrix));
%                     Z = linkage(squareform(-(obj.overlap_matrix - 1),'tovector'), 'single'); % single is default
%                     % dendrogram(Z)
%                     ids = cluster(Z, 'cutoff', threshold, 'Criterion','distance');
%                     
%                     % for those in a same cluster, remove the ones with
%                     % smaller energy
%                     cluster_ids = unique(ids);
%                     to_keep = [];
%                     for i=1:length(cluster_ids)
%                         cluster_id = cluster_ids(i);
%                         in_cluster = (ids==cluster_id);
%                         [min_val, min_id] = min([obj.energies(in_cluster).(obj.deciding_energy_type)]);
%                         ids_in_cluster = find(in_cluster);                        
%                         %subplot_auto_transparent( obj.segments(:,in_cluster),obj.I);
%                         to_keep = [to_keep ids_in_cluster(min_id)];
%                     end
%                     to_remove = ones(numel(ids),1);
%                     to_remove(to_keep) = 0;
%                     to_remove = setdiff(1:length(ids), to_keep);
%                     obj.segments = segm_backup;
                otherwise
                    error('no such threshold type');                                        
            end
            
            if(islogical(to_remove))
                to_remove = find(to_remove);
            end            
            % there was a big bug in the first version, it is was randomizing this step when doing energy filtering!!
            obj.segments = obj.segments(:, setdiff(1:size(obj.segments,2), to_remove));
            %clear obj.segments;
            %obj.segments = segms;
            %clear segms;
            
            if(~isempty(obj.energies))
                obj.energies(to_remove) = [];
            end
            if(~isempty(obj.overlap_matrix) && size(obj.overlap_matrix,1) ~= numel(obj.energies))
                obj.overlap_matrix(to_remove,:) = [];
                obj.overlap_matrix(:,to_remove) = [];
            end
            
            obj.remaining_orig_ids(to_remove) = [];
        end
        
        function [obj] = remove_repeated_segments(obj)
            all_to_remove = [];
            segm_backup = obj.segments;
            segm = reshape(obj.segments, size(obj.I,1), size(obj.I,2), size(obj.segments,2));
            segm = imresize(segm, 0.4, 'nearest');
            obj.segments = reshape(segm, size(segm,1)*size(segm,2), size(obj.segments,2));
            overlap = segm_overlap_mex(obj.segments);              
            overlap = triu(overlap - diag(diag(overlap)));
            
            [ones_i, ones_j] = find(overlap==1);

            to_remove = false(size(overlap,1), 1);
            while(~isempty(ones_i))
                overlap(ones_i(1), :) = 0;
                overlap(:,ones_i(1)) = 0;
                to_remove(ones_i(1)) = true;
                [ones_i, ones_j] = find(overlap>0.98);                
            end

            obj.segments = segm_backup;
            obj.segments(:, to_remove) = [];
            if(~isempty(obj.energies))
                obj.energies(to_remove) = [];
            end
        end
    end
    
    methods(Static)
        function to_remove = cluster_distance(D, threshold, max_n_segms, energies)
            DefaultVal('energies', '[]');
            Z = linkage(squareform(D,'tovector'));
            % dendrogram(Z)
            if(max_n_segms == inf)
                ids = cluster(Z, 'cutoff', threshold, 'Criterion','distance');
            else
                ids = cluster(Z, 'MaxClust', max_n_segms);
            end

            % for those in a same cluster, keep the one with
            % smaller energy
            cluster_ids = unique(ids);
            to_keep = [];
            for i=1:length(cluster_ids)
                cluster_id = cluster_ids(i);
                in_cluster = (ids==cluster_id);
                ids_in_cluster = find(in_cluster);
                            
                if(~isempty(energies))
                    [min_val, min_id] = min(energies(in_cluster));
                else
                    min_id = 1;
                end
                
                to_keep = [to_keep ids_in_cluster(min_id)];
            end
            to_remove = setdiff(1:length(ids), to_keep);
        end
    end

    methods
        function display_segments(obj)
            if(isempty(obj.energies))
                subplot_auto_transparent(obj.segments, obj.I);
            else
                subplot_auto_transparent(obj.segments, obj.I, [obj.energies(:).(obj.deciding_energy_type)]);
            end
        end                
    end
    
    methods(Access=private)
        function [obj] = compute_overlap_matrix(obj)
            assert(islogical(obj.segments))
            %duh = reshape(obj.segments, size(obj.I,1), size(obj.I,2), size(obj.segments,2));
            
            %%% something like this will be much faster
            %for i=1:size(obj.segments,2)
            %    r(:,i) = struct2array(regionprops(duh(:,:,i)));
            %end
            %duh = normalize(r')';
            %obj.overlap_matrix = exp(-100*pdist2(duh',duh'));

            setenv('OMP_NUM_THREADS','8');
            obj.overlap_matrix = double(segm_overlap_mex(obj.segments));
        end
    end
end
