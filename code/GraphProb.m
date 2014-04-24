% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef GraphProb
    % the philosophy is that every time a variable or seed is added, the
    % prob_dgraph is updated.
    % The ids of each type should form a sequence starting in 1.
    properties
        I
        
        mapVars
        n_vars                        
                
        mapAggs
                
        hypConn        
        prob_dgraph
        parametric_dgraph
        upper_bp        
        
        solution
    end
    
    methods (Access=public)     
        function obj = GraphProb(I) 
            if(exist('I', 'var'))
                obj.I = I;
            end
            
            obj.mapVars = containers.Map();
            obj.mapAggs = containers.Map();            
            obj.n_vars = 0;
            obj.hypConn = [];
            obj.prob_dgraph = Dgraph([]);
            obj.parametric_dgraph = Dgraph([]);
            obj.solution = [];
            obj.upper_bp = 50;            
                    
       
        end       
        
        %%%% Basic functionality %%%%
        function [obj] = add_atomic_vars(obj, type, affinities)  
            % assigns graph ids to new variables            
            new_ids =  uint32(obj.n_vars + (1:numel(affinities.nodes)))';    
            obj.n_vars = obj.n_vars + numel(affinities.nodes);        
            n_new = length(new_ids);

            if(obj.mapVars.isKey(type))
                obj.mapVars(type) = [obj.mapVars(type) new_ids];
            else
                obj.mapVars(type) = new_ids;
            end

            % update problem graph
            new_dgraph = affinities.pairwise_to_dgraph();
            obj.prob_dgraph = obj.prob_dgraph.add_size(obj.n_vars);
            obj.parametric_dgraph = obj.parametric_dgraph.add_size(obj.n_vars);
            obj.prob_dgraph = obj.prob_dgraph.set_submatrix(obj.n_vars-n_new+1:obj.n_vars, new_dgraph);
        end
        
        function [obj] = add_aggregate_vars(obj, type_father, affinities, type_children, id_children)
            assert(iscell(id_children));
            [obj] = obj.add_atomic_vars(type_father, affinities);
            internal_ids_children = obj.mapVars(type_children);
            conv_ids_children = cell(size(id_children));
            for i=1:length(id_children)
                conv_ids_children{i} = internal_ids_children(id_children{i});
            end
            
            if(obj.mapAggs.isKey(type_father))
                obj.mapAggs(type_father)= [obj.mapAggs(type_father) conv_ids_children];
            else
                obj.mapAggs(type_father) = conv_ids_children;
            end
        end
        
        function [obj] = add_connections(obj, types, conns)
            % a connection is composed of two sets and one value. 

            for i=1:size(conns,1)         
                conn = conns(i,:);
                
                ids1 = obj.mapVars(types{1});
                vars1 = ids1(conn{1});
                ids2 = obj.mapVars(types{2});
                vars2 = ids2(conn{2});
                val = conn{3};
                
                [vars1, vars2] = obj.make_equal_size(vars1, vars2);
                 %obj.show_regions( types{1}, vars1) 
                 
                % always symmetric
                obj.prob_dgraph = obj.prob_dgraph.add_edges([vars1 vars2; vars2 vars1], val*ones(2*numel(vars1),1));
            end
        end
        
        function [obj] = add_agg_connections(obj,type, val)            
            ids = obj.mapVars(type);
            aggs = obj.mapAggs(type);

            for i=1:length(ids)
                vars1 = ids(i);
                vars2 = aggs{i};            

                if(numel(vars1)>numel(vars2))
                    assert(numel(vars2)==1);
                    vars2 = vars2*ones(numel(vars1),1);
                else
                    assert(numel(vars1)==1);
                    vars1 = vars1*ones(numel(vars2),1);
                end
                obj.prob_dgraph = obj.prob_dgraph.add_edges([vars1 vars2; vars2 vars1], val*ones(2*numel(vars1),1));
            end           
        end
        
        function [obj] = add_hyp_connections(obj, types, sets)
            if(isempty(sets))
                return;
            end
            
            ids1 = cell(size(sets,2),1);
            ids2 = cell(size(sets,2),1);
            for i=1:size(sets,2)
                ids1{i} = obj.mapVars(types{i}{1});
                ids2{i} = obj.mapVars(types{i}{2});
            end
            
            %n_hypothesis_already = size(obj.hypConn,1);
            
            c = {Connection};
            new_hypConn = repmat(c, size(sets,1), size(sets,2));
            
            for i=1:size(sets,1)                
                s = sets(i,:);
                
                counter = 1;
                for j=1:size(s,2)
                    if(isempty(s{j}{1}))
                        continue;
                    end
                    c = Connection;
                    this_s = s{j};
                    
                    % if this fails, it means the number of values is
                    % different from the number of pairs of nodes
                    %assert(numel(this_s{3}) == max(numel(this_s{1}), numel(this_s{2})));
                    c.nodes_a = ids1{j}(this_s{1});
                    c.nodes_b = ids2{j}(this_s{2});
                    c.edge_strength = single(this_s{3});
                    c.parametric_weight = single(this_s{4});
                    %obj.hypConn{n_hypothesis_already+i,counter} = c;
                    new_hypConn{i, counter} = c;
                    counter = counter + 1;
                end
            end
            
            cols_a = size(obj.hypConn,2); 
            cols_b = size(new_hypConn,2);

            if(cols_a == 0)                
            elseif(cols_a>cols_b)
                new_hypConn = [new_hypConn cell(size(new_hypConn,1), cols_a - cols_b)];
            elseif(cols_b>cols_a)
                obj.hypConn = [obj.hypConn cell(size(obj.hypConn,1),cols_b - cols_a)];
            end
            
            obj.hypConn = [obj.hypConn; new_hypConn];
                        
            %obj.show_seed_hyps();           
        end
        
        function [obj] = solve(obj, type, Css, Ct)   
            rows = size(obj.I, 1);
            cols = size(obj.I, 2);
            partitions = cell(1, numel(Css));

            penalty = 0.5*ones(rows,cols);

            for i = 1: numel(Css)
                Cs = Css{i};

                varParas = [rows; cols; 300; 1e-4; 0.3; 0.16];
%                para 0,1 - rows, cols of the given image
%                para 2 - the maximum number of iterations
%                para 3 - the error bound for convergence
%                para 4 - cc for the step-size of augmented Lagrangian method
%                para 5 - the step-size for the graident-projection of p

                [uu, erriter,num,tt] = CMF_GPU(single(penalty), single(Cs), single(Ct), single(varParas));

                partitions{i} = reshape(uu, rows*cols, 1);
            end
            
            obj.solution = ~cell2mat(partitions); % 1 is sink!
        end
        
        function [solution_for_type] = get_type_solution(obj,type)
            ids = obj.mapVars(type);
            solution_for_type = obj.solution(ids,:);
        end
        
        %%%% Variable Set Generation %%%%%
                       
        function [set] = generate_img_frame(obj, type, frame_kind, thickness)
            % frame_kind can be 'all', 'all_but_down', 'horiz', 'vert'.
            %assert(strcmp(type,'pixels'));
            if(nargin ==3)
                thickness = 1;
            end
            pixel_set = frame_pixel_ids(size(obj.I,1), size(obj.I,2), thickness, frame_kind);             
            
            set = [];
            % now convert to actual variable ids
            if(strcmp(type,'pixels'))
                set = pixel_set;
            elseif(strcmp(type, 'superpixels'))
                aggs = obj.mapAggs(type);
                % quite slow, should try improving
                p = pixel_set;
                for j=1:length(aggs)
                    agg = aggs{j};
                    inters = intersect(agg, p);
                    if(~isempty(inters))
                        set = [set j];
                    end
                end
            else
                error('no such type');
            end
        end
        
        function [coords] = generate_rectangle_coords(obj, dims)
            assert(length(dims) == 2); % two numbers required
            if(any(round(dims)~=dims))
                error('expecting integer dims');
            end
            
            if(all(dims == 1))
                coords = [0 0];
            else
                range_x = (0:dims(1)-1)';
                range_y = (0:dims(2)-1)';
                
                [x,y] = cartprod_mex(range_x,range_y);
                
                coords = [(x - (dims(1)-1)/2)  (y - (dims(2)-1)/2)] ;
                
                if(mod(dims(1),2)==0) % if it's even move it right / up
                    coords(:,1) = coords(:,1) + 0.5;
                end
                
                if(mod(dims(2),2)==0) % if it's even move it right / up
                    coords(:,2) = coords(:,2) + 0.5;
                end
            end        
        end
        
        function [sets, obj] = generate_seeds(obj, type, shape_coords, arrangement, n, frame, windows)
            % generates the seeds
            % type is 'pixels' or 'superpixels'
            % shape_coords is the 0-centered coords of each seed, defining its shape
            % if arrangement is 'image_grid', n should be an array [n_rows n_cols]
            % if arrangement is 'ncuts', n is the number of ncut regions            
            persistent orig_I;
            if(isempty(orig_I))
                orig_I = obj.I;
            end            
            
            if(isempty(n))
                sets = [];
                return;
            end
            
            the_frame = [];
            if(exist('frame', 'var') && ~isempty(frame)) 
                [obj.I, x_min, y_min, orig_sz] = obj.crop_img_to_frame(obj.I, frame);
                the_frame = frame;
            else
                x_min = 1; y_min= 1; orig_sz = size(obj.I);
            end
              
            if(strcmp(arrangement, 'grid'))
                assert(numel(n) == 2);
                sets = generate_img_grid(obj,type,shape_coords, n);
            elseif(strcmp(arrangement, 'felzenszwalb'))
                sets = generate_seeds_felzenswalb(obj, type, shape_coords, n);
            elseif(strcmp(arrangement, 'arbelaez'))                
                sets = generate_seeds_arbelaez(obj, type, shape_coords, n);
             elseif(strcmp(arrangement, 'ncuts'))
                sets = generate_seeds_ncuts(obj, type, shape_coords, n);
            elseif(strcmp(arrangement, 'windows'))         
                sets = generate_seeds_prec_windows(obj, windows);
            elseif(strcmp(arrangement, 'saliency'))
                [sets] = generate_seeds_saliency(obj,type, shape_coords,n, the_frame);
            end

            for i=1:numel(sets)
                un_sets = unique(sets{i});
                if(exist('frame', 'var'))
                    [sub_un_sets_1, sub_un_sets_2] = ind2sub([size(obj.I,1) size(obj.I,2)], un_sets);
                    sub_un_sets_1 = sub_un_sets_1 + x_min - 1;
                    sub_un_sets_2 = sub_un_sets_2 + y_min - 1;
                    un_sets = sub2ind(orig_sz, sub_un_sets_1, sub_un_sets_2);
                end
                sets{i} = un_sets;
            end

            %obj.show_seeds(orig_I, sets, true);
        end

        function [newI, x_min, y_min, orig_sz] = crop_img_to_frame(obj, I, frame)            
            orig_sz = [size(I,1) size(I,2)];
            [i,j] = ind2sub(orig_sz,frame);
            x_min = min(i);
            y_min = min(j);
            width = max(i)-min(i);
            height = max(j)-min(j);
            rect = [y_min x_min  height width];
            %theI = obj.I;
            newI = imcrop(I,rect);
        end
        
        function [sets] = generate_img_grid(obj, type, shape_coords, mn)
            if(sum(mn) == 0)
                sets = {};
                return;
            end
                
            shape_dims = [max(shape_coords(:,1)); ...
                                      max(shape_coords(:,2))]; % assumed to be centered coords
            n_rows = size(obj.I,1);
            n_cols = size(obj.I,2);
            
            vert_start = shape_dims(1)+1;
            vert_end = n_rows-shape_dims(1)-1;
            horiz_start = shape_dims(2)+1;
            horiz_end =  n_cols-shape_dims(2)-1;
            
            vert_step = ((vert_end - vert_start) / (mn(1)*2));
            horiz_step = ((horiz_end - horiz_start) / (mn(2)*2));
                      
            range_vert = vert_start:vert_step:vert_end;            
            range_horiz = horiz_start:horiz_step:horiz_end;                
            %assert(range_vert(end) == vert_end);
            %assert(range_horiz(end) == horiz_end);

            range_vert(3:2:end-2) = [];
            range_horiz(3:2:end-2) = [];
            
%             if(length(range_horiz)~=(mn(2)+2)) 
%                 horiz_step = round((horiz_end - horiz_start) / (mn(2) +1));
%                 range_horiz = horiz_start:horiz_step:horiz_end;
%             end            
            
            if(length(range_vert)~=(mn(1)+2)) 
                vert_step = round((vert_end - vert_start) / (mn(1) +1));
                range_vert = vert_start:vert_step:vert_end;
            end
            
            dummy = -1;
            if(isempty(range_vert))
                range_vert = [dummy max(1,floor(n_rows/2)) dummy];
            end
            if(isempty(range_horiz))
                range_horiz = [dummy max(1, floor(n_cols/2)) dummy];
            end
            
            range_vert([1 end]) = [];
            range_horiz([1 end]) = [];                        
            
            [X,Y] = meshgrid(range_vert, range_horiz);

            X = reshape(X, numel(X), 1);
            Y = reshape(Y, numel(Y), 1);
            assert(sqrt(length(X)) == mn(1));
            assert(sqrt(length(Y)) == mn(2));
            centers = round([X Y]);
            coord_sets = cell(length(X), 1);
            sets = coord_sets;
            pixel_sets = sets;
            
            for i=1:length(X)                
                this_center = repmat(centers(i,:), size(shape_coords,1), 1);
                
                c = round(max(1, this_center + shape_coords));
                c(:,1) = min(size(obj.I,1), c(:,1));
                c(:,2) = min(size(obj.I,2), c(:,2));
                pixel_sets{i} = sub2ind([n_rows n_cols], c(:,1), c(:,2));                            
            end
            
            % now convert to actual variable ids
            if(strcmp(type,'pixels'))
                sets = pixel_sets;
            else
                aggs = obj.mapAggs(type);
                % quite slow, should try improving
                for i=1:length(pixel_sets)                    
                    p = pixel_sets{i};
                    for j=1:length(aggs)         
                        agg = aggs{j};
                        inters = intersect(agg, p);
                        if(~isempty(inters))
                            sets{i} = [sets{i} j];
                        end
                    end
                end
            end
            
            %obj.show_seeds(sets);
        end
                 
        function sets = generate_seeds_felzenswalb(obj, type, shape_coords, n)
            k1 = 0.1;
            %k1 = 0.001;
            k2 = 400;
            sets = {};
            
            the_min = prod(n);

            if(sum(n) == 0)
                return;
            end
            
            counter = 0;
            while(1)
                if(size(obj.I,3)==1)
                    obj.I = repmat(obj.I,[1 1 3]);
                end
                
                if(any(size(obj.I)==0))
                    ;
                end
                [L] = vgg_segment_gb(obj.I, k1, k2, 100, true);
                %sc(L, 'rand')
                
                un = unique(L);
                n_segms = numel(un);
                if(the_min > n_segms)
                    k1 = k1 + 0.1*rand();
                    %k2 = k2 + 200*rand();
                else
                    break
                end
                counter= counter+  1;
                if(counter>15)
                    break;
                end
            end
            %sc(L, 'rand')
                        
            n_segms = numel(un);
         
            %%%%%%%%%% selection %%%%%%%%%%%%%
            % get the n biggest
             n_p = zeros(n_segms,1);
             for i=1:n_segms
                 n_p(i) = sum(sum(L==un(i)));
             end
             [val, size_order] = sort(n_p, 'descend');
           
            id_selected = obj.select_around_grid(L, n, size_order);
            id_selected = unique(id_selected);
            
            newL = zeros(size(L,1), size(L,2));
            for i=1:numel(id_selected)
                newL(L==id_selected(i)) = i;
            end
            L = newL;
            
            centroids = regionprops(uint32(L), 'Centroid');
            centroids = round(reshape(struct2array(centroids), 2, numel(centroids)))';
                        
            centroids = centroids(:,[2 1]);
    
            assert(~any(any(isnan(centroids))));

            sets = obj.create_seed_coords(centroids, shape_coords);
%%%% see the selected ones
            %obj.show_selected(L, id_selected);
            
%             for i=1:prod(n)
%                 in_i = find(L==id_selected(i));
%                 %imshow(L==un(order(i)))
%                 
%                 rp = randperm(numel(in_i));
%                 rp = rp(1:size(shape_coords,1));
%                 sets{i} = in_i(rp);
%             end
        end
        
        function sets = generate_seeds_arbelaez(obj, type, shape_coords, n)
            the_min = prod(n);

            if(sum(n) == 0)
                return;
            end
            
            counter = 0;

            var = load([obj.img_name  '_PB.mat'], 'gPb_orient');
            gPb_orient = var.gPb_orient;
            ucm2 = contours2ucm(gPb_orient);

            error('this part missing');
            %sc(L, 'rand')
                        
            n_segms = numel(un);
         
            %%%%%%%%%% selection %%%%%%%%%%%%%
            % get the n biggest
             n_p = zeros(n_segms,1);
             for i=1:n_segms
                 n_p(i) = sum(sum(L==un(i)));
             end
             [val, size_order] = sort(n_p, 'descend');
           
            id_selected = obj.select_around_grid(L, n, size_order);
            
            newL = zeros(size(L,1), size(L,2));
            for i=1:numel(id_selected)
                newL(L==id_selected(i)) = i;
            end
            L = newL;
            
            centroids = regionprops(uint32(L), 'Centroid');
            centroids = round(reshape(struct2array(centroids), 2, numel(centroids)))';
                        
            centroids = centroids(:,[2 1]);
    
            assert(~any(any(isnan(centroids))));
            for i=1:size(centroids,1)
                this_center = centroids(i,:);
                this_center = this_center(ones(size(shape_coords,1), 1), :);
                %this_center = repmat(centers(i,:), size(shape_coords,1), 1);
                
                %coords = this_center + shape_coords;
                coords = round(max(1, this_center + shape_coords));
                coords(:,1) = min(size(obj.I,1), coords(:,1));
                coords(:,2) = min(size(obj.I,2), coords(:,2));                                
                sets{i} = sub2ind([size(obj.I,1) size(obj.I,2)],coords(:,1), coords(:,2));
            end

%%%% see the selected ones
            %obj.show_selected(L, id_selected);            
        end
        
        function sets = generate_seeds_ncuts(obj, type, shape_coords, n)
            nsegms = prod(n);
            segms = imncut(rgb2gray(obj.I), nsegms);      
            L = zeros(size(obj.I,1), size(obj.I,2));
            
            for i=1:nsegms
                L(logical(segms(:,i))) = i;
            end
            
            centroids = regionprops(uint32(L), 'Centroid');
            centroids = round(reshape(struct2array(centroids), 2, numel(centroids)))';
            
            
            centroids = centroids(:,[2 1]);
    
            assert(~any(any(isnan(centroids))));
            for i=1:size(centroids,1)
                this_center = centroids(i,:);
                this_center = this_center(ones(size(shape_coords,1), 1), :);
                %this_center = repmat(centers(i,:), size(shape_coords,1), 1);
                
                %coords = this_center + shape_coords;
                coords = round(max(1, this_center + shape_coords));
                coords(:,1) = min(size(obj.I,1), coords(:,1));
                coords(:,2) = min(size(obj.I,2), coords(:,2));                                
                pixel_sets{i} = sub2ind([size(obj.I,1) size(obj.I,2)],coords(:,1), coords(:,2));
            end
            sets = pixel_sets;

        end
        
        function [sets] = generate_seeds_saliency(obj, type, shape_coords, n, frame)
            persistent sal_struct;
            
            if(isempty(sal_struct))
                sal_struct = gbvs(obj.I);
            end
            
            saliency_map = sal_struct;
            salmap = saliency_map.master_map_resized;
            if(~isempty(frame))
                [salmap, min_x, min_y, orig_sz] = obj.crop_img_to_frame(salmap, frame); 
                %imshow(salmap)
            end            
            
            max_map=  imregionalmax(salmap).*salmap;
            [i, j, val] = find(max_map);
            
            [v, order] = sort(val, 'ascend');
            
            to_keep = order(end-prod(n)+1:end);            
            centroids = [i(to_keep) j(to_keep)];
            sets = obj.create_seed_coords(centroids, shape_coords);
        end
        
        function [sets] = generate_seeds_prec_windows(obj, windows)
            sets = {};
            
            windows = clipboxes(obj.I, windows);
            %showboxes(obj.I, windows);
            for i=1:size(windows,1)
                if(windows(i,3)>windows(i,1))
                    ids2 = floor(windows(i,1):windows(i,3))';
                else
                    ids2 = floor(windows(i,3):windows(i,1))';
                end
                if(windows(i,4)>windows(i,2))
                    ids1 = floor(windows(i,2):windows(i,4))';
                else
                    ids1 = floor(windows(i,4):windows(i,2))';
                end
                [rows,cols] = cartprod_mex(ids1, ids2);
                sets{i} = sub2ind([size(obj.I,1) size(obj.I,2)], rows, cols);
                %max(sets{i})
            end
        end                
        
        function sets = create_seed_coords(obj, centroids, shape_coords)
            for i=1:size(centroids,1)
                this_center = centroids(i,:);
                this_center = this_center(ones(size(shape_coords,1), 1), :);
                %this_center = repmat(centers(i,:), size(shape_coords,1), 1);
                
                %coords = this_center + shape_coords;
                coords = round(max(1, this_center + shape_coords));
                coords(:,1) = min(size(obj.I,1), coords(:,1));
                coords(:,2) = min(size(obj.I,2), coords(:,2));                                
                sets{i} = sub2ind([size(obj.I,1) size(obj.I,2)],coords(:,1), coords(:,2));
           end
        end
        
        function id_selected = select_n_farthest_away(obj, L, n, size_order)
            n = prod(n);
            % compute centroids
            c =  regionprops(L, 'Centroid');
            coords = struct2array(c);
            coords = reshape(coords, 2, numel(coords)/2)';
            
            D = pdist2(coords, coords);
            
            id_selected = size_order([1 2]);
            D(:,size_order([1 2])) = 0;
            for i=3:n                
                % find point farthest away to any of the previously selected
                [val, p] = max(mean(D(id_selected,:)));
                new_id = p;
                id_selected = [id_selected; new_id];
                D(:, new_id) = 0;
            end
        end
        
        function id_selected = select_largest(obj, n, size_order)
            n = prod(n);
            id_selected = size_order(1:n);
        end
       
        function id_selected = select_around_grid(obj, L, n, size_order)
            c =  regionprops(L, 'Centroid');
            coords = struct2array(c);
            coords = reshape(coords, 2, numel(coords)/2)';            

            [sets] = obj.generate_img_grid('pixels', [0 0], n);
            [Y, X] = ind2sub(size(L), cell2mat(sets));
            
            grid_coords = [X Y];
            D = pdist2(grid_coords, coords);
            id_selected = [];
            for i=1:prod(n)
                [val, closest] = min(D(i,:));
                id_selected = [id_selected closest];
                D(:,closest) = inf;
            end
            
            %imshow(L); hold on;
            %plot(X, Y, 'x');
        end
        
        %%%% Interaction and Visualization %%%%  
        
        function show_seed_hyps(obj)
            class_nodes = obj.mapVars('classes');
            
            masks_to_show = zeros(size(obj.I,1), size(obj.I,2), numel(obj.hypConn));
            for i=1:numel(obj.hypConn)
                h = obj.hypConn{i};
                
                if(isempty(intersect(class_nodes, h.nodes_a)))
                    mask = zeros(size(obj.I,1), size(obj.I,2));
                    mask(h.nodes_a) = double(h.edge_strength);
                    masks_to_show(:,:,i) = mask;                    
                end               
            end
            subplot_auto_transparent(masks_to_show, obj.I);
        end
        
        function show_seeds(obj, theI, seeds, on_img)
            if(nargin==2)
                on_img = false;
            end
            
            %imshow(obj.I);
            newI = zeros(size(theI,1), size(theI,2));
            for i=1:numel(seeds)
                newI(seeds{i}) = i;                
            end
            if(on_img)
                imshow(heatmap_overlay(theI , newI));
                %sc(sc(obj.I, 'gray') + sc(newI));%sc(cat(3, obj.I, newI), 'prob');
            else                                
                sc(newI, 'rand');
            end
        end
        
        function show_selected(obj, L, id_selected)
            L2 = zeros(size(L));
            for i=1:numel(id_selected)
                L2(L == id_selected(i)) = i;
            end
            sc(L2, 'rand');
        end
        
        function vars = get_vars_interactive(obj, type)
            % vars are returned in external ids, 1 to x %

            the_pixel_ids = GraphProb.select_rectangle(obj.I);
            
            %duh = rgb2gray(obj.I);
            %duh(conv_pixel_ids) = 0;
            %imshow(duh);
%             
            % get internal ids of pixels
            vars= [];
            in_pixel_ids = obj.mapVars('pixels');
            if(strcmp(type, 'pixels'))
                vars = the_pixel_ids;
            else % it's an aggregate
                type_ids = obj.mapVars(type);                
                aggs = obj.mapAggs(type);
                for i=1:length(type_ids)
                    inters = intersect(aggs{i}, in_pixel_ids(the_pixel_ids));
                    if(~isempty(inters))
                        %duh(aggs{i}) = 0;
                        vars = [vars i]; % external ids
                    end
                end
            end
            %imshow(duh);
            %pause;
        end
        
        
        function show_connectivity_to_vars(obj, type, var_ids)
            map = obj.mapVars(type);
            %if(obj.DEBUG_MODE)
                spy(obj.prob_dgraph.D(map, map));           
            %end
            weights = sum(obj.prob_dgraph.D(map(var_ids), :));
            weights = weights(map);
            
            obj.show_weighted_regions(type, weights)
        end
        
        function show_regions(obj, varargin)
            % arguments are pairs like 'superpixels', []. This show all
            % superpixel regions and doesn't highlight any of them
            optargin = size(varargin,2);
            assert(mod(optargin,2) == 0);
            
            theI = obj.I;            
            
            for h=1:(optargin/2)
                the_color = [255*rand(1,3)];
                type = varargin{1+ 2*(h-1)};
                ids = varargin{2+ 2*(h-1)};
                agg_ids = obj.mapVars(type);
                Aggs = obj.mapAggs(type);

                Iaggs = reshape(obj.cell_of_aggs_to_img_array(Aggs), size(obj.I,1)*size(obj.I,2), 1);
                
                [sel, the_aggs] = intersect(agg_ids, agg_ids(ids));
                
                theI = reshape(theI, size(theI,1)*size(theI,2), 3, 1);
                for i=1:length(the_aggs)
                    theI(Aggs{the_aggs(i)},:) = repmat(the_color, length(Aggs{the_aggs(i)}),1);
                end
                theI = reshape(theI, size(obj.I,1), size(obj.I,2), 3);                
            end
            
            subplot_auto_overlay(Iaggs, theI)
        end
        
        function show_weighted_regions(obj, type, weights)
            Iaggs= zeros(size(obj.I,1) * size(obj.I,2), 1);
            the_color = [255*rand(1,3)];
                        
            agg_ids = obj.mapVars(type);
            the_aggs = obj.mapAggs(type);
            pixel_ids = obj.mapVars('pixels');
            
            for i=1:length(the_aggs)
                Iaggs(pixel_ids(the_aggs{i})) = weights(i);
            end
            Iaggs = reshape(Iaggs, size(obj.I,1), size(obj.I,2));
            
            imshow(Iaggs/max(max(Iaggs)));
        end
        
        function show_results(obj)
            relevant_vars = double(obj.solution(obj.mapVars('pixels'),:));
            subplot_auto_transparent(relevant_vars, obj.I);
        end
        
        function show_results_agg(obj, type)
            % needs to be generalized, it won't work when there's
            % aggregates of aggregates! the elements of the aggregate must
            % be pixels so far.
            relevant_vars = double(obj.solution(obj.mapVars(type), :));
            
            mask2plot = zeros(size(obj.I,1) *size(obj.I,2),size(relevant_vars,2));
            for j=1:size(relevant_vars,2)                
                the_ones = find(relevant_vars(:,j));
                agg = obj.mapAggs(type);
                for i=1:length(the_ones)                                
                    mask2plot(agg{the_ones(i)}, j) = 1;
                end                            
            end
            
            subplot_auto_transparent(mask2plot, obj.I);
        end
    end
    
    methods(Static)        
        function [the_pixel_ids, bbox] = select_rectangle(I)
            f = figure('Visible','on','Position',[360,500,450,285]);
            
            h=subplot(121);
            imshow(I, 'Parent', h);
            
            h2 = uicontrol('Position', [20 20 200 40], 'String', 'Continue', ...
                'Callback', 'uiresume(gcbf)');
            h_rect = imrect(h);
            uiwait();
            pos = getPosition(h_rect); % pos is [x_min y_min width height]      
            
            x = round(pos(1));
            y = round(pos(2));
            width = round(pos(3));
            height = round(pos(4));
            
            bbox = [x y x+width y+height];
            [pixel_coords_x, pixel_coords_y] = cartprod_mex((x:x+width)', (y:y+height)');
            the_pixel_ids = sub2ind([size(I,1) size(I,2)], pixel_coords_y, pixel_coords_x);  
        end
    end
    
    methods (Access=private)
        
        function [vars1, vars2] = make_equal_size(obj, vars1, vars2)    
            if(numel(vars1) == numel(vars2))
                return;
            end
            
            if(numel(vars1)>numel(vars2))
                assert(numel(vars2)==1);
                vars2 = vars2*ones(numel(vars1),1);
            else
                assert(numel(vars1)==1);
                vars1 = vars1*ones(numel(vars2),1);
            end
        end
        
        function Iaggs = cell_of_aggs_to_img_array(obj, aggs)
            Iaggs = zeros(size(obj.I,1), size(obj.I,2));
            for i=1:length(aggs)
                Iaggs(aggs{i}) = i;
            end
        end        
    end
end
