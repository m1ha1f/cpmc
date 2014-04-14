% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

% This class can be used to compute low and mid-level segment features
% Joao Carreira, June 2010

classdef SegmFeatures 
    properties
        features
        t_junctions
        
        I
        original_I 
        
        segments 
        
        boundaries % boundaries of all segments (including internal)
        boundaries_linesegs        
        lineseg_angles
        symmetry_axes
        junc_lut
        
        adj_contours % contours shared by adjacent segments
        adj_contours_ori % local orientation of adjacent contours
        adj_contour_angles
        adj_contour_conn_comp   
        adj_contour_ori_hists
        adj_ls_ids
        others_ls_ids
        
        nb_dgraph % segment neighborhood
        nb_hgraph
        nb_ids         
        
        Intersection
        
        initialized
    end
    
    methods
        function obj = SegmFeatures(I, segments)
            obj.I = I;
            obj.segments = segments;   
            obj.features = obj.create_features_struct();
        end                
       
        function obj = initialize(obj, just_nb, width_to_grow)
            DefaultVal('*width_to_grow', '1');
            %%%% neighborhood and adjacent contours %%%
            if(isempty(obj.nb_ids))
                [obj.nb_hgraph, obj.Intersection] = obj.compute_neigborhood(obj.segments, width_to_grow);
                obj.nb_dgraph = obj.nb_hgraph.pairwise_to_dgraph();
                nb = triu(obj.nb_dgraph.D);
                [nb_ids_i, nb_ids_j] = find(nb);
                obj.nb_ids = [nb_ids_i nb_ids_j];                                    
            end

            DefaultVal('*just_nb', 'false');
            if(just_nb)
                return;
            end                       

            if(isempty(obj.nb_ids))
                return;
            end
            
            obj.boundaries = obj.compute_boundaries(obj.segments);            
            BREAKING_FACTOR = 0.9;     % 0.8
            obj.boundaries_linesegs = obj.linearize_boundaries(BREAKING_FACTOR);
            obj.lineseg_angles = obj.calc_bndr_linesegs_angles();

            obj.adj_contours =  obj.get_adj_contours();
            [obj.adj_ls_ids, obj.others_ls_ids] = obj.get_adj_ls_v2();   

            %obj.symmetry_axes = obj.get_medial_axes();
            obj.symmetry_axes = obj.get_second_moments_axes();
            % always consider reflections across vertical axis
            obj.symmetry_axes = [obj.symmetry_axes (pi/2)*ones(size(obj.symmetry_axes,1),1)];
            obj.junc_lut = makelut(@obj.is_junction, 3);
        end
        
        %%%%%%%%%%%%%%%%%%%%% feature computation %%%%%%%%%%%%%%%%%%%%

        function obj = compute_all_feats(obj)
            if(~obj.initialized)
                error('Run initialize first!');
            end
             
            img_area = (size(obj.I,1) * size(obj.I,2));
            img_diagonal = sqrt((size(obj.I,1)^2 + (size(obj.I,2)^2)));
            n_segms = size(obj.segments,3);
            nb1 = obj.nb_ids(:, 1);     
            nb2 = obj.nb_ids(:, 2);                      
        end
              
        function pairwise_feats = convert_unary_to_pairwise_feats(obj, un_Feats)
            pairwise_feats = un_Feats(:, obj.nb_ids(:,1)) - un_Feats(:, obj.nb_ids(:,2));
        end
        
        function feats = get_feats(obj)
            n_feats = numel(fields(obj.features));
            n_segms = size(obj.segments,3);
            
            feats = zeros(n_segms, n_segms, n_feats);
            feat_names = fieldnames(obj.features);
            for i=1:n_feats
                feats(:,:,i) = obj.features.(feat_names{i});
            end
            
            feats = single(feats);
        end
        
        %%%% Line segment processing %%%%
        
        function len = calc_len_linesegs(obj, linesegs)
            len = vnorm(linesegs(1:2,:) - linesegs(3:4,:));
        end
        
        function angles = calc_bndr_linesegs_angles(obj)
            n = numel(obj.boundaries_linesegs);
            angles = cell(n,1);
            
            for i=1:n
                ls = obj.boundaries_linesegs{i};
                angles{i} = obj.calc_linesegs_angles(ls);
            end
        end
        
        function angles = calc_linesegs_angles(obj, linesegs)
            if(~isempty(linesegs))
                v = linesegs(1:2,:) - linesegs(3:4,:);

                %                      v1n = realsqrt(v1'*v1);
                %                      v1_norm = v1./v1n; % set length = 1
                %                      if v1_norm(2) < 0  % arc cos is unique if ycoord of vector is > 0
                %                        v1_norm = -v1_norm;
                %                      end
                %                     angle1 = acos(v1_norm'*[1 0]');
                angles = atan(v(1,:)./v(2,:));
            else
                angles = [];
            end
        end
        
        function [centers] = centers_of_line_segments(obj, line_segments)
            Pos1 = line_segments(1:2,:);
            Pos2 = line_segments(3:4,:);
            
            V = Pos2 - Pos1;
            
            centers = Pos1 + V/2;
        end
                
        function [ref_line_segments] = reflect_linesegments(obj, line_segments, angle)
        %function [ref_line_segments] = reflect_linesegments(obj, line_segments, line_direction)
            %vec_dir = line_direction(3:4) - line_direction(1:2);
            %lx = vec_dir(1);
            %ly = vec_dir(2);
            %A = (1/((lx*lx) + (ly*ly))) * [(lx*lx) 2*lx*ly; 2*lx*ly (ly*ly - lx*lx)];
            
            A = [cos(2*angle) sin(2*angle); sin(2*angle) -cos(2*angle)];
            ref_line_segments(1:2,:) = A*line_segments(1:2,:);
            ref_line_segments(3:4,:) = A*line_segments(3:4,:);            
        end
        
        function [linesegs_origin_reg] = register_translation_linesegs(obj, linesegs_origin, linesegs_dest)        
            if(isempty(linesegs_origin))
                linesegs_origin_reg = [];
                return;
            end
            
            ep_dest = [linesegs_dest(1:2,:) linesegs_dest(3:4,:)];
            ep_origin = [linesegs_origin(1:2,:) linesegs_origin(3:4,:)];
            
            % 1. compute center of linesegs_dest            
            ep_dest_un = unique(ep_dest', 'rows')';
            center_dest = mean(ep_dest_un,2);    
            
            % 2. compute center of linesegs_origin
            ep_origin_un = unique(ep_origin', 'rows')';
            center_origin = mean(ep_origin_un,2);
            
            % 3. translate linesegs_origin
            ep_origin = ep_origin + repmat((center_dest - center_origin), 1, size(ep_origin,2));
            
            linesegs_origin_reg = [ep_origin(:, 1:size(linesegs_origin,2)); ...
                                                  ep_origin(:, (size(linesegs_origin,2)+1):2*size(linesegs_origin,2))];                                              
        end
              
        function [intersection_point] = find_likely_line_intersection(obj, linesegs1, linesegs2)
            % Finds the most likely intersection point between line
            % segments in linesegs1 and line segments in linesegs2,
            % considering they are lines ( extend to infinity )
            % The elements in each set are assumed to be have high
            % intra-cluster similarity.
            CONST = 10000000000;
            
            dir = linesegs1(3:4,: ) - linesegs1(1:2,:);
            linesegs1(1:2,:) = linesegs1(1:2,:) + CONST * dir;
            linesegs1(3:4,:) = linesegs1(3:4,:) - CONST * dir;
            
            dir = linesegs2(3:4,: ) - linesegs2(1:2,:);
            linesegs2(1:2,:) = linesegs2(1:2,:) + CONST * dir;
            linesegs2(3:4,:) = linesegs2(3:4,:) - CONST * dir;
              
            counter = 1;
            p = zeros(2, size(linesegs1,2) * size(linesegs2,2));
            for i=1:size(linesegs1, 2)                
                for j=1:size(linesegs2,2)
                    % extend line segment to far away                    
                    p(:, counter) = cp2seg(linesegs1(:,i), linesegs2(:,j));
                   
                    
                    %%% visualization %%%
%                     sc(obj.I); hold on;
%                     jl_plot_lines([linesegs1(:,i) linesegs2(:,j)]);
%                     plot(p(2,counter), p(1,counter), 'og');

                    counter = counter + 1;
                end
            end
            intersection_point = mean(p,2);
        end
        
        function [angles, new_junction_centers, leg_region] = get_junction_angles(obj, neighborhood_id, junction_centers)  
            angles = zeros(2, size(junction_centers,2));
            
            segm_ids = obj.nb_ids(neighborhood_id,:);
            segm_id1 = segm_ids(1);
            segm_id2 = segm_ids(2);
            
            angles1 = obj.lineseg_angles{segm_id1};
            angles2 = obj.lineseg_angles{segm_id2};
            
            b_ls1 = obj.boundaries_linesegs{segm_id1};
            b_ls2 = obj.boundaries_linesegs{segm_id2};

            adj_ls_ids_1 = obj.adj_ls_ids{neighborhood_id,1};
            adj_ls_ids_2 = obj.adj_ls_ids{neighborhood_id,2};
            
            all_1 = 1:size(b_ls1,2);
            all_2 = 1:size(b_ls2,2);
                             
            new_junction_centers = zeros(size(junction_centers));
            leg_region = zeros(size(junction_centers,2),1);
            for i=1:size(junction_centers,2)
                d1 = zeros(size(junction_centers,2), size(b_ls1,2));
                for j=1:size(b_ls1,2)
                    d1(i,j) = dist_p2line(b_ls1(:,j), junction_centers(:,i), 'lineseg');
                end
                
                d2 = zeros(size(junction_centers,2), size(b_ls2,2));
                for j=1:size(b_ls2,2)
                    d2(i,j) = dist_p2line(b_ls2(:,j), junction_centers(:,i), 'lineseg');
                end
            
                m_d1 = min(d1);
                m_d2 = min(d2);
                MAX_DIST_LINESEG_FROM_JUNCTION = 4;
                j_ls1 = find(max(MAX_DIST_LINESEG_FROM_JUNCTION, m_d1) > d1(i,:));
                j_ls2 = find(max(MAX_DIST_LINESEG_FROM_JUNCTION, m_d2) > d2(i,:));
                
                if(isempty(j_ls1) || isempty(j_ls2))
                    angles(:,i) = [0;0];
                    continue;
                end
                
                %%% visualization %%%
%                 Img = create_transparent_multiple_colors(obj.segments(:,:,[segm_id1 segm_id2]), obj.I);             
%                 sc(Img);
%                 jl_plot_lines(b_ls1(:, j_ls1));
%                 jl_plot_lines(b_ls2(:, j_ls2));
%                 hold on;
%                 plot(junction_centers(2, i), junction_centers(1,i), 'og', 'Linewidth', 5);
                
                
                %%% find the two dominant angles   
                all_angles = [angles1(j_ls1) angles2(j_ls2)];
                angle_diff = zeros(numel(all_angles), numel(all_angles));
                for m=1:numel(all_angles)
                    for n=1:numel(all_angles)
                        [angle_diff(m,n)] = obj.compute_angle_diffs(all_angles(m), all_angles(n));
                    end
                end
                if(0.2 > max(max(angle_diff)))
                    continue;
                end
                
                Z = linkage(squareform(angle_diff,'tovector'));
                cluster_ids = cluster(Z,'MaxClust',2);
                %[medoid_ids,quality, cluster_ids] = kmedoid(angle_diff, 2, 100);
                                
                a1 = [angles1(j_ls1(cluster_ids(1:numel(j_ls1)) ==1))  angles2(j_ls2(cluster_ids(end-numel(j_ls2)+1:end) == 1))];
                a2 = [angles1(j_ls1(cluster_ids(1:numel(j_ls1)) ==2))  angles2(j_ls2(cluster_ids(end-numel(j_ls2)+1:end) == 2))];
                angle_mean_a1 = obj.compute_angle_mean(a1);
                angle_mean_a2 = obj.compute_angle_mean(a2);
                                
                ls_c1 = [b_ls1(:, j_ls1(cluster_ids(1:numel(j_ls1)) ==1)) b_ls2(:, j_ls2(cluster_ids(end-numel(j_ls2)+1:end) == 1))];                   
                ls_c2 = [b_ls1(:, j_ls1(cluster_ids(1:numel(j_ls1)) ==2)) b_ls2(:, j_ls2(cluster_ids(end-numel(j_ls2)+1:end) == 2))];
                
                [new_junction_centers(:,i)] = obj.find_likely_line_intersection(ls_c1, ls_c2);                
                                  
                %
                % decide which one is the hat and which one is the leg
                %                 
                 
                direction1 = [sin(angle_mean_a1); cos(angle_mean_a1)];  
                direction2 = [sin(angle_mean_a2); cos(angle_mean_a2)];  
                d1_dir1 = zeros(1,size(ls_c1,2));
                d2_dir1 = d1_dir1;
                d1_dir2 = d1_dir1;
                d2_dir2 = d1_dir1;
                
                % get endpoints farther from junction
                counter = 1;
                for j=1:size(ls_c1,2)                    
                    v1 = (ls_c1(1:2,j) - new_junction_centers(:,i));
                    v2 = (ls_c1(3:4,j) - new_junction_centers(:,i));
                    d1_dir1(counter:counter+1) = [(direction1' * v1) (direction1' * v2)];        
                    d1_dir2(counter:counter+1) = [(direction2' * v1) (direction2' * v2)];        
                    counter = counter + 2;                    
                    %%% visualization %%%
                    %figure(2);
                    %hold on;
                    %plot(new_junction_centers(2, 1), new_junction_centers(1,1), 'or', 'Linewidth', 5);
                    %jl_plot_lines([direction1; 0; 0]);
                    %jl_plot_lines([v1; 0 ; 0]);
                    %jl_plot_lines([v2; 0; 0]);
                end
                counter = 1;
                for j=1:size(ls_c2,2)
                    v3 = (ls_c2(1:2,j) - new_junction_centers(:,i));
                    v4 = (ls_c2(3:4,j) - new_junction_centers(:,i));
                    d2_dir1(counter:counter+1) = [(direction1' * v3) (direction1' * v4)];   
                    d2_dir2(counter:counter+1) = [(direction2' * v3) (direction2' * v4)];        
                    counter = counter + 2;
                end        
                
                val1_max_dir1 = abs(max(d1_dir1));
                val1_min_dir1 = abs(min(d1_dir1));
                %val2_max_dir1 = abs(max(d2_dir1));
                %val2_min_dir1 = abs(min(d2_dir1));                                                                
                
                %val1_max_dir2 = abs(max(d1_dir2));
                %val1_min_dir2 = abs(min(d1_dir2));
                val2_max_dir2 = abs(max(d2_dir2));
                val2_min_dir2 = abs(min(d2_dir2));                                                                

                %
                % get angle of leg and base of T
                %
                id_leg = 2;
                id_hat = 1;
                
                leg_in_cluster_two = min([val1_max_dir1 val1_min_dir1])>min([val2_max_dir2 val2_min_dir2]);
                if(leg_in_cluster_two)                                        
                    angles(id_hat,i) = angle_mean_a1; 
                    angles(id_leg,i) = angle_mean_a2;                                                   
                    if(val2_min_dir2>val2_max_dir2)
                    %if(val2_max<=0)
                       angles(id_leg,i) = angles(id_leg,i) + pi;
                    end
                else
                    angles(id_hat,i) = angle_mean_a2;
                    angles(id_leg,i) = angle_mean_a1;                                         
                    if(val1_min_dir1>val1_max_dir1)
                    %if(val1_max<=0)
                        angles(id_leg,i)  = angles(id_leg,i)  + pi;
                    end
                end      
                
                %%%                                                                                        %%%
                %%% Decide to which region the leg belongs, or if to both %%%
                %%%                                                                                        %%%

                % compute cross product between
                %junc = new_junction_centers(:,i);                
                %cross_res = cross([x(end) y(end) 0] - [x(1) y(1) 0], [circle_center(1,k)  circle_center(2,k) 0]  - [x(1) y(1) 0]);
                %disp('estou aqui');
                
                
                one_in_cluster_one = (cluster_ids(1:numel(j_ls1)) ==1);
                one_in_cluster_two = (cluster_ids(1:numel(j_ls1)) ==2);
                two_in_cluster_one = (cluster_ids(end-numel(j_ls2)+1:end) == 1);
                two_in_cluster_two = (cluster_ids(end-numel(j_ls2)+1:end) == 2);
                
                if((leg_in_cluster_two && any(two_in_cluster_two) && ~any(one_in_cluster_two)) || (~leg_in_cluster_two && any(two_in_cluster_one) && ~any(one_in_cluster_one)))
                    leg_region(i) = 2;
                elseif((leg_in_cluster_two && any(one_in_cluster_two) && ~any(two_in_cluster_two)) || (~leg_in_cluster_two && any(one_in_cluster_one) && ~any(two_in_cluster_one)))
                    leg_region(i) = 1;
                else
                    leg_region(i) = 0; % shared
                end
            end                            
        end
        
        function dist_mat = dist_linesegs(l1, l2)
            dist_mat = zeros(size(l1,2), size(l2,2));
            for m=1:size(l1,2)
                for n=1:size(l2,2)
                    %jl_plot_lines(l1(:,m));
                    %jl_plot_lines(l2(:,n));
                    
                    dist_mat(m,n) = cp2seg(l1(:,m) , l2(:,n));
                    
                    %%% visualization %%%
                    %sc(obj.I);
                    %jl_plot_lines(l1(:,m));
                    %jl_plot_lines(l2(:,n));
                    %dist_mat(m,n)
                end
            end
        end
        
        function [dist_mat] = dist_lineseg_1ep(obj, l1, l2)       
            % measure it as the min distance between l1 and some endpoint of l2            
            dist_mat = zeros(size(l1,2), size(l2,2));
            for m=1:size(l1,2)
                for n=1:size(l2,2)                    
                    %jl_plot_lines(l1(:,m));
                    %jl_plot_lines(l2(:,n));        
                    
                    d1 = dist_p2line(l1(:,m), l2(1:2,n), 'lineseg');
                    d2 = dist_p2line(l1(:,m), l2(3:4,n), 'lineseg');
                    dist_mat(m,n) = min(d1,d2);
                                        
                    %%% visualization %%%
                    %sc(obj.I);
                    %jl_plot_lines(l1(:,m));
                    %jl_plot_lines(l2(:,n));                    
                    %dist_mat(m,n)
                end
            end
        end
        
        function [dist_mat] = dist_lineseg_2ep(obj, l1, l2)
             % measure it as the sum of distances between l1 and both endpoints of l2
            dist_mat = zeros(size(l1,2), size(l2,2));
            for m=1:size(l1,2)
                for n=1:size(l2,2)
                    dist_mat(m,n) = dist_p2line(l1(:,m), l2(1:2,n), 'lineseg') + dist_p2line(l1(:,m), l2(3:4,n), 'lineseg');
                end
            end
        end
                
        %%%% Line segments and pixel contour %%%%
                
        function [adj_ls, other_ls] = get_adj_ls_v2(obj)
            n_interfaces = size(obj.adj_contours,1);
            assert(n_interfaces ~= 0);
            adj_ls = cell(n_interfaces,2);
            other_ls = adj_ls;
                                   
            for i=1:n_interfaces
                linesegs1 = obj.boundaries_linesegs{obj.nb_ids(i,1)};
                adj_cont = obj.adj_contours{i,1}';   
                [adj_ls{i,1}, other_ls{i,1}]= obj.assign_linesegs_to_adj_contour(linesegs1, adj_cont);
      
                linesegs2 = obj.boundaries_linesegs{obj.nb_ids(i,2)};
                adj_cont = obj.adj_contours{i,2}';   
                [adj_ls{i,2}, other_ls{i,2}]= obj.assign_linesegs_to_adj_contour(linesegs2, adj_cont);
                
                %%% visualization
                %sc(obj.I);
                %jl_plot_lines(linesegs1(:,adj_ls{i,1}));
                %jl_plot_lines(linesegs2(:,adj_ls{i,2}));
            end
        end
        
        function [adj_ls, other_ls] = assign_linesegs_to_adj_contour(obj, linesegs, adj_contour)           
            if(isempty(adj_contour) || isempty(linesegs) || (size(linesegs,1)~=4))
                adj_ls = [];
                other_ls = [];
                return;
            end
                        
            center_ls= round(obj.centers_of_line_segments(linesegs));            
            ep = [linesegs(1:2,:) linesegs(3:4,:) center_ls];            
            sz = [size(obj.I,1) size(obj.I,2)];
            ep = sub2ind(sz, ep(1,:), ep(2,:));                        
            adj_cont = sub2ind(sz,adj_contour(1,:), adj_contour(2,:));
            
            [inter, id] = intersect(ep,adj_cont);
            map = repmat(1:size(linesegs,2), 1, 3);
            [adj_ls, un] = unique(map(id));
            id = id(un);
            
            % for each adj line segment
            % check that both of its endpoints are close enough to end points 
            % that are known for sure to be in the adj_contour
            to_remove = false(numel(adj_ls), 1);
            for i=1:numel(adj_ls)   
                [x,y] = ind2sub([size(obj.I,1) size(obj.I,2)], inter');
                sure_ep = [x y];
                
                this_ls = linesegs(:,adj_ls(i));

                d1 = min(pdist2(this_ls(1:2)', [x y]));
                d2 = min(pdist2(this_ls(3:4)', [x y]));
                
                to_remove(i) = ((d1+d2)>100);
            end
            
            if(numel(adj_ls)>1)
                adj_ls(to_remove) = [];
            end
                        
            other_ls = setdiff(1:size(linesegs,2),adj_ls);
            
            % visualization
            %jl_plot_lines(linesegs); hold on;            
            %jl_plot_lines(linesegs(:,adj_ls));
            %plot(adj_contour(2,:), adj_contour(1,:), 'og');
        end
                     

        
        %%%%%%%%%%%%%%%%% Visualization %%%%%%%%%%%%%%%%%%

        function show_contour(obj, the_contour)
            sz = [size(obj.I,1) size(obj.I,2)];
            img = zeros(sz);
            img(sub2ind(sz, the_contour(:,1), the_contour(:,2))) = 1;
            sc(img);
        end

        function show_junction_wedges(obj, segms, t_junctions)  
          bg = zeros([size(segms,1) size(segms,2)]);
          new_bg = cat(3, bg, bg);
          new_bg = cat(3, new_bg, bg);          
          
          for i=1:size(segms,3)
              new_segm{i} = cat(3, segms(:,:,i), segms(:,:,i));
              new_segm{i} = cat(3, new_segm{i}, segms(:,:,i));
          end
          
          new_bg(:,:,1) = 0.00;
          new_bg(:,:,2) = 0.00;
          new_bg(:,:,3) = 0.00;

          mults1 = [1 0 0];
          mults2 = [0 0 1];
          for i=1:numel(new_segm)
              if(i==1)
                  mults = mults1;
              elseif(i==2)
                  mults = mults2;
              else
                  mults = rand(3,1);
              end
            new_bg(:,:,1) = new_bg(:,:,1) + new_segm{i}(:,:,1).*mults(1); % .* rand();              
            new_bg(:,:,2) = new_bg(:,:,2) + new_segm{i}(:,:,2).*mults(2); % .* rand();              
            new_bg(:,:,3) = new_bg(:,:,3) + new_segm{i}(:,:,3).*mults(3); % .* rand();              
          end
          new_bg = uint8(255*new_bg);

          alpha_chann = 0.7*ones(size(obj.I,1), size(obj.I,2));   

          Img = immerge(obj.I, new_bg, alpha_chann);
                                                    
          imshow(Img);
          hold on;
          %subplot(1,2,2), plot(cj, ci, 'og', 'Linewidth', 2);
          obj.plot_t_junctions(obj.I, t_junctions, false);
        end

        function show_t_junctions(obj, t_junctions)
            sc(obj.I);
            if(numel(t_junctions)>1)                
                the_t_junctions = struct('coords', [], 'angles', [], 'leg_region', []);
                for i=1:numel(t_junctions)
                    if(~isempty(t_junctions(i).coords))
                        the_t_junctions.coords = [the_t_junctions.coords t_junctions(i).coords];
                        the_t_junctions.angles = [the_t_junctions.angles t_junctions(i).angles];
                        the_t_junctions.leg_region = [the_t_junctions.leg_region; t_junctions(i).leg_region];
                    end
                end
                t_junctions = the_t_junctions;
            end
            obj.plot_t_junctions(obj.I, t_junctions);
        end           
        
        function plot_t_junctions(obj, I, t_junctions, show_image)
            DefaultVal('show_image', 'true');
            
            ci = t_junctions.coords(1,:);
            cj = t_junctions.coords(2,:);
            LEN = 10;
            if(show_image)
                sc(I); hold on;
            end
            for i=1:size(t_junctions.angles,2)
                % plot top of T-junction
                x = [sin( t_junctions.angles(1,i)) cos( t_junctions.angles(1, i )) ];
                top1 = [ci(i) cj(i)] + LEN*x;
                top2 = [ci(i) cj(i)] - LEN*x;
                line([cj(i); top2(2); top1(2)], [ci(i) top2(1) top1(1)],  'Linewidth', 2, 'Color', 'green');
                hold on;
                
                % plot leg of T-junction
                %x = [sin(pi - t_junctions.angles(2, i)) cos(pi - t_junctions.angles(2, i)) ];
                x = [sin(t_junctions.angles(2, i)) cos(t_junctions.angles(2, i)) ];
                top1 = [ci(i) cj(i)] + LEN*x;
                line([cj(i); top1(2)], [ci(i); top1(1)],  'Linewidth', 2, 'Color', 'green');
                hold on;
            end
        end
        
        function show_nb_val(obj, vals1, vals2)
            for i=1:size(obj.nb_ids,1)
                trI = create_transparent_multiple_colors(obj.segments(:,:,obj.nb_ids(i,:)), obj.I);
                sc(trI);
                SvmSegm_plot_labels_at_segms(obj.segments(:,:,obj.nb_ids(i,:)), {num2str(vals1(i)), num2str(vals2(i))});
                pause;
            end
        end   
        
        function visualize_area_theta_circles(obj, relative_area, theta, dist, color)
            DefaultVal('dist', '30');
            DefaultVal('color', '''g'''); % green
            % for visualizing area-angle relationships between neighboring objects
            SZ = 25;
            RADIUS = dist*600;
            LW = 4;
            
            plot(0, 0, 'rO', 'MarkerSize', SZ, 'LineWidth', LW); hold on;
            
            x = RADIUS.*cos(theta);
            y = RADIUS.*sin(theta);
            marker_size = SZ/relative_area;
            m = [color 'O'];
            plot(x, y, m, 'MarkerSize', marker_size, 'LineWidth', LW);
            axis([-500 500 -500 500]);
        end
    
        function show_masks(obj)
            subplot_auto_transparent(obj.segments, obj.I);
        end
    end

    methods(Static)
          
        function [angle_diffs] = compute_angle_diffs(angles1, angles2)
            big = max(angles1, angles2);
            small = min(angles1, angles2);
            angle_diffs_1 =  big - small;
            angle_diffs_2 =  (pi/2 - big) - (-pi/2 - small);
            angle_diffs = min(angle_diffs_1, angle_diffs_2);
        end
        
        function [angle_mean] = compute_angle_mean(angles)
            n = numel(angles);
            if(n == 1)
                angle_mean = angles;
                return;
            end
            new_angles = angles*2;
            new_angle_mean = atan2((1/n) * sum(sin(new_angles)), (1/n) * sum(cos(new_angles)));
            %angle_mean = atan2((1/n) * sum(cos(angles)), (1/n) * sum(sin(angles))); 
            angle_mean = (new_angle_mean/2);
        end
               
        function features = create_features_struct()  
            features = [];
        end
        
        function print_feat_vals(vals)
            feats = obj.create_features_struct();
            feat_names = fieldnames(feats);
            assert(size(vals,1) == numel(feat_names)); 
            for i = 1:numel(feat_names)
               fprintf(['%s: %f\n'], feat_names{i}, vals(i));
            end            
        end
    end
    
    methods
        %%%%%%%%%%%% neighborhood related stuff %%%%%%%%%%%%%%%
        %%%%%%%%%%%% pixel masks version %%%%%%%%%%%%%%%%%%%%
        function [nb_hgraph, Intersection] = compute_neigborhood(obj, segments, width_to_grow)    
            
            % resizing to speed it up
            assert(size(segments,3)>= 1);
            
            res_segments = imresize(segments, 0.4, 'nearest');            
            
            res_boundaries = obj.compute_boundaries(res_segments);            
            bndr_grown = obj.to_bndr_masks(res_boundaries, res_segments, width_to_grow);
            bndr_inters = segm_intersection_mex(reshape(bndr_grown, size(bndr_grown,1)*size(bndr_grown,2), size(bndr_grown,3)));
            bndr_inters = triu(bndr_inters, 1);
                        
            inters = segm_intersection_mex(reshape(res_segments, size(res_segments,1)*size(res_segments,2), size(res_segments,3)));
            Intersection = triu(inters,1);
            
            CONST = 10;
            CONST_BNDR = 20;
            aff_mat = bndr_inters;
            aff_mat(bndr_inters<CONST_BNDR) = 0;
            aff_mat(Intersection>CONST) = 0;
            
            nb_hgraph = obj.create_nn_hgraph(aff_mat);
        end        
        
        function [nb2_ids, intermediary] = compute_2nd_level_neighborhood(obj, above_2nd)
            DefaultVal('*above_2nd', 'false');
            un = unique(obj.nb_ids);
            n_masks = size(obj.segments,3);
            D = sparse2(obj.nb_ids(:,1), obj.nb_ids(:,2), ones(size(obj.nb_ids,1),1), n_masks, n_masks);
            D = Dgraph(D);
            nb2_ids = [];
            intermediary = zeros(numel(un),1);
            
            for i=1:numel(un)
                d = D.get_distances_to_node(un(i));
                
                if(above_2nd)
                    nb2 = find(d>=2);
                else
                    nb2 = find(d==2);
                end
                
                nb2_ids = [nb2_ids; un(i) * ones(numel(nb2),1) nb2];
            end
            
            nb2_ids = unique(sort(nb2_ids,2), 'rows');
            
            % remove overlapping regions
            nb2_ids = nb2_ids(~logical(obj.Intersection(sub2ind([n_masks n_masks], nb2_ids(:,1), nb2_ids(:,2)))), :);
            intermediary = [];
        end
        
        function adj_contours = get_adj_contours(obj)
            bndr_masks = obj.to_bndr_masks(obj.boundaries, obj.segments);
            all_ids = obj.nb_ids;
            %bndr = zeros(size(obj.segments));
            %parfor i=1:size(obj.segments,3)
            %    bndr(:,:,i) = imresize(bndr_grown(:,:,i), [size(obj.segments,1) size(obj.segments,2)], 'nearest');
            %end
            
            HOW_MUCH_TO_GROW = 5;
            bndr_grown = obj.grow_masks(bndr_masks, HOW_MUCH_TO_GROW);

            fat_adj_contours = get_fat_adj_contours(obj, bndr_grown);
                        
            segm_ids = all_ids(:,1);
            adj_contours1 = obj.assign_fat_contours(fat_adj_contours, segm_ids);
            segm_ids = all_ids(:,2);
            adj_contours2 = obj.assign_fat_contours(fat_adj_contours, segm_ids);
            adj_contours = [adj_contours1 adj_contours2];
        end

        function adj_contours = assign_fat_contours(obj, fat_adj_contours, segm_ids)
            adj_contours = cell(size(fat_adj_contours,1),1);
            sz = [size(obj.I,1) size(obj.I,2)];
            %parfor i=1:numel(fat_adj_contours)
            for i=1:numel(fat_adj_contours)
                % intersect real boundaries with these fat contours
                for j=1:size(obj.boundaries{segm_ids(i)})                   
                    bnd = sub2ind(sz, obj.boundaries{segm_ids(i)}{j}(:,1),  obj.boundaries{segm_ids(i)}{j}(:,2));
                    fat = sub2ind(sz, fat_adj_contours{i}(:,1), fat_adj_contours{i}(:,2));
                    [duh, a] = intersect(bnd, fat);
                    
                    a = sort(a);
                    adj_contours{i} = [adj_contours{i}; obj.boundaries{segm_ids(i)}{j}(a,:)];
                    %diff(adj_contours{i})
                end
                
                % visualization
%                 subplot(1,3,1), sc(obj.segments(:,:, obj.nb_hgraph.ids{1}(i,2))*2 + obj.segments(:,:,obj.nb_hgraph.ids{1}(i,1)));
%                 duh = zeros(size(obj.I,1), size(obj.I,2));
%                 duh(sub2ind(size(duh), adj_contours{i}(:,1), adj_contours{i}(:,2))) = 1;
%                 subplot(1,3,2), sc(duh);
%                 duh_bndr = zeros(size(obj.I,1), size(obj.I,2));
%                 duh_bndr(sub2ind(size(duh), obj.boundaries{segm_ids(i)}{1}(:,1), obj.boundaries{segm_ids(i)}{1}(:,2))) = 1;
%                 subplot(1,3,3), sc(duh_bndr);
            end
        end
           
        function fat_adj_contours = get_fat_adj_contours(obj, bndr_grown)            
            ids = obj.nb_ids;
            n_nb = size(ids,1);            
            %fat_adj_contours = zeros(size(bndr_grown,1), size(bndr_grown,2), n_nb);
            
            fat_adj_contours = cell(n_nb,1);    
            %parfor i=1:n_nb                            
            for i=1:n_nb                            
                %fat_adj_contours(:,:,i) = (bndr_grown(:,:,ids(i,1)) + bndr_grown(:,:,ids(i,2))) == 2;
                 [x,y] = find((bndr_grown(:,:,ids(i,1)) + bndr_grown(:,:,ids(i,2))) == 2);
                 fat_adj_contours{i} = [x y];
            end
        end        
        
        function the_endpoints = get_endpoints(obj, x,y)
            % find endpoints
            
            if(size(x,1) < 3)
                the_endpoints = [x y];
                return;
            end
                
            d = pdist2([x y], [x y], 'euclidean');
            endpoints = false(numel(x), 1);            
            for i=1:numel(x)
                di = sort(d(i,:));
                if((di(3) ~= 1) && (di(3) ~= sqrt(2)))
                    endpoints(i) = true;
                end                
            end        
            if(~any(endpoints))
                endpoints(1) = true;
            end
            the_endpoints = [x(endpoints) y(endpoints)];
        end
        
        function bndr_masks = to_bndr_masks(obj, bndr, masks, width_to_grow)
            DefaultVal('width_to_grow', '1');
            % coordinates of the image frame
            nrows = size(masks,1);
            ncols = size(masks,2);

            img_frame = obj.generate_img_frame(nrows, ncols);
            bndr_masks = false(size(masks));
            
            for i=1:numel(bndr)
                tmp = bndr_masks(:,:,i);
                for j=1:numel(bndr{i})
                    tmp(sub2ind(size(tmp), bndr{i}{j}(:,1), bndr{i}{j}(:,2))) = true;
                end
                tmp(img_frame) = false;
                bndr_masks(:,:,i)= tmp;
                %imshow(bndr_masks(:,:,i));
                %pause;
            end
            bndr_masks = obj.grow_masks(bndr_masks, width_to_grow);                        
        end
            
        function img_frame = generate_img_frame(obj,nrows, ncols)
            img_frame = [ones(ncols,1) (1:ncols)'; nrows*ones(ncols,1) (1:ncols)'];
            img_frame = [img_frame; (1:nrows)' ones(nrows,1); (1:nrows)' ncols*ones(nrows,1)];
            img_frame = unique(sub2ind([nrows ncols], img_frame(:,1), img_frame(:,2)));
        end
        
        function grown_masks = grow_masks(obj, masks, width, area_proportional)
            DefaultVal('area_proportional', 'false');
            grown_masks = false(size(masks));
            if(area_proportional)
                areas = sum(reshape(masks, size(masks,1)*size(masks,2), size(masks,3)));
                parfor i=1:numel(areas)                    
                    se = strel('disk', round(width + max(0, 2*(log(areas(i)) - 8)))); % 2* log(areas) - 8 worked so so 
                    grown_masks(:,:,i) = imdilate(masks(:,:,i), se);
                    %sc(5*sc(obj.I) + sc(masks(:,:,i)) + sc(grown_masks(:,:,i)))
                end
            else
                se = strel('disk',width);
                parfor i=1:size(masks,3)
                    grown_masks(:,:,i) = imdilate(masks(:,:,i), se);
                end
            end            
        end
        
        function eroded_masks = erode_masks(obj, masks, width)
            se = strel('disk',width);
            eroded_masks = false(size(masks));
            %parfor i=1:size(masks,3)
            for i=1:size(masks,3)
                eroded_masks(:,:,i) = imerode(masks(:,:,i), se);
            end    
        end
        
        function bndr = compute_boundaries(obj, masks, exclude_img_frame)
            DefaultVal('*exclude_img_frame', 'false');
            % returns the boundary of each segment. If the segment has holes it will return more than one.
            % the larger boundaries are returned first.
            
            bndr = cell(size(masks,3),1);
            sz = [size(obj.I,1) size(obj.I,2)];
       
            parfor i=1:size(masks,3)
            %for i=1:size(masks,3)
                bndr{i} = bwboundaries(masks(:,:,i));
                
                n = cellfun(@numel, bndr{i});
                [n_sorted, ids] = sort(n, 'descend');
                bndr{i} = bndr{i}(ids);
                bndr{i} = {bndr{i}{n_sorted>6}}';
            end
            
            if(exclude_img_frame)
                nrows = size(masks,1);
                ncols = size(masks,2);
                img_frame = obj.generate_img_frame(nrows, ncols);
                for i=1:size(masks,3)
                    % assuming only the first connected component touches
                    % the image frame. The others are holes, i think.

                    if(~isempty(bndr{i}))
                        to_keep = setdiff(sub2ind([size(masks,1) size(masks,2)], bndr{i}{1}(:,1), bndr{i}{1}(:,2)), img_frame);

                        bndr{i}{1} = [];
                        [bndr{i}{1}(:,1), bndr{i}{1}(:,2)] = ind2sub([size(masks,1) size(masks,2)], to_keep);
                    end
                end                
            end
            
            %%%% visualization %%%%
%                 duh = zeros(size(obj.segments(:,:,1)));
%                 for j=1:numel(b)
%                     duh(sub2ind(size(duh), b{j}(:,1), b{j}(:,2))) = 1;
%                 end
%                 %end
%                 sc(duh);
        end        
        
        function bndr_linesegs = linearize_boundaries(obj, breaking_factor)               
            sz = [size(obj.I,1) size(obj.I,2)];
            bndr_linesegs = cell(numel(obj.boundaries),1);
                                  
            for i=1:numel(obj.boundaries)
            %for i=1:numel(obj.boundaries)
                the_lines_adj = [];
                for j=1:numel(obj.boundaries{i})
                    %diff(obj.boundaries{i}{j})
                    %[x,y] = ind2sub(sz, obj.boundaries{i}{j});
                    x = obj.boundaries{i}{j}(:,1);
                    y = obj.boundaries{i}{j}(:,2);
                    thresh = breaking_factor*0.05*log(numel(x))/log(1.1);
                    ln = lineseg({[x y]}, thresh);
                    the_lines_adj = [the_lines_adj seglist2segarray(ln)];
                end         
                the_lines_adj_filtered = obj.filter_truncation_linesegs(the_lines_adj);
                
                bndr_linesegs{i} = the_lines_adj_filtered;                
                
                %%% visualize
                %sc(obj.I); hold on;
                %sc(obj.segments(:,:,i));
                %jl_plot_lines(bndr_linesegs{i} );
            end
        end
        
        function linesegs = filter_truncation_linesegs(obj, linesegs)
            if(isempty(linesegs))
                return;
            end
            lower_coord = size(obj.I,1);
            right_coord = size(obj.I,2);
            linesegs_top = (linesegs(2,:) == 1 & linesegs(4,:) == 1);
            linesegs_bottom = (linesegs(1,:) == lower_coord & linesegs(3,:) == lower_coord);      
            linesegs_left =  (linesegs(1,:) == 1 & linesegs(3,:) == 1);
            linesegs_right = (linesegs(2,:) == right_coord & linesegs(4,:) == right_coord);
            linesegs(:, linesegs_top | linesegs_bottom | linesegs_left | linesegs_right) = [];
        end
        
        function g = create_nn_hgraph(obj, aff_mat)
            aff_mat = triu(aff_mat);
            [ids1, ids2] = find(aff_mat);
            g = Hgraph(size(aff_mat,1));
            val = ones(numel(ids1),1);
            g = g.add_edges([ids1 ids2], val);
        end
        
        function [coords_both_sides] = get_coords_on_both_sides_bndr(obj, nb_id, INTERVAL, SAMPLING_RATE, DIST_FROM_CONTOUR)
            bndr_length = size(obj.adj_contours{nb_id},1);
            n_intervals = numel(1:SAMPLING_RATE:bndr_length) - 1;
            coords_both_sides = zeros(4, n_intervals);
            for k=1:n_intervals
                range = round(max(1, (k*(SAMPLING_RATE) - (INTERVAL/2))):min(bndr_length,(k*(SAMPLING_RATE) + (INTERVAL/2))));
                x = obj.adj_contours{nb_id}(range,1);
                y = obj.adj_contours{nb_id}(range,2);
                
                dir = [x(end) y(end)] - [x(1) y(1)];
                orth_dir = [dir(2); -dir(1)];
                orth_dir = orth_dir/norm(orth_dir);
                center_point = [x(round(INTERVAL/2)); y(round(INTERVAL/2))];
                point_one_side = round(center_point + DIST_FROM_CONTOUR * orth_dir);
                point_other_side = round(center_point - DIST_FROM_CONTOUR * orth_dir);
                point_one_side = max(1, point_one_side);
                point_other_side = max(1, point_other_side);
                point_one_side = min([size(obj.I,1); size(obj.I,2)], point_one_side);
                point_other_side = min([size(obj.I,1); size(obj.I,2)], point_other_side);
                
                coords_both_sides(:,k) = [point_one_side; point_other_side];
                %
                %%% visualization
                %
                %sc(obj.segments(:,:,obj.nb_ids(nb_id,1))); hold on;
                %jl_plot_lines(lineseg_dir);
                %plot(point_one_side(2), point_one_side(1), 'og');
                %plot(point_other_side(2), point_other_side(1), 'or');                
            end
        end
                 
        function show_feature_scores(obj, feature_name)
            the_feats = obj.features.(feature_name);
            %scores = num2cell(the_feats)';
            scores = the_feats';
            obj.nb_ids
            
            for j=1:size(obj.nb_ids,1)
                if(scores(j)>0)
                    titles = {num2str(abs(scores(j))), ''};
                else
                    titles = {'', num2str(abs(scores(j)))};
                end

                figure(1);
                h = subplot_auto_transparent(obj.segments(:,:, obj.nb_ids(j,:)), obj.I, titles);                    
                text(-10, -10, feature_name);
                pause;
            end
        end        
    end
end
