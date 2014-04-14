% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef SimpleSegmentFeatures < SegmFeatures
    properties 
        pb_file
        gPb_thin
        resize_factor
    end
    
    methods
        function obj = SimpleSegmentFeatures(I, segments, pb_file)
            obj = obj@SegmFeatures(I,segments);
            obj.resize_factor = 0.5; % cut values are damaged if you compute them on a different resolution than the segments were computed!
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
            original_masks = obj.segments;
            n_masks = size(original_masks,3);

            % adapt neighborhood for the unary case!                                    
            obj.nb_ids = (1:n_masks)';
            obj.nb_ids = [obj.nb_ids obj.nb_ids+n_masks];                                                           
            
            resh_masks = reshape(original_masks, size(original_masks,1) * size(original_masks,2), size(original_masks,3));

            no_separate = true;
            S = SegmentProcessor([], resh_masks, obj.original_I, obj.pb_file,  0, no_separate, obj.resize_factor);
            %
            %%%% Cut features %%%%
            %
            if(isfield(obj.features, 'cut_values'))
                get_all = true;
                S = S.compute_energies(get_all);

                D_cut = zeros(8, size(original_masks,3));
                D_cut(1,:) = [S.energies(:).cut_ratio];
                D_cut(2,:) = [S.energies(:).cut];
                D_cut(3,:) = [S.energies(:).normalized_cut];
                D_cut(4,:) = [S.energies(:).unbalanced_normalized_cut];
                tmp_var = [S.energies(:).fraction_of_healthy_boundary];
                D_cut(5:end,:) = reshape(tmp_var, 4,size(tmp_var,2)/4);
                
                obj.features.cut_values = D_cut;
            end
            
            %
            %%%% Coarse Region Shape and Location features %%%%%%
            %
            
            if(isfield(obj.features, 'coarse_shape_location'))
                D_crsl = zeros(19, size(original_masks,3)); 

                % absolute quantities
                masks = reshape(S.segments, size(S.I,1), size(S.I,2), size(S.segments,2));
                for i=1:size(masks,3)
                    props = regionprops(masks(:,:,i), 'Area', 'Centroid', 'BoundingBox', 'MajorAxisLength', ...
                                                        'MinorAxisLength', 'Eccentricity', 'Orientation', 'ConvexArea', 'EulerNumber', 'EquivDiameter', 'Solidity', 'Extent', 'Perimeter');
                    if(~isempty(props))
                        props = props(1);  % might be more than one, describe the first
                        D_crsl(1:17,i) = struct2array(props(1)); 

                        % convexity
                        D_crsl(18, i) = props.Area / props.ConvexArea; % solidity is the same as this one ( repeated but kept for backwards compatibility )
                        % absolute distance to center of image
                        D_crsl(19, i) = norm([size(obj.I,1) size(obj.I,2)]./2 - props.Centroid);
                    end
                end    
                
                obj.features.coarse_shape_location = D_crsl;            
            end
                        
            % visualization
            %
%              feats = fieldnames(obj.features);
%              for i=1:numel(feats)
%                [val, id] = sort(obj.features.(feats{i})(1,:), 'descend');
%                subplot_auto_transparent(obj.segments(:,:,id), obj.I, val);
%              end
        end                                   
        
        %%%%%%%%%%%%%%%% Features %%%%%%%%%%%%%%%%%%%%%%
    end
    methods (Static)
        function features = create_features_struct()
            %features.cut_values = [];
            features.coarse_shape_location= [];
        end
        
        function print_feat_vals(vals)
            feats = SimpleSegmentFeatures.create_features_struct();
            feat_names = fieldnames(feats);
            assert(size(vals,1) == numel(feat_names));
            for i = 1:numel(feat_names)
                fprintf(['%s: %f\n'], feat_names{i}, vals(i));
            end
        end    
    end       
end
