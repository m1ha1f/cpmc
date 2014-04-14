% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef LongRangeSegmenter < Segmenter
    % unary potentials are different for different pixels
    properties     
        F
        D
    end
    
    methods
        function obj = LongRangeSegmenter(I, img_name)
            obj = obj@Segmenter(I,img_name);
            
            obj.RECT_DIMS = [15 15]; % 15
            obj.grid_dims = [6 6];

            obj.upper_bp = 300; % 400
            obj.SEED_FRAME_WEIGHT = 1000; % 1000 original
            obj.MAX_ENERGY = 250; % 250
            obj.filter_segments = true;
            obj.resize_factor = 0.5; %%%% important %%%%%%%%            
            obj.DONT_BREAK_DISCONNECTED = false;
                        
            % use FH segmentation algorithm to generate superpixels, then
            % place square seed in centers close to regular grid.
            obj.arrangement = 'felzenszwalb';     
            
            % use just a regular grid
            %obj.arrangement = 'grid';              
        end

        
        function obj = set_params(obj, par)
            if(nargin==1)
                
                % params
                obj.params.DEBUG_MODE = false;
                obj.params.CONSTANT = 1000; % 1000
                
                 if(strcmp(obj.params.local_dist_type, 'color'))
                    if(isempty(obj.params.CONTRAST_SENSITIVE_SIGMA))
                        obj.params.CONTRAST_SENSITIVE_SIGMA = 0.10; % 0.08
                        obj.MAX_ENERGY = 400;
                    end
                    obj.params.PAIRWISE_OFFSET = 2/obj.params.CONSTANT; % 2 /
                elseif(strcmp(obj.params.local_dist_type, 'pb'))
                    if(isempty(obj.params.CONTRAST_SENSITIVE_SIGMA))
                        obj.params.CONTRAST_SENSITIVE_SIGMA = 3;
                    end
                    %obj.MAX_ENERGY = 100;
                    obj.params.PAIRWISE_OFFSET = 4/obj.params.CONSTANT; % 6 /
                 else
                    obj.params.PAIRWISE_OFFSET = 2/obj.params.CONSTANT;
                    obj.params.CONTRAST_SENSITIVE_SIGMA = 0.05;
                 end
                 obj.params.CONTRAST_SENSITIVE_WEIGHT = 1; % 1
                obj.params.n_neighbors =  4;
            else
                obj.params = par;
            end
        end          
        
        function [hyp_conns, types] = generate_interactive_hyp(obj)
            [hyp_conns, types] = obj.generate_interactive_growth_hyps_unary();
            %[hyp_conns, types] = obj.generate_interactive_external_growth_hyps_frame_unary();
            %[hyp_conns, types] = obj.generate_interactive_frame();
            %[hyp_conns, types] = obj.generate_interactive_internal();
        end
        
        function obj = add_hyp_conns(obj)        
            [hyp_conns, types] = obj.generate_growth_hyps_unary('internal');
            obj.P = obj.P.add_hyp_connections(types, hyp_conns); 
            [hyp_conns, types] = obj.generate_growth_hyps_unary('external');            
            obj.P = obj.P.add_hyp_connections(types, hyp_conns); 
            %[hyp_conns, types] = obj.generate_subframes_growth_hyps_unary(false);            
            %obj.P = obj.P.add_hyp_connections(types, hyp_conns); 
        end
        
        % just call the superclass methods
        
        %function [hyp_conns, types] = generate_interactive_external_growth_hyps_unary(obj)
        %    [hyp_conns, types] = generate_interactive_external_growth_hyps_unary@Segmenter(obj);
        %end
        
        function [hyp_conns, types] = generate_growth_hyps_unary(obj, internal_external)
            [hyp_conns, types] = generate_growth_hyps_unary@Segmenter(obj, internal_external);
        end        
               
        function [hyp_conns, types] = generate_subframes_growth_hyps_unary(obj, iter)
            obj.subframe_dims = [size(obj.I,1)-2 size(obj.I,2)-2];
            [hyp_conns, types] = generate_subframes_growth_hyps_unary@Segmenter(obj, iter);
        end 
    end
end
