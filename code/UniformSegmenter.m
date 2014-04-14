% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef UniformSegmenter < Segmenter
    % uniform unary potentials (it's the same value for all pixels, except   
    % seeds which have infinite )
    properties      
    end
    
    methods
        function obj = UniformSegmenter(I, img_name)
            obj = obj@Segmenter(I,img_name);

            obj.grid_dims = [6 6];
            obj.RECT_DIMS = [40 40]; % [25 25]
            obj.SEED_FRAME_WEIGHT = 1000; % 500
            obj.filter_segments = true;
            obj.DONT_BREAK_DISCONNECTED = false;            
            obj.upper_bp = 150; % 15, 15 is better than 100
            obj.resize_factor = 0.5;         
            obj.arrangement = 'grid'; 
        end
        
        function obj = set_params(obj, par)
            if(nargin==1)
                
                obj.params.CONSTANT = 1000;
                if(strcmp(obj.params.local_dist_type, 'color'))
                    if(isempty(obj.params.CONTRAST_SENSITIVE_SIGMA))
                        obj.params.CONTRAST_SENSITIVE_SIGMA = 0.07;
                    end
                    obj.params.PAIRWISE_OFFSET = 1/obj.params.CONSTANT;
                    %obj.RECT_DIMS = [10 10];
                elseif(strcmp(obj.params.local_dist_type, 'pb'))
                    %obj.params.CONTRAST_SENSITIVE_SIGMA = 0.05; %(0.2)
                    if(isempty(obj.params.CONTRAST_SENSITIVE_SIGMA))
                        obj.params.CONTRAST_SENSITIVE_SIGMA = 1.5; % 1.5
                    end
                    obj.params.PAIRWISE_OFFSET = 1/obj.params.CONSTANT; % 5
                else
                    obj.params.CONTRAST_SENSITIVE_SIGMA = 0.02;
                    obj.params.PAIRWISE_OFFSET = 1/obj.params.CONSTANT;
                end
                obj.params.CONTRAST_SENSITIVE_WEIGHT =1;
                obj.params.n_neighbors =  4;
            else
                obj.params = par;
            end            
        end

         
        function obj = add_hyp_conns(obj)
            [hyp_conns, types] = obj.generate_growth_hyps('internal');
            obj.P = obj.P.add_hyp_connections(types, hyp_conns);
            
            obj.RECT_DIMS = round(obj.RECT_DIMS/2); % these can be smaller
            [hyp_conns, types] = obj.generate_growth_hyps('external'); 
            obj.P = obj.P.add_hyp_connections(types, hyp_conns);   
            
            % just frame
            [hyp_conns, types] = obj.generate_subframes_growth_hyps();
            obj.P = obj.P.add_hyp_connections(types, hyp_conns);
        end
        
        % just call the superclass methods
        
        function [hyp_conns, types] = generate_growth_hyps(obj, internal_external)
            [hyp_conns, types] = generate_growth_hyps@Segmenter(obj, internal_external);
        end
        
        function [hyp_conns, types] = generate_subframes_growth_hyps(obj)
            obj.subframe_dims = [size(obj.I,1)-2 size(obj.I,2)-2];
            [hyp_conns, types] = generate_subframes_growth_hyps@Segmenter(obj);
        end
        
        %%%%%%%%%%%%%%% interactive %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         function [hyp_conns, types] = generate_interactive_hyp(obj)
            [hyp_conns, types] = obj.generate_interactive_growth_hyps();
            %[hyp_conns, types] = obj.generate_interactive_internal();
         end        
         
        function [hyp_conns, types] = generate_interactive_internal_growth_frame(obj)
            [hyp_conns, types] = generate_interactive_internal_growth_frame@Segmenter(obj);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    end
end
