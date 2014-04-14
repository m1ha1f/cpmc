% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef GridOfFramesSegmenter < Segmenter
    % uniform unary potentials (it's the same value for all pixels, except   
    % seeds which have infinite )
    properties              
    end

    methods
        function obj = GridOfFramesSegmenter(I, img_name)
            obj = obj@Segmenter(I,img_name);            
            obj.subframe_dims = [150 150];
            obj.frame_grid_dims = [5 5];
            obj.RECT_DIMS = [5 5];
            obj.grid_dims = [1 1]; 
            obj.upper_bp = 150; % 100
            obj.filter_segments = true;
            obj.SEED_FRAME_WEIGHT = 100;
            obj.arrangement = 'windows';
            obj.resize_factor = 0.5; % 0.5
            obj.DONT_BREAK_DISCONNECTED = false;
            obj.MAX_ENERGY = 100; % 75
        end

        function obj = set_params(obj, par)
            if(nargin==1)                
                %if(~exist(obj, 'params'))
                %    obj.params = JustSegmGen.create_pars();
                %    obj.params.local_dist_type = 'pb';
                %end
                
                % params
                obj.params.DEBUG_MODE = false;
                obj.params.CONSTANT = 1000; % 1000
                
                if(strcmp(obj.params.local_dist_type, 'color'))
                    if(isempty(obj.params.CONTRAST_SENSITIVE_SIGMA))
                        obj.params.CONTRAST_SENSITIVE_SIGMA = 0.07; % 0.06                    
                    end
                    obj.params.PAIRWISE_OFFSET = 2/obj.params.CONSTANT; % 2 /     
                elseif(strcmp(obj.params.local_dist_type, 'pb'))
                    if(isempty(obj.params.CONTRAST_SENSITIVE_SIGMA))
                        obj.params.CONTRAST_SENSITIVE_SIGMA = 1.2;
                    end
                    obj.params.PAIRWISE_OFFSET = 2/obj.params.CONSTANT; % 2 /
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
           %[hyp_conns, types] = obj.generate_interactive_subframe(); % not implemented
           %[hyp_conns, types] = obj.generate_interactive_subframe_plus_internal();
           [hyp_conns, types] = obj.generate_interactive_subframe_plus_internal_unary();
        end

        function obj = add_hyp_conns(obj)            
            
            obj = obj.compute_windows();
            
            %%%%% color %%%%%
            
              
            
            %%%% these are the ones in use %%%%
            %%%%% simple %%%%%
            if(~isempty(obj.windows))
                %iterative = false;
                %[hyp_conns, types] = obj.generate_subframes_growth_hyps_unary(iterative);
                %obj.P = obj.P.add_hyp_connections(types, hyp_conns);
            
                %%%% no color, internal seed %%%%%
                %[hyp_conns, types] = obj.generate_subframe_plus_internal_growth_hyps();
                %obj.P = obj.P.add_hyp_connections(types, hyp_conns);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% these are the three we use more %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                [hyp_conns, types] = obj.generate_subframes_growth_hyps();
                obj.P = obj.P.add_hyp_connections(types, hyp_conns);

                %%%% color, internal seed %%%
                [hyp_conns, types] = obj.generate_subframe_plus_internal_growth_hyps_unary();
                obj.P = obj.P.add_hyp_connections(types, hyp_conns);            

                %%%% color, iterative (it searches for good foreground pixels for seed )%%%
                %ok for wiry stuff %            
                %iterative = true;            
                %[hyp_conns, types] = obj.generate_subframes_growth_hyps_unary(iterative);
                %obj.P = obj.P.add_hyp_connections(types, hyp_conns);            
            end
        end

        % just call the superclass methods
        function [hyp_conns, types] = generate_interactive_subframe_plus_internal(obj) 
            [hyp_conns, types] = generate_interactive_subframe_plus_internal@Segmenter(obj);
        end        

        function [hyp_conns, types] = generate_interactive_subframe_plus_internal_unary(obj)
            [hyp_conns, types] = generate_interactive_subframe_plus_internal_unary@Segmenter(obj);
        end

        function [hyp_conns, types] = generate_interactive_subframe(obj)
            [hyp_conns, types] = generate_interactive_subframe@Segmenter(obj);
        end

        function [hyp_conns, types] = generate_subframes_growth_hyps(obj)
            [hyp_conns, types] = generate_subframes_growth_hyps@Segmenter(obj);
        end        

        function [hyp_conns, types] = generate_subframe_plus_internal_growth_hyps_unary(obj)
            [hyp_conns, types] = generate_subframe_plus_internal_growth_hyps_unary@Segmenter(obj);
        end
        
        function [hyp_conns, types] = generate_subframes_growth_hyps_unary(obj, iterative)
            [hyp_conns, types] = generate_subframes_growth_hyps_unary@Segmenter(obj, iterative);
        end
        
    end
end