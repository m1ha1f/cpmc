% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef Connection 
    properties
        nodes_a
        nodes_b    
        edge_strength        
        is_parametric
        parametric_weight
    end    
end
    