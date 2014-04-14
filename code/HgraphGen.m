% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef HgraphGen
 % Property data is private to the class
  properties (SetAccess = public, GetAccess = public)
    
  end % properties
  
  methods
    
    % Construct an object
    function obj = HgraphGen()  
    end     
  end  % methods
  
  
  methods (Abstract)
    hgraph = apply(annot,inf_rep) 
  end
  
  
end % class

      