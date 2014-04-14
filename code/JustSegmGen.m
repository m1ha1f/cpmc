% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef JustSegmGen < HgraphGen
  properties
    extractor % the collection of hyperedge extractor functionalities
    local_dist_type
    type_gen
    
    DEBUG_MODE
    DO_CACHING
    CONTRAST_SENSITIVE_WEIGHT % fill this!
    CONTRAST_SENSITIVE_SIGMA  % contrast sensitive potential sigma parameter
    PAIRWISE_OFFSET  % weight of potts potential    
  end
  
  methods
    % Construct an object
    function obj = JustSegmGen(params)   
        if(nargin==0)
            params = JustSegmGen.create_pars();
        end
        
        obj = obj.set_pars(params);
        obj.extractor = HgraphExtractor(obj.DEBUG_MODE, obj.DO_CACHING);
    end
    
    
    function obj = set_pars(obj,params)                      
      obj.CONTRAST_SENSITIVE_WEIGHT = params.CONTRAST_SENSITIVE_WEIGHT;
      obj.CONTRAST_SENSITIVE_SIGMA = params.CONTRAST_SENSITIVE_SIGMA; % bigger is smoother, smaller is sharper      
      obj.PAIRWISE_OFFSET = params.PAIRWISE_OFFSET;
      obj.local_dist_type = params.local_dist_type;      
      obj.DEBUG_MODE = params.DEBUG_MODE;
      obj.DO_CACHING = params.DO_CACHING;
      %obj.type_gen = params.type_gen;
    end
    
    function [hgraph] = apply(obj,Img, tmp_filename, bbox) 
        if(exist('bbox', 'var'))
            nrows = bbox(3) - bbox(1) + 1;
            ncols = bbox(4) - bbox(2) + 1;
        else  
          bbox = [];
          nrows = size(Img,1);
          ncols = size(Img,2);
        end
        
      I_sz = [nrows ncols];
      n_nodes = prod(I_sz);
      hgraph = Hgraph(n_nodes);
            
      if(strcmp(obj.local_dist_type,'color'))
        [ids, val] = obj.extractor.pairwise_diff_aff(Img,obj.CONTRAST_SENSITIVE_WEIGHT,obj.PAIRWISE_OFFSET,obj.CONTRAST_SENSITIVE_SIGMA, bbox);
      elseif(strcmp(obj.local_dist_type,'pb'))
        [ids, val] = obj.extractor.pairwise_pb_aff(Img,obj.CONTRAST_SENSITIVE_WEIGHT,obj.PAIRWISE_OFFSET,obj.CONTRAST_SENSITIVE_SIGMA, tmp_filename, bbox);   
      elseif(strcmp(obj.local_dist_type,'pb_fat'))
        [ids, val] = obj.extractor.pairwise_pb_aff(Img,obj.CONTRAST_SENSITIVE_WEIGHT,obj.PAIRWISE_OFFSET,obj.CONTRAST_SENSITIVE_SIGMA, tmp_filename, bbox, true);           
      end
      
      if(obj.DEBUG_MODE)
        %hist(val,20);
        %pause;
      end
      
      if (~isempty(val>0))
        hgraph = hgraph.add_edges(ids, val);
      end      
    end       
  end  % methods

  methods (Static)
    function params = create_pars()
        params.CONTRAST_SENSITIVE_WEIGHT = 1;
        params.CONTRAST_SENSITIVE_SIGMA = 0.1; % bigger is smoother, smaller is sharper
        params.PAIRWISE_OFFSET = 4/1000;
        params.local_dist_type = 'pb';
        params.type_gen = 'local';
        params.DEBUG_MODE = false;
        params.DO_CACHING = false;
    end    
  end
end % class

      
