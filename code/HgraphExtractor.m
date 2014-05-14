% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

 classdef HgraphExtractor < handle
  properties (SetAccess = private, GetAccess = public)    
    DEBUG_MODE
    DO_CACHING
  end
  
  methods
    function obj = HgraphExtractor(DEBUG_MODE, DO_CACHING)
      if(nargin == 0)
        obj.DEBUG_MODE = true;
        obj.DO_CACHING = false;
      else
        obj.DEBUG_MODE = DEBUG_MODE;
        obj.DO_CACHING = DO_CACHING;
      end
    end

    
    function [ids, val] = pairwise_diff_aff(obj, I, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA, bbox)
      
      a = 0.5 / ((SIGMA * SIGMA));

      if(exist('bbox', 'var') && ~isempty(bbox))
          bbox_range_1 = bbox(2):bbox(4);
          bbox_range_2 = bbox(1):bbox(3);
          
          nrows = numel(bbox_range_1);
          ncols = numel(bbox_range_2);
          
          I = I(bbox_range_1, bbox_range_2, :);
      end
      
      nrows = size(I,1);
      ncols = size(I,2);

      ids = create_4_or_8_neighborhood([nrows ncols], 4);

      I = double(I);            
      I = I/255.0;
      
      I = histeq(rgb2gray(I));

%      [edges,ephase] = quadeg(I);   
%      I = ephase;
%       canny_edges = edge(rgb2gray(I), 'canny');
%       edges_thin = edges.* canny_edges;
%       val = intens_pixel_diff_mex(edges_thin(:,:,1), uint32(ids(:,1)),uint32(ids(:,2)));         
%    
           
      if(size(I,3) == 3)        
        the_val = rgb_pixel_diff_mex(I(:,:,1), I(:,:,2), I(:,:,3),uint32(ids(:,1)),uint32(ids(:,2)));  % should check out this code, see if there's some problem
        the_val = the_val/3;
      elseif(size(I,3) == 1)
        the_val = intens_pixel_diff_mex(I(:,:,1), uint32(ids(:,1)),uint32(ids(:,2)));         
      else        
        error('not ready for that!');
      end


      %mat = sparse(ids(:,1), ids(:,2), the_val, size(I,1)*size(I,2), size(I,1)*size(I,2));
      %newI = sum(mat);
      %newI = reshape(newI, size(I,1), size(I,2));
      %sc(full(newI));
      
      DIV_FACTOR = 5;      
      the_val = the_val/(DIV_FACTOR);
      val = (CONTRAST_SENSITIVE_WEIGHT*exp(-(the_val)*a)) + POTTS_WEIGHT + 0.007;
    end
    
    function [ids, val] = pairwise_pb_aff(obj, I, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA, tmp_filename, bbox, fat_pb)       
      a = 0.5 / (SIGMA * SIGMA);

      
      if(exist('bbox', 'var') && ~isempty(bbox))
          bbox_range_1 = bbox(2):bbox(4);
          bbox_range_2 = bbox(1):bbox(3);
          
          nrows = numel(bbox_range_1);
          ncols = numel(bbox_range_2);
      else      
        nrows = size(I,1);
        ncols = size(I,2);
      end      
      
      ids = create_4_or_8_neighborhood([nrows ncols], 4);
            
      %if(~isempty(tmp_filename) && exist(tmp_filename, 'file'))
      if(exist('fat_pb', 'var') && fat_pb)
          load(tmp_filename, 'gPb_orient');
          gPb = mean(gPb_orient,3);
          DIV_FACTOR = 0.3;   % smaller is sharper edges, larger is smoother edges ( 5 smooth, 0.5 very sharp )
      else
          load(tmp_filename, 'gPb_thin');
          gPb = gPb_thin;
          DIV_FACTOR = 200;   % smaller is sharper edges, larger is smoother edges ( 5 smooth, 0.5 very sharp )
      end
      
      %else
      %    [gPb_orient, gPb_thin, textons] = globalPb_new_nofiles(I, [], 0.6);
      %    if(~isempty(tmp_filename))
      %        save(tmp_filename, 'gPb_thin', 'gPb_orient', 'textons');
      %    end
      %    %gPb_thin = double(gPb_thin);
      %end
      
      if(any(any(0 > gPb)))
          gPb = gPb+abs(min(min(gPb)));
      end
%       
%       DEBUG = true;
%       if(DEBUG)
%           hist(reshape(pb_prob, prod(size(pb_prob)), 1));
%       end
%       
      if(exist('bbox', 'var')&& ~isempty(bbox))
          gPb = gPb(bbox_range_1, bbox_range_2);
          I = I(bbox_range_1, bbox_range_2, :);
      end      
      
      if(ncols ~=size(gPb,2))
          gPb = imresize(gPb, 'OutputSize', [nrows ncols]);
      end
           
      val = intens_pixel_diff_mex(double(gPb), uint32(ids(:,1)),uint32(ids(:,2)));
      
      %mat = sparse(ids(:,1), ids(:,2), val, size(I,1)*size(I,2), size(I,1)*size(I,2));
      %newI = sum(mat);
      %newI = reshape(newI, size(I,1), size(I,2));
      %sc(full(newI));
      %val = intens_pixel_max_mex(double(gPb), uint32(ids(:,1)),uint32(ids(:,2)));
      
      %maxo = pi./maxo;
      %maxo(maxo==inf) = 0;
      %val = ic_mine(gPb_thin, maxo, uint32(ids(:,1)),uint32(ids(:,2))); 
            
      %DIV_FACTOR = 1.3;
      %val = val/DIV_FACTOR;
      new_val = normVal(val, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
    end 

    function [leftTranspose, rightTranspose, top, bottom] = pairwise_pb_aff_matrices(obj, I, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA, tmp_filename, bbox, fat_pb)       
      if(exist('bbox', 'var') && ~isempty(bbox))
          bbox_range_1 = bbox(2):bbox(4);
          bbox_range_2 = bbox(1):bbox(3);
          
          nrows = numel(bbox_range_1);
          ncols = numel(bbox_range_2);
      else      
        nrows = size(I,1);
        ncols = size(I,2);
      end      
                  
      if(exist('fat_pb', 'var') && fat_pb)
          load(tmp_filename, 'gPb_orient');
          gPb = mean(gPb_orient,3);
          DIV_FACTOR = 0.3;   % smaller is sharper edges, larger is smoother edges ( 5 smooth, 0.5 very sharp )
      else
          load(tmp_filename, 'gPb_thin');
          gPb = gPb_thin;
          DIV_FACTOR = 200;   % smaller is sharper edges, larger is smoother edges ( 5 smooth, 0.5 very sharp )
      end
      
      if(any(any(0 > gPb)))
          gPb = gPb+abs(min(min(gPb)));
      end
       
      if(exist('bbox', 'var')&& ~isempty(bbox))
          gPb = gPb(bbox_range_1, bbox_range_2);
          I = I(bbox_range_1, bbox_range_2, :);
      end      
      
      if(ncols ~=size(gPb,2))
          gPb = imresize(gPb, 'OutputSize', [nrows ncols]);
      end

      leftTranspose = gPb(1: nrows, 2:ncols) - gPb(1:nrows, 1:ncols-1);
      leftTranspose = abs(leftTranspose);
      leftTranspose = obj.normVal(leftTranspose, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
      leftTranspose = [zeros(nrows, 1) leftTranspose];
      leftTranspose = leftTranspose';

      rightTranspose = gPb(1:nrows, 1:ncols-1) - gPb(1:nrows, 2:ncols);
      rightTranspose = abs(rightTranspose);
      rightTranspose = obj.normVal(rightTranspose, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
      rightTranspose = [rightTranspose zeros(nrows, 1)];
      rightTranspose = rightTranspose';

      top = gPb(2:nrows, :) - gPb(1:nrows-1, :);
      top = abs(top);
      top = obj.normVal(top, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
      top = [zeros(1, ncols); top];

      bottom = gPb(1:nrows-1, :) - gPb(2:nrows, :);
      bottom = abs(bottom);
      bottom = obj.normVal(bottom, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA);
      bottom = [bottom; zeros(1, ncols)];
    end       
    
  end % methods
  
  methods (Access=private)
    function new_val = normVal(obj, val, CONTRAST_SENSITIVE_WEIGHT, POTTS_WEIGHT, SIGMA)
      a = 0.5 / (SIGMA * SIGMA);
      val = single(val);
      new_val = (CONTRAST_SENSITIVE_WEIGHT*exp(-(val)*a)) + POTTS_WEIGHT+ 0.007;           
      duh = hist(new_val);
      if(duh(2)>duh(1))
          a = a + 0.07;
          new_val = (CONTRAST_SENSITIVE_WEIGHT*exp(-(val)*a)) + POTTS_WEIGHT+ 0.007;  
      end
    end
  end
end
