% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

%function pixel_ids = frame_pixel_ids(nrows, ncols, width, custom)
% creates a frame whose external border has nrows and ncols, and width = 'width'
function pixel_ids = frame_pixel_ids(nrows, ncols, width, custom)
  assert(width>=1);
  
  [x,y] = cartprod_mex((1:width)', (1:ncols)');
  boundary_pixels_horiz_top = [x y];
  [x,y] = cartprod_mex((nrows-width+1:nrows)', (1:ncols)');
  boundary_pixels_horiz_bottom = [x y];
  [x,y] = cartprod_mex((1:nrows)', (1:width)');
  boundary_pixels_vert_left = [x y];
  [x,y] = cartprod_mex((1:nrows)', (ncols-width+1:ncols)');
  boundary_pixels_vert_right = [x y];
  
  if(strcmp(custom, 'horiz'))
    seeds = [boundary_pixels_horiz_bottom; boundary_pixels_horiz_top];
  elseif(strcmp(custom, 'down'))
    seeds = [boundary_pixels_horiz_bottom];
  elseif(strcmp(custom, 'up'))    
    seeds = [boundary_pixels_horiz_top];
  elseif(strcmp(custom, 'vert'))
    seeds = [boundary_pixels_vert_left; boundary_pixels_vert_right];
  elseif(strcmp(custom, 'all_but_down'))
    seeds = [boundary_pixels_horiz_top; boundary_pixels_vert_left;...
             boundary_pixels_vert_right];
  elseif(strcmp(custom, 'all')) % whole frame      
    seeds = [boundary_pixels_horiz_top; boundary_pixels_horiz_bottom; ...
             boundary_pixels_vert_left; boundary_pixels_vert_right];             
  else
    error('no such option');
  end

  pixel_ids = sub2ind([nrows ncols], seeds(:,1), seeds(:,2));
  pixel_ids = unique(pixel_ids);