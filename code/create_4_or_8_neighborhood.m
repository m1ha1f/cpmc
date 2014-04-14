% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

% Retrieves a neighborhood matrix, either 4 or 8 connected
% Input: 
%   - the_size is the size of the matrix [m n]
%   - n_neighbors is either 4 or 8
function pairs = create_4_or_8_neighborhood(the_size, n_neighbors)
  padd_nrows = the_size(1)+2; % padded quantities
  padd_ncols = the_size(2)+2;
  sq_side = padd_nrows*padd_ncols;
  
  assert(n_neighbors == 4 || n_neighbors == 8);

  E = 1; W = 2; N = 3; S = 4; SE = 5; SW = 6; NE = 7; NW = 8;
  
  off(E) = padd_nrows;       % east_offset  
  off(W) = -padd_nrows;      % west_offset
  off(N) = -1;               % north_offset
  off(S) = 1;                % south_offset
  off(SE) = padd_nrows + 1;   % southeast_offset
  off(SW) = - padd_nrows + 1; % southwest_offset
  off(NE) = padd_nrows - 1;   % northeast_offset
  off(NW) = -padd_nrows - 1;  % northwest_offset

  if(n_neighbors == 8)
    the_off = off([NW N NE W]);
  else
    the_off = off([E S]);
  end
  
  [xy_img_i xy_img_j] = cartprod_mex((2:padd_nrows-1)', (2:padd_ncols-1)');
  linear_indices = sub2ind([padd_nrows padd_ncols],xy_img_i, xy_img_j);
  
  pairs = [];
  
  for i=1:length(the_off)
    ind = linear_indices + the_off(i);
    these_ind = ind(find(linear_indices));
    new_pairs = [linear_indices these_ind];
    new_pairs(these_ind<1,:) = [];
    new_pairs(these_ind>sq_side,:) = [];
    
    pairs = [pairs; new_pairs];
  end
    
  first_row = sub2ind([padd_nrows padd_ncols], ones(1,padd_ncols), 1:padd_ncols);
  last_row = sub2ind([padd_nrows padd_ncols], padd_nrows*ones(1,padd_ncols), 1:padd_ncols);
  first_col = sub2ind([padd_nrows padd_ncols], 1:padd_nrows, ones(1,padd_nrows));
  last_col = sub2ind([padd_nrows padd_ncols], 1:padd_nrows, padd_ncols*ones(1,padd_nrows));  
  forbidden_ids = unique([first_row last_row first_col last_col]);

  while(1) % remove borders
    [c, i2] = intersect(pairs(:,2),forbidden_ids);
    if(isempty(i2))
      break;
    end
    pairs(i2, :) = [];
  end
 
  old_values = unique(pairs); % this sorts
  new_values = 1:length(old_values);
  
  map = zeros(size(pairs,1),1);
  map(old_values) = new_values;
  
  pairs = map(pairs);

