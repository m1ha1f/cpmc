% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

%function h = subplot_auto_transparent(segments, I, titles)
function [Imgs] = subplot_auto_transparent(segments, I, titles, voc_cmap_ids, grid_type)
  segments(segments==inf) = 10000;
  
  border_side = 0.05;
  border_top = 0.05;
  if(length(size(I))==2)
    I = repmat(I, [1 1 3]);
  end
  if(issparse(segments))
    segments = full(segments);
  end
  
  if(size(segments,1) ~= size(I,1))
    n_imgs = size(segments,2);
  else
    n_imgs = size(segments,3);
    segments = reshape(segments, size(segments,1) * size(segments,2), n_imgs);
  end

  n_rows = round(sqrt(n_imgs));
  n_cols = ceil(n_imgs/n_rows);

  assert(n_cols*n_rows >= n_imgs);

  counter = 1;
  
  if(exist('voc_cmap_ids', 'var') && ~isempty(voc_cmap_ids))      
      cmap = VOClabelcolormap();
      assert(size(segments,2) == numel(voc_cmap_ids));
      
      for i=1:size(segments,2)
          bground{i} = zeros(size(I,1),size(I,2),3);
          bground{i}(:,:,1) = cmap(voc_cmap_ids(i),1);
          bground{i}(:,:,2) = cmap(voc_cmap_ids(i),2);
          bground{i}(:,:,3) = cmap(voc_cmap_ids(i),3);
          bground{i} = uint8(255*bground{i});
      end
  else      
      bground = zeros(size(I,1),size(I,2),3);
      bground(:,:,2) = 255;
  end
  
  Imgs = cell(n_imgs,1);
  for i=1:n_imgs        
    if(size(segments,1) ~= size(I,1))
        alpha_chann = reshape(segments(:,i), size(I,1), size(I,2))*0.5;    
    else 
        alpha_chann = segments(:,:,i)*0.5;
    end
    
    %sc(sc(I).*sc(alpha_chann))
    if(exist('voc_cmap_ids', 'var') && ~isempty(voc_cmap_ids))
        Imgs{i} = immerge(I, bground{i}, alpha_chann);
    else
        Imgs{i} = immerge(I, bground, alpha_chann);
    end

    counter = counter + 1;
  end  

  %montage_new(Imgs, titles, 'Size', [n_rows n_cols], 'Border', 0.1);
  if(~exist('grid_type', 'var'))
      grid_type = [];
  end
  
  if(exist('titles', 'var') && ~isempty(titles))
    if(~iscell(titles)) % if not a cell assumes they're numbers
        for i=1:length(titles)
            new_titles{i} = sprintf('%f', titles(i));
        end
        titles = new_titles;
    end
    if(~exist('grid_type', 'var'))
        montage_new(Imgs, titles, 'Border',  [border_top border_side]);
    else
        montage_new(Imgs, titles, 'Border', [border_top border_side], 'Size', grid_type);
    end
  else
    montage_new(Imgs, [], 'Border',  [border_top border_side], 'Size', grid_type);
  end
%   hold on;
%   for i=1:n_imgs
%     if(nargin==3)
%       h(counter) = subplot(n_rows,n_cols,i); title(titles{i});
%     end
%   end