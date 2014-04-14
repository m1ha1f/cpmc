% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Calculate the F-score of the evaluated method            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% best_s is the filename of the best figure-ground segmentation (needs to
% be binary!!!) Only considers the segment with label 0

%%% P and R are ignored, the score is in F
function [P R F]=myCalcCandScoreFigureGroundAll(Segmaps,HumanSeg, type, care)
  Fmax=0;
  F = zeros(size(Segmaps,3),1);
  P = F;
  R = F;

  if(strcmp(type, 'overlap_and_bbox_edges') || strcmp(type, 'overlap_action'))
      % do a priori computation
      var = regionprops(HumanSeg, 'BoundingBox');
      bb_human = var.BoundingBox;
      
      bb_segms = zeros(size(Segmaps,3),4);
      for i=1:size(Segmaps,3)
        var = regionprops(Segmaps(:,:,i), 'BoundingBox');
        if(~isempty([var.BoundingBox]))
          bb_segms(i,:) = var.BoundingBox;
	else
          bb_segms(i,:) = [-1 -1 1 1];
	end
      end
  end
  
  gt_img_fraction = sum(sum(HumanSeg)) /numel(HumanSeg);
  for i=1:size(Segmaps,3)
     FragMap = Segmaps(:,:,i);

     if(strcmp(type, 'fmeasure'))
        [p r f]=CalcPRPixel(HumanSeg,FragMap, care);%Calculate the F-score
     elseif(strcmp(type, 'overlap'))
        [p r f]=CalcOverlap(HumanSeg,FragMap, care);%Calculate the overlap, puts it in f
     elseif(strcmp(type, 'overlap_and_bbox_edges') || strcmp(type, 'overlap_action'))
         [p r f] = CalcOverlapBBoxEdges(HumanSeg, FragMap, bb_human, bb_segms(i,:), care);
     elseif(strcmp(type, 'overlap_times_area_fraction'))
         [p r f]=CalcOverlap(HumanSeg,FragMap, care);%Calculate the overlap, puts it in f
         % now multiple by area fraction
         f = f*gt_img_fraction;
     end
     
     F(i)=f;
     P(i)=p;
     R(i)=r;
   end;%Go over all segmentations in the Dir
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (ADAPTED FROM WEIZMANN DATABASE CODE)             Calculate the F-score      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p r f]=CalcPRPixel(GT,mask, care)
    if(exist('care', 'var'))
        error('not ready yet, code this!');
    end
    
    % CARREIRA if
    if(size(GT,3) == 3)
      GT = GT(:,:,1);
    end
    
    sum_GT_and_mask = sum(GT(:)&mask(:));
    if (sum_GT_and_mask==0)
        p=0;r=0;f=0;
        return;
    end;
    r=sum_GT_and_mask./sum(GT(:));
    c=sum(mask(:))-sum_GT_and_mask;
    p=sum_GT_and_mask./(sum_GT_and_mask+c);
    f=(r*p)/(0.5*(r+p));
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Overlap %%%%%%%%%%%%%%%%%%%%%%%

function [p, r, f] = CalcOverlap(GT, mask, care)
    p  = -1;
    r = -1;
    %f = sum(flat(GT & mask & care)) / sum(flat((GT | mask) & care));
    f=overlap_care(logical(GT),logical(mask),logical(care));
    %x-f
    %assert(x==f);

end

function [p, r, f] = CalcOverlapBBoxEdges(GT, mask, bb_gt, bb_mask, care)
    p  = -1;
    r = -1;

    width_gt = bb_gt(3);
    height_gt = bb_gt(4);
    width_mask = bb_mask(3);
    height_mask = bb_mask(4);

    inters = rectint(bb_gt, bb_mask);
    area_bb_gt = width_gt*height_gt;
    area_mask = width_mask*height_mask;
    
    reun = (area_bb_gt+area_mask - inters);
    
    f = inters/reun;
    
%     min_x_gt = bb_gt(1);
%     max_x_gt = bb_gt(1)+bb_gt(3);
%     min_y_gt = bb_gt(2);
%     max_y_gt = bb_gt(2) + bb_gt(4);
%     
%     

%     min_x_mask = bb_mask(1);
%     max_x_mask = bb_mask(1) + bb_mask(3);
%     min_y_mask = bb_mask(2);
%     max_y_mask  = bb_mask(2) + bb_mask(4);
%     
%     diffs(1) = abs(min_x_mask - min_x_gt) / width_mask;
%     diffs(2) = abs(max_x_mask - max_x_gt ) / width_mask;
%     diffs(3) = abs(min_y_mask - min_y_gt ) / height_mask;
%     diffs(4) = abs(max_y_mask - max_y_gt ) / height_mask;
%     
%     f = mean(diffs);
    
    %%% visualization %%%
    %showboxes(GT*255, [bb_gt(1) bb_gt(2) bb_gt(1) + bb_gt(3) bb_gt(2) + bb_gt(4)])
    %showboxes(mask(:,:,1)*255, [bb_mask(1,1) bb_mask(1, 2) bb_mask(1,1) + bb_mask(1,3) bb_mask(1,2) + bb_mask(1,4)])
end


function y=flat(x)
  y=x(:);
end
