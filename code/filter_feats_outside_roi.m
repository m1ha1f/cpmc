% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

% function [D1,F1,idx] = filter_feats_outside_roi(D1,F1,Roi,type)
%
% Removes features with center outside region of interest.
% INPUT:
%  F1  - Frames of the descriptors
%  D1  - Descriptors
%  Roi - binary image. White means inside region of interest. Invert
%  contrast to filter out features inside the region of interest.
%  type - 'points' or 'linesegmentpairs', depending on kind of feature
%
% OUTPUT:
%  D1, F1 - filtered descriptors and frames
%  idx    - indices of the features kept.

function [D1,F1,idx] = filter_feats_outside_roi(D1,F1,Roi,type)
  if strcmp(type, 'points')
    Pos = F1(1:2,:);
    Pos = single(round(Pos));
    sz = size(Roi);

    Ind = sub2ind([sz(1) sz(2)],Pos(2,:),Pos(1,:));
    In = (Roi(Ind) == 1);

    F1 = F1(:,In);
    if(~isempty(D1))
      D1 = D1(:,In);
    end
    idx = find(In);
  elseif strcmp(type, 'linesegmentpairs')
    Pos1 = F1(1:2,:);
    Pos1 = single(round(Pos1));
    Pos2 = F1(3:4,:);
    Pos2 = single(round(Pos2));    
    
    sz = size(Roi);

    Ind1 = sub2ind([sz(1) sz(2)],Pos1(1,:),Pos1(2,:));
    Ind2 = sub2ind([sz(1) sz(2)],Pos2(1,:),Pos2(2,:));    
    In = ((Roi(Ind1) == 255) & (Roi(Ind2) == 255));

    F1 = F1(:,In);
    D1 = D1(:,In);
    idx = find(In > eps);
  else 
    error('No such type of feature!');
  end