function [ucm] = contours2ucm(pb_oriented, fmt)
% Creates Hierarchical Regions from oriented contours
%
% syntax:
%   [ucm] = contours2ucm(pb_oriented, fmt)
%
% description:
%   Given an oriented contour signal pb_oriented, this code computes first
%   an over-segmentation using the Oriented Watershed Transform (OWT) 
%   and then an Ultrametric Contour Map (UCM) by considering 
%   the mean pb value on the boundary between regions as dissimilarity.
%
% arguments:
%   pb_oriented: Oriented Probability of Boundary in uint8
%   fmt:         Output format. 'imageSize' (default) or 'doubleSize'
%
% output:
%   ucm:    Ultrametric Contour Map in uint8
%
% Pablo Arbelaez <arbelaez@eecs.berkeley.edu>
% April 2009

if nargin<2, fmt = 'imageSize'; end;

if ~strcmp(fmt,'imageSize') && ~strcmp(fmt,'doubleSize'),
    error('possible values for fmt are: imageSize and doubleSize');
end

% Oriented Watershed Transform:
% create finest partition and transfer contour strength
[ows] = oriented_watershed_transform(pb_oriented);

% prepare pb for ucm
ows2 = double(super_contour_4c(ows));
ows2 = clean_watersheds(ows2);
labels2 = bwlabel(ows2 == 0, 8);
labels = labels2(2:2:end, 2:2:end) - 1; % labels begin at 0 in mex file.
ows2(end+1, :) = ows2(end, :);
ows2(:, end+1) = ows2(:, end);

% Ultrametric Contour Map with mean pb.
super_ucm = uint8(ucm_mean_pb(ows2, labels));

% output
if strcmp(fmt,'doubleSize'),
    ucm = super_ucm;
else
    ucm = super_ucm(3:2:end, 3:2:end);
end

%%

function ows = oriented_watershed_transform(pb_oriented)

pb = max(pb_oriented,[],3);
ws = watershed(pb);
ws_bw = (ws == 0);

contours = fit_contour(double(ws_bw));
angles = zeros(numel(contours.edge_x_coords), 1);

for e = 1 : numel(contours.edge_x_coords)
    if contours.is_completion(e), continue; end
    v1 = contours.vertices(contours.edges(e, 1), :);
    v2 = contours.vertices(contours.edges(e, 2), :);

    if v1(2) == v2(2),
        ang = 90;
    else
        ang = atan((v1(1)-v2(1)) / (v1(2)-v2(2)));
    end
    angles(e) = ang*180/pi;
end

orient = zeros(numel(contours.edge_x_coords), 1);
orient((angles<-78.75) | (angles>=78.75)) = 1;
orient((angles<78.75) & (angles>=56.25)) = 2;
orient((angles<56.25) & (angles>=33.75)) = 3;
orient((angles<33.75) & (angles>=11.25)) = 4;
orient((angles<11.25) & (angles>=-11.25)) =5;
orient((angles<-11.25) & (angles>=-33.75)) = 6;
orient((angles<-33.75) & (angles>=-56.25)) = 7;
orient((angles<-56.25) & (angles>=-78.75)) = 8;

ows = zeros(size(ws_bw));
for e = 1 : numel(contours.edge_x_coords)
    if contours.is_completion(e), continue; end
    for p = 1 : numel(contours.edge_x_coords{e}),
        ows(contours.edge_x_coords{e}(p), contours.edge_y_coords{e}(p)) = ...
            max(pb_oriented(contours.edge_x_coords{e}(p), contours.edge_y_coords{e}(p), orient(e)), ows(contours.edge_x_coords{e}(p), contours.edge_y_coords{e}(p)));
    end
    v1=contours.vertices(contours.edges(e,1),:);
    v2=contours.vertices(contours.edges(e,2),:);
    ows(v1(1),v1(2))=max( pb_oriented(v1(1),v1(2), orient(e)),ows(v1(1),v1(2)));
    ows(v2(1),v2(2))=max( pb_oriented(v2(1),v2(2), orient(e)),ows(v2(1),v2(2)));
end
ows=uint8(ows);
%%

function pb2 = super_contour_4c(pb)
% transfer pixel boundaries to the dual lattice in Khalimsky space.
V = min(pb(1:end-1,:), pb(2:end,:));
H = min(pb(:,1:end-1), pb(:,2:end));

[tx, ty] = size(pb);
pb2 = zeros(2*tx, 2*ty);
pb2(1:2:end, 1:2:end) = pb;
pb2(1:2:end, 2:2:end-2) = H;
pb2(2:2:end-2, 1:2:end) = V;
pb2(end,:) = pb2(end-1, :);
pb2(:,end) = max(pb2(:,end), pb2(:,end-1));

%%

function [ws_clean] = clean_watersheds(ws)
% remove artifacts created by non-thin watersheds (2x2 blocks) that produce
% isolated pixels in super_contour

ws_clean = ws;

c = bwmorph(ws_clean == 0, 'clean', inf);

artifacts = ( c==0 & ws_clean==0 );
R = regionprops(bwlabel(artifacts), 'PixelList');

for r = 1 : numel(R),
    xc = R(r).PixelList(1,2);
    yc = R(r).PixelList(1,1);
    
    vec = [ max(ws_clean(xc-2, yc-1), ws_clean(xc-1, yc-2)) ...
            max(ws_clean(xc+2, yc-1), ws_clean(xc+1, yc-2)) ... 
            max(ws_clean(xc+2, yc+1), ws_clean(xc+1, yc+2)) ...
            max(ws_clean(xc-2, yc+1), ws_clean(xc-1, yc+2)) ];
    
    [nd,id] = min(vec);
    switch id,
        case 1,
            if ws_clean(xc-2, yc-1) < ws_clean(xc-1, yc-2),
               ws_clean(xc, yc-1) = 0;
               ws_clean(xc-1, yc) = vec(1);
            else
               ws_clean(xc, yc-1) = vec(1);
               ws_clean(xc-1, yc) = 0;
               
            end
            ws_clean(xc-1, yc-1) = vec(1);
        case 2,
           if ws_clean(xc+2, yc-1) < ws_clean(xc+1, yc-2),
               ws_clean(xc, yc-1) = 0;
               ws_clean(xc+1, yc) = vec(2);
           else
               ws_clean(xc, yc-1) = vec(2);
               ws_clean(xc+1, yc) = 0;
            end
            ws_clean(xc+1, yc-1) = vec(2);
            
        case 3,
            if ws_clean(xc+2, yc+1) < ws_clean(xc+1, yc+2), 
               ws_clean(xc, yc+1) = 0;
               ws_clean(xc+1, yc) = vec(3);
            else
                ws_clean(xc, yc+1) = vec(3);
                ws_clean(xc+1, yc) = 0;
            end
            ws_clean(xc+1, yc+1) = vec(3);
        case 4, 
            if ws_clean(xc-2, yc+1) < ws_clean(xc-1, yc+2), 
               ws_clean(xc, yc+1) = 0;
               ws_clean(xc-1, yc) = vec(4);
            else
               ws_clean(xc, yc+1) = vec(4);
               ws_clean(xc-1, yc) = 0;
            end
            ws_clean(xc-1, yc+1) = vec(4);
    end 
end
%%


