function [gPb_thin, gPb, maxo] = globalPbnofiles(I, outFile, rsz)
% syntax:
%   [gPb_thin, gPb, maxo] = globalPb(I, rsz)
%
% description:
%   compute Globalized Probability of Boundary of a color image.
%
% arguments:
%   imgFile : jpg format
%   outFile:  bmp format (optional)
%   rsz:      resizing factor in (0,1], to speed-up eigenvector computation
%
% output:
%   gPb:      globalized probability of boundary.
%   gPb_thin: thinned contour image.
%   maxo:     orientation of maximal response
%
% Pablo Arbelaez <arbelaez@eecs.berkeley.edu>
% April 2008
if nargin<3, rsz = 1; end
if nargin<2, outFile = ''; end

if ((rsz<=0) || (rsz>1)),
    error('resizing factor rsz out of range (0,1]');
end

total_time = 0;
tic;

im = double(I) / 255;
[tx, ty, nchan] = size(im);
orig_sz = [tx, ty];

% default feature weights
if nchan == 3,
    weights = [ 0   0    0.0028    0.0041    0.0042    0.0047    0.0033    0.0033    0.0035    0.0025    0.0025    0.0137    0.0139];
else
    weights = [ 0   0    0.0054         0         0         0         0         0         0    0.0048    0.0049    0.0264    0.0090];
end

%% mPb
[mPb, mPb_rsz, bg1, bg2, bg3, cga1, cga2, cga3, cgb1, cgb2, cgb3, tg1, tg2, tg3] = multiscalePb(im, rsz);
t = toc;
total_time = total_time + t;

%% sPb
tic;
twist = floor(sum(sum(sum(I))));
twist
rand('twister', twist);
outFile2 = strcat(outFile, [sprintf('%s', rand) '_pbs.mat']);
[sPb] = spectralPb(mPb_rsz, orig_sz, outFile2);
delete(outFile2);
t = toc;
fprintf('Spectral Pb:%g\n', t);
total_time = total_time + t;

%% gPb
tic;
gPb_all = zeros(size(tg1));
for o = 1 : size(gPb_all, 3),
    l1 = weights(1)*bg1(:, :, o);
    l2 = weights(2)*bg2(:, :, o);
    l3 = weights(3)*bg3(:, :, o);

    a1 = weights(4)*cga1(:, :, o);
    a2 = weights(5)*cga2(:, :, o);
    a3 = weights(6)*cga3(:, :, o);

    b1 = weights(7)*cgb1(:, :, o);
    b2 = weights(8)*cgb2(:, :, o);
    b3 = weights(9)*cgb3(:, :, o);

    t1 = weights(10)*tg1(:, :, o);
    t2 = weights(11)*tg2(:, :, o);
    t3 = weights(12)*tg3(:, :, o);

    sc = weights(13)*sPb(:, :, o);

    gPb_all(:, :, o) = l1 + a1 + b1 + t1 + l2 + a2 + b2 + t2 + l3 + a3 + b3 + t3 + sc;
end

%% outputs
[gPb, maxo] = max(gPb_all, [], 3);

gPb_thin = gPb .* (mPb>0.05);
gPb_thin = gPb_thin .* bwmorph(gPb_thin, 'skel', inf);

maxo = maxo .* (gPb_thin>0);

gPb_thin = normalize_output(gPb_thin);
gPb = normalize_output(gPb);

if ~strcmp(outFile,'') && ~isempty(outFile)
    imwrite(gPb_thin, outFile); 
end

t=toc;
total_time = total_time + t;
fprintf('Total Time for Global Pb:%g s.\n', total_time);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [pb_norm] = normalize_output(pb)

[tx, ty] = size(pb);
pb_norm = max(0, min(1, 1.2*pb));

%contrast change with sigmoid 
beta = [-2.6433; 10.7998];
pb_norm = pb_norm(:);
x = [ones(size(pb_norm)) pb_norm]';
pb_norm = 1 ./ (1 + (exp(-x'*beta)));
pb_norm = (pb_norm-0.0667) / 0.9333;
pb_norm = reshape(pb_norm, [tx ty]);



