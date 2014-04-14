% Global Pb Example
%
% Notes:
%  - Threshold for optimal F-measure (for rsz = 0.5) is 0.13
%  - Calls pre-compiled executable files:
%     - local cues by Michael Maire (mex files)
%     - Intervening Contour by Charless Fowlkes (segment program)
%    You can obtain the source code for these files from:
%    http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/ 
%
% Pablo Arbelaez.

% usage example:
addpath('lib')

imgFile='101087.jpg';
outFile='101087.bmp';
rsz=0.5;
tic;
[gPb_thin, gPb, maxo] = globalPb(imgFile, outFile, rsz);
toc;
