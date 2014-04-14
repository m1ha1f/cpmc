% 
% 32bits and 64bits mex files and Makefiles for chi2_mex are provided
% the same for segm_overlap_mex.c
% You need to set the right path to matlab in the Makefiles, in order to
% recompile them.
% 
% These files have makefiles because they use multiple cores with omp
%

% for the other files 
mex -O code/cartprod_mex.c -o code/cartprod_mex
cd ./code/
!make
cd ..
mex -O code/int_hist.c -o code/int_hist
mex -O code/intens_pixel_diff_mex.c -o code/intens_pixel_diff_mex
cd ./external_code/paraFmex/
make_pseudo()
cd ../..

% these two files contributed by andreas mueller
mex -cxx -O external_code/my_phog_desc_mex.cpp -o external_code/my_phog_desc_mex % requires boost development files
mex -O external_code/overlap_care.c -o external_code/overlap_care