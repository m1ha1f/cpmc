function debug_paraFmex()
  img_side = 20;
  grid_spacing = 1;
  noise_level = 0.1;
  seed_side = 1;
  cut_factor = 1000000/(img_side^2);
  sq_side  = [10]; % fixed
  
%   [Img, gtruth] = rect_images(1, 'Gaussian', noise_level, img_side, sq_side);   
%   Img = Img{1};
%   sz = size(Img);
% 

  Img = zeros(4,4);
  Img(1,2) = 0.05;
  Img(2,3) = 0.2;
  Img(2,4) = 0.12;
  Img(1,3) = 0.4;
  Img(3,4) = 0.03;
  Img = imnoise(Img, 'Gaussian', 0.1);
  
  GS_correct = SharonGridSegm(Img);
  
  t_normal = tic();
  [segm_corr,cost_corr] = apply(GS_correct, cut_factor);
  time_unoptimized = toc(t_normal)
  
  disp('end');



  
  function [segm,cost] = apply(theGS, the_cut_factor)
    dg = theGS.normal_dgraph;
    theGS = theGS.set_seed_type('Square', seed_side);
    theGS.Grid.spacing = grid_spacing;
    theGS.DEBUG_MODE = true;
    theGS.DISC_FACTOR_CUT = the_cut_factor;

    [segm, cost] = theGS.compute_best();
  end
end