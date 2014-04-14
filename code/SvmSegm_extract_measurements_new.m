% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function SvmSegm_extract_measurements_new(type, img_names, exp_dir, MAX_DIM)
  % Extracts measurements on the whole image
  
  %type = 'dense_sift_4_scales';
  %type = 'dense_color_sift_3_scales';  
  %type = 'mean_intensity';
  %type = 'weighted_centroid';
  %type = 'mean_color';
  %type = 'phog';
  %type = 'dense_self_similarity_4_scales';
  
  DefaultVal('*MAX_DIM', 'inf');
  if(strcmp(type, 'sift'))
      type_func = @run_exp_database_do_sift;
  elseif(strcmp(type, 'dense_sift_4_scales'))
      type_func = @run_exp_database_do_dense_sift_4_scales;
  elseif(strcmp(type, 'dense_color_sift_3_scales'))
      type_func = @run_exp_database_do_dense_color_sift_3_scales;
  elseif(strcmp(type, 'dense_self_similarity_4_scales'))
      type_func = @run_exp_database_do_dense_self_similarity_4_scales;
  elseif(strcmp(type, 'mean_intensity') || strcmp(type, 'mean_color') || strcmp(type, 'weighted_centroid'))
      type_func = @run_exp_database_do_simple;
  elseif(strcmp(type, 'phog'))
      type_func = @run_exp_database_do_phog;
  elseif(strcmp(type, 'shape_context'))
      type_func = @run_exp_database_do_sc;
      error('not ready');
  else
      error('no such type defined');
  end
    
  our_folder = [exp_dir  'MyMeasurements/' type];  
  if(~isdir(our_folder))
      mkdir(our_folder);
  end
  
  %parfor (j=1:length(img_names),8)
  for j=1:length(img_names)
      filename = [our_folder '/' img_names{j} '.mat'];
      if(exist(filename,'file'))
         continue;
      end

      t = tic();
      img_file = [exp_dir '/JPEGImages/' img_names{j}];
      I = imread([img_file '.jpg'], 'jpg');

      max_dim = max(size(I));
      if(max_dim>MAX_DIM)
        I = imresize(I, MAX_DIM/max_dim);
      end

      pb_path = [exp_dir 'PB/'  img_names{j} '_PB.mat'];
      [F,D] = type_func(I, type,pb_path, exp_dir);
            
      my_save_function(F,D, filename);
      %toc(t)
  end 
end

function my_save_function(F, D,filename)
    save(filename, 'F', 'D');
end

function [F,D] = run_exp_database_do_phog(I,type, pb_path, exp_dir)
  if(size(I,3)>1)
    I = rgb2gray(I); 
  end
  nbins = 20;
  roi = [1 size(I,1) 1 size(I,2)]'; %[ytop,ybottom,xleft,xright];
  time_image=tic();
  n_levels = 2;
  F = roi;
  [D] = anna_phog(double(I), nbins, 180, n_levels, roi);
  toc(time_image)  
end

function [F,D] = run_exp_database_do_sift(I, type, pb_path, exp_dir)  
  if(size(I,3)>1) % sift doesn't use color
    I = rgb2gray(I); 
  end
  
  time_image=tic();
  [F,D] = vl_sift(single(I));
  toc(time_image)  
end

function [F,D] = run_exp_database_do_dense_sift(I, type, pb_path, exp_dir)  
  if(size(I,3)>1) % sift doesn't use color
    I = rgb2gray(I); 
  end
  
  time_image=tic();
  [F,D] = vl_dsift(single(I), 'Step', 3);
  toc(time_image)  
end

function [F,D] = run_exp_database_do_dense_sift_4_scales(I,type, pb_path, exp_dir)
  % it's really just using 3 scales!!
  
  if(size(I,3)>1) % sift doesn't use color
    I = rgb2gray(I); 
  end
  
  spacing = 4;
  time_image=tic();
  
  square_edge = [16 24 36 54 ];
  bin_size = square_edge ./ 4;
   
  D = [];
  F = [];
  for i=1:3    
    [F1,D1] = vl_dsift(single(I), 'Step', spacing, 'Size', bin_size(i));
    D = [D D1];
    F = [F F1];
  end

  toc(time_image)  
end

function [F,D] = run_exp_database_do_dense_color_sift_3_scales(I, type, pb_path, exp_dir) 
  spacing = 5;
  time_image=tic();
  
  scales = [1 3 5.3]; 
  
  D = [];
  F = [];
    
  [F,D] = csift(I, spacing, scales);

  toc(time_image)  
end

function [F,D] = run_exp_database_do_dense_self_similarity_4_scales(I, type, pb_path, exp_dir)  
    parms.size=5;
    parms.numRadiiIntervals=6;
    parms.numThetaIntervals=12;
    parms.varNoise=25*3*36;
    parms.autoVarRadius=1;
    parms.saliencyThresh=0; % I usually disable saliency checking
    parms.nChannels=size(I,3);

    radius=(parms.size-1)/2; % the radius of the patch
    
    spacing = 7;
    
    square_edge = [47 31 23 15];    
    bin_size = square_edge ./ 4;

    last_id = 0;
    D = [];
    F = [];
    
    for i=1:length(square_edge)  
        scale = base_sigma_from_sift_support_side(square_edge(i));
        parms.coRelWindowRadius=(square_edge(i)/2 - parms.size/2) ;
        marg=radius+parms.coRelWindowRadius;
        
        [X,Y]=meshgrid([marg+1:spacing:size(I,2)-marg], ...
            [marg+1:5:size(I,1)-marg]);
        
        X=X(:)';
        Y=Y(:)';
        
        F_XY = [X; Y; scale*ones(1,length(X))];
        
        n_points = size(F_XY,2);
        
        [D1]=ssimDescriptor(double(I) ,parms ,uint32(X),uint32(Y));
        D = [D D1];
        F = [F F_XY];
        last_id = last_id+n_points;
    end        
end

function [F,D] = run_exp_database_do_sc(I, type, pb_path, exp_dir)
    error('not ready');
    % 1. load pb
    load(pb_path);
    % 2. extract shape contexts along edges 
    [Bsampx, Bsampy] = find(gPb_thin > 0.1);
    Bsamp = [Bsampx';  Bsampy'];
    Tsamp = zeros(1,size(Bsamp,2));    
    nbins_theta = 12;
    nbins_r = 8;
    r_inner = 10;
    r_outer = 10;
    out_vec = Tsamp;
    D = sc_compute(Bsamp,Tsamp,[],nbins_theta,nbins_r,r_inner,r_outer, out_vec);
end

function [F,D] = run_exp_database_do_simple(I, type, pb_path, exp_dir)   
  time_image=tic(); 
  F = [];
  if(strcmp(type, 'mean_intensity'))    
    I = rgb2gray(I); 
    D = mean(mean(I));
  elseif(strcmp(type, 'weighted_centroid'))    
    I = rgb2gray(I); 
    LabelImg = ones(size(I));
    D = regionprops(LabelImg, I, 'WeightedCentroid');
    D = D.WeightedCentroid';
  elseif(strcmp(type, 'mean_color'))
    if(size(I,3)==1)
      D(1) = mean(mean(I));
      D(2) = 0;
      D(3) = 0;
    else
      D(1) = mean(mean(I(:,:,1)));
      D(2) = mean(mean(I(:,:,2)));
      D(3) = mean(mean(I(:,:,3)));
    end
  end
  
  toc(time_image)  
end
    
