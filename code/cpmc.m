% function [masks, scores] = cpmc(exp_dir, img_name, diversify_const, ranker_file, segm_pars)
%
% Implementation of the CPMC algorithm. It computes a ranked set of figure-ground segmentation, with associated scores.
% 
% Inputs: 
%   - exp_dir is the folder where the data is (images, ranker models,
%   codebooks, etc.)
%   - img_name is the name of the image without the ending, assuming a jpg file
%   - diversify_const is the the diversification parameter (between 0 and 1).
%       zero will give a random ranking, one will be a very redundant
%       ranking. Default value is 0.75.
%   - ranker_file is the filename of the ranker model (defaults to the
%   version trained on the training set of the VOC2010 segmentation
%   dataset). Please switch to the version trained on trainval for testing 
%   on imgs not in the train/validation sets.
%   - segm_pars is a struct with default value defined below.
%   Check the code for details.
%
% Outputs: 
%   - masks: the figure-ground segmentations in logical matrices, sorted by score
%   - scores: scores for the sorted list of masks
% 
% If you use this code please cite the following references:
%
% @inproceedings{carreira_cvpr10,
%  author = {J. Carreira and C. Sminchisescu},
%  title = {{Constrained Parametric Min-Cuts for Automatic Object Segmentation}},
%  booktitle = {IEEE International Conference on Computer Vision and Pattern Recognition},
%  year = {2010},
%  month = {June}, 
%  pdf = {http://sminchisescu.ins.uni-bonn.de/papers/cs-cvpr10.pdf 1}
% }
%
% @misc{cpmc-release1,
%  author = "J. Carreira and C. Sminchisescu",
%  title = "Constrained Parametric Min-Cuts for Automatic Object Segmentation, Release 1",
%  howpublished = "http://sminchisescu.uni-bonn.de/code/cpmc-release1/"
% }
%
% (C) Joao Carreira 2010
% 
function [masks, scores] = cpmc(exp_dir, img_name, diversify_const, segm_pars)
    DefaultVal('*diversify_const', '0.75');
    
    %  do switch to '''attention_model_fewfeats_lambda_10.00_trainval.mat'''
    %  for imgs not in the VOC2010 trainval image set.
    DefaultVal('*ranker_file', '''attention_model_fewfeats_lambda_10.00_train.mat'''); 
    DefaultVal('*segm_pars', '[]');
    
    if(isempty(segm_pars))
        segm_pars.pb_folder = [exp_dir './PB/'];
        segm_pars.name = 'dummy_masks';

        % UniformSegmenter uses a uniform unary term. LongRangeSegmenter
        % uses a color-based unary term. GridOfFramesSegmenter is defined
        % inside subframes (rectangular regions of interest), and it's good for smaller objects.       
        % Each will generate and solve different energy functions.
        segm_pars.segm_methods = {'UniformSegmenter', 'LongRangeSegmenter', 'GridOfFramesSegmenter'};        
        
        % can set a limit on the number of segments for each kind of
        % Segmenter. The joint set will be filtered in the end.
        segm_pars.max_n_segms = [inf inf inf]; 
        segm_pars.min_n_pixels = 200; % 1000
        segm_pars.sigma = {1, 2, 0.8};
        %segm_pars.sigma = {1.2, 2.2, 1}; % larger is smoother (fewer segments)
        
        % how much to resize the image (there doesn't seem to be any numerical negative 
        % impact on the VOC segmentation dataset to resize to half).
        segm_pars.resize_factor= 0.5;
        
        % does a morphological operation to remove little wires in the
        % image. In general a good idea. But might negatively affect wiry structures
        % like bikes a little bit.
        segm_pars.morph_open = true;
        
        % Do the fast filtering step (always recommended).
        segm_pars.filter_segments = {true, true, true};
        
        % Dimensions of the grid, for each segmenter type.
        %
        % Be careful with the third set. In GridOfFramesSegmenter, a grid
        % is set up inside each subframe . [1 1] sets a single
        % foreground seed, [2 2] would set 4 foreground seeds. If you have 40 subframes like currently, it's
        % really better not to have more than [1 1].
        % On the other two [5 5] should be fine. We used [6 6] on
        % the VOC2010 challenge.
        segm_pars.grid_dims = {[5 5], [5 5], [1 1]};  
        
        % The maximum energy for the initial filtering step, for each
        % method. The proper value depends on the sigma value, which is
        % unfortunate. If  you don't to spend time on it, just leave the
        % default value. If you don't want to filter based on the energy uncomment the following line        
        %the_segm_pars.max_energy = {inf, inf, inf};        
                
        %%%% these are parameters related to the GridOfFramesSegmenter
        % Two options are implemented: a regular grid of subframes, and
        % using a bounding box detector. Using the detector gets 1% better
        % covering.        
        segm_pars.window_gen_parms.kind = 'grid_sampler';        
        segm_pars.windows_folder = [exp_dir 'WindowsOfInterest/grid_sampler'];                
        % If you want to use instead the detector you'll have to install
        % the latent svm code from http://people.cs.uchicago.edu/~pff/latent/, including the star
        % cascade. Then you can uncomment the  two bottom lines
        %segm_pars.window_gen_parms.kind = 'sliding_window_detection';
        %segm_pars.windows_folder = [exp_dir 'WindowsOfInterest/sliding_window_detection']; 

        segm_pars.window_gen_parms.det_classes = []; % this needs to be here (bug)

        segm_pars.randomize_N = 1000; % maximum number of segments to pass to the clustering step, for each type of segmenter. For speed considerations.
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% 1. compute masks %%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    t = tic();
    disp('Starting computation of figure-ground segmentations');
    masks = cpmc_masks(exp_dir, img_name, segm_pars);
    time_segm = toc(t);
    fprintf('Time computing figure-ground segmentations: %f\n', time_segm);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% 2. compute segment ranking features %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % load ranker model (trained on training set of voc segmentation dataset)   
    % substitute train by trainval for images outside the segmentation validation set
    % of VOC2010.
    load([exp_dir '/MySegmentRankers/' ranker_file]);
    
    t = tic();
    disp('Starting feature extraction.');
    D = get_features(exp_dir, img_name, segm_pars.name, masks, segment_measurements);
    t_feats = toc(t);
    fprintf('Time getting all features: %f\n', t_feats);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% 3. rank masks %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
    % normalize the features
    D = cellfun(@double, D, 'UniformOutput', false);
    D = D(:);
    
    if(iscell(scaling))
        for i=1:numel(scaling)
            D{i} =  scale_data(D{i}, scaling_type, scaling{i});
        end
        D_scaled = cell2mat(D);
    else
        D_scaled = cell2mat(D);
        D_scaled = scale_data(D_scaled, scaling_type, scaling);
    end
    
    dims = cellfun(@(a) size(a,1), D);
    Feats = D_scaled;
    [Feats] = SvmSegm_mix_linear_and_random_features(Feats, w.kernels, scaling_type, w.scaling_RF, dims, [], w.N_OUT_DIMS, w.obj, w.merge_kernels);
    pred = predict_regressor(Feats, w.weights);
   
    small_masks = imresize(masks, 0.25, 'nearest');
    small_masks = reshape(small_masks, size(small_masks,1) * size(small_masks,2), size(small_masks,3));
    overlap_mat = single(segm_overlap_mex(small_masks));
    scores = diversify_ranking(pred, overlap_mat, diversify_const);   
    
    [sorted_scores, ids] = sort(scores, 'descend');
    masks = masks(:,:,ids);
    scores = sorted_scores;
end

function D = get_features(exp_dir, img_name, mask_type, masks, feat_types)
    segm_feats = {'SimpleSegmentFeatures', 'GestaltSegmentFeatures'};
    sf_in = [];
    
    t1 = tic();
    for i=1:numel(segm_feats)
        sf_in = [sf_in find(strcmp(segm_feats{i}, feat_types))];
    end    
    
    D_1 = SvmSegm_compute_segm_feats_nofiles(exp_dir, img_name, masks, segm_feats(sf_in));    
    t_segm_feats1 = toc(t1);
    fprintf('Time getting first set of features: %f\n', t_segm_feats1);
    
    t2 = tic();
    recog_feats_in = setdiff(1:numel(feat_types), sf_in);
    delete_local_feats = false;
    D_2 = SvmSegm_extract_all_class_features_nofiles(exp_dir, img_name, mask_type, feat_types(recog_feats_in), delete_local_feats);    
    t_segm_feats2 = toc(t2);
    fprintf('Time getting second set of features: %f\n', t_segm_feats2);
    D = [D_1 D_2];
    
    % reorder
    D = D([sf_in recog_feats_in]);
end
