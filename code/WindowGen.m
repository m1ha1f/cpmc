% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

% This class generates windows that can be used as additional seeds by cpmc
% 
classdef WindowGen
    properties
        kind
        detection_classes            
        target_I_size
        toolboxes_path
        I
        load_save_path
        threshs
        DEBUG
        MAX_N
    end
    
    methods        
        function obj = WindowGen(I, kind, target_I_size, load_save_path)
            obj.MAX_N = 40;
            obj.DEBUG = false;
            paths = SvmSegm_paths();
            obj.toolboxes_path = [paths.repos_path 'external_code/'];
                        
            % is load_save_path is defined, then tries to load the windows from
            % file. if file doesn't exist, computes them and saves them to
            % file.
            % if it's not set then just computes windows
            obj.kind = kind;
            
            assert(~isempty(I));
             
            obj.detection_classes = {'aeroplane', 'bottle', 'chair', 'horse', 'sheep', 'bicycle', 'bus', 'cow', 'motorbike', 'sofa', 'bird', 'car',  'diningtable', 'person', 'train', 'boat', 'cat', 'dog', 'pottedplant', 'tvmonitor'};
            
            %obj.detection_classes = {'person'};
            obj.target_I_size = target_I_size;
            obj.I = I;
            
            if(exist('load_save_path', 'var') && ~isempty(load_save_path))
                obj.load_save_path = load_save_path;
            else
                obj.load_save_path = [];
            end
        end
        
        function windows = apply(obj, detection_classes)
            if(~isempty(obj.load_save_path))
                if(exist(obj.load_save_path, 'file'))
                    var = load(obj.load_save_path);
                    windows = var.windows;
                    if(obj.DEBUG)
                        showboxes(imresize(obj.I, obj.target_I_size(1)/size(obj.I,1)), windows);
                        pause;
                    end
                    return;
                end
            end
            
            if(exist('detection_classes', 'var') && ~isempty(detection_classes))
                assert(iscell(detection_classes));
                obj.detection_classes = detection_classes;
            end
            
            switch(obj.kind)
                case 'sliding_window_detection'
                    windows = obj.run_detectors();
                case 'objectness_window_sampler'
                    windows = obj.run_objectness();
                case 'grid_sampler'
                    windows = obj.run_grid();
                otherwise
                    error('not ready for that!');
            end            
            
            if(~isempty(windows))
                rsz_factor = obj.target_I_size(1)/size(obj.I,1);
                windows = windows*rsz_factor;
                windows((abs(windows(:,1) - windows(:,3)) <8), :) = [];
                windows((abs(windows(:,2) - windows(:,4)) <8), :) = [];
            end
            
            if(~isempty(obj.load_save_path))
                save(obj.load_save_path, 'windows');
            end
            
            if(obj.DEBUG)
                showboxes(imresize(obj.I, rsz_factor), windows);
                pause;
            end
        end
        
        function windows = run_detectors(obj)
            topbbox = [];
            
            t = tic();
            parfor i=1:numel(obj.detection_classes)                
            %for i=1:numel(obj.detection_classes)                
                if(strcmp(obj.detection_classes{i}, 'face'))                    
                    model = [obj.toolboxes_path '/FaceDetect/haarcascade_frontalface_alt2.xml'];
                    if(size(obj.I,3) ~=1)
                        Igray = double(rgb2gray(obj.I));
                    else
                        Igray = double(obj.I);
                    end
                    
                    bbox = FaceDetect(model,  Igray);
                    if(numel(bbox)==1 && bbox == -1)
                        continue;
                    else
                        bbox(:,3) = bbox(:,1) + bbox(:,3);
                        bbox(:,4) = bbox(:,2) + bbox(:,4);
                        bbox = [bbox zeros(size(bbox,1),1)];
                    end
                else
                    var = load([obj.toolboxes_path '/felz_detector_4/VOC2009/' obj.detection_classes{i} '_final.mat']);
                    model = var.model;
                                                          
                    % slow version
                   %boxes = detect(obj.I, model, obj.threshs(i));
                   
                   csc_model = cascade_model(model, '2009', 5, -1);

                  % compute the feature pyramid for image im
                  pyra = featpyramid(obj.I, csc_model);

                  % get cascade detections
                  [dets, bbox] = cascade_detect(pyra, csc_model, csc_model.thresh);
                  %bbox = process(obj.I, model);
                   
%                     if(isempty(boxes))
%                         continue;
%                     end
%                     bbox = getboxes(model, boxes);
                end
                             
                %bbox = clipboxes(obj.I, bbox);          
                %bbox = nms(bbox, 0.85); % 0.85
                
                topbbox = [topbbox; bbox];
                classes = [];
            end
            t_detectors = toc(t)
            
            if(isempty(topbbox))
                windows = [];
                return;
            end
            
            windows = obj.postprocess_boxes(topbbox, obj.MAX_N);
        end
        
        function windows = run_objectness(obj)
            topbbox = runObjectness(obj.I, 10000);
            [val, id] = sort(topbbox(:,5), 'descend');
            topbbox = topbbox(1:1000,:);
            %similarity = pdist2(topbbox(:,1:4), topbbox(:,1:4), 'L1');
            %sy = exp(-0.001*similarity);
            %[new_scores, sorted_Indexes] =  diversify_ranking(topbbox(:,5), sy, 0.7)
            %topbbox(:,5) = new_scores;
            
            %topbbox(:,5) = rand(size(topbbox,1),1);
            %showboxes(obj.I, topbbox(id(1:20),:))
            %ids = nms(topbbox, 0.99);
            %topbbox = topbbox(ids,:);
            windows = obj.postprocess_boxes(topbbox, obj.MAX_N, 1/8);
        end
        
        function windows = run_grid(obj, type, DIM)
            DefaultVal('*type', '''grid''');
            DefaultVal('*DIM', '0.6');
            
            %obj.MAX_N = ;
            
            gp = GraphProb(obj.I);
            [sets] = gp.generate_seeds('pixels', [0 0], type, [floor(sqrt(obj.MAX_N)) floor(sqrt(obj.MAX_N))]);
            if(size(sets,1) > size(sets,2))
                sets = sets';
            end
            
            [x, y] = ind2sub([size(obj.I,1) size(obj.I,2)], cell2mat(sets)');
            centers = [x y];
            %imshow(obj.I); hold on; plot(y, x, 'o');
            
            DIM = DIM * max([size(obj.I,1), size(obj.I,2)]);
            windows = [];
            for i=1:size(centers,1)                
                window = [(centers(i,2) - 0.5*DIM) (centers(i,1) - 0.5*DIM) (centers(i,2) + 0.5*DIM) (centers(i,1) + 0.5*DIM)];
                windows = [windows; window];
            end
            windows = clipboxes(obj.I, windows);
        end
            
        function windows = postprocess_boxes(obj, topbbox, MAX_N, MAX_FRACTION_OF_IMAGE)
            % remove big boys
            DefaultVal('*MAX_FRACTION_OF_IMAGE', '1/4');
            npixels = (topbbox(:,4) - topbbox(:,2)) .* (topbbox(:,3) - topbbox(:,1));
            topbbox_too_big = (npixels > (MAX_FRACTION_OF_IMAGE*(size(obj.I,1)*size(obj.I,2))));
            topbbox(topbbox_too_big, :) = [];                                                       
            
            [sorted_scores, sorted_ids] = sort(topbbox(:,5), 'descend');            
                         
            % keep up to the maximum number
            topbbox = topbbox(sorted_ids(1:min(size(topbbox,1), MAX_N)), :);            
            topbbox = obj.enlarge_bbox(topbbox, 1.1, size(obj.I));
            
            % for the small ones, enlarge them a bit, and add them to the pool
            npixels = (topbbox(:,4) - topbbox(:,2)) .* (topbbox(:,3) - topbbox(:,1)); % less than half of the image
            small_ones = npixels < ((2/5)*(size(obj.I,1)*size(obj.I,2)));
            %small_ones = 1:size(topbbox(:,1));
            if(any(small_ones))
                topbbox_large = obj.enlarge_bbox(topbbox(small_ones,:),1.1, size(obj.I));
                topbbox_large = clipboxes(obj.I, topbbox_large);             
                topbbox = [topbbox_large];            
            end
            %showboxes(obj.I, topbbox);
            windows = topbbox;
        end
        
        function new_bbox = enlarge_bbox(obj, bbox, factor, sz_I)
            % bbox format ist ([tlx tly brx bry])
            width = bbox(:,3) - bbox(:,1);
            height = bbox(:,4) - bbox(:,2);
            
            bbox_centers = [(bbox(:,1) + (width/2)) (bbox(:,2) + (height/2))];
            new_width = width * factor;
            new_height = height * factor;
            new_leftop = [(bbox_centers(:,1) - (new_width/2)) (bbox_centers(:,2) - (new_height/2))];
            new_right_bottom = [(bbox_centers(:,1) + (new_width/2)) (bbox_centers(:,2) + (new_height/2))];
            x = new_right_bottom(:,1);
            y = new_right_bottom(:,2);
            
            %new_leftop(new_leftop<1) = 1;
            %x(x>sz_I(2)) = sz_I(2);
            %y(y>sz_I(1)) = sz_I(1);
            new_bbox = [new_leftop x y bbox(:,5)];
        end
        
    end    
end
