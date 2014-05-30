function test()
	addpath('../code');
    addpath('../external_code');
    addpath('../external_code/npp');
    addpath('../external_code/paraFmex');
    addpath('../external_code/vlfeats/toolbox/kmeans/');
    addpath('../external_code/vlfeats/toolbox/mex/mexa64/');
    addpath('../external_code/vlfeats/toolbox/mex/mexglx/');
	
	%img_name = '2010_000238'; % airplane and people   
    %img_name = '2007_009084'; % dogs, motorbike, chairs, people    
    %img_name = '2010_002868'; % buses   
    %img_name = '2010_003781'; % cat, bottle, potted plants
    
    imgs = dir('/home/mihai/Downloads/VOCdevkit/VOC2012/JPEGImages/');
    
    K = 1;
    
    total_overlap = 0;
    n_segments = 0;

    l = 1/1000;
    u = 15;
    
    lambda_range =  [logspace(log10(l), log10(u), 20)];
    lambdas_mismatch = 0;
    
    bad_imgs = [];
    bad_lambdas = [];
    
    ov = [];
    od = [];
    
    for i = 1:numel(imgs)
        img = imgs(i);
        idx = strfind(img.name, '.jpg');
        if isempty(idx) || idx(end) ~= numel(img.name)-3
            continue;
        end
        
%         I = imread(['/home/mihai/Downloads/VOCdevkit/VOC2012/JPEGImages/' img.name]);
        I = imread(['/home/mihai/Downloads/VOCdevkit/VOC2012/JPEGImages/' '2007_000123.jpg']);
        I = single(I)/255;

        width = 15;
        height = 15;
        prow = randi(size(I, 1) - height, 1);
        pcol = randi(size(I, 2) - width, 1);
        [foregroundids, foregroundSeeds, backgroundids, backgroundSeeds] = get_seeds(I, width, height, prow, pcol);
        
        for j = 1:numel(lambda_range)
            l = lambda_range(j);
            u = lambda_range(j);
            [cuts_nppi, lambdas_nppi] = test_nppiGraphCut(I, K, foregroundids, foregroundSeeds, backgroundids, backgroundSeeds, l, u);
            [cuts_hochbaum, lambdas_hochbaum] = test_hochbaum(I, K, foregroundids, foregroundSeeds, backgroundids, backgroundSeeds, l, u);
        
            assert(all(size(cuts_nppi) == size(cuts_hochbaum)));
            if lambdas_nppi ~= round(lambdas_hochbaum*1000)
                lambdas_mismatch = lambdas_mismatch + 1;
                lambdas_hochbaum
                lambdas_nppi
                jsalkfjslkdjf;
            end

            for k = 1:size(cuts_nppi, 2)
                n_segments = n_segments + 1
                overlap = sum(cuts_nppi(:, k) & cuts_hochbaum(:, k)) / sum(cuts_nppi(:, k) | cuts_hochbaum(:, k));
                ov = [ov overlap];
                total_overlap = (total_overlap*(n_segments - 1) + overlap)/n_segments
                lambdas_mismatch
                
                if overlap < 0.2
                    bad_imgs = [bad_imgs img.name];
                    bad_lambdas = [bad_lambdas l];
                end
                
                mismatch = cuts_nppi(:, k) ~= cuts_hochbaum(:, k);
                mismatch_ids = find(mismatch);
                [mx my] = ind2sub([size(I, 1) size(I, 2)], mismatch_ids);
                d = mean(pdist([mx my]))/sqrt((size(I, 1)-1)^2 + (size(I, 2)-1)^2);
                od = [od d];
            end
        end
        
        if n_segments > 2000
                break;
            end
    end
    
    hist(ov);
    bad_imgs
    bad_lambdas 
    
    avg_dist = mean(od)
end