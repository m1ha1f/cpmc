% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function [new_scores, sortedIndexes] = diversify_ranking(scores, overlap, lambda)
    % [new_scores, sortedIndexes] = diversify_ranking(scores, overlap, lambda)
    %
    % follows the Maximal Marginal Relevance algorithm
    %
    % parameters:
    %   scores: vector of scores for each "region/segment" to be ranked
    %   overlap: matrix with overlap between regions
    %   lambda: lambda*original_score+(1-lambda)*overlap
    %
    % return:
    %   new_scores: new scores for each of the input elements
    %   sortedIndexes: indexes of input segments, sorted by new_scores
    
    if(size(scores,2) > size(scores,1))
        scores = scores';
    end
    
    new_scores = zeros(size(scores));
    sortedIndexes = new_scores;
    [s, id] = max(scores);
    new_scores(id) = lambda*s;
    sortedIndexes(1) = id;
    
%     t = tic();
%     R = 1:numel(scores);
%     S = id;
%     R_S = setdiff(R, S);
%     %each iteration selects the MMR id
%     for h=1:numel(R)-1   
%         temp_scores = lambda* scores(R_S) -  (1-lambda) * max(overlap(R_S, S), [], 2);                
%                 
%         [s, id] =  max(temp_scores);
%         new_scores(R_S(id)) = s;
%         S = [S R_S(id)];
%         R_S(id) = [];
%     end
%     new_scores_old = new_scores;
%     toc(t)
    
 %   t = tic();
    S = false(numel(scores),1);    
    S(id) = true;
    R_S = ~S;
    
    for h=2:numel(S)
        %t1 = tic()        
        %temp_scores = lambda* scores(R_S) -  (1-lambda) * max(overlap(R_S, S), [], 2);
        
        scores(id) = [];
        temp_scores = lambda* scores -  (1-lambda) * max(overlap(R_S, S), [], 2);
        %t1 = toc(t1)
        
        [s, id] = max(temp_scores);        
        
        options = find(R_S);                
        real_id = options(id);
        
        new_scores(real_id) = s;
        sortedIndexes(h) = real_id;
        
        S(real_id) = true;
        R_S(real_id) = false;
    end
    
%    toc(t)
end
   
