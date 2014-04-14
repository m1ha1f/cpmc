% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function FeatsChunks = chunkify(Feats, n_chunks)
    n_examples = size(Feats,2);
    if n_examples < n_chunks
        FeatsChunks{1} = Feats;
        return
    end
    ex_per_chunk = ceil(n_examples/n_chunks)*ones(n_chunks,1);
    ex_per_chunk(end) = ex_per_chunk(end) - mod(n_chunks*ex_per_chunk(1),n_examples);
    previous_id = 0;
    for i=1:n_chunks
        FeatsChunks{i} = Feats(:, previous_id+1:previous_id+ex_per_chunk(i));
        previous_id = previous_id+ex_per_chunk(i);
    end
end