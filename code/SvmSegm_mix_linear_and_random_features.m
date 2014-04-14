% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function [Feats, rf_obj, RF_scaling] = SvmSegm_mix_linear_and_random_features(Feats, kernels, scaling_type, RF_scaling, dims, gamma, N_OUT_DIMS, rf_obj, merge_kernels)
    DefaultVal('merge_kernels', '[]');
    
    counter = 1;
    kernel_range = [];
       
    for i=1:numel(dims)
        curr_range = counter:counter+dims(i)-1;
        kernel_range{i} = curr_range;
        counter = counter + dims(i);
    end    
     
    if(~isempty(merge_kernels))
        for i=1:numel(merge_kernels)
            kernel_range{min(merge_kernels{i})} = [kernel_range{merge_kernels{i}}];
            to_clear = setdiff(merge_kernels{i}, min(merge_kernels{i}));
            for j=1:numel(to_clear)
                kernel_range{to_clear(j)} = [];
                kernels{to_clear(j)} = [];
            end
            
            merge_kernels{i} = [];
        end
        kernel_range(cellfun(@isempty, kernel_range)) = [];
        kernels(cellfun(@isempty, kernels)) = [];
    end
               
    for i=1:numel(kernel_range)
        kernel_type = kernels{i};

        if(~strcmp(kernel_type, 'linear'))
            if(~exist('rf_obj', 'var') || ~(numel(rf_obj) >= i))               
                rf_obj{i} = InitExplicitKernel( kernel_type, gamma{i}, numel(kernel_range{i}), N_OUT_DIMS{i});
            end
            new_feats{i} = rf_featurize(rf_obj{i}, Feats(kernel_range{i},:)')';
            
            %if(~isempty(RF_scaling))
            %    new_feats{i} = scale_data(new_feats{i},scaling_type, RF_scaling{i});
            %else
            %    [new_feats{i}, RF_scaling{i}] = scale_data(new_feats{i},scaling_type);
            %end            
            RF_scaling = []; % not being used right now (but helps a little!)
        else
            rf_obj{i} = [];
            new_feats{i} = Feats(kernel_range{i}, :);
        end
    end
    
    if(~exist('rf_obj', 'var'))
        rf_obj = [];
    end
    
    Feats = [];
    Feats = cell2mat(new_feats');
end