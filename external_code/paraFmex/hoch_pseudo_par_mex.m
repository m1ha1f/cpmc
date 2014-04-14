%[cuts, lambdas] = hoch_pseudo_par_mex(N, lambda_edges, lambda_weights, s, t, l, u)  
function [cuts, lambdas] = hoch_pseudo_par_mex(N, lambda_edges, lambda_weights, lambda_offsets, s, t, l, u, real_n_lambdas)   
  % replace very big values (including infinity) and discretize, just
  
  DISC_FACTOR = 1000; % original is 1000
  original_l = 1/DISC_FACTOR;
  original_u = 20;
  N_LAMBDAS = 20;% 20 is fine
  
  if(exist('real_n_lambdas', 'var'))
      N_LAMBDAS = real_n_lambdas;
  end
  
  if(l==-1)
      l = original_l;
  end
  
  if(u==-1)
      u = original_u;
  end
  
  assert(l>=original_l);

  BIG_VALUE = 214750000; % i don't think it can be bigger (beyond certain value goes crazy)  
  
  N = N*DISC_FACTOR;
  %N(N>BIG_VALUE) = BIG_VALUE;
  N(t,:) = 0; % don't want edges out of sink
  
  lambda_weights = (lambda_weights) * DISC_FACTOR;
  lambda_offsets = lambda_offsets * DISC_FACTOR;
  lambda_offsets(lambda_offsets>BIG_VALUE) = BIG_VALUE;
  lambda_weights(lambda_weights>BIG_VALUE) = BIG_VALUE;
  
  [i,j, vals] = find(N);
  subs = [i j];
  
  new_value = BIG_VALUE;
  [i,j] = find(N >= new_value);
  others = subs(vals<new_value,:);
  other_i = others(:,1);
  other_j = others(:,2);
  other_values = vals(vals < new_value);
  N = sparse2([i; other_i],[j; other_j],[new_value*ones(length(i),1); other_values], size(N,1), size(N,2));

  pars.s = s; 
  pars.t = t; 
  pars.lambda_edges = lambda_edges; 
  pars.lambda_weights = lambda_weights;
  pars.lambda_offsets = lambda_offsets;
  pars.lambda_range =  [logspace(log10(l), log10(u), N_LAMBDAS)];
  %pars.lambda_range =  [linspace(l, u, N_LAMBDAS)];
  pars.lambda_range(pars.lambda_range > (u+eps)) = [];
  pars.lambda_range(l > pars.lambda_range) = [];
  pars.lambda_range = [0 pars.lambda_range];
  
  [n_nodes, n_arcs, n_params, s, t, normal_edges, special_edges] = aux_func(N, pars);
  %[ids, lambdas] = hoch_pseudo_par(n_nodes, n_arcs, n_params, s, t, normal_edges, special_edges);
  out = hoch_pseudo_par(n_nodes, n_arcs, n_params, s, t, normal_edges, special_edges);
    
  ids = out(:,1);
  x = out(:,2);
  un = unique(x);
  s_number = min(un);
  t_number = max(un);
  un(un==s_number) = [];
  un(un==t_number) = [];
  n_bp = length(un); % minus source and set
  
  cuts = true(size(N,1), n_bp);  
  cuts(s,:) = false(1,n_bp);
  
  if(isempty(un))
      cuts = [];
      lambdas = [];
      return;
  end
  
  cuts(ids(x==un(1)),:) = 0;
  for i=2:n_bp
    cuts(ids(x==un(i)),i:end) = 0;
  end

  % we only want before-extremal values (the all-but-sink solution isn't
  % interesting)
  if(all(cuts(:,end) == [zeros(size(N,1)-1,1); 1]))
    cuts(:,end) = [];
  end
  
  lambda_ids = un-1;
  if(lambda_ids(1) == 0)
      lambda_ids(1) = [];
      initial_lambda = 0; %pars.lambda_range(lambda_ids(1))-eps;
  else
      initial_lambda = [];
  end
  lambdas = [pars.lambda_range(sort(lambda_ids, 'ascend')) initial_lambda];
  if(length(lambdas) > size(cuts,2)) % hack
      lambdas(end) = [];
  end
end

function  [n_nodes, n_arcs, n_params, s, t, normal_edges, special_edges] = aux_func(N, pars)
    % should either have source or sink parametric edges, didn't think
    % about both!
    s = pars.s; 
    t = pars.t; 
    
    lambda_edges = pars.lambda_edges;
    lambda_weights = pars.lambda_weights;
    
    %assert(all(lambda_edges(:,1) == s) || all(lambda_edges(:,2) == t));
    
    %n_non_lambda_edges = size(N.vals,1);
    %n_non_lambda_edges = numel(find(N));
    
    n_non_lambda_edges = numel(nonzeros(N));
    n_lambda_edges = size(lambda_edges,1);

    n_nodes = size(N,1);
    n_arcs = n_lambda_edges + n_non_lambda_edges;

    %%%%%%
    n_params = size(pars.lambda_range,2);
    %%%%%%
    
    % a < source > < to-node > < capacity1 > < capacity2 > < capacity3 > ... < capacityk >  
    % a < from-node > < sink > < capacity1 > < capacity2 > < capacity3 > ... < capacityk >

    lambda_weights = lambda_weights(:, ones(1,length( pars.lambda_range)));        
        
    s_lambda = (lambda_edges(:,1) == s);
    source_range = sort(pars.lambda_range, 2, 'ascend');   
    source_lambda_weights = lambda_weights(s_lambda,:);    
    for i=1:length(source_range)
        source_lambda_weights(:,i) = source_lambda_weights(:,i)*source_range(i);
    end
        
    source_lambda_weights = source_lambda_weights + pars.lambda_offsets(s_lambda, ones(1, size(source_lambda_weights,2)));
    %source_lambda_weights = source_lambda_weights + repmat(pars.lambda_offsets(s_lambda), 1, size(source_lambda_weights,2));
    
    t_lambda = (lambda_edges(:,2) == t);
    sink_range = sort(pars.lambda_range, 2, 'descend');    
    sink_lambda_weights = lambda_weights(t_lambda,:);    
    for i=1:length(sink_range)
        sink_lambda_weights(:,i) = sink_lambda_weights(:,i)*sink_range(i);
    end
        
    sink_lambda_weights = sink_lambda_weights + pars.lambda_offsets(t_lambda, ones(1, size(sink_lambda_weights,2)));
    %sink_lambda_weights = sink_lambda_weights + repmat( pars.lambda_offsets(t_lambda), 1, size(sink_lambda_weights,2));
     
    source_lambda_edges_rows = full([ lambda_edges(s_lambda,:) source_lambda_weights]);
    sink_lambda_edges_rows = full([ lambda_edges(t_lambda,:)  sink_lambda_weights]);

    duh_s = N(s,:);
    [duh, source_edges, weights] = find(duh_s);      
    source_edges = source_edges';
    weights = weights';
    if(~isempty(source_edges))
        %weights = duh_s.vals;
        source_edges_rows = [s*ones(length(source_edges),1) source_edges weights(:,ones(1,length(pars.lambda_range)))];    
    else
        source_edges_rows = [];
    end
    
    duh_t = N(:,t);    
    [sink_edges, duh, weights] = find(duh_t);
    
    if(~isempty(sink_edges))        
        %weights = duh_t.vals;
        sink_edges_rows = [sink_edges t*ones(length(sink_edges),1) weights(:,ones(1,length(pars.lambda_range)))];    
    else
        sink_edges_rows = [];
    end
    
    %%%%%%%%
    %special_edges = round([source_lambda_edges_rows; source_edges_rows; sink_lambda_edges_rows; sink_edges_rows]);
    special_edges = [source_lambda_edges_rows; source_edges_rows; sink_lambda_edges_rows; sink_edges_rows];
    %%%%%%%%
    
    N = N';
    N(:,s) = 0;
    N = N';
    N(:,t) = 0;
    
    % a < from-node > < to-node > < capacity >       
    [i, j, vals] = find(N);
    assert(all(i - j ~= 0)); % don't want self edges nowhere!!!    
    intermed_edges = [i j];
    new_vals = vals(intermed_edges(:,2) ~= t);
    intermed_edges(intermed_edges(:,2) == t, :) = [];
    %%%%%%
    normal_edges = full([intermed_edges new_vals]);
    %%%%%%
    
    
    assert((size(special_edges,1) + size(normal_edges,1)) == n_arcs);
end

