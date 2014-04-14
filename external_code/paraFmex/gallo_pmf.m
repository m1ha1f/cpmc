function [cuts, lambdas, min_possible_lambda] = gallo_pmf(this_dgraph, lambda_edges, lambda_slopes, s, t, lambda_min, lambda_max, first)  
  %
  %
  % lambda_min and lambda_max are not used as absolute values but instead
  % are instead multiplied by min_possible_lambdas
  
  %the_time = tic();
  assert(~isempty(lambda_edges));
  assert(length(t) == 1);
  assert((lambda_edges(1,1) == s) || (lambda_edges(1,2) == t));

  DISC_FACTOR = 1000; % need to multiply by big number, the algorithm rounds the values..
    
  %BIG_VALUE = 214750000;
  BIG_VALUE = 1000000;
  
  % the algorithm will multiply all capacities by m so that their sum is
  % 2^62, m will determine the precision of the algorithm, or the smallest
  % lambda it can search for.
  
  %  assert(BIG_VALUE<=100000000); % be careful with this!!
  
  this_dgraph = this_dgraph.discretize(DISC_FACTOR);  
  this_dgraph = this_dgraph.replace_bigger_or_equal_than(BIG_VALUE, BIG_VALUE);  
  lambda_slopes = round(lambda_slopes); % * DISC_FACTOR);
  lambda_slopes(lambda_slopes>BIG_VALUE) = BIG_VALUE;
  
  lambda_offsets = this_dgraph.get_edges_weights(lambda_edges);  
  this_dgraph = this_dgraph.set_edges_weights(lambda_edges,0);
  
  [i j] = find(this_dgraph.D);
  normal_weights = this_dgraph.get_weights();

  % remove self edges if they exist  
  self_edges = find(i==j);
  i(self_edges) = [];
  j(self_edges) = [];
  normal_weights(self_edges) = [];
  
  the_edges = [i j; lambda_edges];
  
  % assumes that lambda edges have no offset weight!!!    
  
  assert(BIG_VALUE >= max(normal_weights));
  %cap = [normal_weights; zeros(length(lambda_slopes),1)];
  
  lambda_offsets(:) = 0;
  cap = full([normal_weights; lambda_offsets]);

  % only lambda_edges have slope
  slope = full([zeros(length(i), 1); lambda_slopes]);

  n_edges = length(cap);
  n_nodes = this_dgraph.n_nodes;

  tt = tic();
  if(nargin > 6)
    assert(lambda_max <= intmax('int32'));
    assert(lambda_min >= intmin('int32')); % unsure about this assertion
    assert(lambda_min <= lambda_max);
    % assert(length(unique(the_edges, 'rows')) == size(the_edges,1));
    [y, lambdas, min_possible_lambda] = paraFmex([the_edges cap slope], s, t, n_edges, n_nodes, lambda_min, lambda_max);
  else
    [y, lambdas, min_possible_lambda] = paraFmex([the_edges cap slope], s, t, n_edges, n_nodes);
    %y = paraFmex([the_edges cap slope], s, t, n_edges, n_nodes);
  end
  toc(tt);
  
  x = y(:,2);
  ids = y(:,1);

  MAX_N_BP = 50;
  n_bp = max(x);
  if(n_bp>MAX_N_BP)
    n_bp = MAX_N_BP;
    if(exist('first', 'var') && first)
      disp(['too many breakpoints, ignoring past the ' int2str(MAX_N_BP) 'th']);
      lambdas(MAX_N_BP+1:end) = [];
    else
      disp(['too many breakpoints, ignoring all but the last ' int2str(MAX_N_BP)]);
      lambdas(1:end-MAX_N_BP) = [];
    end
  end

  cuts = zeros(this_dgraph.n_nodes, n_bp);
  cuts(ids(x==1),1) = 1;
  for i=2:n_bp
    cuts(ids(x==i),i) = 1;
    cuts(:,i) = cuts(:,i) + cuts(:,i-1); % they're nested
  end

  cuts(t,:) = zeros(size(cuts,2),1);
  
  % switch (source to 0 sink to 1)
  %  cuts(:,1) = [];
  %  cuts(:,end) = [];
  %  lambdas(end) = [];
  cuts(logical(cuts)) = -1;
  cuts(cuts == 0) = 1;
  cuts(cuts==-1) = 0;

  lambdas = lambdas';
  if(length(lambdas) ~= size(cuts,2))
    disp('the thing ocurred');
    lambdas = 0;
    cuts(s,:) = 0;
    cuts(t,:) = 1;    
  end
  %assert(length(lambdas) == size(cuts,2));

%%%%%%%%% for debugging  
%   for i=1:size(lambda_edges,1)
%     two_nodes = find(this_dgraph.D(lambda_edges(i,2),:))
%     cuts(two_nodes,1)
%     if(cuts(lambda_edges(i,2),1)==0)
%       assert(all(cuts(two_nodes,1) == 0));
%     end
%   end
  %toc(the_time);