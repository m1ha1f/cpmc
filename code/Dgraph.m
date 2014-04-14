% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

% directed graph class
%
classdef Dgraph 
  properties (SetAccess = public, GetAccess = public)
    D % the directed graph matrix
    n_nodes
  end % properties  
  methods
    % Constructor
    
    function obj = Dgraph(input)
		% Dgraph(input)
		%
		% input: one of the following
		%	a Dgraph object - creates graph with no edges and obj.n_nodes nodes
		%   D - adjacency matrix
      if(length(input) == 1) % if it's n_nodes
        obj.n_nodes = input;
        obj.D = sparse2(obj.n_nodes, obj.n_nodes);
        %obj.D = sparse(obj.n_nodes, obj.n_nodes);
      else % if it's the adjacency matrix
        assert(size(input,1) == size(input,2));
        obj.D = input;
        obj.n_nodes = size(input,1);
      end
    end        

    %%%%%%%%%%%%%%%%%%%% plotting code %%%%%%%%%%%%%%%%%%%%%%    

    function h = plot(obj, names)
      % adds a return edge for all edges, if they don't have it already
      % this is to make it easy to update
      
      %obj = obj.set_bidirectional_with_epsilon();
      assert(100>obj.n_nodes); % don't visualize if it's too big!
      
      %type = 'Equilibrium';
      type = 'Hierarchical';
      %type = 'Radial';
      if(nargin == 2)
        b = biograph(obj.D, names, 'ShowWeights','on', 'LayoutType', type, 'LayoutScale', 1.5);
      elseif(nargin == 1)
        b = biograph(obj.D, [], 'ShowWeights','on', 'LayoutType', type, 'LayoutScale', 1.5);
      end
      
      for i=1:length(obj.n_nodes)
        b.Nodes(i).UserData = i; % id
      end
      
      for i=1:length(b.Edges) 
        nodeIDS = {};
        for j=1:length(b.Nodes)
          nodeIDS{j} = b.Nodes(j).ID;
        end

        numeric_ids = regexp(b.Edges(i).ID, '\s* -> \s*', 'split');
        for k=1:length(nodeIDS)
          if(strcmp(numeric_ids(1), nodeIDS{k}))
            id1 = k;
          end
          if(strcmp(numeric_ids(2), nodeIDS{k}))
            id2 = k;
          end
        end       
        
        b.Edges(i).UserData = [id1 id2];
        
        if(b.Edges(i).Weight == eps) % don't show the epsilons (show 0)
          b.Edges(i).Weight = 0;
        end
      end
      
      h = view(b);
    end

    function h = plot_path(obj,h, path)
      assert(length(path)>1);
      path_color = [0 1 0];
      set(h.Nodes(path), 'Color', path_color);
      for i=1:length(path)-1
        edge_id = getedgesbynodeid(h, {h.Nodes(path(i)).ID},{h.Nodes(path(i+1)).ID});
        set(edge_id, 'LineColor', path_color);
      end      
    end
    
    function h = plot_color_nodes(obj, h, nodes_to_color, color)
      if(nargin == 3)
        color = [0.0 0.0 1.0];
      end

      set(h.Nodes(nodes_to_color), 'Color', color);     
    end

    function h = plot_reset_colors(obj,h)
      default_node_color = [1.0000    1.0000    0.7000];
      default_edge_color = [0.5000 0.5000 0.5000];
      
      for i=1:length(h.Nodes)        
        set(h.Nodes(i), 'Color', default_node_color);
      end
      for i=1:length(h.Edges)
        set(h.Edges(i), 'LineColor', default_edge_color);
      end
    end
        
    function plot_laplacian_image(obj, img_sz)
      d = sum(obj.D);
      I = full(reshape(d, img_sz(1), img_sz(2)));
      imshow(I);
    end
    
    function h = update_plot_edges(obj, h)      
      for i=1:length(h.Edges)
        id1 = h.Edges(i).UserData(1);
        id2 = h.Edges(i).UserData(2);
              
        weight = obj.D(id1, id2);
        set(h.Edges(i), 'Weight', weight);
      end      
      h = view(h);
    end

    function h = update_plot_ids(obj, h, names)
      for i=1:length(names)
        set(h.Nodes(i), 'ID', names{i});
      end
    end    

    %%%%%%%%%%%%%%%%%%%%%%%%% End of plotting code %%%%%%%%%%%%%%%%%%%%%%%%
        
    %
    % sort of overloaded matlab operators
    %

    function b = plus(obj, a)
      if(isempty(a))
        b = obj;
        return;
      end
      
      assert(strcmp(class(a), 'Dgraph'));
      new_side = max(obj.n_nodes, a.n_nodes);
      b = Dgraph(new_side);
      
      if(obj.n_nodes ~= new_side) % this one is the smallest
        obj = obj.add_size(new_side);
      elseif (a.n_nodes ~= new_side) % else if this one is the smallest
        a = a.add_size(new_side);
      end % else they're the same size
      
      b.D = obj.D + a.D;
    end
 
    % functionality
    
    function obj = remove_nodes(obj, nodes)
      remaining_nodes = setdiff(1:obj.n_nodes, nodes);

      newD = obj.D(remaining_nodes, remaining_nodes);
      obj.D = newD;

      obj.n_nodes = size(obj.D,1);
    end
     
    function obj = add_edges(obj, node_ids, values)      
      new_d = sparse2(node_ids(:,1), node_ids(:,2), values, ...
                      obj.n_nodes, obj.n_nodes);
     
%      new_d = sparse(node_ids(:,1), node_ids(:,2), values, ...
%                      obj.n_nodes, obj.n_nodes);     
                  
      if(isempty(obj.D))
        obj.D = new_d;
      else
        obj.D = obj.D + new_d;
      end
      %toc(t)
    end
        
    function obj = add_size(obj, new_side)
      % i'm assuming it's always square, therefore only one parameter
      assert(new_side > obj.n_nodes);
      obj.n_nodes = new_side;     


      %%% this used to work, now must be sparse
      %n_to_add = new_side - obj.n_nodes;      
      %obj.D = [obj.D zeros(new_side - n_to_add, n_to_add)];
      %obj.D = [obj.D; zeros(n_to_add, new_side)];

      [i,j,s] = find(obj.D);
      obj.D = sparse2(i,j, s, new_side, new_side);      
      %obj.D = sparse(i,j, s, new_side, new_side);      
    end

    function [obj, permutation] = permute_a_b(obj, ids, target_ids)    
      %puts elements in positions 'ids' in positions 'target_ids'
      %and arranges the elements in target_ids to the empty positions
      permutation = 1:obj.n_nodes;

      dislodged = setdiff(target_ids, ids); % target_ids not in ids
      left_positions = setdiff(ids,target_ids); % ids not in target_ids
      permutation(target_ids) = ids;
      permutation(left_positions) = dislodged;
            
      obj = obj.permute(permutation);
    end

    function obj = permute(obj, perm)
      % permutes rows and columns of the matrix representation of the
      % digraph.
      
      temp_D = obj.D;
      temp_D = obj.D(:,perm);
      obj.D = temp_D;
      temp_D = obj.D(perm,:);
      
      obj.D = temp_D;
    end
    
    function array = convert_to_array(obj)
      array = obj.D;
    end
    
    function [this_dgraph, new_ids_this, new_ids_other] = combine_with(this_dgraph, other_dgraph, map_ids_this, map_ids_other)
      % map_ids_this: ids in this digraph that map to ids map_ids_other in
      % digraph other_dgraph
      
      % there are no multivalued mapping (every
      % node in this_dgraph maps to at most one node in other_dgraph)
      assert(length(unique(map_ids_this)) == length(map_ids_this));
      assert(length(unique(map_ids_other)) == length(map_ids_other));      
      
      % 1. Identify unique nodes in the digraphs
      ids_this = 1:this_dgraph.n_nodes;
      unique_this = setdiff(ids_this, map_ids_this);
      number_unique_this = length(unique_this);      
      
      ids_other = 1:other_dgraph.n_nodes;
      unique_other = setdiff(ids_other, map_ids_other);
      number_unique_other = length(unique_other);
      
      % 2. Grow both digraphs
      ids_grown = 1:(number_unique_this + length(ids_other));
      other_dgraph = other_dgraph.add_size(length(ids_grown));    
      this_dgraph = this_dgraph.add_size(length(ids_grown));      
      
      assert((number_unique_this + length(ids_other)) == ...
             (number_unique_other + length(ids_this)));
           
      % 3. Create a complete map
      final_map_ids_other = zeros(length(ids_grown),1);
      nulls = (length(ids_other)+1):ids_grown(end);
      final_map_ids_other(map_ids_this) = [map_ids_other]; 
      final_map_ids_other(unique_this) = nulls;
      final_map_ids_other(final_map_ids_other==0) = unique_other;
      
      % 4. permute other_dgraph
      other_dgraph = other_dgraph.permute(final_map_ids_other);
      
      % 5. sum together
      this_dgraph = this_dgraph + other_dgraph;
      
      % 6. return useful info
      new_ids_this = ids_this; % they're the same
      [dummy, new_ids_other] = sort(final_map_ids_other);  
      new_ids_other = new_ids_other(1:length(ids_other))';
    end

    function obj = discretize(obj, disc_factor)
      obj.D = round(obj.D*disc_factor);     
    end
    
    function obj = replace_bigger_or_equal_than(obj,threshold, new_value)        
        obj.D(obj.D>=threshold) = new_value;
%       N = obj.D;
%       spN = sptensor(N);      
%       vals = spN.vals;
%       subs = spN.subs;
%       [i] = find(spN >= threshold);
%       if(length(i)>0)
%         j = i(:,2);
%         i = i(:,1);
%         others = subs(vals<threshold,:);
%         other_i = others(:,1);
%         other_j = others(:,2);
%         other_values = vals(vals < threshold);
%         obj.D = sparse2([i; other_i],[j; other_j],[new_value*ones(length(i),1); other_values], size(N,1), size(N,2)); 
%         %obj.D = sparse([i; other_i],[j; other_j],[new_value*ones(length(i),1); other_values], size(N,1), size(N,2)); 
%       end
    end
    
    function obj = zero_diagonal(obj)
      N = obj.D;
      N(diag(N)) = 0;
      obj.D = N;
    end
    
    function truth = is_symmetric(obj)
      if(isequal(obj.D,obj.D'))
        truth = true;
      else
        truth = false;
      end
    end
    
    function obj = symmetrize(obj)     
      % symmetrizes a dgraph. If there's already bidirectional edges, and
      % they have different values I don't know what the precise output is,
      % you think about it
      error('should remove sptensor');
      N = obj.D;
      spN = sptensor(N);      
      vals = spN.vals;
      subs = spN.subs;
      subs = [subs; subs(:,2) subs(:,1)];
      vals = [vals; vals];
      
      [subs, rows] = unique(subs, 'rows');
      vals = vals(rows);
      
      obj.D = sparse2(subs(:,1),subs(:,2),vals, size(N,1), size(N,2));   
      %obj.D = sparse(subs(:,1),subs(:,2),vals, size(N,1), size(N,2));    
    end    
    
    function obj = symmetrize_except_inf(obj)    
      % same as previous function but leaves Infinite edges alone
      
      
      error('should remove sptensor');
      N = obj.D;
      spN = sptensor(N);      
      vals = spN.vals;
      subs = spN.subs;
      
      subs_without_inf = subs(vals~=Inf,:);
      vals_without_inf = vals(vals~=Inf);
      
      subs = [subs; subs_without_inf(:,2) subs_without_inf(:,1)];
      vals = [vals; vals_without_inf];
      
      [subs, rows] = unique(subs, 'rows');
      vals = vals(rows);
      
      obj.D = sparse2(subs(:,1),subs(:,2),vals, size(N,1), size(N,2));       
      %obj.D = sparse(subs(:,1),subs(:,2),vals, size(N,1), size(N,2));    
    end    

    function obj = set_bidirectional_with_epsilon(obj, eps_substitute)
      % assures that all connections become bidirectional, with epsilon
      % weights where there was none.
      if(nargin == 1)
        eps_substitute = eps;
      end
      N = obj.D;
      spN = sptensor(N);      
      vals = spN.vals;
      subs = spN.subs;
      subs = [subs; subs(:,2) subs(:,1)];
      vals = [vals; eps_substitute*ones(length(vals),1)];
      
      [subs, i, j] = unique(subs, 'rows', 'first');
      vals = vals(i);
      
      obj.D = sparse2(subs(:,1),subs(:,2),vals, size(N,1), size(N,2));       
      %obj.D = sparse(subs(:,1),subs(:,2),vals, size(N,1), size(N,2));       
    end      
    
    function weights = get_edges_weights(obj, the_edges)
      un1 = unique(the_edges(:,1));
      if(numel(un1) == 1)
        weights = obj.D(un1, the_edges(:,2))'; % faster
        return;
      else
        un2 = unique(the_edges(:,2));
        if(numel(un2) == 1)
          weights = obj.D(the_edges(:,1), un2);
          return;
        end
      end
      
      [I J s]=find(obj.D);
      [tf loc]= ismember([the_edges(:,1) the_edges(:,2)],[I J],'rows');
      weights = zeros(size(tf));
      weights(tf) = s(loc(tf));

      
%       un1 = unique(the_edges(:,1));
%       if(numel(un1) == 1)
%           weights = obj.D(un1, the_edges(:,2))'; % faster
%           return;
%       else
%           un2 = unique(the_edges(:,2));
%           if(numel(un2) == 1)
%             weights = obj.D(the_edges(:,1), un2);
%             return;
%           end
%       end
%                   
%       linear_ids= sub2ind(size(obj.D), the_edges(:,1), the_edges(:,2));
%       weights = obj.D(linear_ids);
    end
    
    function weights = get_weights(obj)
        
      error('using sptensor, should remove!');
      the_d = sptensor(obj.D);
      
      weights = the_d.vals;
    end

    function obj = set_submatrix(obj, range, other_dgraph)          
      assert(length(range) == other_dgraph.n_nodes); % only considering "square" graphs
      if(all(size(other_dgraph.D) == size(obj.D)))
          obj.D = other_dgraph.D;
          return;
      end
           
      [i,j,weights2] = find(other_dgraph.D);
      ids2 = [i j];
      
      if(isempty(ids2))
          return;
      end
      
      error('using sptensor, should remove!');
      altered_ids2 = ids2 + min(range)-1;
      
      t = tic();
      D1 = sptensor(obj.D);
      D1([altered_ids2]) = weights2;
      [coords] = find(D1);
      vals = D1.vals;         
      obj.D = sparse2(coords(:,1), coords(:,2),vals, size(D1,1), size(D1,2));
      toc(t)
    end
    
    function obj = set_edges_weights(obj, the_edges, weights)
        
      error('should remove sptensor');
      D = sptensor(obj.D);
      D(the_edges) = weights;
      [coords] = find(D);
      vals = D.vals;
      obj.D = sparse2(coords(:,1), coords(:,2),vals, size(D,1), size(D,2));
      %obj.D = sparse(coords(:,1), coords(:,2),vals, size(D,1), size(D,2));
    end
    
    function lambda_offsets = compute_lambda_offsets(obj,s,t,lambda_nodes)
      D = obj.D;
      disp('offsets are a little buggy');
      lambda_offsets.s_side = D(s,lambda_nodes.t_side);
      lambda_offsets.t_side = D(lambda_nodes.s_side,t);
    end
    
    function lambda_nodes = compute_lambda_nodes(obj, s, t, lambda_edges)
      % assumes the normal situation with lambda_nodes being connected only
      % to source or sink, if not it won't work
      if(size(lambda_edges,2)>size(lambda_edges,1))
        lambda_edges = lambda_edges';
      end
      if(any(lambda_edges(:,1) == s) || any(lambda_edges(:,1) == t))
        lambda_nodes = lambda_edges(:,2);
      else
        lambda_nodes = lambda_edges(:,1);
      end
    end

    function [obj, ids_map_orig_to_new, contraction_node] = contract(obj, nodeset)
      % contracts a set of nodes
      % new node gets the node with highest id in the graph
      % ids_map_orig_to_new gives the new ids of each node
      if(size(nodeset,2) < size(nodeset,1))
        % nodeset should be a row
        nodeset = nodeset';
      end
      
      assert(length(unique(nodeset))==length(nodeset));
      
      G = obj.D;
      ids = 1:length(G);

      obj.n_nodes = obj.n_nodes - length(nodeset) + 1;
      other_nodes = setdiff(ids, nodeset);

      duh = G';
      new_row = full(sum(duh(:,nodeset),2))';
      %new_row = full(sum(G(nodeset,:),1));
      new_row(nodeset) = [];
      new_row = [new_row 0]; % add a space for itself
      new_col = full(sum(G(:,nodeset),2));
      new_col(nodeset) = []; 
      new_col = [new_col; 0]; % add a space for itself
      
      %t = tic();
      newG = G(:, [other_nodes 1]);
      newG = newG';
      newG = newG(:, [other_nodes 1]);
      [i, j, vals] = find(newG);
      coords = [i j];
      %toc(t);      
      
      G = sparse2(coords(:,1), coords(:,2),vals, numel(other_nodes)+1, numel(other_nodes)+1);
      G(end,:) = new_row;
      G(:,end) = new_col;
      
      ids_map_orig_to_new = zeros(length(ids),1);
      ids_map_orig_to_new(nodeset) = size(G,1);
      ids_to_change = other_nodes';
      
      %nodeset_replicated = repmat(nodeset,length(ids_to_change), 1);
      %ids_to_change_replicated = repmat(ids_to_change, 1, length(nodeset));
      %n_ids_in_nodeset_before_id = sum(nodeset_replicated<ids_to_change_replicated,2);
      %ids_map_orig_to_new(ids_to_change) = ids_to_change - n_ids_in_nodeset_before_id;
      
      ids_map_orig_to_new(ids_to_change) = 1:length(ids_to_change);
      obj.D = G;     
      contraction_node = length(obj.D);
    end
    
    function [new_obj, new_lambda_edges, new_lambda_weights, new_lambda_offsets, ids_map_orig_to_new, ids_map_new_to_orig, contraction_node] = contract_with_lambda(obj,nodeset, lambda_edges, lambda_weights, ids_map_orig_to_new)
      if(length(nodeset) == 1)         
        new_obj = obj;
        new_lambda_edges = lambda_edges;
        new_lambda_weights = lambda_weights;
        contraction_node = nodeset;
        ids_map_new_to_orig = [];
        error('not ready for this');
        ids_map_orig_to_new = 1000; % see later
        return;
      end
      
      % the contraction
      [new_obj, ids_map_orig_to_new2, contraction_node] = obj.contract(nodeset);
            
      ids_map_orig_to_new = ids_map_orig_to_new2(ids_map_orig_to_new);
      assert(contraction_node == new_obj.n_nodes); % we also require this for the simplification three lines below!            
      [ids_sorted,ids_map_new_to_orig] = sort(ids_map_orig_to_new);
      ids_map_new_to_orig(end-length(nodeset)+1:end) = []; % (using the simplification allowed by the previous assertion
      ids_map_new_to_orig = [ids_map_new_to_orig; -1]; % -1 means it's mapping to more than one element
      
      % Take care of lambda edges      
      new_lambda_edges = ids_map_orig_to_new(lambda_edges);  
      if(size(new_lambda_edges,2) == 1)
        new_lambda_edges = [];
        new_lambda_weights = [];
        return;
      end
      
      [new_lambda_edge_dest, j] = sort(new_lambda_edges(:,2));
      new_lambda_weights = lambda_weights(j);
      
      new_lambda_weights = accumarray(new_lambda_edge_dest, new_lambda_weights);
      the_zeros = (new_lambda_weights==0);
      new_lambda_weights(the_zeros) = [];
      new_lambda_edge_dest = unique(new_lambda_edge_dest);
      new_lambda_edges = [new_lambda_edges(1,1)*ones(length(new_lambda_edge_dest),1) new_lambda_edge_dest];
      if(new_lambda_edges(end,1) == new_lambda_edges(end,2))
        new_lambda_edges(end,:) = [];
        new_lambda_weights(end) = [];
      end
            
      % in our current problem, we don't lambda offsets, change it later if
      % you need more generality
      new_lambda_offsets = zeros(size(new_lambda_weights));
      
      % remove all lambda edges into the contracted node (he should be
      % connected to the sink with infinity weights)
      to_remove = (new_lambda_edges(:,2)==contraction_node);
      new_lambda_edges(to_remove,:) = [];
      new_lambda_weights(to_remove,:) = [];
      new_lambda_offsets(to_remove,:) = [];
      assert(length(new_lambda_weights) == size(new_lambda_edges,1));    
      assert(length(nodeset) + size(new_obj.D,1) == size(obj.D,1) + 1);
    end

    function [new_obj, new_lambda_edges, new_lambda_weights, new_lambda_offsets, ids_map_orig_to_new, ids_map_new_to_orig, contraction_node] = contract_to_source_with_lambda(obj, nodeset, lambda_edges, lambda_weights, lambda_offsets, ids_map_orig_to_new)
      %%%% need to pass source node inside nodeset!! %%%%%%%%%%%%%
      % assumes all lambda edges originate from the same node
      assert(length(unique(lambda_edges(:,1))) == 1);

      if(length(nodeset) == 1)
        new_obj = obj;
        new_lambda_edges = lambda_edges;
        new_lambda_weights = lambda_weights;
        contraction_node = nodeset;
        ids_map_new_to_orig = [];
        return;
      end
      
      % the contraction
      [new_obj, ids_map_orig_to_new2, contraction_node] = obj.contract(nodeset);
      
      ids_map_orig_to_new = ids_map_orig_to_new2(ids_map_orig_to_new);
      assert(contraction_node == new_obj.n_nodes); % we also require this for the simplification three lines below!
            
      [ids_sorted,ids_map_new_to_orig] = sort(ids_map_orig_to_new);
      ids_map_new_to_orig(end-length(nodeset)+1:end) = []; % (using the simplification allowed by the previous assertion
      ids_map_new_to_orig = [ids_map_new_to_orig; -1]; % -1 means it's mapping to more than one element
      
      
      % Take care of lambda edges
      new_lambda_edges = ids_map_orig_to_new2(lambda_edges);  
      if(size(new_lambda_edges,2) == 1)
        new_lambda_edges = [];
        new_lambda_weights = [];
        return;
      end
      
      [new_lambda_edge_dest, j] = sort(new_lambda_edges(:,2));
      new_lambda_weights = lambda_weights(j);
      
      new_lambda_weights = accumarray(new_lambda_edge_dest, new_lambda_weights);
      the_zeros = (new_lambda_weights==0);
      new_lambda_weights(the_zeros) = [];
      new_lambda_edge_dest = unique(new_lambda_edge_dest);
      new_lambda_edges = [new_lambda_edges(1,1)*ones(length(new_lambda_edge_dest),1) new_lambda_edge_dest];
      if(new_lambda_edges(end,1) == new_lambda_edges(end,2))
        new_lambda_edges(end,:) = [];
        new_lambda_weights(end) = [];
      end
      
      assert(length(new_lambda_weights) == size(new_lambda_edges,1));
      
      % take care of lambda_offsets
      new_lambda_offsets = lambda_offsets;
      new_lambda_offsets(ids_map_orig_to_new(1:end-1)==contraction_node) = [];
      ids = find(new_obj.D(contraction_node,:));
      new_lambda_offsets(ids) = new_obj.D(contraction_node,ids);
      new_obj.D(contraction_node,:) =0;
      
      if(~isempty(new_lambda_edges))
        assert(new_lambda_edges(1,1) == contraction_node);
      end
    end
     
    function [new_obj, new_lambda_edges, new_lambda_weights, new_lambda_offsets, ids_map_orig_to_new, ids_map_new_to_orig, contraction_node] = contract_to_sink_remove_lambdas(obj, nodeset, lambda_edges, lambda_weights, lambda_offsets, ids_map_orig_to_new)
      assert(length(unique(lambda_edges(:,1))) == 1);

      if(length(nodeset) == 1)
        new_obj = obj;
        new_lambda_edges = lambda_edges;
        new_lambda_weights = lambda_weights;
        contraction_node = nodeset;
        ids_map_new_to_orig = [];
        return;
      end
      
      % the contraction
      [new_obj, ids_map_orig_to_new2, contraction_node] = obj.contract(nodeset);
      
      ids_map_orig_to_new = ids_map_orig_to_new2(ids_map_orig_to_new);
      assert(contraction_node == new_obj.n_nodes); % we also require this for the simplification three lines below!
            
      [ids_sorted,ids_map_new_to_orig] = sort(ids_map_orig_to_new);
      ids_map_new_to_orig(end-length(nodeset)+1:end) = []; % (using the simplification allowed by the previous assertion
      ids_map_new_to_orig = [ids_map_new_to_orig; -1]; % -1 means it's mapping to more than one element
      
      
      % Take care of lambda edges
      new_lambda_edges = ids_map_orig_to_new2(lambda_edges);  
      if(size(new_lambda_edges,2) == 1)
        new_lambda_edges = [];
        new_lambda_weights = [];
        return;
      end
      
      [new_lambda_edge_dest, j] = sort(new_lambda_edges(:,2));
      new_lambda_weights = lambda_weights(j);
      new_lambda_offsets = lambda_offsets(j);
      
      new_lambda_weights = accumarray(new_lambda_edge_dest, new_lambda_weights);      
      new_lambda_offsets =  accumarray(new_lambda_edge_dest, new_lambda_offsets);
      
      the_zeros = (new_lambda_weights==0);
      new_lambda_weights(the_zeros) = [];
      new_lambda_offsets(the_zeros) = [];
      new_lambda_edge_dest = unique(new_lambda_edge_dest);
      new_lambda_edges = [new_lambda_edges(1,1)*ones(length(new_lambda_edge_dest),1) new_lambda_edge_dest];
      if(new_lambda_edges(end,1) == new_lambda_edges(end,2))
        new_lambda_edges(end,:) = [];
        new_lambda_weights(end) = [];
      end
      
      assert(length(new_lambda_weights) == size(new_lambda_edges,1));
      
      % take care of lambda_offsets
%       new_lambda_offsets(ids_map_orig_to_new(1:end-1)==ids_map_orig_to_new(source)) = [];
%       ids = find(new_obj.D(contraction_node,:));
%       new_lambda_offsets(ids) = new_obj.D(contraction_node,ids);
%       new_obj.D(contraction_node,:) =0;
             
      % now remove lambda edges between source and contracted node
      ids_to_remove = find(new_lambda_edges(:,2) == contraction_node);
      new_lambda_edges(ids_to_remove,:) = [];
      new_lambda_weights(ids_to_remove,:) = [];
      new_lambda_offsets(ids_to_remove,:) = [];
      new_obj.D(contraction_node,:) = 0; % remove all edges out of here (it's gonna be the sink!)
    end

    function [linkage_value] = linkage(obj, nodeset1, nodeset2)
      % returns the sum of weights of edges from nodeset1 to nodeset2
      % these sets can overlap  
      
      if((length(nodeset1) == 1) && (length(nodeset2) == 1) && (nodeset1 == nodeset2))
        linkage_value = 0;
        return;
      end      
      
      all_conn = obj.D(nodeset1, nodeset2);
      linkage_value = sum(sum(all_conn));      
    end
        
    function truth = is_hamiltonian(obj, path)
      if(length(path) ~= obj.n_nodes)
        truth = false;
      else
        truth = is_hamiltonian_of_subgraph(obj, path);
      end
    end
    
    function truth = is_hamiltonian_of_subgraph(obj, path)
      % is true if path has no repeated nodes, and it's a path ( there are
      % links between the consecutive nodes )
      if(length(unique(path)) ~= length(path))
        truth = false;
        return;
      end
      
      for i=1:length(path)-1
        if(obj.D(path(i),path(i+1)) <= 0)
          truth = false;
          return;
        end
      end
      
      truth = true;
    end

    function distances = get_distances_to_node(obj, node)      
      D = obj.D;
      [d pred] = dijkstra_sp(D,node);
      distances = d;
    end
  end

  %
  %
  % cut functionality
  %
  %
  methods  
    function [cut, flow] = min_st_cut(obj, s, t, alg)
      % 1 in cut is sink, 0 is source
      
      if(nargin == 3)        
        alg = 'push_relabel'; % by default
      end
      
      inf_replacement = double(intmax('int32'))*0.1;
      h = tic;     
      if(strcmp(alg, 'kolmogorov'))
        %disp(['Running Kolmogorov''s Min Cut... The cut is correct, ' ...
        %      'but sometimes the flow is wrong: it seems like it''s ' ...
        %      'subtracting the capacities on the direction t--s']);
        [flow, cut] = obj.kolmogorov_min_cut(s,t);       
      elseif(strcmp(alg, 'push_relabel')) 
        % there's a problem with the precision
        MULT_FACTOR = 1000;
        obj.D = obj.D*MULT_FACTOR;
        obj = obj.replace_bigger_or_equal_than(inf_replacement,inf_replacement);
        
        [flow, cut] = max_flow(obj.D, s, t);
        
        flow = flow/MULT_FACTOR;
        cut(cut==-1) = 2; % uniform output over different implementations
        cut = cut - 1;
      else
        error('no such algorithm implemented, only push_relabel and kolmogorov');
      end      
      assert(all(cut) >= 0);
      toc(h)
    end
    
    function [cuts, lambdas] = parametric_min_st_cut(obj, lambda_edges, lambda_weights, s, t, l, u, alg)
        assert(~isempty(t));
        assert(~isempty(s));
                
        if(nargin == 8)
            the_alg = alg;
        else
            the_alg = 'gallo';
        end
        
        if((strcmp(the_alg, 'gallo')))
            [cuts, lambdas] = gallo_pmf(obj, lambda_edges, lambda_weights, s, t, l, u, true);
        elseif(strcmp(the_alg, 'hochbaum'))
            %the_time =tic();
            lambda_offsets = obj.get_edges_weights(lambda_edges); 
            %toc(the_time)
            [cuts,lambdas] = hoch_pseudo_par_mex(obj.D, lambda_edges, lambda_weights, lambda_offsets, s, t, l, u);
        elseif(strcmp(the_alg, 'bin_search'))
            assert('fix this');
            cuts = bin_search_par_min_st_cut(obj, s,t, parametric_pars, visualiz);
        else
            error('no such algorithm, only gallo and hochbaum available');
        end
    end

    function [cut, flow] = min_global_cut_undirected(obj)
      if(size(obj.D,1) < 2)
        cut = [];
        flow = [];
        return;
      end
      
      % uses hao & orlin's algorithm for undirected graphs, an executable
      % called ho.
      
      pars = [];
      BIG_VALUE = 1000000;      

      obj2 = obj.replace_bigger_or_equal_than(Inf, BIG_VALUE);
            
      filename = network2dimacs(obj2.D,'global_min_cut', pars);

      [cut, flow] = global_min_cut(filename, obj2.n_nodes, 'mincutlib')

      flow = flow
    end
    
    function [cut, flow] = min_global_cut(obj,source, sinks)
      if(size(obj.D,1) < 2)
        cut = [];
        flow = [];
        return;
      end
      
      % uses hao & orlin's algorithm for directed graphs, an executable
      % called gdmincut, or gmindcut_with_source
      % (must be connected, e.g. any nodes should be accessible from any other)
      % source is optional
      pars = [];  
      BIG_VALUE = 1000000;   

      obj2 = obj.replace_bigger_or_equal_than(Inf, BIG_VALUE);
      
      if(nargin==3) 
        pars.s = source;
        pars.t = sinks;       
        
        % verify if sinks form a hamiltonian path on the subgraph induced
        % by the nodes of the path
        assert(obj.is_hamiltonian_of_subgraph(sinks));
        
        filename = network2dimacs(obj2.D,'global_min_cut_dir_with_source_and_sinks', pars);

        [cut, flow] = global_min_cut_with_source_and_sinks(filename, obj2.n_nodes);
              
%         names{1} = '1';
%         names{2} = '2';
%         names{3} = '3';
%         names{4} = '4';
%         names{5} = 't-NCUT 1 2';
%         names{6} = 't-NCUT 1 3';
%         names{7} = 't-NCUT 2 3';
%         names{8} = 't-NCUT 1 4';      
%         names{9} = 't-NCUT 2 4';      
%         names{10} = 't-NCUT 3 4';
%         names{11} = 's';
%         obj2.plot_cut_with_names(cut,names);
%         title(sprintf('flow: %f', flow));
%         disp('one more');
      elseif(nargin==2)
        pars.s = source;
        
        [cut, flow] = gmindcut_mex(obj2.D, source);        
      elseif(nargin==1)
        filename = network2dimacs(obj2.D,'global_min_cut_dir', pars);        
        [cut, flow] = global_min_cut(filename, obj2.n_nodes);
        % let the smaller segment on the source side
        if(sum(cut==1) < sum(cut==0))
          cut = ~cut;
        end
      end
      flow = flow
      %value1 = Partitioner.measure_st_cut(obj2, cut)
      %cut_optimal = [1 0 1 1 0 1 0 1 0 1 1];
      %value2 = Partitioner.measure_st_cut(obj2, cut_optimal)  
      %obj2.plot_cut_with_names(cut_optimal,names);
    end          

    function [cut, flow, the_lambda] = bin_search_par_min_global_cut_s_t(obj, s, t, l, u, lambda_nodes, synth_nodes, source_or_sink, max_or_min, visualiz)
      counter = 0;
      if(max_or_min == 'min')
        the_lambda = inf;
      else
        the_lambda = 0;
      end
      
      previous_lambda = 0;
      n_iters = 0;
      %while (u-l) >= (1/(n*(n-1)))
      cut = zeros(size(obj.D,1),1);
      flow = -inf;
 
      if(strcmp(source_or_sink, 'source'))
        lambda_edges = [lambda_nodes.pivot*ones(length(lambda_nodes.s_side),1) lambda_nodes.s_side];
      else
        lambda_edges = [lambda_nodes.t_side lambda_nodes.pivot*ones(length(lambda_nodes.t_side),1)];
      end
      lambda_nodes.pivot = [];
      
      lambda_weights = obj.get_edges_weights(lambda_edges);
      
      for i=1:10 % doing no more than 1 iterations
        if(isempty(lambda_edges))
          return;
        end
        
        counter = counter + 1;

        new_lambda = (u+l)/2; % g is the mid value, u is the upper and l the lower
        new_lambda = new_lambda

        if(new_lambda == previous_lambda)
          break;
        else 
          previous_lambda = new_lambda;
          n_iters = n_iters+1
        end

        new_weights = lambda_weights*new_lambda;
        
        [obj] = obj.set_edges_weights(lambda_edges, new_weights(:,1));
        if(ncols(new_weights) == 2)
          [obj] = obj.set_edges_weights([lambda_edges(:,2) lambda_edges(:,1)], new_weights(:,2));
        end
        
        source = s;
        sinks = t;
        
        %obj.sharon_ncut_check(source, sinks);
        %[this_cut, flow] = obj.min_global_cut(source, sinks);   
        [this_cut, flow] = obj.exhaustive_min_global_cut(source, sinks);
        flow = flow;
        %this_cut
        
        %obj.plot_cut(this_cut,lambda_edges(1,1));
        
        if(nargin>=8)
          obj.show_intermediate_cut(this_cut,synth_nodes,visualiz);
        end
        
        if(strcmp(max_or_min, 'min'))
          if(any(this_cut(t)) && new_lambda<the_lambda)
            sum(this_cut(struct2array(lambda_nodes)))
            the_lambda = new_lambda;
            cut = this_cut;

            u = new_lambda;
          else
            l=new_lambda;
          end         
        else
          if(any(this_cut(lambda_nodes)))
            sum(this_cut(t))
            the_lambda = new_lambda;
            cut = this_cut;

            l = new_lambda; % upper = mid
          else % else lower = mid
            u=new_lambda;
          end
        end
      end      
    end
    
    function [cut, flow] = exhaustive_min_global_cut(obj, s, t)
      % computes the minimum directed cut among all s and t
      cut = zeros(obj.n_nodes,1);
      flow = -1;
      for i=1:length(s)
        for j=1:length(t)
          [this_cut,this_flow] = obj.min_st_cut(s(i),t(j));
          if(this_flow > flow)
            cut = this_cut;
            flow = this_flow;
          end
        end
      end      
    end

    function segmentation = recursive_spectral_bipartitioning(obj, criteria, n_segments)  
      % it's not returning the exact number of segments yet, just going
      % into the depth level where this number would be attained
      function segm = bipartition(obj, criteria, max_depth, current_depth)
        [cut] = obj.spectral_cuts(criteria);

        if(current_depth < max_depth)      
          ids_1 = find(cut==1);
          ids_0 = find(cut==0);
          the_dgraph_1 = obj.remove_nodes(ids_0);
          the_dgraph_0 = obj.remove_nodes(ids_1);
          segm_1 = bipartition(the_dgraph_1, criteria, max_depth, current_depth+1);    
          segm_0 = bipartition(the_dgraph_0, criteria, max_depth, current_depth+1);

          un0 = unique(segm_0);
          max_un0 = max(un0);

          assert(length(segm_1) + length(segm_0) == length(cut));
          segm_1 = segm_1 + max_un0 + 1;
          segm = zeros(length(cut),1);
          segm(ids_1) = segm_1;
          segm(ids_0) = segm_0;
        else
          segm = cut;
        end
      end
      max_depth = sqrt(n_segments);
      segmentation = bipartition(obj, criteria, max_depth, 1); 
    end 
    
    function [cut] = spectral_cuts(obj, type, par_name, value)
      %assert(obj.is_symmetric()); 
      threshold = -0.0000000000000005;
      npixels = 1000;
      use_threshold = false;
      use_npixels = false;
      cut = zeros(size(obj.D,1),1);
      
      if(nargin() == 2)
        use_threshold = true;
      elseif(nargin() == 4)
        if(strcmp(par_name, 'Threshold'))
          % threshold to apply
          threshold = value;
        elseif(strcmp(par_name, 'NPixels'))          
          % number of pixels to include in the segment
          use_npixels = true;
          npixels = value;
        elseif(strcmp(par_name, 'NSegments'))
          n_segments = value;
        else
          error('no such parameter');
        end
      else
        error('wrong number of parameters');
      end
      W = obj.D;
      
      np = size(W,1);
      
      reg = 0.000;
      
      %assert(reg~=0);
      d = sum(abs(W),2);
      
      dr = reg*ones(size(d,1),1);

      D = spdiags(d,0,np,np);
      Dr =  D + spdiags(dr,0,np,np);
      WplusD = W + D;

%       Wr = W + spdiags(reg*ones(length(W),1),0,np,np);
%       WDr = W + Dr;
%       WD = W + D;
%       L = D - W; % the unnormalized laplacian     
%       Lr = Dr - W;
      
      options.issym = 1;
      options.disp = 0; % put zero for no diagnostic
      options.tol = 0.001;
      %options.maxit = 20;
      
      if(strcmp(type, 'perona'))
        [E, V] = eigs(W,1,'LM', options);
      elseif(strcmp(type, 'shi'))
        [E, V] = eigs(D - W, D, 2,'SM', options); %,options);
        E = E(:,end-1); % the second smallest
      elseif(strcmp(type, 'sharon_ncut'))
        SIGMA = 0.00001;
        [E, V] = eigs(D - W, W,2,'SM', options); %,options);
        %[E, V] = eigs(D - W, W,2,SIGMA, options); %,options);
        [id_i, id_j] = find(V==max(max(V)));
        E = E(:,id_j);
      elseif(strcmp(type, 'yu'))
        dih = 1./(sqrt(d)+eps);
        Wd = spmd_nonameconflict(W + Dr,dih,dih);
        [E,V] = eigs(Wd,n_segments,'LM',options);
        cut = getbinsol(E);
      else
        error('no such type of spectral cut');
      end
      
      if(use_threshold)
        cut = (E<=-threshold);
        %cut = (E >=-threshold);
      elseif(use_npixels)
        %imshow(reshape(E, 157, 243));
        [sorted_vals, ids] = sort(E,'descend');
        cut(ids(1:npixels)) = 1;
      end
    end
    
    function energy = measure_partition_energy(obj, partition,criteria, bground_label)
      % partition should be a vector with integer labels, one for each node
      % of the directed graph the_dgraph      
      % if bground_label is provided, the cost of the segment with that
      % label won't enter the energy
      assert(all(partition>=0));
      
      if(nargin==3)
        bground_label = -1;
      end
      
      assert(obj.n_nodes == size(partition,1));
      
      labels = unique(partition);
           
      n_segments = length(labels);
      energy = 0;

      % pixels ids for each segment
      segments_to_consider = [];
      for i=1:n_segments
        [k{i}] = find(partition==labels(i))';
        if(labels(i) ~= bground_label)
          segments_to_consider = [segments_to_consider i];
        end
      end
      
      % links between all segments
      assert(obj.is_symmetric()); % assumes symmetric weights
      for i=1:n_segments
        for j=i:n_segments                    
          L(i,j) = obj.linkage(k{i}, k{j});
          if(i==j)
            L(i,j) = L(i,j)/2; % cost functions assume undirected graph
          end
        end
      end
      L = L + triu(L,1)'; % symmetrize
      
      degree = sum(L,2); % all edges that have at least one node in the segment
      internal = diag(L);
      crossing = L - diag(internal);
          
      if(strcmp(criteria, 'cut'))
        % computes the sum of weights of all edges connecting different segments
        if(length(crossing) == 1) % if there's no edge crossing outside
          energy = NaN;
        else
          energy = sum(sum(crossing)); 
        end
      elseif (strcmp(criteria, 'sharon_ncut')) 

        % for each segment computes the ratio of weights of edges going
        % out and the weights of edges inside        
        if(length(crossing) == 1)
          energy = NaN;
        else
          for i=1:length(segments_to_consider)
            %energy = energy + ((sum(crossing(segments_to_consider(i),:))+sum(crossing(:,segments_to_consider(i))))/internal(segments_to_consider(i)));
            energy = energy + sum(crossing(:,segments_to_consider(i)))/internal(segments_to_consider(i));
          end        
        end
      elseif (strcmp(criteria, 'shi_cut'))
        if(length(crossing) == 1) % can't have all pixels in the same segment
          energy = NaN;
        else
          for i=1:n_segments
            energy = energy + (sum(crossing(i,:))+sum(crossing(:,i)))/(degree(i));
          end        
        end
      elseif (strcmp(criteria, 'max_density'))
        assert(length(unique(partition)) < 3); % up to binary segmentations only!
        partition(bground_label) = [];
        n_inside = length(partition);
        energy = -internal(segments_to_consider)/n_inside;
      else
        error('no such criteria defined');
      end
    end    
  end

  methods(Access = private)
    function [flow, cut] = kolmogorov_min_cut(obj,s,t)                
      % in the cut, 0 is source 1 is sink
      % it needs s and t to be the last nodes, bummer
      if( (s~=obj.n_nodes-1) || (t~=obj.n_nodes) )
        [permuted_obj, permutation] = obj.permute_a_b([s t], [obj.n_nodes-1 obj.n_nodes]);
      else
        permuted_obj.D = obj.D;
        permutation = 1:length(obj.D);
      end     

      [flow, perm_cut] = maxflowkolm(permuted_obj.D);

      % repermute it back
      perm_cut = [perm_cut; 0; 1];      
      cut = perm_cut;
      cut(permutation) = perm_cut; 
    end

    function [cut, flow, the_lambda] = bin_search_par_min_st_cut(obj, s,t, parametric_pars, visualiz)      
      % - s is source, t sink, net_update_func is a function that updates
      %               the lambda edges.
      % - u is the upper bound, l the lower bound
      % - lambda_edges are the pairs of nodes whose edges are parametric
      % - synth_nodes is just in case you want to see the interesting
      %               nodes only, then you can remove synth_nodes from the cut
      % - max_or_min is either 'max' or 'min'            
      
      DISC_FACTOR = 128; % could improve discretization
      obj = obj.discretize(DISC_FACTOR);
      
      the_lambda = -inf;
      counter = 0;
      previous_lambda = 0;
      n_iters = 0;
      %while (u-l) >= (1/(n*(n-1)))
      cut = zeros(size(obj.D,1),1);
      cut(t) = 1;
      
      lambda_edges = parametric_pars.lambda_edges;
      lambda_offsets = parametric_pars.lambda_offsets;
      lambda_weights = parametric_pars.lambda_weights;
      l = parametric_pars.l;
      u = parametric_pars.u;
      max_or_min = parametric_pars.max_or_min;
      
      % lambda_edges should be connected to source right now
      assert(lambda_edges(1,1) == s);
      
      for i=1:20 % doing no more than 20 iterations
        counter = counter + 1;

        new_lambda = (u+l)/2; % g is the mid value, u is the upper and l the lower        

        if(new_lambda == previous_lambda)
          break;
        else 
          previous_lambda = new_lambda;
          n_iters = n_iters+1
        end

        obj = obj.update_lambda_layer(lambda_edges, lambda_weights, lambda_offsets, new_lambda);
        [this_cut, flow] = obj.min_st_cut(s,t);
        
        sum(this_cut)
         
          % see if any lambda node is connected to source
        the_criteria = any((~this_cut(lambda_edges(:,2)))>0);
        
        if(the_criteria)          
          the_lambda = new_lambda;
          cut = this_cut;
          if(strcmp(max_or_min, 'max'))
            l = new_lambda; % upper = mid
          else
            u = new_lambda;
          end
        else % else lower = mid
          if(strcmp(max_or_min, 'max'))
            u=new_lambda;
          else
            l = new_lambda;
          end
        end
      end
      
      flow = flow/DISC_FACTOR;
      the_lambda = the_lambda/DISC_FACTOR;
    end
  end

  % 
  % other
  % 
  methods  
    function show_intermediate_cut(obj, cut, synth_nodes, visualiz)
      if(strcmp(visualiz.type, 'segm'))
        plottable_cut = zeros(size(visualiz.I,1)*size(visualiz.I,2),1);
        plottable_cut(visualiz.original_ids) = cut;

        plottable_cut(struct2array(visualiz.original_synth_nodes)) = [];
        imshow(reshape(plottable_cut,[size(visualiz.I,1) size(visualiz.I,2)])); % put by hand, just for now
      elseif(strcmp(visualiz.type, 'outlier_rejection'))        
        %assert('fix me');
        if(isempty(synth_nodes.lambda_nodes.t_side))
          show_surviving_max_dense_points(visualiz.I_tmpl,visualiz.I_test,...
              visualiz.F_tmpl,visualiz.F_test, find(plottable_cut(visualiz.ids_tmpl)), find(plottable_cut(visualiz.ids_test)));
        else
          plottable_cut = ~plottable_cut;
          show_surviving_max_dense_points(visualiz.I_tmpl,visualiz.I_test,...
              visualiz.F_tmpl,visualiz.F_test, find(plottable_cut(visualiz.ids_tmpl)), find(plottable_cut(visualiz.ids_test)));          
        end
      end

      %title(['iteration: ' int2str(i)]);
      drawnow();
    end    
        
    function [the_dgraph] = update_lambda_layer(the_dgraph,lambda_edges, lambda_weights, lambda_offsets, new_lambda)
      % might have to change this when having lambda simultaneously in
      % source and sink layer
      the_dgraph =the_dgraph.set_edges_weights(lambda_edges, lambda_offsets + new_lambda*lambda_weights);
    end    
          
    function integrity_check(obj)      
      assert(all(all(obj.D))>=0);
      assert(~any(any(isnan(obj.D))));
%      assert(~any(any(isinf(obj.D))));
      %assert(all(all(obj.D)) < 
    end    
%     
%     function sharon_ncut_check(obj, source, sinks)
%       % check wether it's a sharon ncut max flow problem
%       obj.D(sinks,source)
%     end


	function answer = has_cycles(obj) % return true if the graph has any cycles
		answer = has_cycles_(obj.D)
	end

  end  % methods
  
  methods(Static)
	  
	  
	  function L = DAGTopologicalSorting_(M)
		  % M is a logical (directed) adjacency matrix
		  % M(i,j) -> means edge from i to j
		  
		  % use algorithm Kahn, A. B. (1962), "Topological sorting of large networks",
		  % Communications of the ACM 5 (11): 558--562, doi:10.1145/368996.369025
		  % http://en.wikipedia.org/wiki/Topological_sorting
		  
		  L = zeros(size(M,1),1); Lpos = 1;
		  
		  S = find(sum(M)==0); % no incoming edges
		  Sst = 1; Send = length(S)+1; S(size(M,1)) = 0;
		  
		  while  Sst <= length(S)
			  n = S(Sst); Sst = Sst + 1; % remove a node n from S
			  L(Lpos) = n; Lpos = Lpos+1; % insert n into L
			  
			  
			  ms = find(M(n,:)); % for each node m with an edge e from n to m do
			  M(n, ms) = false; % remove edge e from the graph
			  
			  for m = ms
				% (Adrian: this line can be made faster by keepeng a vector
				% with the indegree of each vertex, and just counting there)
				if isempty(find(M(:, m),1)) % if m has no other incoming edges then
					S(Send) = m; Send = Send+1; % insert m into S
				end
			  end
		  end
	  end
	  
	  function answer = has_cycles_(D)
		  % run Floyd-Warshall algorithm
		  M = false(size(D));
		  
		  M(D~=0) = true;
		  
          n_nodes = size(D,1);
		  for k=1:n_nodes
			  for i=1:n_nodes
				  for j=1:n_nodes
					  M(i,j) = M(i,j) | (M(i,k) & M(k,j));
				  end
			  end
		  end
		  
		  % if we got something else than zeros on the diagonal, the graph has cycles
		  answer = ~isequal(M(sub2ind(size(M),1:n_nodes,1:n_nodes)), zeros(1, n_nodes));
	  end
	  
	  function [param_pars] = get_parametric_pars(obj)
		  param_pars = struct('lambda_edges',[],'lambda_weights',[],...
			  'max_or_min', 'max', ...
			  'l', 0, 'u', 100000000);
	  end
  end
end 
% class