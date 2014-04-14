% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

classdef Hgraph
 % Property data is private to the class
  properties (SetAccess = public, GetAccess = public)
    arities
    values
    ids
    maintain
    img_size
    inf_replacement    
     % fraction of non smooth pixels before it gets max clique energy
  end
  properties (SetAccess = private, GetAccess = public)
    special_nodes
    nodes
  end % properties
  
  methods
    % Construct an object
    function obj = Hgraph(n_nodes)
      obj.ids = {};
      obj.values = {};
      obj.arities = [];
      obj.maintain = [];
      obj.nodes = (1:n_nodes)';
      obj.maintain = [];
      obj.special_nodes = [];  
    end
         
    % returns pairwise terms as a directed graph,
    % with unidirectional arrows from s to the others, 
    % and from the others to t
    function the_dgraph = pairwise_to_st_dgraph(obj, s, t)
      the_dgraph = obj.pairwise_to_dgraph();
      if(~isempty(the_dgraph))
        the_dgraph(:,s) = 0;
        the_dgraph(t,:) = 0;
      end
    end

    function the_dgraph = pairwise_to_dgraph(obj)
      % returns pairwise terms as a directed graph, with arrows in both
      % directions everywhere
      slot = obj.get_pairwise_slot();
      if(isempty(slot))
        the_dgraph = Dgraph(length(obj.nodes));
        return;
      end
      the_ids = obj.ids{slot};
      the_ids = [the_ids; the_ids(:,2) the_ids(:,1)]; % bidirectional
      the_values = [obj.values{slot}; obj.values{slot}];
      n_nodes = length(obj.nodes);
      the_dgraph = Dgraph(n_nodes);
      the_dgraph = the_dgraph.add_edges(the_ids, the_values);
    end

    function [the_dgraph, the_synth_nodes] = highorder_potts_to_st_dgraph(obj,s,t, HO_POTTS_TRUNC_FRACTION)
      % returns high-order terms (arity>2) as a directed graph, 
      % with arrows from s to new synthetic nodes, and from other new
      % synthetic nodes to t
      
      the_synth_nodes = obj.create_synth_node_struct();

      ho_slots = get_highorder_slots(obj);

      % will add 2*n synthetic nodes, then will remove the unneeded
      n_synth = obj.number_of_hyperedges() * 2;
      synth_nodes = (length(obj.nodes) + 1):(length(obj.nodes) + n_synth);

      the_dgraph = Dgraph(n_synth + length(obj.nodes));

      current_synth_node = length(obj.nodes) + 1;

      n_to_preallocate = 10000000;
      ids_to_save = zeros(n_to_preallocate,2);
      values_to_save = zeros(n_to_preallocate,1);
      cur_line = 1;
      for i=1:length(ho_slots) % one arity at a time
        slot = ho_slots(i);
        the_values = obj.values{slot};
        the_ids = obj.ids{slot};
        n_hyperedges = size(the_ids,1);
        n_elem_ids = obj.arities(slot);

        for j=1:n_hyperedges % add one hyperedge at a time (should add all hyperedges of same arity at the same time)
          %%%%%% first take care of node ids
          if(~isempty(obj.special_nodes))
            pres1 = find(obj.special_nodes(1) == the_ids(j,:));
            pres2 = find(obj.special_nodes(2) == the_ids(j,:));
          else
            pres1 = [];
            pres2 = [];
          end
            
          
          if(~isempty(pres1))
            source_or_sink = 'source';
          elseif(~isempty(pres2))
            source_or_sink = 'sink';            
          end

          assert((~isempty(pres1) + ~isempty(pres2)) ~=2);

          if(isempty(pres1) && isempty(pres2)) % universal (smoothing term)
            the_synth_nodes.s_side_smooth = [the_synth_nodes.s_side_smooth; current_synth_node];
            the_synth_nodes.t_side_smooth = [the_synth_nodes.t_side_smooth; current_synth_node+1];
            [final_ids, final_values, current_synth_node] = ...
                    convert_hyp_edge_to_pn_smooth_potts(obj, the_ids(j,:),the_values(j), ...
                                     current_synth_node,n_elem_ids, s,t, HO_POTTS_TRUNC_FRACTION);
          else % class specific
            if(strcmp(source_or_sink,'source'))
              the_synth_nodes.s_side_spec = [the_synth_nodes.s_side_spec; current_synth_node]
            elseif(strcmp(source_or_sink,'sink'))              
              the_synth_nodes.t_side_spec = [the_synth_nodes.t_side_spec; current_synth_node]
            end

            [final_ids, final_values, current_synth_node] = ...
                      convert_hyp_edge_to_pn_class_specific_potts(obj, ...
                the_ids(j,:), the_values(j), current_synth_node,n_elem_ids,s,t, source_or_sink, HO_POTTS_TRUNC_FRACTION);
          end
          
          next_cur_line = cur_line + size(final_ids,1) - 1;
          ids_to_save(cur_line:next_cur_line,:) = final_ids;
          values_to_save(cur_line:next_cur_line,:) = final_values;
          cur_line = next_cur_line+1;
          %the_dgraph = the_dgraph.add_edges(final_ids, final_values);
        end
      end

      ids_to_save(cur_line:end,:) = [];
      values_to_save(cur_line:end) = [];
      the_dgraph = the_dgraph.add_edges(ids_to_save, values_to_save);
      
      % remove unneeded synthetic nodes (should i ? )
      nodes_to_remove = (current_synth_node):max(synth_nodes);
      the_dgraph = the_dgraph.remove_nodes( nodes_to_remove );
      for i=1:length(nodes_to_remove)
       synth_nodes(synth_nodes == nodes_to_remove(i))  = [];
      end
    end

    function [the_dgraph, s,t, synth_nodes, lambda_edges, lambda_weights, id_map] = density_to_st_dgraph(obj,source_or_sink)
      % source_or_sink can be 'source' or 'sink'
      if(nargin==1)
        source_or_sink = 'source';
      else
        assert(strcmp(source_or_sink, 'source') || strcmp(source_or_sink, 'sink'));
      end
      
      n_synth = obj.number_of_links();           
       
      n_nodes = length(obj.nodes);
      synth_nodes = (n_nodes + 1):(n_nodes + n_synth);
      s = length(obj.nodes) + n_synth + 1;
      t = length(obj.nodes) + n_synth + 2;
      
      id_map = 1:n_nodes; % it's the identity map (if you change this implementation change this too!
      
      the_dgraph = Dgraph(n_synth + length(obj.nodes) + 2);

      % now the construction
      % fill directed edges from nodes Ni to nodes xi
      final_subs = [];
      n_hyperedges_so_far = 0;
      
      for n=1:length(obj.arities)
        ids_arity_n = obj.ids{n};
        n_hp_this_arity = size(obj.ids{n},1);
        ids = n_nodes + ((n_hyperedges_so_far+1):(n_hyperedges_so_far+n_hp_this_arity));
        n_hyperedges_so_far = n_hyperedges_so_far + n_hp_this_arity;
        ids_col_1 = repmat(ids',obj.arities(n),1);
        ids_col_2 = reshape(ids_arity_n,obj.arities(n)*n_hp_this_arity,1);
        if(strcmp(source_or_sink, 'sink'))
          subs = [ids_col_1 ids_col_2];
        elseif(strcmp(source_or_sink, 'source'))
          subs = [ids_col_2 ids_col_1 ];
        end
        
        final_subs = [final_subs; subs];
      end
      the_dgraph = the_dgraph.add_edges(final_subs, inf);

      % weight of each hyperedge
      in_w = cell2mat(obj.values')';

      %
      % Fill directed edge capacities
      %
      
      normal_nodes = 1:n_nodes;
      if(strcmp(source_or_sink, 'source'))             
        col_sink = [zeros(1,n_nodes) in_w 0 0]';
        the_dgraph.D(:,t) = col_sink;
        lambda_edges = [s*ones(n_nodes,1) normal_nodes'];
      elseif(strcmp(source_or_sink, 'sink'))
        row_source = [zeros(1,n_nodes) in_w 0 0];
        the_dgraph.D(s,:) = row_source;
        lambda_edges = [normal_nodes' t*ones(n_nodes,1)];        
      end

      lambda_weights = ones(1,n_nodes)';
    end
    
    function [the_dgraph, s,t, synth_nodes, id_map] = sharons_ncut_to_st_dgraph(obj)  
      n_special_nodes = length(obj.special_nodes);
      assert(length(obj.special_nodes)==2);
      
      tmp1 = obj.special_nodes(1);
      tmp2 = obj.special_nodes(2);
      obj.special_nodes(1) = tmp2;
      obj.special_nodes(2) = tmp1;      
      
      assert((length(obj.arities) == 1) && (obj.arities(1) == 2)); % have to do analysis for the high-order case

      initial_lambda = 1;
 
      n_nodes = length(obj.nodes);      
      the_ids = obj.ids{1};

      n_synth = obj.number_of_links();
      synth_nodes = (n_nodes + 1):(n_nodes + n_synth);
      s = length(obj.nodes) + n_synth + 1;
      t = length(obj.nodes) + n_synth + 2;

      id_map = 1:n_nodes; % it's the identity map (if you change this implementation change this too!

      the_dgraph = Dgraph(n_synth + length(obj.nodes) + 2);

      % now the construction

      %
      % fill directed edges from nodes Ni to nodes xi
      %
      final_subs = [];
      n_hyperedges_so_far = 0;
      synth_with_special_node_2 = [];     
      for n=1:length(obj.arities) % this part is ready for high-order          
        ids_arity_n = obj.ids{n};        
        n_hp_this_arity = size(ids_arity_n,1);
        
        synth_ids = n_nodes + ((n_hyperedges_so_far+1):(n_hyperedges_so_far+n_hp_this_arity));
        
        % special nodes should be handled specially, they should be
        % separated, therefore one should be grounded to source, and the
        % other to sink        
        synth_with_special_node_2 = [synth_with_special_node_2; synth_ids(ids_arity_n==obj.special_nodes(2))];   
        
        n_hyperedges_so_far = n_hyperedges_so_far + n_hp_this_arity;
        ids_col_1 = repmat(synth_ids',obj.arities(n),1);
        ids_col_2 = reshape(ids_arity_n,obj.arities(n)*n_hp_this_arity,1);
        subs = [ids_col_1 ids_col_2];
        final_subs = [final_subs; subs]; 
      end
      
      the_dgraph = the_dgraph.add_edges(final_subs, inf*ones(size(final_subs,1),1));

      %
      % fill directed edges between nodes xi
      %
      the_ids = [the_ids; the_ids(:,2) the_ids(:,1)]; % bidirectional
      the_values = [obj.values{1}; obj.values{1}]; % assumes only one arity
      the_dgraph = the_dgraph.add_edges(the_ids, the_values);

      %
      % fill directed edges from nodes source to nodes Ni (except Ni
      % representing infinity connections between an xi and special_node(2)
      % !!! 
      % this would create an infinite path between source and sink
      %            
      in_w = cell2mat(obj.values')';

      ids_with_sn2 = synth_with_special_node_2-n_nodes;
      must_change = (in_w(ids_with_sn2)==inf);
      in_w(ids_with_sn2(must_change)) = 0;


      row_source = [zeros(1,n_nodes) initial_lambda*in_w  0 0];
      the_dgraph(s,:) = row_source;      

      the_dgraph(s,obj.special_nodes(1)) = inf; % slightly different from hochbaum's paper       
      
      % fill directed edges from xi to sink
      the_dgraph(:,t) = zeros(1,length(row_source));
      the_dgraph(obj.special_nodes(2), t) = inf; % slightly different from hochbaum's paper        
    end    
    
    function [the_dgraph, s, synth_nodes, lambda_edges, lambda_slopes, id_map] = sharons_ncut_to_s_dgraph(obj)  
      n_special_nodes = length(obj.special_nodes);
      assert(n_special_nodes == 0);      
      assert((length(obj.arities) == 1) && (obj.arities(1) == 2)); % have to do analysis for the high-order case

      n_nodes = length(obj.nodes);      
      the_ids = obj.ids{1};

      n_synth = obj.number_of_links();
      synth_nodes = (n_nodes + 1):(n_nodes + n_synth);
      s = length(obj.nodes) + n_synth + 1;

      id_map = 1:n_nodes; % it's the identity map (if you change this implementation change this too!

      the_dgraph = Dgraph(n_synth + length(obj.nodes) + 1);

      % now the construction

      %
      % fill directed edges from nodes Ni to nodes xi
      %
      final_subs = [];
      n_hyperedges_so_far = 0;
      for n=1:length(obj.arities) % this part is ready for high-order          
        ids_arity_n = obj.ids{n};        
        n_hp_this_arity = size(ids_arity_n,1);
        
        synth_ids = n_nodes + ((n_hyperedges_so_far+1):(n_hyperedges_so_far+n_hp_this_arity));
      
        n_hyperedges_so_far = n_hyperedges_so_far + n_hp_this_arity;
        ids_col_1 = repmat(synth_ids',obj.arities(n),1);
        ids_col_2 = reshape(ids_arity_n,obj.arities(n)*n_hp_this_arity,1);
        subs = [ids_col_1 ids_col_2];
        final_subs = [final_subs; subs]; 
      end
      
      the_dgraph = the_dgraph.add_edges(final_subs, inf*ones(size(final_subs,1),1));

      %
      % fill directed edges between nodes xi
      %
      the_ids = [the_ids; the_ids(:,2) the_ids(:,1)]; % bidirectional
      the_values = [obj.values{1}; obj.values{1}]; % assumes only one arity
      the_dgraph = the_dgraph.add_edges(the_ids, the_values);

      %
      % fill lambda_edges from nodes source to nodes Ni (except Ni
      % representing infinity connections between an xi and special_node(2)
      % !!! 
      % this would create an infinite path between source and sink
      %            
      in_w = cell2mat(obj.values')';
  
      lambda_slopes = in_w';
      lambda_edges = [s*ones(length(in_w), 1) [(n_nodes+1):(n_nodes)+length(in_w)]'];
      %row_source = [zeros(1,n_nodes) initial_lambda*in_w  0];
      %the_dgraph(s,:) = row_source;            
    end

    function [the_dgraph, s, synth_nodes, lambda_edges, lambda_slopes, id_map] = ratior_to_s_dgraph(obj)
      the_dgraph = obj.pairwise_to_dgraph();
      n_nodes = the_dgraph.n_nodes;
      degree = sum(the_dgraph.D)';
      
      % add the source node
      the_dgraph = the_dgraph.add_size(the_dgraph.n_nodes + 1);
      source_node = n_nodes + 1;
      synth_nodes = source_node; % just one
      s = source_node;
      
      % add the lambda connections

      lambda_edges = [source_node*ones(n_nodes,1) (1:n_nodes)'];
      %lambda_slopes = degree; % assumes 4-connected
      lambda_slopes = ones(n_nodes,1);
      
      id_map = 1:n_nodes;
    end
    
    function [the_dgraph, s, t, synth_nodes, lambda_edges, lambda_slopes, id_map] = ratior_to_st_dgraph(obj)
      the_dgraph = obj.pairwise_to_dgraph();
      n_nodes = the_dgraph.n_nodes;
      %degree = sum(the_dgraph.D)';
      
      % add the source node
      the_dgraph = the_dgraph.add_size(the_dgraph.n_nodes + 2);
      source_node = n_nodes + 1;
      sink_node = n_nodes + 2;
      synth_nodes = [source_node sink_node];        

      s = source_node;
      t = sink_node;
      % add the lambda connections

      lambda_edges = [source_node*ones(n_nodes,1) (1:n_nodes)'];
      %lambda_slopes = degree; % assumes 4-connected
      lambda_slopes = ones(n_nodes,1);
      
      id_map = 1:n_nodes;
    end
    
    % test if state is valid
    function test_integrity(obj)
      % number of cells
      assert((size(obj.arities,1) + size(obj.ids,1) + size(obj.values,1)) == 3*size(obj.arities,1));
      % number of values and number of ids
      for i=1:length(obj.arities)
        the_ids = obj.ids{i};
        the_values = obj.values{i};
        assert(length(the_values) == size(the_ids,1));
      end

      % should test if special nodes come first in row
    end

    function node_ids = get_topological_connections(obj,from_node)
      id_2_arity = get_pairwise_slot(obj);
      the_ids = obj.ids{id_2_arity};
      the_vals = obj.values{id_2_arity};
      
      matches = (the_ids(:,1) == from_node) & (the_vals == obj.inf_replacement);
      node_ids = the_ids(matches,:);
    end

    function [obj, sn_ids] = add_special_nodes(obj, n_special_nodes)
      sn_ids = obj.nodes(end) + (1:n_special_nodes)';
      obj.nodes = [obj.nodes; sn_ids];
      obj.special_nodes = [obj.special_nodes; sn_ids];
    end
    
    function obj = add_edges(obj, ids,values)
      if(isempty(values))
        return;
      end
      
      assert(size(values,2) == 1);
      assert(max(max(ids)) <= length(obj.nodes));
      null = (values == 0);
      values(null) = [];
      ids(null,:) = [];      
      ar = size(ids,2); % arity
      n_edges = size(ids,1); % number of edges
      slot = find(obj.arities==ar);

      if(isempty(slot)) % create new slot, and put them there
        last_slot = length(obj.arities);
        obj.ids{last_slot+1} = ids;
        obj.values{last_slot+1} = values;
        obj.arities(last_slot+1) = ar;
      else
        present_ids = obj.ids{slot};
        obj.ids{slot} = [present_ids; ids];
        present_values = obj.values{slot};
        obj.values{slot} = [present_values; values];
      end
    end
    
    % just high-order hyperedges (excludes order two)
    function n = number_of_hyperedges(obj)
      ho_slots = obj.get_highorder_slots();
      n = sum(cellfun(@length, obj.values(ho_slots)));
    end
    
    % all links, any order
    function n = number_of_links(obj)
      n = sum(cellfun(@length, obj.values));
    end

    function obj = replace_infinity(obj, new_inf)
      for i=1:length(obj.arities)
        v = obj.values{i};
        v(v == inf) = new_inf;
        obj.values{i} = v;
      end
    end
    
    function obj = replace_bigger_or_equal_than(obj, threshold, new_value)
      for i=1:length(obj.arities)
        v = obj.values{i};
        v(v >= threshold) = new_value;
        obj.values{i} = v;
      end
    end
    
    function obj = discretize(obj, disc_factor)
      for i=1:length(obj.arities)
        v = obj.values{i};
        v = round(v*disc_factor);
        obj.values{i} = v;
      end
    end    
  end  % methods
  
  methods (Access=private)
    function slot = get_pairwise_slot(obj)
      slot = find(obj.arities == 2);
    end

    function slots = get_highorder_slots(obj)
      slots = find(obj.arities > 2);
    end

    % it must contain s or t
    function [st_side_ids, final_values, current_synth_node] = ...
            convert_hyp_edge_to_pn_class_specific_potts(obj,the_ids, the_value, ...
                           current_synth_node,n_elem_ids,s,t, source_or_sink,HO_POTTS_TRUNC_FRACTION)
      synth_node = current_synth_node;
      current_synth_node = current_synth_node + 1;

      the_ids_col = the_ids';

      if(strcmp(source_or_sink, 'source'))
        the_ids_col(the_ids_col==s) = [];
        st_side_ids = [synth_node*ones(n_elem_ids-1,1) the_ids_col];
        st_side_ids = [s synth_node; st_side_ids];
      elseif(strcmp(source_or_sink, 'sink'))
        the_ids_col(the_ids_col==t) = [];
        st_side_ids = [the_ids_col synth_node*ones(n_elem_ids-1,1)];        
        st_side_ids = [synth_node t; st_side_ids];
      end

      %%%%%% now take care of values
      final_values = the_value*ones(size(st_side_ids,1),1);
      
      if(HO_POTTS_TRUNC_FRACTION > 0)
        assert(HO_POTTS_TRUNC_FRACTION < 0.4); % otherwise it would defy reasonable expectations
        % robust Pn potential. Energy goes linearly up with number of 
        % nonsmooth pixels, untill fraction is bigger than HO_POTTS_TRUNC_FRACTION
        final_values(2:end) = final_values(2:end)/(HO_POTTS_TRUNC_FRACTION*n_elem_ids-1);
      end                                   
    end

    function [final_ids, final_values, current_synth_node] = ...
            convert_hyp_edge_to_pn_smooth_potts(obj,the_ids, the_value, ...
                                          current_synth_node,n_elem_ids, s,t,HO_POTTS_TRUNC_FRACTION)
      s_side_synth_node = current_synth_node;
      current_synth_node = current_synth_node + 1;
      t_side_synth_node = current_synth_node;
      current_synth_node = current_synth_node + 1;

      the_ids_col = the_ids';
      s_side_ids = [s_side_synth_node*ones(n_elem_ids,1) the_ids_col];
      t_side_ids = [the_ids_col t_side_synth_node*ones(n_elem_ids,1)];
      s_ids = [s s_side_synth_node];
      t_ids = [t_side_synth_node t];
      final_ids = [s_side_ids; t_side_ids; s_ids; t_ids];

      %%%%%% now take care of values
      final_values = the_value*ones(size(final_ids,1),1);
      if(HO_POTTS_TRUNC_FRACTION > 0)
        assert(HO_POTTS_TRUNC_FRACTION < 0.4); % otherwise it would defy reasonable expectations
        % robust Pn potential. Energy goes linearly up with number of 
        % nonsmooth pixels, untill fraction is bigger than HO_POTTS_TRUNC_FRACTION
        final_values(1:(end-2)) = final_values(1:(end-2))/(HO_POTTS_TRUNC_FRACTION*n_elem_ids);
      else 
        disp('seeing if there are numerical problems, change me later');
        final_values(1:(end-2)) = final_values(1:(end-2));
      end
    end

    function synth_nodes = create_synth_node_struct(obj)
      synth_nodes.s_side_smooth = [];
      synth_nodes.t_side_smooth = [];      
      synth_nodes.s_side_spec = [];
      synth_nodes.t_side_spec = [];      
      synth_nodes.s_side_density = [];
      synth_nodes.t_side_density = [];
      synth_nodes.s_side_aux = [];
      synth_nodes.t_side_aux = [];      
    end

  end
end % class