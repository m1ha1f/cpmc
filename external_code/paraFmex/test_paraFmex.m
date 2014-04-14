function test_paraFmex(I)  
  img_sz = size(I);
  pars = JustSegmGen.create_pars();
  hgraph_gen = JustSegmGen(pars);  
  the_hgraph = hgraph_gen.apply(I);
  the_hgraph.test_integrity();     
    
 [best_partit, min_flow_so_far, this_dgraph, extra_nodes] = ...
   exhaustive_1_times_n_search(the_hgraph, I);        

  best_partit(extra_nodes) = [];

  segmI_to_show = reshape(best_partit, img_sz(1), img_sz(2));
  imshow(segmI_to_show);
end

% a2) searches only over all pairs of t, rest as above
function [best_partit, min_flow_so_far, this_dgraph,extra_nodes] = ...
    exhaustive_1_times_n_search(the_hgraph, I)
  % let's test every combination of an edge in the source and an edge in
  % the sink  

  sharons_conv = SharonsNcutConverter(the_hgraph, 'segm');  
  [original_dgraph, s_pivot, t, sharons_synth_nodes, lambda_edges, lambda_slopes] = sharons_conv.apply(); 
  
  pars.s = s_pivot;
  
  new_n_nodes = original_dgraph.n_nodes + 1;
  original_dgraph = original_dgraph.add_size(new_n_nodes);
  t = new_n_nodes;

  synth_nodes = [s_pivot sharons_synth_nodes t];
  
  min_flow_so_far = inf;
   
  sink_edge_nodes =  setdiff(1:original_dgraph.n_nodes, synth_nodes);
  
  n_cuts_computed = 0;
  for j=1:length(sink_edge_nodes)
    the_sink_pixel = sink_edge_nodes(j);

    this_dgraph = original_dgraph;
    this_dgraph.D(the_sink_pixel,t) = inf;

    [cuts, lambdas] = gallo_pmf(this_dgraph, lambda_edges, lambda_slopes, s_pivot, t);

    n_cuts_computed = n_cuts_computed+1;
  end

  best_partit = cuts(:,1);
  extra_nodes = [sharons_synth_nodes s_pivot t the_hgraph.special_nodes'];
end
 