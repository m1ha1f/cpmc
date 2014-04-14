function benchmark_parafmex()
  load benchmark.mat;

  lambda_min = 1;
  lambda_max = 10;
  disc_factor = 1;
  
  time = tic();
  [cuts, lambda, min_possible_lambda] = gallo_pmf(the_dgraph, lambda_edges, lambda_weights, s, t, disc_factor, lambda_min, lambda_max);
  toc(time);
  sum(cuts)
  lambda
  min_possible_lambda