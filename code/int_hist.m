% Copyright (C) 2010 Joao Carreira
%
% This code is part of the extended implementation of the paper:
% 
% J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
% 

function h = int_hist(x, n)
% INT_HIST(x, n) is a histogram of all integer values 1..n in x.
% If n is not given, max(x) is used.

% Hans Olsson's one-liner from matlab faq
%h = full(sum(sparse(1:length(x(:)),x(:),1)));

h = hist(x(:),1:max(x(:))); % andreas suggestion

if nargin == 2
  if n > length(h)
    % pad with zeros
    h = [h zeros(1,n-length(h))];
  end
end
