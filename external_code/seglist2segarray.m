function [the_lines] = seglist2segarray(seglist)
  the_lines = [];
  for i=1:length(seglist)
    nodes = seglist{i}';
    j = 1:(size(nodes,2)-1);
    new_lines = [nodes(:,j); nodes(:,j+1)];
    the_lines = [the_lines new_lines];
  end
  