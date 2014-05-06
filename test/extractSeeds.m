function [seeds] = extractSeeds(I, px, py, width, height)
	assert(px >= 1 & py >= 1 & width >= 1 & height >= 1);
	assert(px + height -1 <= size(I, 1));
	assert(py + width -1 <= size(I, 2));

	seeds = I(px:px+height-1, py:py+width-1, :);
    seeds = reshape(seeds, width*height, 3);
end