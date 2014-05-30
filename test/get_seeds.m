function [foregroundids, foregroundSeeds, backgroundids, backgroundSeeds] = getSeeds(I, width, height, prow, pcol)
	rows = size(I, 1);
    cols = size(I, 2);

    foregroundSeeds = extractSeeds(I, prow, pcol, width, height);
    
    [ix, iy] = meshgrid(1:height, 1:width);
    foregroundids = [ix(:)+prow-1 iy(:)+pcol-1];
    foregroundids = sub2ind([rows cols], foregroundids(:, 1), foregroundids(:, 2));

    backgroundSeeds = [];
%     backgroundSeeds = [backgroundSeeds; extractSeeds(I, 1, 1, cols, 1)];
    backgroundSeeds = [backgroundSeeds; extractSeeds(I, 1, 1, 1, rows)];
    backgroundSeeds = [backgroundSeeds; extractSeeds(I, 1, cols, 1, rows)];
%     backgroundSeeds = [backgroundSeeds; extractSeeds(I, rows, 1, cols, 1)];

    backgroundids = sub2ind([rows cols], 1:rows, ones(1, rows));
    backgroundids = [backgroundids sub2ind([rows cols], 1:rows, ones(1, rows)*cols)];
end