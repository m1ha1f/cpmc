function [C] = computeCapacity(I, seeds)
	rows = size(I, 1);
	cols = size(I, 2);

    I = reshape(I, rows*cols, 3);
	C = ones(rows, cols);
        
    for s  = seeds'
        Cnow = sum(abs(I - repmat(s', rows*cols, 1)), 2);
        Cnow = reshape(Cnow, rows, cols);
        C = min(C, Cnow);
%         C = C + Cnow;
    end

    C = (C+0.5).^-1;
    
%     C = C / numel(seeds);
end