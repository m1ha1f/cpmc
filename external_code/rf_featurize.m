function F = rf_featurize(obj, X, Napp)
%rf_featurize get the features corresponding to the inputs X
% obj is a random feature object initialized by InitExplicitKernel.
% Don't need to specify Napp, but if specified, it would only extract so
% many features.
% (C) Fuxin Li and Catalin Ionescu 2010

[N D] = size(X);

if ~exist('Napp','var')
    if ~isfield(obj,'period') && isfield(obj,'omega')
        Napp = size(obj.omega,2);
    else
        Napp = obj.Napp;
    end
end
if isfield(obj,'omega') && Napp > size(obj.omega,2) && ~isfield(obj,'period')
    disp(['Warning: selected number of random features ' num2str(Napp) 'more than built-in number of random features ' num2str(size(obj.omega,2)) '.']);
    disp(['Changing the number of random features to ' num2str(size(obj.omega,2)) '.']);
    disp('You can increase the built-in number in rf_init()');
    Napp = size(obj.omega,2);
end

if D ~= obj.dim
  error('Dimension mismatch!');
end

switch obj.name
    case 'gaussian'
        F = sqrt(2) * (cos( X * obj.omega(:,1:Napp) + repmat(obj.beta'*2*pi,N,1)));        
    case 'laplace'
        F = sqrt(2) * (cos( X * obj.omega(:,1:Napp) + repmat(obj.beta'*2*pi,N,1)));
        % Linear is just replicate
    case 'linear'
        F = X;
    otherwise
        error('Unknown kernel approximation scheme');
end
end
