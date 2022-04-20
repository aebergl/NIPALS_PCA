function [T] = NIPALS_PCA_Predict(X,PCAmodel,varargin)
% NIPALS_PCA PCA calculation using NIPALS algorithm, handles missing values
%
% USAGE:
%


% Input checking
if nargin < 2
    error('NIPALS_PCA_Predict requires at least a matrix and PCA model as input');
end

if ~ismatrix(X) || ~isnumeric(X) || min(size(X)) < 2
    error('The first input needs to be a NxM matrix with N & M > 1');
end

if ~isstruct(PCAmodel) || ~isfield(PCAmodel,'P') || ~isfield(PCAmodel,'x_mean') || ~isfield(PCAmodel,'NumComp')
    error('The second input needs to be a PCA model structure');
end

%Check Dimensions
if size(PCAmodel.P,1) ~= size(X,2)
    error('The number of variables must agree');
end

% Check for missing values
MVX = isnan(X); %Logical matrix with true for NaNs
if ~any(MVX,'all')
    MVX = [];  %No missing values,  MVX=[];
else
    X(MVX) = 0;
    MVX = ~MVX; % Create the "Hole" matrix with 0 where there are NaNs
end
% Remove mean
if ~isempty(PCAmodel.x_mean)
    X = bsxfun(@minus, X, PCAmodel.x_mean); %same as X = X - (ones(N,1) * x_mean);
    if ~isempty(MVX)
        X(~MVX) = 0; %replace MV's with zeros
    end

end

% Scale data
if ~isempty(PCAmodel.x_weight)
    X = bsxfun(@times, X, PCAmodel.x_weight); %same as X = X .* (ones(N,1) * x_weight);
    if ~isempty(MVX)
        X(~MVX) = 0; %replace MV's with zeros
    end

end
[N,K] = size(X);
T = zeros(N,PCAmodel.NumComp);


for i=1:PCAmodel.NumComp
    if isempty(MVX)
        t = X * PCAmodel.P(:,i);
    else
        t = X * PCAmodel.W(:,i) ./ (MVX * PC.Model.W(:,i).^2);
    end
    T(:,i) = t;
    X = X - t * PCAmodel.P(:,i)';
    if ~isempty(MVX)
        X(MVX) = 0;
    end
end


% %Calculate DModX
% if isempty(MV)
%     [S1]= sqrt(((sum((X.*X)'))) / (K-NumComp));
% else
%     [S1]= sqrt(((sum((X.*X)'))) ./ (MV.mv_row-NumComp)');
% end
% RES.DmodX =  S1/PCAmodel.Model.S0(end);
% RES.fcrit = PCAmodel.Model.fcrit(end);


