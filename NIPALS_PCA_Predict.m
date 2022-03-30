function [T,DModX] = NIPALS_PCA_Predict(X,PCAmodel,varargin)
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

if size(PCAmodel.P,2) ~= size(X,2)
    
end