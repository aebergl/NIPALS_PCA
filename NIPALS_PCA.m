function [T,P,PCAmodel,X] = NIPALS_PCA(X,varargin)
% NIPALS PCA


% Input checking
if nargin < 1
    error('NIPALS_PCA requires at least a matrix as input');
end
if ~ismatrix(X) || ~isnumeric(X) || min(size(X)) < 2
    error('The first input needs to be a NxM matrix with N & M > 1');
end

if nargin > 1 && isscalar(varargin{1}) && isnumeric(varargin{1})
    if varargin{1} >= 1
        NumComp = varargin{1};
        Options = parseArguments(varargin{2:end});
        Options.NumComp = NumComp;
    elseif varargin{1} > 0 && varargin{1} < 1
        ExplVarStop = varargin{1} * 100;
        Options = parseArguments(varargin{2:end});
        Options.ExplVarStop = ExplVarStop;
    else
        error('The second argumnet need to be a integer for NumComp or 0-1 for explained variation stop');
    end
else
    Options = parseArguments(varargin{:});
end

% Check for Missing values
MVX=[];
[N,K] = size(X);
if Options.MVCheck
    
    if Options.Verbose
        fprintf('\n')
        fprintf('Checking for missing values....\n')
    end
    
    MVX = isnan(X); %Logical matrix with true for NaNs
    if ~any(MVX,'all')
        MVX = [];  %No missing values, returns MV=[];
    else
        % check for colums with too many missing values
        mv_col = sum(MVX,1);
        col_rem = find(mv_col / N * 100 > Options.MVTolCol);
        if length(col_rem) >= 1
            X(:,col_rem) = []; %Remove column with too many mv's
            MVX(:,col_rem) = [];
            mv_col(col_rem) = [];
            disp (' ')
            fprintf('%u columns removed with more than %u %% missing values\n',length(col_rem),Options.MVTolCol);
        end
        % check for rows with too many missing values
        mv_row = sum(MVX,2);
        row_rem = find(mv_row / K * 100 > Options.MVTolRow);
        if length(row_rem) >= 1
            X(row_rem,:) = []; %Remove rows with to many mv's
            MVX(row_rem,:) = [];
            mv_row(row_rem) = [];
            disp (' ')
            fprintf('%u rows removed with more than %u %% missing values\n',length(row_rem),Options.MVTolRow);
            mv_col = sum(MVX); %Checks again
        end
        
        %Replace NaN with zeros
        X(MVX) = 0;
        MVX = ~MVX; % Create the "Hole" matrix with 0 where there are NaNs
    end
    [N,K] = size(X);
end

min_N_K = min([N,K]);
% Remove mean
if Options.CentreX == 1
    if Options.Verbose
        fprintf('\n')
        fprintf('Removing mean....\n')
    end
    if isempty(MVX) %No missing values
        PCAmodel.x_mean = sum(X,1,'includenan','default') ./ N;
    else %Missing values
        PCAmodel.x_mean = sum(X,1) ./ (N - mv_col);
    end
    X = bsxfun(@minus, X, PCAmodel.x_mean); %same as X = X - (ones(N,1) * x_mean);
    if ~isempty(MVX)
        X(~MVX) = 0; %replace MV's with zeros
    end
end

% Scale each variable to unit variance
if Options.ScaleX == 1
    if Options.Verbose
        fprintf('\n')
        fprintf('Normalizing to unit variance....\n')
    end
    if isempty(MVX)% No MV
        PCAmodel.x_weight = sqrt(sum(X.^2,1)*(1/(N-1)));
    else
        PCAmodel.x_weight = sqrt(sum(X.^2,1) ./ ((N - mv_col) - 1));
    end
    indx = (PCAmodel.x_weight < 1e-8); %Rplace
    PCAmodel.x_weight = 1 ./ PCAmodel.x_weight;
    PCAmodel.x_weight(indx) = 0;
    X = bsxfun(@times, X, PCAmodel.x_weight); %same as X = X .* (ones(N,1) * x_weight);
    if ~isempty(MVX)
        X(~MVX) = 0;%replace MV's with zeros
    end
end

if Options.MVAverage % Use zeros instaed of Missing Values
    MVX = [];
end

% Decide/estimate PCA components to calculate
if nargin == 1
    Options.NumComp = 10;
elseif Options.NumComp == 0
    Options.NumComp = Options.MaxComp;
end

% Check if NumComp is smaller then min_N_K
if Options.NumComp > min_N_K
    Options.NumComp = min_N_K;
end

% Initiate PCA calculations
if ~strcmpi('none',Options.Verbose)
    fprintf('\n')
    fprintf('Initiate PCA calculations....\n')
end

if nargout > 0
    T = zeros(N,Options.NumComp);
end
if nargout > 1
    P = zeros(K,Options.NumComp);
end
if nargout > 2
    PCAmodel.ExplVar = zeros(Options.NumComp,1);
    PCAmodel.ExplVarCum = zeros(Options.NumComp,1);
    PCAmodel.Eig = ones(Options.NumComp,1);
    PCAmodel.nIter = ones(Options.NumComp,1);
    PCAmodel.MaxOrth = ones(Options.NumComp,1);
    PCAmodel.StopCrit = cell(Options.NumComp,1);
end

ssx_orig = sum(X.^2,'all'); % Original sum of squares
if ssx_orig < 1e-10
    error('X has close to zero variance');
end

if strcmp('Comp',Options.Verbose)
    fprintf('#Comp Iter    Eig   %%Var   %%VarCum ConvValue Orthogonality\n')
elseif strcmp('Iteration',Options.Verbose)
    fprintf('#Comp\tIter\tConv\n')
end

OneMoreComp = true;
CurrentComp = false;
MaxOrth = 1;

while OneMoreComp
    CurrentComp = CurrentComp + 1;
    nit=0;
    p0 = ones(K,1) ./ sqrt(K);
    switch Options.Tstart
        case 'Ones'
            t0 = X * p0; % Ensures no sign flipping and that changes in number of variables or samples give the same sign
        case 'MaxVar'
            [~,indx] = max(sum(X.^2)); %Start with x-variable which has the largest variance
            t0 = X(:,indx);
        case 'Random'
            t0 =  X(:,randi(K));
        case 'First'
            t0 =  X(:,1);
    end
    
    
    Converged = false;
    MaxIterStop = false;
    ConvergenceValueStop = false;
    ConvergenceRatioStop = false;
    MaxOrthStop = false;
    nIncreasedConv = 0;
    while ~Converged
        if isempty(MVX)
            p = X' * t0; % same as p = ((t0' * X))'; no speed difference in MATLAB
        else
            %p = X' * t0 ./ (t0.^2' * MVX)'; % Adjust for missing values
            p = (t0' * X ./ (t0.^2' * MVX))';
        end
        
        p = p / norm(p);
        
        if isempty(MVX)
            t = X * p;
        else
            t = X * p ./ (MVX * p.^2); % Adjust for missing values
            %t = X * p ./ (MV.MVX * p.^2);
        end
        
        %Check Convergence
        Conv_Value_t = sum((t-t0).^2);
        Conv_Value_p = sum((p-p0).^2);
        
        if CurrentComp > 1
            % Check orthogonality to previous components
            MaxOrth = max(abs(sum(t.*T(:,1:CurrentComp-1),1)));
        end
        if nit > 0
            if Conv_Value_t > 0
                ConvergenceRatio = ConvergnceValue_old/Conv_Value_t;
            end
        else
            ConvergenceRatio = NaN;
        end
        
        if strcmp('Iteration',Options.Verbose)
            fprintf('%u\t%u\t%g\t%g\t%g\n',CurrentComp,nit,Conv_Value_t,ConvergenceRatio,MaxOrth)
        end
        
        if ConvergenceRatio <= 1 && Conv_Value_t < Options.ConvValue
            nIncreasedConv = nIncreasedConv + 1;
        end
        
        % Check individual stopping criteria
        if nit >= Options.MaxIter
            MaxIterStop = true;
        end
        
        if Conv_Value_t < Options.ConvValue && Conv_Value_p < Options.ConvValue
            ConvergenceValueStop = true;
        end
        if nIncreasedConv > 0
            ConvergenceRatioStop = true;
        end
        
        if MaxIterStop || (ConvergenceValueStop && ConvergenceRatioStop) || Conv_Value_t == 0
            Converged = true;
        else
            t0 = t;
            p0 = p;
            nit = nit + 1;
            ConvergnceValue_old = Conv_Value_t;
        end
    end
    
    X = X - (t*p'); % Remove compaonent and start over
    if ~isempty(MVX)
        X(~MVX) = 0;
    end
    
    if nargout > 0
        T(:,CurrentComp) = t;
    end
    if nargout > 1
        P(:,CurrentComp) = p;
    end
    ExplVarCum = (ssx_orig - sum(X.^2,'all')) / ssx_orig * 100; % Calculate cumulative explained variation
    
    if CurrentComp == 1
        ExplVar = ExplVarCum;
    else
        ExplVar = ExplVarCum - ExplVarCum_PrevComp;
    end
    ExplVarCum_PrevComp = ExplVarCum;
    Eig = ExplVar * min_N_K / 100; %Calculates the eigen value

    if nargout > 2
        PCAmodel.ExplVarCum(CurrentComp) = ExplVarCum;
        
        PCAmodel.ExplVar(CurrentComp) = ExplVar;
        
        PCAmodel.Eig(CurrentComp) = Eig;
        PCAmodel.nIter(CurrentComp) = nit;
        PCAmodel.MaxOrth(CurrentComp) = MaxOrth;
    end
    fprintf('%4u %5u %8.1f %6.2f %6.2f %8.2g %8.2g\n',CurrentComp,nit,Eig,ExplVar,ExplVarCum,Conv_Value_t,MaxOrth)
    if CurrentComp >= Options.NumComp || ExplVarCum >= Options.ExplVarStop
        OneMoreComp = false;
    end
end
if nargout > 0
    T = T(:,1:CurrentComp);
end
if nargout > 1
    P = P(:,1:CurrentComp);
end
if nargout > 2
    PCAmodel.ExplVarCum  = PCAmodel.ExplVarCum(1:CurrentComp);
    PCAmodel.ExplVar = PCAmodel.ExplVar(1:CurrentComp);
    PCAmodel.Eig  = PCAmodel.Eig(1:CurrentComp);
    PCAmodel.nIter = PCAmodel.nIter(1:CurrentComp);
    PCAmodel.MaxOrth = PCAmodel.MaxOrth(1:CurrentComp);
end


end

function options = parseArguments(varargin)
options = inputParser;

expectedVerbose = {'Comp', 'Iteration', 'None'};
expectedStopCriteria = {'Bottom', 'ConvValue'};
expectedTstart = {'Ones', 'MaxVar','Random','First'};
addParameter(options,'NumComp', 0, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(options,'StopCriteria', 'Bottom', @(x) any(validatestring(x,expectedStopCriteria)));
addParameter(options,'Tstart', 'P', @(x) any(validatestring(x,expectedTstart)));
addParameter(options,'ScaleX', false, @islogical);
addParameter(options,'CentreX', true, @islogical);
addParameter(options,'MVCheck', true, @islogical);
addParameter(options,'MVTolCol', 20, @(x) isnumeric(x) && isscalar(x) && x>0 && x<100);
addParameter(options,'MVTolRow', 20, @(x) isnumeric(x) && isscalar(x) && x>0 && x<100);
addParameter(options,'MVAverage', false, @islogical);
addParameter(options,'ConvValue', 1e-14, @(x) isnumeric(x) && isscalar(x));
addParameter(options,'MaxOrtho', 1e-8, @(x) isnumeric(x) && isscalar(x));
addParameter(options,'ExplVarStop', 100, @(x) isnumeric(x) && isscalar(x) && x>0 && x<100);
addParameter(options,'MaxComp', 100, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(options,'MaxIter', 500, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(options,'Verbose', 'Comp', @(x) any(validatestring(x,expectedVerbose)));
parse(options,varargin{:});
options = options.Results;
end



