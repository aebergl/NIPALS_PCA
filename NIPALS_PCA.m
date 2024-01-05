 function [PCAmodel,X] = NIPALS_PCA(X,varargin)
% NIPALS_PCA PCA calculation using NIPALS algorithm, handles missing values
%
% USAGE:
%
% PCA_model = NIPALS_PCA(X) Calculates 10 components and returns the PCA_model structure
% PCA_model = NIPALS_PCA(X, A) If A=integer A components are calculated
% PCA_model = NIPALS_PCA(X, A) If 0<A<1 calculates enough components to explain A*100 percent of the variation in X
% T         = NIPALS_PCA(X, A,'Tonly',true) Only outputs T scores
% PCA_model = PCA_model(...,Name,Value) specifies options using one or more Name,Value pair arguments.
% [PCA_model, E] = PCA_model(X,A) also returns the residual matrix E
%
% INPUTS:
% * X is a NxK matrix which may contain missing values (MV)
%
% OUTPUTS:
% * PCA_Model: PCA_model structure with the following fields
%    col_rem: [] vector with columns that have been removed with too many MVs
%    row_rem: [] vector with rows that have been removed with too many MVs
%     x_mean: [1×K] vector with column mean
%   x_weight: [1xK] vector with scaling value for each column [] if no scaling
%    NumComp: A number of PCA components
%          T: [N×A] matrix with scores
%          P: [KxA] matrix with loadings
%    ExplVar: [A×1] vector with explained variance for each component in percent
% ExplVarCum: [A×1] vector with cumulative explained variance for each component in percent
%        Eig: [A×1] vector with eigenvalue for each component in percent
%      nIter: [A×1] vector with number of iteration used for each component
%    MaxOrth: [A×1] vectore with the ortogonality to previous components
%   ssx_orig: Original sum of squares
%      DmodX: [N×A] matrix with distance to model for each sample
%
% * E:        [NxK] residual matrix
%
% OTHER PARAMETERS passed as parameter-value pairs, defaults in []
% 'NumComp':     Integer,  number of components to calculate [10]
% 'AddComp':     PCA model, adds components to an existing PCA model, X is E from PCA model
% 'FullConv':    true/false stops after full convergence (or maxiter) [true]
% 'Tstart':      Type of starting vector for NIPALS algorithm [RowSum]
% 'Tonly':       true/false only OUTPUT T [false]
% 'CentreX':     true/false for removing the mean from each column in X [true]
% 'ScaleX':      true/false for scaling each column in X to unit variance [false]
% 'MVCheck':     true/false for check for missing valeus in X [true]
% 'MVTolCol':    Amount of column missing value tolerence in percent [20]
% 'MVTolRow':    Amount of row missing value tolerence in percent [20]
% 'MVAverage':   true/false replave MVs with column mean [false]
% 'ConvValue':   convergence criteria for t and p [1e-14]
% 'ExplVarStop': double, explained variance used to decide NumComp []
% 'MaxComp':     integer, maximum number of components to calculate [100]
% 'MaxIter':     integer, maximun number of iteration [5000]
% 'Verbose':     'Component'/'Iteration'/'None' [Component]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% by Anders Berglund, 2021 aebergl@gmail.com                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Input checking
if nargin < 1
    error('NIPALS_PCA requires at least a matrix as input');
end

if ~ismatrix(X) || ~isnumeric(X) || min(size(X)) < 2
    error('The first input needs to be a NxM matrix with N & M > 1');
end

if nargin > 2 && any(strcmpi('AddComp',varargin)) % Add more components to an existing PCA model
    if varargin{1} >= 1
        NumComp = varargin{1};
        options = parseArguments(varargin{2:end});
        options.NumComp = NumComp;
    elseif varargin{1} > 0 && varargin{1} < 1
        ExplVarStop = varargin{1} * 100;
        options = parseArguments(varargin{2:end});
        options.ExplVarStop = ExplVarStop;
    else
        error('The second argumnet need to be a integer for NumComp or 0-1 for explained variation stop');
    end
    PCAmodel=options.AddComp;
    options.AddComp = true;
    options.CentreX = false;
    options.ScaleX = false;
elseif nargin > 1 && isscalar(varargin{1}) && isnumeric(varargin{1}) % Second argument used for deciding number of components
    if varargin{1} >= 1
        NumComp = varargin{1};
        options = parseArguments(varargin{2:end});
        options.NumComp = NumComp;
    elseif varargin{1} > 0 && varargin{1} < 1
        ExplVarStop = varargin{1} * 100;
        options = parseArguments(varargin{2:end});
        options.ExplVarStop = ExplVarStop;
    else
        error('The second argumnet need to be a integer for NumComp or 0-1 for explained variation stop');
    end
else
    options = parseArguments(varargin{:});
end

% Check for Missing Values
MVX=[];
[N,K] = size(X);
PCAmodel.col_rem = [];
PCAmodel.row_rem = [];

if options.MVCheck
    if options.Verbose
        fprintf('\n')
        fprintf('Checking for missing values....\n')
    end
    MVX = ~isfinite(X); %Logical matrix with true for NaNs. Updated to include -inf and inf, use isfinite instead of isnan
    if ~any(MVX,'all')
        MVX = [];  %No missing values,  MV=[];
    else
        % check for colums with too many missing values
        mv_col = sum(MVX,1);
        PCAmodel.col_rem = mv_col / N * 100 > options.MVTolCol; % Keep indx for which columns were removed
        if any(PCAmodel.col_rem)
            X(:,PCAmodel.col_rem) = []; %Remove column with too many mv's
            MVX(:,PCAmodel.col_rem) = [];
            mv_col(PCAmodel.col_rem) = [];
            disp (' ')
            fprintf('%u columns removed with more than %u %% missing values\n',sum(PCAmodel.col_rem),options.MVTolCol);
        end
        % check for rows with too many missing values
        mv_row = sum(MVX,2);
        PCAmodel.row_rem = mv_row / K * 100 > options.MVTolRow; % Keep indx for which rows were removed
        if any(PCAmodel.row_rem)
            X(PCAmodel.row_rem,:) = []; %Remove rows with to many mv's
            MVX(PCAmodel.row_rem,:) = [];
            mv_row(PCAmodel.row_rem) = [];
            disp (' ')
            fprintf('%u rows removed with more than %u %% missing values\n',sum(PCAmodel.row_rem),options.MVTolRow);
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
PCAmodel.x_mean = [];
if options.CentreX == 1
    if options.Verbose
        fprintf('\n')
        fprintf('Removing mean....\n')
    end
    if isempty(MVX) %No missing values
        PCAmodel.x_mean = sum(X,1) ./ N;
    else %Missing values
        PCAmodel.x_mean = sum(X,1) ./ (N - mv_col); % Make sure it works with older version without NaN support for mean
    end
    X = bsxfun(@minus, X, PCAmodel.x_mean); %same as X = X - (ones(N,1) * x_mean);
    if ~isempty(MVX)
        X(~MVX) = 0; %replace MV's with zeros
    end
end

% Scale each variable to unit variance
PCAmodel.x_weight = [];
if options.ScaleX == 1
    if options.Verbose
        fprintf('\n')
        fprintf('Normalizing to unit variance....\n')
    end
    if isempty(MVX)% No MV
        PCAmodel.x_weight = sqrt(sum(X.^2,1)*(1/(N-1)));
    else
        PCAmodel.x_weight = sqrt(sum(X.^2,1) ./ ((N - mv_col) - 1)); % Make sure it works with older version without NaN support for std
    end
    indx = (PCAmodel.x_weight < 1e-8); % No divide by zero and also makes sure that variables with very small stDev do not get inflated
    PCAmodel.x_weight = 1 ./ PCAmodel.x_weight;
    PCAmodel.x_weight(indx) = 0;
    X = bsxfun(@times, X, PCAmodel.x_weight); %same as X = X .* (ones(N,1) * x_weight);
    if ~isempty(MVX)
        X(~MVX) = 0;%replace MV's with zeros
    end
end

if options.MVAverage % Use zeros instead of MVs, same as replacing MVs with column mean
    MVX = [];
end

% Decide/estimate number of PCA components to calculate
if nargin == 1
    options.NumComp = 10;
elseif options.NumComp == 0
    options.NumComp = options.MaxComp;
end

% Check if NumComp is smaller then min_N_K
if options.NumComp > min_N_K
    options.NumComp = min_N_K;
end

% Initiate PCA calculations
if options.Verbose
    fprintf('\n')
    fprintf('Initiate PCA calculations....\n')
end

% Create PCA model structure

if options.AddComp % Adding to the old PCA models
    CurrentComp =  PCAmodel.NumComp;
    ExplVarCum_PrevComp = PCAmodel.ExplVarCum(end);
    PCAmodel.NumComp = PCAmodel.NumComp + options.NumComp;
    PCAmodel.T = [PCAmodel.T zeros(N,options.NumComp)];
    PCAmodel.P = zeros(K,options.NumComp);
    PCAmodel.ExplVar = zeros(options.NumComp,1);
    PCAmodel.ExplVarCum = zeros(options.NumComp,1);
    PCAmodel.Eig = ones(options.NumComp,1);
    PCAmodel.nIter = ones(options.NumComp,1);
    PCAmodel.MaxOrth = ones(options.NumComp,1);
else % Create a new PCA model scructure
    CurrentComp = 0;
    PCAmodel.ssx_orig = sum(X.^2,'all'); % Original sum of squares
    PCAmodel.T = zeros(N,options.NumComp);
    if ~options.Tonly
        PCAmodel.NumComp = options.NumComp;
        PCAmodel.P = zeros(K,options.NumComp);
        PCAmodel.ExplVar = zeros(options.NumComp,1);
        PCAmodel.ExplVarCum = zeros(options.NumComp,1);
        PCAmodel.Eig = ones(options.NumComp,1);
        PCAmodel.nIter = ones(options.NumComp,1);
        PCAmodel.MaxOrth = ones(options.NumComp,1);
    end
end

% make sure that we have some variation to work with
if PCAmodel.ssx_orig < 1e-10
    error('X has close to zero variance');
end

if options.Verbose
    fprintf('\n')
    fprintf('Number of rows: %i\t Number of columns: %i\n',N,K);
    HeadingNames = {'#Comp','Iter','EigVal','%Var','%VarCum','ConvValue','Orthogonality'};
    fprintf('%5s %5s %8s %6s %6s %9s %13s\n',HeadingNames{:})
end
if strcmp('Iteration',options.Verbose)
    fprintf('\n')
    fprintf('#Comp\tIter\tConv\n')
end

OneMoreComp = true;
MaxOrth = NaN;
% Start PCA calculations, one PCA component at the time
while OneMoreComp
    CurrentComp = CurrentComp + 1;
    p0 = ones(K,1) ./ sqrt(K);
    % How to pick starting vector
    if ismatrix(options.Tstart) && isnumeric(options.Tstart) % Starting vector was provided by user
        if size(options.Tstart,2) >= CurrentComp
            t0 = options.Tstart(:,CurrentComp);
        else
            options.Tstart = 'RowSum';
        end
    else
        switch lower(options.Tstart)
            case 'rowsum'
                t0 = X * p0; % Ensures no sign flipping and that changes in number of variables or samples give the same sign
            case 'maxvar'
                [~,indx] = max(sum(X.^2)); %Start with x-variable with has the largest variance
                t0 = X(:,indx);
            case 'random' % Why?
                t0 =  X(:,randi(K));
            case 'first'
                t0 =  X(:,1); % Not recomended but sometimes used
        end
    end
    nit=0;
    Converged = false;
    MaxIterStop = false;
    ConvergenceValueStop = false;
    ConvergenceRatioStop = false;
    %MaxOrthStop = false;
    nIncreaseFlips = 0;
    IncreaseTrend = 0;

    while ~Converged  % Iterate until convergence
        % Calculate Loadings
        if isempty(MVX)
            p = X' * t0; % same as p = ((t0' * X))'; no speed difference in MATLAB
        else
            p = (t0' * X ./ (t0.^2' * MVX))';% Adjust for missing values
        end
        p = p / norm(p); % Normalise p to unit length

        % Calculate Scores
        if isempty(MVX)
            t = X * p;
        else
            t = X * p ./ (MVX * p.^2); % Adjust for missing values
        end

        %Check Convergence for both t and p
        Conv_Value_t = sum((t-t0).^2);
        Conv_Value_p = sum((p-p0).^2);

        if CurrentComp > 1
            % Check orthogonality to previous components, returns the
            % largest value
            MaxOrth = max(abs(sum(t .* PCAmodel.T(:,1:CurrentComp-1),1)));
        end

        if nit > 0
            if Conv_Value_t > 0 % Compares the rate of convergence
                ConvergenceRatio = ConvergenceValue_old/Conv_Value_t;
            end
        else
            ConvergenceRatio = NaN;
        end

        if strcmp('Iteration',options.Verbose)
            fprintf('%u\t%u\t%g\t%g\t%g\n',CurrentComp,nit,Conv_Value_t,ConvergenceRatio,MaxOrth)
        end

        if ConvergenceRatio <= 1 && IncreaseTrend == 0 % If the true bottom is reached the Convergence Ratio flips between >1 to <1
            nIncreaseFlips = nIncreaseFlips + 1;
            IncreaseTrend = 1;
        end

        if ConvergenceRatio > 1
            IncreaseTrend = 0;
        end

        % Check individual stopping criteria
        if nit >= options.MaxIter
            MaxIterStop = true;
        end

        if Conv_Value_t < options.ConvValue && Conv_Value_p < options.ConvValue
            ConvergenceValueStop = true;
        end

        if nIncreaseFlips >= 3
            ConvergenceRatioStop = true;
        end


        if MaxIterStop || (ConvergenceValueStop && (ConvergenceRatioStop || ~options.FullConv)) || Conv_Value_t == 0
            Converged = true;
        else
            t0 = t;
            p0 = p;
            nit = nit + 1;
            ConvergenceValue_old = Conv_Value_t;
        end
    end

    X = X - (t*p'); % Remove component and start over

    if ~isempty(MVX)
        X(~MVX) = 0;
    end

    % Calculate sum of squares once and only once
    sumSquaresRow = sum(X.^2,2)';
    sumSquares = sum(sumSquaresRow);

    % Calculate DModX
    if isempty(MVX)
        [S1]= sqrt(sumSquaresRow / (K-CurrentComp));
    else
        [S1]= sqrt(sumSquaresRow ./ (mv_row-CurrentComp)');
    end
    v=sqrt(N/(N-CurrentComp-1));
    S1=S1 * v;
    if isempty(MVX)
        df = ((N-CurrentComp-1)*(K-CurrentComp));
    else
        df=((N-CurrentComp-1)*(K-CurrentComp))-sum(K-mv_row);
    end
    S0=sqrt(sumSquares / df);

    % Calculate explained variation
    ExplVarCum = (PCAmodel.ssx_orig - sumSquares) / PCAmodel.ssx_orig * 100; % Calculate cumulative explained variation

    if CurrentComp == 1
        ExplVar = ExplVarCum;
    else
        ExplVar = ExplVarCum - ExplVarCum_PrevComp;
    end
    ExplVarCum_PrevComp = ExplVarCum;
    Eig = ExplVar * min_N_K / 100; %Calculates the eigen value

    % Save PCA component results
    PCAmodel.T(:,CurrentComp) = t;
    if ~options.Tonly
        PCAmodel.DmodX(:,CurrentComp) =  S1/S0;

        PCAmodel.P(:,CurrentComp) = p;
        PCAmodel.ExplVarCum(CurrentComp) = ExplVarCum;
        PCAmodel.ExplVar(CurrentComp) = ExplVar;
        PCAmodel.Eig(CurrentComp) = Eig;
        PCAmodel.nIter(CurrentComp) = nit;
        PCAmodel.MaxOrth(CurrentComp) = MaxOrth;
    end
    if options.Verbose
        fprintf('%5u %5u %8.2f %6.2f %7.2f %9.2g %13.2g\n',CurrentComp,nit,Eig,ExplVar,ExplVarCum,Conv_Value_t,MaxOrth)
    end
    % Are we done with all components?
    if CurrentComp >= options.NumComp || ExplVarCum >= options.ExplVarStop
        OneMoreComp = false;
    end
end
% Only keep the PCA components we calculated
if ~options.Tonly
    PCAmodel.NumComp = CurrentComp;
    PCAmodel.T = PCAmodel.T(:,1:CurrentComp);
    PCAmodel.P = PCAmodel.P(:,1:CurrentComp);
    PCAmodel.DmodX = PCAmodel.DmodX(:,1:CurrentComp);
    PCAmodel.ExplVarCum  = PCAmodel.ExplVarCum(1:CurrentComp);
    PCAmodel.ExplVar = PCAmodel.ExplVar(1:CurrentComp);
    PCAmodel.Eig  = PCAmodel.Eig(1:CurrentComp);
    PCAmodel.nIter = PCAmodel.nIter(1:CurrentComp);
    PCAmodel.MaxOrth = PCAmodel.MaxOrth(1:CurrentComp);
else
    PCAmodel  = PCAmodel.T(:,1:CurrentComp);
end

end

function options = parseArguments(varargin)
options = inputParser;

expectedVerbose = {'Component', 'Iteration', 'None'};
expectedTstart = {'RowSum', 'MaxVar','Random','First'};

addParameter(options,'NumComp', 0, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(options,'AddComp', [],  @(x) isstruct(x) );
addParameter(options,'Tstart', 'RowSum', @(x) ismatrix(x) || any(validatestring(x,expectedTstart)));
addParameter(options,'ScaleX', false,  @(x) islogical(x) || x==1 || x==0);
addParameter(options,'CentreX', true, @(x) islogical(x) || x==1 || x==0);
addParameter(options,'MVCheck', true, @(x) islogical(x) || x==1 || x==0);
addParameter(options,'Tonly', false, @(x) islogical(x) || x==1 || x==0);
addParameter(options,'MVTolCol', 20, @(x) isnumeric(x) && isscalar(x) && x>0 && x<100);
addParameter(options,'MVTolRow', 20, @(x) isnumeric(x) && isscalar(x) && x>0 && x<100);
addParameter(options,'MVAverage', false, @(x) islogical(x) || x==1 || x==0);
addParameter(options,'FullConv', true, @(x) islogical(x) || x==1 || x==0);
addParameter(options,'ConvValue', 1e-14, @(x) isnumeric(x) && isscalar(x));
addParameter(options,'ExplVarStop', 100, @(x) isnumeric(x) && isscalar(x) && x>0 && x<100);
addParameter(options,'MaxComp', 100, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(options,'MaxIter', 5000, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(options,'Verbose', 'Component', @(x) islogical(x) || any(validatestring(x,expectedVerbose)));
parse(options,varargin{:});
options = options.Results;

if strcmpi('none',options.Verbose)
    options.Verbose = false;
end

end



