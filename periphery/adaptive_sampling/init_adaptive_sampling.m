function sampling_struct = init_adaptive_sampling(x0,stimRange,expType,options)
% function sampling_struct = init_adaptive_sampling(x0,stimRange,expType,options)
% initialization of the struct for the adaptive sampling method. 
%

sampling_struct = struct;

if ~exist('x0','var') || isempty(x0)
    sampling_struct.last_stim = .1;
else
    sampling_struct.last_stim = x0;
end

if ~exist('stimRange','var') || isempty(stimRange)
    stimRange = [eps,.5];
elseif length(stimRange)>2
    sampling_struct.possibleX = stimRange;
    stimRange = [min(stimRange),max(stimRange)];
else 
    sampling_struct.possibleX = exp(linspace(log(stimRange(1)),log(stimRange(2)),81));
end
sampling_struct.stimRange = stimRange;

    %% options
if ~exist('options','var') || isempty(options),options             = struct;  end
if exist('expType','var') && ~isempty(expType),options.expType     = expType; end
if ~isfield(options,'expType'),                options.expType     = '2AFC';  end
if ~isfield(options,'sigmoidName'),            options.sigmoidName = 'norm';  end

if strcmp(options.expType,'2AFC'),             options.expType     = 'nAFC';
                                               options.expN        = 2;       end
                                           
if any(strcmpi(options.sigmoidName,{'Weibull','logn','weibull'}))
    options.logspace = 1;
else
    options.logspace=0;
end
options.useGPU           = 0;
options.verbose          = -5;                                           
options.widthalpha       = .05;
options.sigmoidHandle    = getSigmoidHandle(options);
options.derivativeHandle = getDerivativeHandle(options);

options.priors{1}= @(x) 1;
options.priors{2}= @(x) 1;
options.priors{3}= @(x) betapdf(x,1,20);
options.priors{4}= @(x) betapdf(x,1,20);
options.priors{5}= @(x) 1;


options.poolMaxGap    = inf;        % max gap between two trials of a block
options.poolMaxLength = inf;        % maximal blocklength 
options.poolxTol      = 50*eps(stimRange(2));             % maximal difference to elements pooled in a block



sampling_struct.data       = [];

sampling_struct.options = options;

%% first grid
if options.logspace
    sampling_struct.X1D{1} = linspace(log(stimRange(1)),log(stimRange(2)),20);
    sampling_struct.X1D{2} = linspace(eps(log(stimRange(2))),log(stimRange(2))-log(stimRange(1)),20);
else
    sampling_struct.X1D{1} = linspace(stimRange(1),stimRange(2),20);
    sampling_struct.X1D{2} = linspace(eps(stimRange(2)),stimRange(2)-stimRange(1),20);
end
sampling_struct.X1D{3} = linspace(0,.5,10);
if strcmp(sampling_struct.options.expType,'YesNo')
    sampling_struct.X1D{4} = linspace(0,.5,10);
elseif strcmp(sampling_struct.options.expType,'nAFC')
    sampling_struct.X1D{4} = 1/sampling_struct.options.expN;
else
    sampling_struct.X1D{4} = 0;
end
sampling_struct.X1D{5} = 0;

X1D = sampling_struct.X1D;
sampling_struct.posterior  = ones(length(X1D{1}),length(X1D{2}),length(X1D{3}),length(X1D{4}),length(X1D{5}));
sampling_struct.posterior  = bsxfun(@times,sampling_struct.posterior, sampling_struct.options.priors{3}(reshape(sampling_struct.X1D{3},1,1,[])));

sampling_struct.weights    = getWeights(X1D);
normalizer                 = sampling_struct.weights.*sampling_struct.posterior;
sampling_struct.posterior  = sampling_struct.posterior./sum(normalizer(:));
