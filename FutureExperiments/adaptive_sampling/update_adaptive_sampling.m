function sampling_struct = update_adaptive_sampling(sampling_struct,data)
% function sampling_struct = update_adaptive_sampling
% this function updates the sampling struct to include the new collected
% data data.
%
% This updates the grid for the posterior, whenever this is needed, which
% might take (a little) time...

optionsNew = sampling_struct.options;
optionsNew.priors = [];

%% calculate new posterior on old grid
sampling_struct.posterior = sampling_struct.posterior.*exp(...
    logLikelihood(data,optionsNew,sampling_struct.X1D{1},sampling_struct.X1D{2},sampling_struct.X1D{3},sampling_struct.X1D{4},sampling_struct.X1D{5}));


normalizer                = sampling_struct.weights.*sampling_struct.posterior;
sampling_struct.posterior = sampling_struct.posterior./sum(normalizer(:));

%% update struct fields
sampling_struct.data = [sampling_struct.data;data];
sampling_struct.last_stim = data(end,1);


%% Done if grid is fine-> test grid

%calculate marginals and test whether one is near enough to 0 to drop a
%part of the grid
recalculate = false;
for i = 1:4
    if length(sampling_struct.X1D{i})>2
        % actually not density, but the probability mass is calculated
        switch i
            case 1
                marginal = sum(sum(sum(sampling_struct.posterior.*sampling_struct.weights,2),3),4);
            case 2
                marginal = sum(sum(sum(sampling_struct.posterior.*sampling_struct.weights,1),3),4);
            case 3
                marginal = sum(sum(sum(sampling_struct.posterior.*sampling_struct.weights,1),2),4);
            case 4
                marginal = sum(sum(sum(sampling_struct.posterior.*sampling_struct.weights,1),2),3);
        end
        retest = true;
        i1     = 0; % shift for start of grid
        i2     = 0; % shift for end of grid
        while retest
            retest = false;
            if marginal(2+i1)<(1/1000)
                retest = true;
                i1     = i1+1;
            end
            if marginal(end-i2-1)<(1/1000)
                retest = true;
                i2     = i2+1;
            end
        end
        if i1>0 ||i2>0
            x1 = sampling_struct.X1D{i}(i1+1);
            x2 = sampling_struct.X1D{i}(end-i2);
            sampling_struct.X1D{i} = linspace(x1,x2,length(sampling_struct.X1D{i}));
            recalculate = true;
        end
    end
end


%% recalculate if necessary
if recalculate
    % pool data
    sampling_struct.data = poolData(sampling_struct.data,sampling_struct.options);
    % calculate posterior on new grid
    sampling_struct.posterior = exp(...
        logLikelihood(sampling_struct.data,sampling_struct.options,sampling_struct.X1D{1},sampling_struct.X1D{2},sampling_struct.X1D{3},sampling_struct.X1D{4},sampling_struct.X1D{5}));
    sampling_struct.weights   = getWeights(sampling_struct.X1D);
    normalizer                = sampling_struct.weights.*sampling_struct.posterior;
    sampling_struct.posterior = sampling_struct.posterior./sum(normalizer(:));
end