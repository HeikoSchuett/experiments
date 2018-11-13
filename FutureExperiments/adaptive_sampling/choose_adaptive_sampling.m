function x0 = choose_adaptive_sampling(sampling_struct,type)
% function x = choose_adaptive_sampling(sampling_struct)
% this function implements the choice of the next stimulus level based on
% the posterior

if ~exist('type','var') || isempty(type)
    type = 1;
end

if type <=4 % choose among the possible X values
%% init
x = sampling_struct.possibleX;

optionsNew = sampling_struct.options;
optionsNew.priors = [];


%% calculate posterior predictive for the possible stimulus levels

%reshaping
alpha   = reshape(sampling_struct.X1D{1}   ,[],1);
beta    = reshape(sampling_struct.X1D{2}    ,1,[]);
lambda  = reshape(sampling_struct.X1D{3}  ,1,1,[]);
gamma   = reshape(sampling_struct.X1D{4}   ,1,1,1,[]);


scale = bsxfun(@(g,l) 1-g-l,gamma,lambda);
p = nan(length(x),1);
criterion = nan(length(x),1);
for i = 1:length(x)
    xi    = x(i);
    psi   = bsxfun(@(a,b) sampling_struct.options.sigmoidHandle(xi, a, b), alpha,beta);
    psi   = bsxfun(@plus,bsxfun(@times, psi,scale), gamma);% probability of success predicted
%    psi(psi==1) = psi-eps(1); % correction to everything is possible 
%    psi(psi==0) = eps(1);     % replaced by entropy contribution of
%    impossible parts = 0;
    psiW  = psi.*sampling_struct.posterior.*sampling_struct.weights;
    p(i)  = sum(psiW(:));
    %% calculate the expected gain
    switch type
        case 1 % overall entropy needs correction, 
            %others do not because averaging over lambda and/or gamma
            %always yields values away from 0 or 1.
            entropy1   = psi.*sampling_struct.posterior.*sampling_struct.weights;
            %normalizer = entropy1.*sampling_struct.weights;
            entropy1   = entropy1./sum(entropy1(:));
            entropy1   = -log2(entropy1).*entropy1;
            entropy1(isnan(entropy1)) = 0;  % correct for probability 0 results
            entropy1s  = sum(entropy1(:));
            entropy2   = (1-psi).*sampling_struct.posterior.*sampling_struct.weights;
            %normalizer = entropy2.*sampling_struct.weights;
            entropy2   = entropy2./sum(entropy2(:));
            entropy2   = -log2(entropy2).*entropy2;
            entropy2(isnan(entropy2)) = 0;  % correct for probability 0 results
            entropy2s  = sum(entropy2(:));
            criterion(i) = -entropy1s*p(i)-entropy2s*(1-p(i));
        case 2 % entropy threshold & width
            entropy1   = psi.*sampling_struct.posterior.*sampling_struct.weights;
            normalizer = entropy1.*sampling_struct.weights;
            entropy1   = entropy1./sum(normalizer(:));
            entropy1   = sum(sum(entropy1,3),4); % This is the averaging over lambda,gamma
            entropy1   = -log2(entropy1).*entropy1;
            entropy1s  = sum(entropy1(:));
            entropy2   = (1-psi).*sampling_struct.posterior.*sampling_struct.weights;
            normalizer = entropy2.*sampling_struct.weights;
            entropy2   = sum(sum(entropy2,3),4); % This is the averaging over lambda,gamma
            entropy2   = entropy2./sum(normalizer(:));
            entropy2   = -log2(entropy2).*entropy2;
            entropy2s  = sum(entropy2(:));
            criterion(i) = -entropy1s*p(i)-entropy2s*(1-p(i));
        case 3 % entropy threshold
            entropy1   = psi.*sampling_struct.posterior.*sampling_struct.weights;
            normalizer = entropy1.*sampling_struct.weights;
            entropy1   = entropy1./sum(normalizer(:));
            entropy1   = sum(sum(sum(entropy1,3),4),2); % This is the averaging over lambda,gamma
            entropy1   = -log2(entropy1).*entropy1;
            entropy1s  = sum(entropy1(:));
            entropy2   = (1-psi).*sampling_struct.posterior.*sampling_struct.weights;
            normalizer = entropy2.*sampling_struct.weights;
            entropy2   = sum(sum(sum(entropy2,3),4),2); % This is the averaging over lambda,gamma
            entropy2   = entropy2./sum(normalizer(:));
            entropy2   = -log2(entropy2).*entropy2;
            entropy2s  = sum(entropy2(:));
            criterion(i) = -entropy1s*p(i)-entropy2s*(1-p(i));
        case 4 % entropy width
            entropy1   = psi.*sampling_struct.posterior.*sampling_struct.weights;
            normalizer = entropy1.*sampling_struct.weights;
            entropy1   = entropy1./sum(normalizer(:));
            entropy1   = sum(sum(sum(entropy1,3),4),1); % This is the averaging over lambda,gamma
            entropy1   = -log2(entropy1).*entropy1;
            entropy1s  = sum(entropy1(:));
            entropy2   = (1-psi).*sampling_struct.posterior.*sampling_struct.weights;
            normalizer = entropy2.*sampling_struct.weights;
            entropy2   = sum(sum(sum(entropy2,3),4),1); % This is the averaging over lambda,gamma
            entropy2   = entropy2./sum(normalizer(:));
            entropy2   = -log2(entropy2).*entropy2;
            entropy2s  = sum(entropy2(:));
            criterion(i) = -entropy1s*p(i)-entropy2s*(1-p(i));
            
        otherwise
            error('unrecognized type for the generation of the next stimulus level');
    end
end

[~,idx] = max(criterion);

x0 = x(idx);
else
    switch type
        case 5
            x0 = minimize(sampling_struct.last_stim,'calc_entropy',-250,sampling_struct);
        otherwise
            error('unrecognized type for the generation of the next stimulus level');
    end
end
