function Handle = getDerivativeHandle(options)


if isstruct(options)
    %actually is checked in psignifit as well
    if ~isfield(options,'widthalpha'), options.widthalpha = .05; end
    alpha   = options.widthalpha;
    sigmoid = options.sigmoidName;
elseif ischar(options)
    sigmoid = options;
    alpha   = .05;
end

if ischar(sigmoid)
    switch sigmoid
        case {'norm','gauss'}   % cumulative normal distribution
            C      = my_norminv(1-alpha,0,1) - my_norminv(alpha,0,1);
            %Handle = @(X, m, width) my_normcdf(X, m, width ./ C);
            Handle = @(X,m,width) normpdf(X,m,width./C);
        case 'logistic'         % logistic function
            %Handle = @(X, m, width) 1 ./ (1 + exp(-2 * log(1./alpha - 1) ./ width .* (X-m)));
            Handle = @(X, m, width) exp(-2 * log(1./alpha - 1) ./ width .* (X-m))./((1 + exp(-2 * log(1./alpha - 1) ./ width .* (X-m))).^2);
        case 'gumbel'           % gumbel
            % note that gumbel and reversed gumbel definitions are sometimesswapped
            % and sometimes called extreme value distributions
            C      = log(-log(alpha)) - log(-log(1-alpha));
            %Handle = @(X, m, width) 1 - exp(-exp(C ./ width .* (X-m) + log(-log(.5))));
            Handle = @(X, m, width) exp(-exp( C ./ width .* (X-m) + log(-log(.5)))).*log(2).*exp(C./width.*(X-m)).*(C./width);
        case 'rgumbel'          % reversed gumbel 
            % note that gumbel and reversed gumbel definitions are sometimesswapped
            % and sometimes called extreme value distributions
            C      = log(-log(1-alpha)) - log(-log(alpha));
            %Handle = @(X, m, width) exp(-exp( C ./ width .* (X-m) + log(-log(.5))));
            Handle = @(X, m, width) -exp(-exp( C ./ width .* (X-m) + log(-log(.5)))).*log(2).*exp(C./width.*(X-m)).*(C./width);
        case 'logn'             % cumulative lognormal distribution
            C      = my_norminv(1-alpha,0,1) - my_norminv(alpha,0,1);
            Handle = @(X, m, width) my_normcdf(log(X), m, width ./ C);
        case {'Weibull','weibull'} % Weibull
            C      = log(-log(alpha)) - log(-log(1-alpha));
            %Handle = @(X, m, width) 1 - exp(-exp(C./ width .* (log(X)-m) + log(-log(.5))));
            Handle = @(X, m, width) exp(-exp( C ./ width .* (log(X)-m) + log(-log(.5)))).*log(2).*exp(C./width.*(log(X)-m)).*(C./width)./X;
        case {'tdist','student','heavytail'} 
                                % student T distribution with 1 df
                                %-> heavy tail distribution
            C      = (my_t1icdf(1-alpha) - my_t1icdf(alpha));
            Handle = @(X, m, width) my_t1cdf(C.*(X-m)./ width);
        otherwise
            error('unknown sigmoid function');
    end
elseif isa(sigmoid,'function_handle')
  %Handle = @(X,m,width) sigmoid(X,m,width);
    Handle = sigmoid;  %???
else
    error('sigmoid must be either handle or name of a sigmoid');
end