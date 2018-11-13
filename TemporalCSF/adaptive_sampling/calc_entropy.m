function [entropy,derivative] = calc_entropy(x,sampling_struct)
%function [entropy,derivative] = calc_entropy(x,sampling_struct)
% this function calculates the expected entropy from a posterior and a
% stimulus level tested.
% derivative requires that the derivative of the sigmoid is supplied in the
% sampling struct


alpha   = reshape(sampling_struct.X1D{1}   ,[],1);
beta    = reshape(sampling_struct.X1D{2}    ,1,[]);
lambda  = reshape(sampling_struct.X1D{3}  ,1,1,[]);
gamma   = reshape(sampling_struct.X1D{4}   ,1,1,1,[]);
scale = bsxfun(@(g,l) 1-g-l,gamma,lambda);

entropy = NaN(size(x));

if nargout >= 2
    derivative = NaN(size(x));
end

pri = sampling_struct.posterior.*sampling_struct.weights;
p   = NaN(size(x));

for i = 1:length(x)
    %% calc posterior
    xi    = x(i);
    psi   = bsxfun(@(a,b) sampling_struct.options.sigmoidHandle(xi, a, b), alpha,beta);
    psi   = bsxfun(@plus,bsxfun(@times, psi,scale), gamma);% probability of success predicted
    %    psi(psi==1) = psi-eps(1); % correction to everything is possible
    %    psi(psi==0) = eps(1);     % replaced by entropy contribution of
    %    impossible parts = 0;
    psiW  = psi.*sampling_struct.posterior.*sampling_struct.weights;
    p(i)  = sum(psiW(:));
    %% calculate Entropy
    posterior1   = psi.*sampling_struct.posterior;
    posterior1   = posterior1.*sampling_struct.weights;
    normalizer1s = sum(posterior1(:));
    posterior1   = posterior1./normalizer1s;
    entropy1     = log(posterior1).*posterior1;
    entropy1(isnan(entropy1)) = 0;  % correct for probability 0 results
    entropy1s    = sum(entropy1(:));
    posterior2   = (1-psi).*sampling_struct.posterior;
    posterior2   = posterior2.*sampling_struct.weights;
    normalizer2s = sum(posterior2(:));
    posterior2   = posterior2./normalizer2s;
    entropy2     = log(posterior2).*posterior2;
    entropy2(isnan(entropy2)) = 0;  % correct for probability 0 results
    entropy2s    = sum(entropy2(:));
    entropy(i)   = -entropy1s*p(i)-entropy2s*(1-p(i));
    %% calculate derivative
    if nargout >= 2
        Dpsi      = bsxfun(@(a,b) sampling_struct.options.derivativeHandle(xi, a, b), alpha,beta);
        Dpsi      = bsxfun(@times, Dpsi,scale);
        DpsiNorm  = Dpsi.*pri./normalizer1s;
        DpsiSum   = sum(DpsiNorm(:));
        %Dentropy1 = DpsiSum.* log(posterior1) ...
        %    + DpsiSum - posterior1.*DpsiNorm...
        %    - posterior1 .* log(posterior1).*DpsiNorm;
        Dentropy1 = (DpsiNorm - posterior1.*DpsiSum).*(log(posterior1)+1);
        DpsiNorm  = -Dpsi.*pri./normalizer2s;
        DpsiSum   = sum(DpsiNorm(:));
        %Dentropy2 = DpsiSum.* log(posterior2) ...
        %    + DpsiSum - posterior2.*DpsiNorm...
        %    - posterior2 .* log(posterior2).*DpsiNorm;
        Dentropy2 = (DpsiNorm - posterior2.*DpsiSum).*(log(posterior2)+1);
        Dentropy1(~isfinite(Dentropy1)) = 0;  % correct for probability 0 results
        Dentropy2(~isfinite(Dentropy2)) = 0;  % correct for probability 0 results
        derivative(i) = -p(i).*sum(Dentropy1(:))-(1-p(i)).*sum(Dentropy2(:))...
            -sum(pri(:).*Dpsi(:)).*(entropy1s-entropy2s);
    end
    
end