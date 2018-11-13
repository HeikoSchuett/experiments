function sampling_struct = adaptive_sampling(sampling_struct,data)
% function adaptive_sampling
% this function calculates optimal stimulus levels to measure a
% psychometric function. It implements a variant of the psi method.
% main difference is that we use the posterior predictive to evaluate the
% probability for an correct answer not the mean estimate. 


