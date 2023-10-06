function [pars, apps] = buildFreqParsSurrogate(freq_test, pars, tol_freq, tol_pars, getSample, getWeights)
% This function creates a barycentric struct from frequency data, using
%   greedy sampling! This function makes calls like 
%       buildFreqSurrogate(freq_test, tol_freq, ...);
%   where samples are taken at the parameter value pars.t(i, :). The vector
%   freq_test is also used for computation of the approximation error.
% see 'buildFreqSurrogate' for more info.
% 
% The additional inputs are:
%   .       pars: the locations of the sample points in parameter space;
%                 should be in sparse grid format.
%   .   tol_pars: tolerance on approximation error that drives
%                 pars-adaptivity.
%   . getWeights: function that computes weights for interpolation over
%                 pars-space, with signature
%                     getWeights(par, pars_supp, getWeightsArgs)
%                 where par is the parameter value where the weights have
%                 to be computed and pars_supp is the sparse grid of
%                 current support points.
% 
% The outputs 'barys' and 'Qs' are cells of entry-wise outputs of
%   'buildRationalFromData'.
    
    buildArgs = {getSample, freq_test, tol_freq};

    T = size(pars.t, 1); apps = cell(1, T);
    % get initial samples and surrogates
    for i = 1:T
        par = pars.t(i, :);
        fprintf("Adding pars-sample #%d (out of %d) at [ ", i, T); fprintf("%.4e ", par); fprintf("]\n");
        apps{i} = buildFreqSurrogateWrapper(par, buildArgs);
    end
    
    [pars, apps] = templateSampleAdaptive(apps, pars, tol_pars, @buildFreqSurrogateWrapper, buildArgs, getWeights, freq_test);
end
