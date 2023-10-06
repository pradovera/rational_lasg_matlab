clear; clc; close all
rng(42)
addpath('./auxiliary')
addpath('./sparse_grid')
addpath('./build_surrogate')

% define frequency and parameter ranges
freq_range = [1 10];
pars_range = {[0 1] [0 .5]};

tol = 2.5e-3; tol_freq = 1e-3; % define tolerances

n_freq_test = 101; % number of test frequencies
freq_test = linspace(freq_range(1), freq_range(2), n_freq_test);
freq_0 = mean(freq_range); % arbitrary freq point

% define initial pars sparse grid
pars = []; pars.t = [mean(pars_range{1}) mean(pars_range{2})];
pars.dt = [.5*(pars_range{1}(2)-pars_range{1}(1)) .5*(pars_range{2}(2)-pars_range{2}(1))];
pars.used = 0; pars.level = [-1 -1]; pars.hierarchy_hash = pars.level; pars.hierarchy = {1};

%%% train impedance model
[parsImpedance, appsImpedance] = buildFreqParsSurrogate(freq_test, pars, tol_freq, tol, ...
                                                        @getSampleImpedance_1p2, ...
                                                        @evaluatePiecewiseLinearInterpolant);
n_impedance = size(parsImpedance.t, 1);
[~, ~, vanderInvImpedance] = evaluatePiecewiseLinearInterpolant([0 0], parsImpedance, 1);

% prepare data structures for force surrogate
getForceWrapperArgs = {@getSampleForce_1p2, freq_0};
appsForce = cell(size(appsImpedance));
n_force_keep = sum(parsImpedance.used);
parsForce = parsImpedance;
parsForce.t = parsForce.t(parsForce.used == 1, :);
parsForce.level = parsForce.level(parsForce.used == 1, :);
parsForce.used = zeros(1, n_force_keep);
parsForce = pointsBuildHierarchy(parsForce);
jGood = 1; jBad = n_force_keep + 1;
for j = 1:n_impedance
    if parsImpedance.used(j)
        jEff = jGood; jGood = jGood + 1;
    else
        jEff = jBad; jBad = jBad + 1;
    end
    % build constant rational function for force at all explored par points
    appsForce{jEff} = buildConstSurrogateWrapper(parsImpedance.t(j, :), getForceWrapperArgs);
    if ~parsImpedance.used(j)
        appsForce{jEff}.t_ref = parsImpedance.t(j, :);
    end
end

%%% train force surrogate
[parsForce, appsForce] = templateSampleAdaptive(appsForce, parsForce, tol, ...
                                                @buildConstSurrogateWrapper, getForceWrapperArgs, ...
                                                @evaluatePiecewiseLinearInterpolant, freq_0);
n_force = size(parsForce.t, 1);
[~, ~, vanderInvForce] = evaluatePiecewiseLinearInterpolant([0 0], parsForce, 1);

%%% export results
out_filename = "output_1p2.mat";
save(out_filename, "tol", "tol_freq", "freq_test", "parsImpedance", "appsImpedance", ...
     "vanderInvImpedance", "parsForce", "appsForce", "vanderInvForce")
