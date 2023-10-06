clear; clc; close all
rng(42)
addpath('./auxiliary')
load("output_1p2.mat")

freq_range = [freq_test(1) freq_test(end)]; freq_0 = mean(freq_range);
pars_range = {[parsImpedance.t(1,1)-parsImpedance.dt(1), parsImpedance.t(1,1)+parsImpedance.dt(1)] ...
              [parsImpedance.t(1,2)-parsImpedance.dt(2), parsImpedance.t(1,2)+parsImpedance.dt(2)]};

%%% plot sample points
figure()
for i = 1:size(parsForce.t, 1)
    if i <= size(parsImpedance.t, 1)
        subplot(131)
        plot(appsImpedance{i}.supp, parsImpedance.t(i, 1), 'x')
        hold all
        subplot(132)
        plot(appsImpedance{i}.supp, parsImpedance.t(i, 2), 'x')
        hold all
        symbol = 'x';
    else
        symbol = 'o';
    end
    subplot(133)
    plot(parsForce.t(i, 1), parsForce.t(i, 2), symbol)
    hold all
end
subplot(131)
grid on; xlabel("freq"); ylabel("par_1")
subplot(132)
grid on; xlabel("freq"); ylabel("par_2")
subplot(133)
grid on; xlabel("par_1"); ylabel("par_2")
drawnow

%%% define grid for surrogate evaluation
n_pars_1_app = 50; n_pars_2_app = 50;
pars_1_app = linspace(pars_range{1}(1), pars_range{1}(2), n_pars_1_app);
pars_2_app = linspace(pars_range{2}(1), pars_range{2}(2), n_pars_2_app);

%%% evaluate surrogate
appF = zeros(101, n_pars_1_app, n_pars_2_app);
for j1 = 1:n_pars_1_app
    for j2 = 1:n_pars_2_app
        weightsForce = evaluatePiecewiseLinearInterpolant([pars_1_app(j1) pars_2_app(j2)], parsForce, 1, "", vanderInvForce);
        appF(:, j1, j2) = evaluateMultiBarycentric(freq_0, appsForce, weightsForce);
    end
end

%%% plot force norm
[Pars_1_app, Pars_2_app] = meshgrid(pars_1_app, pars_2_app);
normF = reshape(sum(appF.^2, 1), [n_pars_1_app, n_pars_2_app]).^.5;
figure()
surf(Pars_1_app, Pars_2_app, normF)
xlabel("app F"); ylabel("par")
drawnow
