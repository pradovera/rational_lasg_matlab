clear; clc; close all
rng(42)
addpath('./auxiliary')
load("output_1p1.mat")

freq_range = [freq_test(1) freq_test(end)]; freq_0 = mean(freq_range);
pars_range = [parsImpedance.t(1,1)-parsImpedance.dt, parsImpedance.t(1,1)+parsImpedance.dt];

%%% plot sample points
figure()
for i = 1:size(parsForce.t, 1)
    if i <= size(parsImpedance.t, 1)
        subplot(1,5,1:4)
        plot(appsImpedance{i}.supp, parsImpedance.t(i, 1), 'x')
        hold all
    end
    subplot(155)
    plot(appsForce{i}.supp, parsForce.t(i, 1), 'x')
    hold all
end
subplot(1,5,1:4)
grid on; xlabel("freq"); ylabel("par")
subplot(155)
grid on; ylabel("par")
drawnow

%%% define grid for error computation
n_freq_err = 100; n_pars_err = 50;
freq_err = linspace(freq_range(1), freq_range(2), n_freq_err+2);
pars_err = linspace(pars_range(1), pars_range(2), n_pars_err+2);
freq_err = freq_err(2:end-1); pars_err = pars_err(2:end-1);

%%% compute high-fidelity reference
exF = zeros(101, n_pars_err); exZ = zeros(10201, n_freq_err, n_pars_err);
for j = 1:n_pars_err
    exF(:, j) = getSampleForce_1p1(pars_err(j));
    exZ(:, :, j) = getSampleImpedance_1p1(freq_err, pars_err(j));
end

%%% evaluate surrogate
appF = zeros(size(exF)); appZ = zeros(size(exZ));
weightsImpedance = evaluatePiecewiseLinearInterpolant(pars_err.', parsImpedance, 1, "", vanderInvImpedance);
weightsForce = evaluatePiecewiseLinearInterpolant(pars_err.', parsForce, 1, "", vanderInvForce);
for j = 1:n_pars_err
    fprintf("test sample %d\n", j)
    appF(:, j) = evaluateMultiBarycentric(freq_0, appsForce, weightsForce(j, :));
    appZ(:, :, j) = evaluateMultiBarycentric(freq_err, appsImpedance, weightsImpedance(j, :));
end

%%% measure error
discrF = zeros(size(exF)); errF = zeros(1, n_pars_err);
discrZ = zeros(size(exZ)); errZ = zeros(1, n_pars_err);
for j = 1:n_pars_err
    discrF(:, j) = appF(:, j) - exF(:, j);
    discrZ(:, :, j) = appZ(:, :, j) - exZ(:, :, j);
    refF_norm = sum(abs(exF(:, j)).^2,1);
    refF_norm(refF_norm < 1) = 1;
    errF(j) = (sum(abs(discrF(:, j)).^2,1) ./ refF_norm).^.5;
    refZ_norm = sum(abs(exZ(:, :, j)).^2, 1);
    refZ_norm(refZ_norm < 1) = 1;
    errZ(j) = mean(sum(abs(discrZ(:, :, j)).^2, 1) ./ refZ_norm, 2).^.5;
end

%%% plot error
figure()
subplot(131)
plot(appF', pars_err)
xlabel("app F"); ylabel("par")
subplot(132)
plot(exF', pars_err)
xlabel("ex F"); ylabel("par")
subplot(133)
plot(discrF', pars_err)
xlabel("abs err F"); ylabel("par")
drawnow
figure()
subplot(121)
semilogx(errF, pars_err)
xlabel("rel err (RMS) F"); ylabel("par")
subplot(122)
semilogx(errZ, pars_err)
xlabel("rel err (RMS) Z"); ylabel("par")
drawnow
