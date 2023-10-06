function f = getSampleForce_1p1(par)
% get high-fidelity value of force at some parameter
    N = 100;
    f = zeros(N+1,1);
    f(end) = 25*cos(3*pi*par/(par+1))*2*par/(par+1); % parametric dependence (scaling)
end