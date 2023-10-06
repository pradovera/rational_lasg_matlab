function f = getSampleForce_1p2(par)
% get high-fidelity value of force at some parameter
    N = 100;
    f = zeros(N+1,1);
    f(end) = 25*cos(3*pi*par(1)/(par(1)+1))*2*par(2)/(par(2)+1); % parametric dependence (scaling)
end