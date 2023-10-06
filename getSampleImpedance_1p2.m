function Z = getSampleImpedance_1p2(freq, par)
% get high-fidelity value of impedance at some frequency and parameter
    N = 100;
    Z = (2*eye(N+1)-diag(ones(N,1),1)-diag(ones(N,1),-1))/N; %discrete Laplacian
    Z(end,end-1) = -1j/N; Z(end,end) = 1j/N; % damped Neumann BC
    Z = Z*(10+par(1)+20*par(2)/(par(2)+1)); % parametric dependence (scaling)
    Z(1,1) = 1; Z(1,2) = 0; % Dirichlet BC

    S = numel(freq); freq = reshape(freq,1,S);
    Z = repmat(reshape(Z,[],1),1,S);
    % add frequency dependence
    Z(end-N-1, :) = Z(end-N-1, :) .* freq;
    Z(end, :) = Z(end-N-1, :) .* freq;
end