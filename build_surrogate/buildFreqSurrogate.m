function bary = buildFreqSurrogate(freq_test_0, tol_freq, getSampleEff)
% This function creates a barycentric struct from frequency data, using
%   greedy sampling!
% 
% The inputs are:
%   . freq_test: vector containing a grid of test frequencies (independent
%                variable). The function might select an arbitrary subset
%                of them as sample points (which are also support points).
%   .  tol_freq: the desired tolerance on the relative approximation
%                error.
%   . getSample: function to be called to obtain a new sample. The
%                signature is 'u = getSample(z, getArgs)', where the input
%                'z' is a scalar and the output 'u' is a column vector.
% See 'evaluateBarycentric' for more info on the barycentric format.

    % initialize training points and remove them from testing points
    freq_train = [freq_test_0(1, 1); freq_test_0(1, end)];
    freq_test = freq_test_0(1, 2:end-1);

    % get training samples and their QR decomposition
    data = [getSampleEff(freq_train(1, 1)) getSampleEff(freq_train(2, 1))];
    [Q, R] = qr(data, 0);

    break_zloop = false; S = 2;
    % each iteration adds a new sample point
    while true
        % build surrogate from current samples
        bary = buildMinimalRationalInterpolant(freq_train, freq_train, R);
        bary.data = data;
        if numel(freq_test) <= 0 || break_zloop; break; end
        
        % find worst-approximated point (in a certain sense) with current surrogate
        % this can be done just by looking at the surrogate denominator
        Qnorm = abs(bary.q * (freq_test - freq_train).^-1);
        [~, next_idx] = min(Qnorm);

        % move new sample point from testing points to training points
        freq_train = [freq_train; freq_test(next_idx)]; freq_test(next_idx) = [];
        fprintf("\tAdding sample #%d at freq = %.4e.\n", numel(freq_train), freq_train(end));

        % get new training sample and update the QR decomposition by Gram-Schmidt
        data = [data getSampleEff(freq_train(end))];
        S = S + 1;
        next_reduced_basis = data(:, end);
        next_reduced_R = zeros(size(R, 1), 1);
        for j = 1:2 % orthogonalizing twice is enough for numerical stability
            reduced_component = Q' * next_reduced_basis;
            next_reduced_basis = next_reduced_basis - Q * reduced_component;
            next_reduced_R = next_reduced_R + reduced_component;
        end
        next_norm = norm(next_reduced_basis);
        Q = [Q next_reduced_basis / next_norm];
        R = [R next_reduced_R; zeros(1, size(R, 2)) next_norm];
        
        % check relative approximation error at new training sample
        error = norm(evaluateBarycentric(freq_train(end), bary) - data(:, end)) / max([next_norm 1]);
        fprintf("\t\tComputed test error = %.4e.\n",error);
        if error < tol_freq; break_zloop = true; end
    end
end
