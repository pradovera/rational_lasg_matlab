function bary = buildMinimalRationalInterpolant(freq, supp, R)
% This function creates a barycentric struct from frequency data. If 'freq'
%   and 'supp' have the same length, the rational surrogate interpolates
%   the data. If 'freq' has more entries than 'supp', the surrogate
%   approximates the data in a LS sense. 'freq' cannot have less entries
%   than 'supp'.
%
% Note that bary is low-dimensional, since it does not approximate the data
%   directly. Instead, it approximates the coefficients of data on the
%   basis given by data! I.e., evaluating 'bary' at a point x yields an
%   S-dimensional vector, with S the number of sample frequencies, that is,
%   the number of entries of 'freq' and the number of columns of 'data'.
%   Thin means that, for extrapolation at a new point x, one should compute
%   'data*bary(x)', and not just 'bary(x)'!
% 
% The inputs are:
%   .      freq: the locations of the sample points (independent variable).
%   .      supp: the desired support points of the rational surrogate.
%   .         R: the upper triangular matrix appearing in the QR
%                decomposition of the data.
% 
% The output 'bary' has 5 properties:
%   . supp: supp (column vector).
%   .    p: the coefficients of the barycentric numerator (matrix where
%           each column corresponds to a different support point).
%   .    q: the coefficients of the barycentric denominator (row vector
%           where each entry corresponds to a different support point).
% The output 'Q' is the orthogonal matrix appearing in the QR
% decomposition of data.
% See 'evaluateBarycentric' for more info on the barycentric format.

    N = numel(supp) - 1; S = numel(freq);
    
    % assign support points
    bary.supp = supp(:);

    %%% build denominator
    % build cauchy matrix and its variants
    cauchy = freq - bary.supp.';
    % clean cauchy matrix (inf values might be present)
    ibad = []; jbad = [];
    for j = 1:N+1
        [dj, ij] = min(abs(cauchy(:, j)));
        if dj < 1e-8; ibad = [ibad ij]; jbad = [jbad j]; end
    end
    cauchy = cauchy.^-1; cauchy(ibad, :) = 0;
    cauchyR = cauchy; cauchy2 = cauchy' * cauchy;
    cauchy(:, jbad) = 0; cauchy2(jbad, :) = 0;
    for i = 1:numel(ibad)
        cauchy(ibad(i), jbad(i)) = 1;
        cauchyR(ibad(i), jbad(i)) = 1;
        cauchy2(jbad(i), jbad(i)) = 1;
    end
    % build vector that extracts dominant coefficient of barycentric form
    c_sum = (ones(1, size(cauchy, 2)) / cauchy2) * cauchy';

    % stack modified data info
    Rstack = zeros(S, N + 1);
    for j = 1:S; Rstack(j, :) = (c_sum .* R(j, :)) * cauchyR; end

    % get denominator as minimizer of a quadratic form
    [~,gramhs,gramv] = svd(Rstack);
    [gramhs,idx] = sort(abs(diag(gramhs)));
    bary.q = gramv(:, idx(1)).';
    % print stability info
    if S > 1
        fprintf("\tMRI of size %d from %d samples.\tMinimal singular value: %.4e. Spectral gap: %.4e.\n",...
                N + 1, S, gramhs(1) / gramhs(end), (gramhs(end) - gramhs(1)) / (gramhs(2) - gramhs(1)));
    else
        fprintf("\tMRI of size %d from 1 sample.\n", N + 1);
    end

    % build numerator by interpolation or LS approximation
    bary.p = (conj(cauchy) .* (cauchyR * bary.q.')) / cauchy2.';
end
