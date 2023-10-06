function [u_int, pw_weights, pw_weights_generator] = evaluatePiecewiseLinearInterpolant(par, pars_supp, u_supp, pw_weights, pw_weights_generator)
% This function evaluates a piecewise linear interpolant at points "t". The
%   support values "u_supp" are prescribed at the support points "t_supp".
%   The interpolation weights "pw_weights" are either given by the user or
%   computed by Vandermonde inversion.
    if nargin < 4 || (nargin < 5 && (isstring(pw_weights) || ischar(pw_weights) || nargout > 2))
        % compute pw_weights (and maybe pw_weights_generator) from scratch
        vander = getHierarchicalValuesSelf(pars_supp);
        if nargout > 2
            pw_weights_generator = inv(vander);
            pw_weights = getHierarchicalValues(pars_supp, par) * pw_weights_generator;
        else
            pw_weights = getHierarchicalValues(pars_supp, par) / vander;
        end
    elseif isstring(pw_weights) || ischar(pw_weights)
        % compute pw_weights from pw_weights_generator
        pw_weights = getHierarchicalValues(pars_supp, par) * pw_weights_generator;
    end
    u_int = pw_weights * u_supp;
end

function hv = getHierarchicalValuesSelf(pars)
    % set up hierarchical function weights over t space (for efficiency, bound checks are skipped!)
    [T, d] = size(pars.t);
    size_IJV = ceil(T * min(1 + log(T), 2));
    IJV = zeros(size_IJV, 3);
    % populate diagonal and first column
    IJV(1:2 * T - 1, :) = [(1:T)' (1:T)' ones(T, 1); ...
                           (2:T)' ones(T - 1, 1) ones(T - 1, 1)];
    used_IJV = 2 * T - 1;
    baseline = pars.level(1, 1);
    for i = 2:size(pars.hierarchy_hash, 1)
        % populate columns with indices ts.hierarchy{i}
        hash_i = pars.hierarchy_hash(i, :); hier_i = pars.hierarchy{i};
        width = .5 .^ hash_i .* pars.dt; N_i = numel(hier_i);
        % find descendants
        goodh = find(all(pars.hierarchy_hash >= hash_i, 2));
        for j_j = 1:numel(goodh)
            j = goodh(j_j);
            if j == i; continue; end
            hier_j = pars.hierarchy{j}; dtk = ones(N_i, numel(hier_j));
            for k = 1:d
                if hash_i(k) > baseline
                    dtk = dtk .* max(1 - abs(pars.t(hier_i, k) - pars.t(hier_j, k).') / width(k), 0);
                end
            end
            IJ = find(dtk > 1e-15);
            if IJ
                % add nonzero elements
                new_used_IJV = used_IJV + numel(IJ);
                if new_used_IJV > size_IJV % increase size of IJV
                    IJV = [IJV; zeros(T, 3)];
                    size_IJV = size_IJV + T;
                end
                I = mod(IJ - 1, N_i) + 1;
                J = floor((IJ - 1) / N_i) + 1;
                IJV(used_IJV + 1 : new_used_IJV, :) = [reshape(hier_j(J), [], 1), ...
                                                       reshape(hier_i(I), [], 1), ...
                                                       reshape(dtk(IJ), [], 1)];
                used_IJV = new_used_IJV;
            end
        end
    end
    hv = sparse(IJV(1:used_IJV, 1), IJV(1:used_IJV, 2), IJV(1:used_IJV, 3), T, T);
end

function hv = getHierarchicalValues(pars, par)
    % set up hierarchical function weights over t space (for efficiency, bound checks are skipped!)
    [T, d] = size(pars.t); N = size(par, 1); hv = ones(N, T);
    if all(pars.level(1, :) == 0); baseline = 0; else; baseline = -1; end
    for k = 1:d
        goodk = pars.level(:, k) > baseline;
        widthk = .5.^pars.level(goodk, k).'*pars.dt(k);
        hv(:, goodk) = hv(:, goodk) .* max(1 - bsxfun(@rdivide, abs(par(:, k) - pars.t(goodk, k).'), widthk), 0);
    end
end
