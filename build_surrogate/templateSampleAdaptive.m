function [pars, apps] = templateSampleAdaptive(apps, pars, tol_pars, buildFreqSurrogateWrapper, buildArgs, getWeights, getErrorArgs)
% It is allowed to provide apps and pars of different length. In such case,
%   it is assumed that the additional elements of apps are exactly those
%   that would be computed by running the for loop on the refined elements
%   of pars. In particular, exactly T_max extra elements of apps are given,
%   with T_max being the value computed at the first 'while' iteration.
    idx_refine = find(~pars.used);
    while numel(idx_refine) > 0
        pars_old = pars; T = size(pars.t, 1); Tapps = numel(apps);
        [pars, T_new] = pointsRefine(pars, idx_refine);
        idx_refine = [];
        % compute interpolation weights at all new points
        wsi = getWeights(pars.t(T + 1 : end, :), pars_old, 1);
        generate_samples = (T == numel(apps));
        if generate_samples
            apps{T + T_new} = {zeros(0)};
        elseif T + T_new ~= Tapps
            error("Extra values of apps are inconsistent with refinement of pars.")
        end
        % get new samples and surrogates
        for i = T + 1:T + T_new
            par_test = pars.t(i, :);
            if generate_samples
                fprintf("Adding pars-sample #%d (%d out of %d) at [ ", i, i-T, T_new); fprintf("%.4e ", par_test); fprintf("]\n");
                apps{i} = buildFreqSurrogateWrapper(par_test, buildArgs);
            else
                % reorder precomputed samples so that they match pars(T + 1 : end)
    		    % this *should* be unnecessary, but it's cheap and safe
                if ~all(abs(apps{i}.t_ref - par_test) < 1e-10)
                    % look for the right app!
                    j_correct = 0;
                    for j = i + 1:T + T_new
                        if all(abs(apps{j}.t_ref - par_test) < 1e-10)
                            j_correct = j; break
                        end
                    end
                    if j_correct
                        app = apps{j};
                        apps{j} = apps{i};
                        apps{i} = app;
                    else
                        error("Extra values of apps are inconsistent with refinement of pars.")
                    end
                end
                try
                    apps{i} = rmfield(apps{i}, "t_ref");
                catch
                end
                fprintf("Reusing pars-sample #%d (%d out of %d) at [ ", i, i-T, T_new); fprintf("%.4e ", par_test); fprintf("]\n");
            end

            % compute error between surrogate and actual value at par_test
            ws = [wsi(i - T, :) zeros(1, T_new)];
            approx_error = getError(apps, ws, i, getErrorArgs);
            fprintf(2, "\tApproximation error = %.4e ", approx_error);

            if approx_error > tol_pars
                fprintf(" !!!!!");
                idx_refine = [idx_refine i];
            end
            fprintf("\n");
        end
    end
end
