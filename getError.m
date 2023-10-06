function err_norm = getError(apps, ws, i, freq_test)
% get hybrid error between current surrogate and reference surrogate
    ws0 = [zeros(1, i - 1) 1 zeros(1, numel(ws) - i)];
    ref = evaluateMultiBarycentric(freq_test, apps, ws);
    app = evaluateMultiBarycentric(freq_test, apps, ws0);
    ref_norm = sum(abs(ref).^2,1);
    ref_norm(ref_norm < 1.) = 1.;
    err_norm = sum(abs(app - ref).^2,1);
    err_norm = mean(err_norm./ref_norm).^.5;
end
