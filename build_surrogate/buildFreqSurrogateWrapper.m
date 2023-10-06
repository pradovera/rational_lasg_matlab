function bary = buildFreqSurrogateWrapper(par, args)
    % bargs = {getSample, freq_test, tol, getArgs, compress_fun, err_compr, evalError, Smax}
    getSampleEff = @(f) args{1}(f, par);
    bary = buildFreqSurrogate(args{2}, args{3}, getSampleEff);
end