function bary = buildConstSurrogateWrapper(par, args)
    bary.data = args{1}(par); % take sample
    bary.supp = args{2};
    bary.q = 1; bary.p = 1; % constant rational function
end