function r = evaluateBarycentric(x, bary)
% This function evaluates a barycentric struct at arbitrary points.
% 
% The inputs are:
%   .        x: the locations of the evaluation points (independent
%               variable).
%   .     bary: the struct containing the barycentric function.
% 
% The output 'r' is a matrix with as many rows as 'bary.p', and as many
%   columns as x has entries. The value of bary at some scalar value z is:
%                           S     bary.p(:,j)
%                         \sum ------------------
%                          j=1  z - bary.supp(j)
%             bary(z) = ---------------------------
%                           S      bary.q(j)
%                         \sum ------------------
%                          j=1  z - bary.supp(j)
%   where S is the number of support points 'bary.supp' and the number of
%   columns of the numerator 'bary.p' and the number of entries of the
%   denominator 'bary.q'.

    if ~ismatrix(x) || size(x, 1) > 1; x = reshape(x, 1, []); end

    % build cauchy matrix
    cauchy = x - bary.supp;
    % clean cauchy matrix (inf values might be present)
    ibad = []; jbad = []; ibadz = []; jbadz = []; i = 1;
    while i <= numel(bary.supp)
        [di, ij] = min(abs(cauchy(i, :)));
        if di < 1e-8
            if bary.q(i) == 0 % remove small value in cauchy
                ibadz = [ibadz i]; jbadz = [jbadz ij];
                cauchy(i, ij) = 1;
            else  % remove column values in cauchy and set small to 1
                ibad = [ibad i]; jbad = [jbad ij];
                cauchy(:, ij) = 1;
            end
        else
            i = i + 1;
        end
    end
    
    % cauchy matrix is the same for numerator and denominator
    cauchy = cauchy.^-1; cauchy(:, jbad) = 0;
    for i = 1:numel(ibad); cauchy(ibad(i), jbad(i)) = 1; end
    for i = 1:numel(ibadz); cauchy(ibadz(i), jbadz(i)) = 0; end
    r = bary.data*((bary.p * cauchy) ./ (bary.q * cauchy));
end
