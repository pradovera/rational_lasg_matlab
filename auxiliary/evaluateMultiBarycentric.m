function r = evaluateMultiBarycentric(x, barys, ws)
% This function evaluates a barycentric struct at arbitrary points.
% 
% The inputs are:
%   .        x: the locations of the evaluation points (independent
%       `       variable).
%   .    barys: the cell containing barycentric structs.
%   .       ws: the vector containing the t-weights.
% 
% The output 'r' is a matrix with as many rows as 'barys{1}.p', and as many
%   columns as x has entries. The value of bary at some scalar value z is:
%                               S_i    barys{i}.p(:,j)
%                              \sum ----------------------
%                         T     j=1  z - barys{i}.supp(j)
%          barys(z,t) = \sum ------------------------------- ws(i)
%                        i=1    S_i     barys{i}.q(j)
%                              \sum ----------------------
%                               j=1  z - barys{i}.supp(j)
%   where S_i is the number of support points 'barys{i}.supp' and the number
%   of columns of the numerator 'barys{i}.p' and the number of entries of
%   the denominator 'barys{i}.q', and T is the number of elements of barys.

    if ~ismatrix(x) || size(x, 1) > 1; x = reshape(x, 1, []); end

    r = 0;
    goodi = find(ws);
    for i = goodi
        r = r + ws(i) * evaluateBarycentric(x, barys{i});
    end
end
