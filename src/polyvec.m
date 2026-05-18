function v = polyvec(P,s,delta)
% POLYVEC
% -------------------------------------------------------------------------
% Converts a symbolic polynomial into a coefficient column vector expressed
% in the ordered monomial basis:
%
%   {s^(delta-1), s^(delta-2), ..., s, 1}.
%
% If the polynomial degree is lower than delta - 1, leading zeros are added
% so that the resulting vector has length delta.
%
% Inputs:
%   P     : symbolic polynomial
%   s     : symbolic variable
%   delta : target vector length
%
% Output:
%   v     : coefficient column vector
%
% Author:
%   Carlos Alberto Pereyda-Pierre
% -------------------------------------------------------------------------

    P = expand(simplify(P));

    [numP, denP] = numden(P);

    if simplify(denP) ~= 1
        error('The input expression is not a polynomial.')
    end

    c = sym2poly(expand(numP));
    c = sym(c(:));

    if length(c) < delta
        c = [sym(zeros(delta-length(c),1)); c];
    end

    v = c;

end