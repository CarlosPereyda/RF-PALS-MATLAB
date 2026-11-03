function v = polyvec(P,s,delta)
%POLYVEC  Coefficient vector of a polynomial in the basis {s^(delta-1),...,s,1}.
%
%   v = POLYVEC(P,s,delta) returns the delta-by-1 column vector of
%   coefficients of the polynomial P(s), from s^(delta-1) down to the
%   constant term. Polynomials of lower degree are padded with leading
%   zeros. Coefficients may be numbers or symbolic parameters.
%
%   Example
%     syms s
%     polyvec(s + 1, s, 3)      % returns [0; 1; 1]
%     polyvec(s/2 + 3, s, 2)    % returns [1/2; 3]
%
%   See also RF2LS, COEFFS.
%
%   Author: Carlos Alberto Pereyda-Pierre. License: MIT.

    P = expand(sym(P));
    [~,dP] = numden(P);
    if has(dP,s)
        error('polyvec:notPolynomial', 'The input is not a polynomial in %s.', char(s))
    end

    if isAlways(P == 0)
        v = sym(zeros(delta,1));
        return
    end

    c = coeffs(P,s,'All');               % leading coefficient first
    c = c(:);
    if numel(c) > delta
        error('polyvec:degreeTooHigh', ...
              'The polynomial has degree %d, larger than delta-1 = %d.', ...
              numel(c)-1, delta-1)
    end
    v = [sym(zeros(delta-numel(c),1)); c];
end
