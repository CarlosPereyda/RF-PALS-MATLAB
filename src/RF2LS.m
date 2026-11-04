function [M,n,P,den,numPF,lambda,Fpf,Df,N,D] = RF2LS(Nin,Din,s)
%RF2LS  Linear system of a partial fraction decomposition, built column by column.
%
%   [M,n,P,den,numPF,lambda,Fpf] = RF2LS(N,D,s) takes a proper rational
%   function F(s) = N(s)/D(s) and returns the square linear system
%
%                           M * lambda = n
%
%   whose solution lambda contains the partial fraction coefficients of F.
%
%   Each unknown coefficient has an ASSOCIATED POLYNOMIAL: D(s) divided by
%   the denominator of its fraction, multiplied by s when the unknown is the
%   coefficient of s in a quadratic numerator (B*s + C). The column of M for
%   that unknown is the coefficient vector of its associated polynomial in
%   the ordered monomial basis {s^(delta-1), ..., s, 1}, delta = deg D.
%   The associated polynomials form a basis of the polynomials of degree
%   less than delta, so M is always invertible: M is the change-of-basis
%   matrix from the partial fraction basis to the monomial basis, and
%   lambda holds the coordinates of N in the partial fraction basis.
%
%   Inputs
%     N, D : symbolic polynomials in s with deg N < deg D. Coefficients may
%            be numbers or symbolic parameters (e.g. syms a b c d).
%     s    : symbolic variable.
%
%   Outputs
%     M      : delta-by-delta coefficient matrix (one column per unknown)
%     n      : coefficient vector of the (reduced) numerator
%     P      : associated polynomials, one per column of M
%     den    : denominator of each partial fraction
%     numPF  : numerator pattern of each fraction (1 or s)
%     lambda : partial fraction coefficients, lambda = M\n
%     Fpf    : sum of lambda(k)*numPF(k)/den(k), the decomposition of F
%     Df     : factored denominator
%     N, D   : reduced numerator and monic denominator actually used
%
%   Notes
%     * Common factors of N and D are cancelled first, and D is made monic.
%     * D is factored with FACTOR, i.e. over the rationals (or over the
%       field of the symbolic parameters). A quadratic factor that is
%       irreducible over the rationals but has real roots, such as s^2 - 2,
%       is kept as a quadratic block; the result is still a valid
%       decomposition, only not fully split over the reals. Factors of
%       degree greater than 2 are not supported.
%     * The unknowns are ordered factor by factor: linear factors first,
%       then quadratic ones, each group by increasing constant term (for
%       symbolic coefficients the order returned by FACTOR is kept). Within
%       each factor the powers increase, and for a quadratic factor the
%       coefficient of s (B) comes before the constant (C). For example,
%       s*(s+1)^2 gives the unknowns of 1/s, 1/(s+1), 1/(s+1)^2.
%
%   Example
%     syms s
%     [M,n,~,~,~,lambda,Fpf] = RF2LS(3*s + 1, (s + 1)*(s^2 + 4), s)
%     % M = [1 1 0; 0 1 1; 4 0 1],  n = [0; 3; 1],
%     % lambda = [-2/5; 2/5; 13/5]
%
%   See also POLYVEC, FACTOR, PARTFRAC.
%
%   Author: Carlos Alberto Pereyda-Pierre. License: MIT.

    %% 1) Reduction to lowest terms with a monic denominator
    [N,D] = numden(simplify(Nin/Din));
    N = expand(N);
    D = expand(D);

    cD  = coeffs(D,s,'All');          % leading coefficient first
    lcD = cD(1);
    N = expand(N/lcD);
    D = expand(D/lcD);

    degN = polynomialDegree(N,s);
    delta = polynomialDegree(D,s);
    if delta < 1
        error('RF2LS:constantDenominator', ...
              'The denominator must have degree at least 1.')
    end
    if degN >= delta
        error('RF2LS:improper', ...
              ['RF2LS requires a proper rational function, deg N < deg D ', ...
               '(after cancelling common factors). Divide first.'])
    end
    Df = factor(D);

    %% 2) Distinct factors and their multiplicities
    fac = Df(:);
    uniqueFac = sym([]);
    mult = [];
    for i = 1:numel(fac)
        f = expand(fac(i));
        if isempty(symvar(f)) || ~has(f,s)
            continue                    % skip constant factors, if any
        end
        j = find(arrayfun(@(g) isequal(expand(f - g), sym(0)), uniqueFac), 1);
        if isempty(j)
            uniqueFac(end+1,1) = f;    %#ok<AGROW>
            mult(end+1,1) = 1;         %#ok<AGROW>
        else
            mult(j) = mult(j) + 1;
        end
    end

    % Deterministic order: linear factors first, then quadratic ones, each
    % group by increasing constant term (when the coefficients are numeric).
    c0 = subs(uniqueFac, s, 0);
    if isempty(symvar(c0))
        degs = arrayfun(@(f) polynomialDegree(f,s), uniqueFac);
        [~,ord] = sortrows([degs(:), double(c0(:))]);
        uniqueFac = uniqueFac(ord);
        mult = mult(ord);
    end

    %% 3) Associated polynomials, one per unknown (column by column)
    P = sym([]);  den = sym([]);  numPF = sym([]);
    col = 0;
    for i = 1:numel(uniqueFac)
        f = uniqueFac(i);
        m = mult(i);
        switch polynomialDegree(f,s)
            case 1                      % A_k/(s+q)^k,  k = 1..m
                for k = 1:m
                    col = col + 1;
                    den(col,1)   = expand(f^k);
                    numPF(col,1) = sym(1);
                    P(col,1)     = exactQuotient(D,den(col),s);
                end
            case 2                      % (B_k s + C_k)/(s^2+bs+c)^k, k = 1..m
                for k = 1:m
                    Q = exactQuotient(D,expand(f^k),s);
                    col = col + 1;      % column of B_k: s*D/f^k
                    den(col,1)   = expand(f^k);
                    numPF(col,1) = s;
                    P(col,1)     = expand(s*Q);
                    col = col + 1;      % column of C_k: D/f^k
                    den(col,1)   = expand(f^k);
                    numPF(col,1) = sym(1);
                    P(col,1)     = Q;
                end
            otherwise
                error('RF2LS:highDegreeFactor', ...
                      ['FACTOR returned a factor of degree greater than 2: %s. ', ...
                       'Only linear and quadratic factors are supported.'], ...
                      char(f))
        end
    end

    if col ~= delta                     % sanity check: the system is square
        error('RF2LS:notSquare', ...
              'Internal error: %d columns for degree %d.', col, delta)
    end

    %% 4) Matrix assembly: column k is the coefficient vector of P(k)
    M = sym(zeros(delta,delta));
    for k = 1:delta
        M(:,k) = polyvec(P(k),s,delta);
    end
    n = polyvec(N,s,delta);

    %% 5) Solution and decomposition
    lambda = simplify(M\n);
    Fpf = sym(0);
    for k = 1:delta
        Fpf = Fpf + lambda(k)*numPF(k)/den(k);
    end
end

function Q = exactQuotient(D,g,s)
%EXACTQUOTIENT  Quotient D/g, checking that the division is exact.
    [Q,R] = quorem(D,g,s);
    if ~isequal(simplify(expand(R)), sym(0))
        error('RF2LS:inexactDivision', 'Non-exact division of D by %s.', char(g))
    end
    Q = expand(Q);
end
