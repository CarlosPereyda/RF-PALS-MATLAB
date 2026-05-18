function [F,D,N,P,Q,lambda,Fpf,Df] = PALS2RF(A,n,s)
% PALS2RF
% -------------------------------------------------------------------------
% Reconstructs an associated rational function from a Polynomially
% Associated Linear System (PALS).
%
% Given a matrix A and a right-hand-side vector n, this function interprets
% the columns of A as coefficient vectors of generated polynomials P_r(s).
% The denominator D(s) is reconstructed as the least common multiple of
% these generated polynomials. The numerator N(s) is obtained from n.
%
% Inputs:
%   A : PALS coefficient matrix
%   n : numerator coefficient vector
%   s : symbolic variable
%
% Outputs:
%   F      : reconstructed rational function N(s)/D(s)
%   D      : reconstructed denominator polynomial
%   N      : reconstructed numerator polynomial
%   P      : generated polynomials associated with the columns of A
%   Q      : partial denominators D(s)/P_r(s)
%   lambda : solution vector of the PALS
%   Fpf    : reconstructed partial-fraction expression
%   Df     : factorized denominator
%
% Author:
%   Carlos Alberto Pereyda-Pierre
% -------------------------------------------------------------------------

    A = sym(A);
    n = sym(n);

    [~,numCols] = size(A);

    %% 1) Reconstruct generated polynomials P_r(s)

    P = sym(zeros(numCols,1));

    for r = 1:numCols
        P(r) = poly2sym(A(:,r).',s);
        P(r) = expand(P(r));
    end

    %% 2) Reconstruct numerator polynomial N(s)

    N = poly2sym(n(:).',s);
    N = expand(N);

    %% 3) Reconstruct denominator D(s) using LCM

    D = P(1);

    for r = 2:numCols
        D = lcm(D,P(r));
        D = expand(D);
    end

    D = expand(D);
    Df = factor(D);

    %% 4) Compute partial denominators Q_r(s)

    Q = sym(zeros(numCols,1));

    for r = 1:numCols
        Q(r) = simplify(D/P(r));
    end

    %% 5) Solve the PALS

    lambda = simplify(A\n);

    %% 6) Reconstruct rational function and partial-fraction expression

    F = simplify(N/D);

    Fpf = sym(0);

    for r = 1:numCols
        Fpf = Fpf + lambda(r)/Q(r);
    end

    Fpf = simplify(Fpf);

end