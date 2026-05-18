function [A,n,P,den,numPF,lambda,Fpf,Df,N,D] = RF2PALS(Nin,Din,s)
% RF2PALS
% -------------------------------------------------------------------------
% Constructs the Polynomially Associated Linear System (PALS) induced by a
% proper rational function.
%
% Given a proper rational function F(s) = N(s)/D(s), this function performs
% a canonical reduction, factors the denominator, generates the associated
% polynomials, assembles the PALS matrix, and computes the corresponding
% partial-fraction coefficients.
%
% Inputs:
%   Nin : numerator polynomial
%   Din : denominator polynomial
%   s   : symbolic variable
%
% Outputs:
%   A      : PALS coefficient matrix
%   n      : numerator coefficient vector
%   P      : generated polynomials associated with the PALS columns
%   den    : partial-fraction denominators
%   numPF  : partial-fraction numerators
%   lambda : solution vector of the associated PALS
%   Fpf    : reconstructed partial-fraction expression
%   Df     : factorized denominator
%   N      : reduced numerator polynomial
%   D      : reduced monic denominator polynomial
%
% Author:
%   Carlos Alberto Pereyda-Pierre
% -------------------------------------------------------------------------

    %% 1) Canonical reduction

    F = simplify(Nin/Din);
    [N,D] = numden(F);

    N = expand(N);
    D = expand(D);

    % Normalize the denominator to monic form
    cD = sym2poly(D);
    lcD = cD(1);

    N = expand(simplify(N/lcD));
    D = expand(simplify(D/lcD));

    degN = polynomialDegree(N,s);
    degD = polynomialDegree(D,s);

    if degN >= degD
        error('RF2PALS requires deg(N) < deg(D) after simplification.')
    end

    delta = degD;
    Df = factor(D);

    %% 2) Exact factor extraction

    fac = factor(D);
    fac = fac(:);

    uniqueFac = sym([]);
    mult = [];

    for i = 1:length(fac)

        f = expand(fac(i));
        found = false;

        for j = 1:length(uniqueFac)
            if simplify(f - uniqueFac(j)) == 0
                mult(j) = mult(j) + 1;
                found = true;
                break
            end
        end

        if ~found
            uniqueFac(end+1,1) = f;
            mult(end+1,1) = 1;
        end
    end

    %% 3) PALS column construction

    P = sym([]);
    den = sym([]);
    numPF = sym([]);

    col = 1;

    for i = 1:length(uniqueFac)

        f = expand(uniqueFac(i));
        m = mult(i);

        degf = polynomialDegree(f,s);

        if degf == 1

            % Linear factor contribution
            for k = 1:m

                den(col,1) = expand(f^k);
                numPF(col,1) = sym(1);

                [Q,R] = quorem(D,den(col),s);

                if simplify(R) ~= 0
                    error('Non-exact division detected in D/den.')
                end

                P(col,1) = expand(numPF(col)*Q);

                col = col + 1;
            end

        elseif degf == 2

            % Irreducible quadratic factor contribution
            for k = 1:m

                den(col,1) = expand(f^k);
                numPF(col,1) = s;

                [Q,R] = quorem(D,den(col),s);

                if simplify(R) ~= 0
                    error('Non-exact division detected in D/den.')
                end

                P(col,1) = expand(numPF(col)*Q);

                col = col + 1;

                den(col,1) = expand(f^k);
                numPF(col,1) = sym(1);

                P(col,1) = expand(numPF(col)*Q);

                col = col + 1;
            end

        else
            error(['This implementation supports only linear factors ', ...
                   'and irreducible quadratic factors.'])
        end
    end

    %% 4) Matrix assembly

    A = sym(zeros(delta,length(P)));

    for k = 1:length(P)
        A(:,k) = polyvec(P(k),s,delta);
    end

    n = polyvec(N,s,delta);

    %% 5) PALS solution

    lambda = simplify(A\n);

    %% 6) Partial-fraction reconstruction

    Fpf = sym(0);

    for k = 1:length(lambda)
        Fpf = Fpf + lambda(k)*numPF(k)/den(k);
    end

    Fpf = simplify(Fpf);

end