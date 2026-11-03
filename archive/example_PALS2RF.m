% example_PALS2RF
% -------------------------------------------------------------------------
% Illustrative example of the PALS-to-RF reconstruction.
% -------------------------------------------------------------------------

clear; clc;
syms s

% PALS matrix
A = [1 1 1;
     5 3 2;
     6 0 0];

% Right-hand-side vector
n = [4;
     13;
     12];

% PALS-to-RF reconstruction
[F,D,N,P,Q,lambda,Fpf,Df] = PALS2RF(A,n,s);

disp('Reconstructed rational function:')
pretty(F)

disp('Partial-fraction coefficients:')
disp(lambda)

disp('Reconstructed partial-fraction expression:')
pretty(Fpf)