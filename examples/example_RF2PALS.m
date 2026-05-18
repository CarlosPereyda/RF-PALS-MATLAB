% example_RF2PALS
% -------------------------------------------------------------------------
% Illustrative example of the RF-to-PALS construction.
% -------------------------------------------------------------------------

clear; clc;
syms s

% Proper rational function
N = 3*s + 1;
D = (s + 1)*(s^2 + 4);

% RF-to-PALS construction
[A,n,P,den,numPF,lambda,Fpf,Df,Nred,Dred] = RF2PALS(N,D,s);

disp('PALS matrix A:')
disp(A)

disp('Right-hand-side vector n:')
disp(n)

disp('Partial-fraction coefficients:')
disp(lambda)

disp('Reconstructed partial-fraction expression:')
pretty(Fpf)