% EXAMPLE_RF2LS  Worked examples of the column-by-column construction.
%
% Reproduces the examples of the accompanying article: distinct linear
% factors, repeated linear factors, an irreducible quadratic factor, and
% the personalized exam item with symbolic parameters a, b, c, d.
% Run from the repository root:  addpath src; run examples/example_RF2LS

clear; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
syms s

cases = {
    'Distinct linear factors',    5*s + 7,  (s + 1)*(s + 2)
    'Repeated linear factors',    4*s + 7,  s*(s + 1)^2
    'Irreducible quadratic factor', 3*s + 1, (s + 1)*(s^2 + 4)
};

for k = 1:size(cases,1)
    [M,n,P,~,~,lambda,Fpf] = RF2LS(cases{k,2}, cases{k,3}, s);
    fprintf('\n=== %s:  F(s) = %s\n', cases{k,1}, char(cases{k,2}/cases{k,3}));
    disp('Associated polynomials (one per column):'), disp(P.')
    disp('M ='), disp(M)
    disp('n ='), disp(n)
    disp('lambda = M\n ='), disp(lambda)
    disp('Decomposition:'), pretty(Fpf)
    % Check against MATLAB's own partial fraction routine
    assert(isAlways(simplify(Fpf - cases{k,2}/cases{k,3}) == 0))
end

% Personalized exam item: split the student ID into blocks a, b, c, d
syms a b c d
[M,n,~,~,~,lambda] = RF2LS(a*s + b, (s + c)^2*(s + d), s);
fprintf('\n=== Personalized item (as + b)/((s + c)^2 (s + d))\n');
disp('M ='), disp(M)
disp('lambda (A1, A2, B) ='), disp(simplify(lambda))

% Numerical instance: ID 22339876 -> a = 22, b = 33, c = 98, d = 76
lambda_num = subs(lambda, [a b c d], [22 33 98 76]);
disp('With a = 22, b = 33, c = 98, d = 76:'), disp(lambda_num)
% Expected: [149/44; 193/2; -149/44]

% Expected matrices for the three cases above (same as in the article):
%   [1 1; 2 1],  [1 1 0; 2 1 1; 1 0 0],  [1 1 0; 0 1 1; 4 0 1]
