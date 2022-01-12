% run_tests.m   Scripts to reproduce the data used in [1, Fig. 8.1--8.4].
%
% References:
%   [1] M. Croci, M. Fasi, N. J. Higham, T. Mary, M. Mikaitis.
%       Stochastic Rounding: Implementation, Error Analysis, and
%       Applications. Tech. Report 2021.17, Manchester Institute for
%       Mathematical Sciences, The University of Manchester, UK.
%       October 2022. Revised January 2022.

% Add dependencies
addpath('./deps/chop');
addpath('./deps/cpfloat/mex');
currdir = pwd();
cd('./deps/cpfloat/mex');
cpfloat_compile_nomake;
cd(currdir);

% Run tests that produce Figure 8.1
test_sum

% Run tests that produce Figure 8.2
test_matvec

% Run tests that produce Figure 8.3
for testcase = [0, 1]
    ODE_tests
end

% Run tests that produce Figure 8.4
format = 'bfloat16';
for N = [32 512 2048 8192]
    unit_circle_ODE
end

format = 'fp16';
for N = [32 512 16384 65536]
    unit_circle_ODE
end
