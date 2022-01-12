% ODE_tests.m   Test accuracy of ODE solvers in different rounding modes.
%  This code uses the function chop at https://github.com/higham/chop.
%
% References:
%   [1] N. J. Higham, S. Pranesh. Simulating Low Precision Floating-Point
%       Arithmetic. SIAM J. Sci. Comput., 41(5), pp. 585–602. October 2019.
%       http://dx.doi.org/10.1137/19M1251308
%
%   [2] M. Fasi, M. Mikaitis. Algorithms for stochastically rounded
%       elementary arithmetic operations in IEEE 754 floating-point
%       arithmetic. IEEE Trans. Emerg. Topics Comput., 9(3): 1451–1466.
%       July 2021.
%       http://dx.doi.org/10.1109/TETC.2021.3069165
%
%   [3] M. Croci, M. Fasi, N. J. Higham, T. Mary, M. Mikaitis.
%       Stochastic Rounding: Implementation, Error Analysis, and
%       Applications. Tech. Report 2021.17, Manchester Institute for
%       Mathematical Sciences, The University of Manchester, UK.
%       October 2022. Revised January 2022.

% Clear chop options and reset the PRNG seed.
clear options
rng(1)

% Set the number of times to repeat the SR experiments.
rep = 10;

% Choose a testcase (only 0 and 1 are supported)
if ~exist('testcase', 'var')
    testcase = 1;
end

% Set up initial ODE conditions.
a = 0;
if (testcase == 0)
    b = 1.0;
    y0 = 0.015625;
elseif (testcase == 1)
    b = 0.015625;
    y0 = 1.0;
else
    error('This value of testcase is not supported.');
end

% Exact solution to the exponential decay ODE.
if (testcase == 0)
    yexact = exp(-b)*y0;
else
    yexact = exp(-b/20)*y0;
end

% Decay function.
if (testcase == 0)
    decay_ODE_acc = @(y, options)(-y);
    decay_ODE_inacc = @(y, options)(-chop(y, options));
else
    decay_ODE_acc = @(y, options)(-y/20);
    decay_ODE_inacc = @(y, options)(-chop(y/20, options));
end

nrange = round(10.^linspace(1, 6, 16));
m = length(nrange);

% Solution in binary64.
for j = 1:m
    n = nrange(j);
    x_dp = a;
    h_dp = (b-a)/n;
    y_dp = y0;

    for i=1:n
        y_dp = Euler(true, decay_ODE_acc, h_dp, y_dp, 0, []);
    end
    efp(j, 1) = abs(y_dp - yexact);
end

options.round = 1; % RN

% All chop formats.
for k = 1:6
    switch k
      case 1, options.format = 'b'; ...
            options.subnormal = 1; sr = 0;
      case 2, options.format = 'b'; ...
            options.subnormal = 1; sr = 1;
      case 3, options.format = 'h'; ...
            options.subnormal = 1; sr = 0;
      case 4, options.format = 'h'; ...
            options.subnormal = 1; sr = 1;
      case 5, options.format = 's'; ...
            options.subnormal = 1; sr = 0;
      case 6, options.format = 's'; ...
            options.subnormal = 1; sr = 1;
    end

    fprintf('k = %1.0f, prec = %s, subnormal = %1.0f\n', ...
            k, options.format, options.subnormal)
    chop([],options)

    a = chop(a); b = chop(b); y0 = chop(y0);

    for j = 1:m
        n = nrange(j);
        h_fp = chop((b-a)/n);
        y_fp_init = chop(y0);

        if sr
            repeat = rep;
        else
            repeat = 1;
        end

        avg_err = 0;
        max_err(j,k+1) = 0;
        min_err(j,k+1) = Inf;

        for l=1:repeat
            y_fp = y_fp_init;
            for i=1:n
                y_fp = Euler(false, decay_ODE_inacc, h_fp, y_fp, ...
                             sr, options);
            end
            err = abs(y_fp - yexact);
            if err > max_err(j, k+1)
                max_err(j, k+1) = err;
            end
            if err < min_err(j, k+1)
                min_err(j, k+1) = err;
            end
            avg_err = avg_err + err;
        end
        efp(j, k+1) = avg_err/repeat;
    end
end

fileName = sprintf('euler%d.dat', testcase);
fileID = fopen(fileName, 'w');
fprintf(fileID, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n', 'fp64', ...
        'bf16-rn','bf16-sr-avg', 'bf16-sr-worst', 'bf16-sr-best', ...
        'fp16-rn','fp16-sr-avg', 'fp16-sr-worst', 'fp16-sr-best', ...
        'fp32-rn', 'fp32-sr-avg', 'fp32-sr-worst', 'fp32-sr-best', 'n');
fprintf(fileID, '%e %e %e %e %e %e %e %e %e %e %e %e %e %d\n', ...
        [efp(:, 1:3), max_err(:,3), min_err(:,3), ...
         efp(:, 4:5), max_err(:,5), min_err(:,5), ...
         efp(:, 6:7), max_err(:,7), min_err(:,7), nrange']');
fclose(fileID);

function f = Euler(accurate, decay_ODE, h, y, SR, options)
    if (accurate)
        f = y + h*decay_ODE(y, options);
    elseif (SR)
        temp = options.round;
        options.round = 5;
        f = chop(y + ...
                 chop(h*chop(decay_ODE(y, options), options), ...
                      options), options);
        options.round = temp;
    else
        f = chop(y + ...
                 chop(h*chop(decay_ODE(y, options), options), ...
                      options), options);
    end
end
