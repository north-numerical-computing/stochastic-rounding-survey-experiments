% test_sum.m   Test accuracy of partial sums of harmonic series.
%  This script evaluates the harmonic sum \sum_{i=1}^{\infty} 1/i
%  in different arithmetics and with different rounding modes.
%  This code uses the function chop at https://github.com/higham/chop.
%
% References:
%   [1] M. Fasi, M. Mikaitis. Algorithms for stochastically rounded
%       elementary arithmetic operations in IEEE 754 floating-point
%       arithmetic. IEEE Trans. Emerg. Topics Comput., 9(3), pp. 1451–1466.
%       July 2021.
%       http://dx.doi.org/10.1109/TETC.2021.3069165
%
%   [2] M. Croci, M. Fasi, N. J. Higham, T. Mary, M. Mikaitis.
%       Stochastic Rounding: Implementation, Error Analysis, and
%       Applications. Tech. Report 2021.17, Manchester Institute for
%       Mathematical Sciences, The University of Manchester, UK.
%       October 2022. Revised January 2022.

rng(4)

formats = ['h', 'b'];
t = [11, 8]; % Number of precision bits.

for global_index=1:2
    options.format = formats(global_index);

    % nvals = [10 50 100 5e2 1000 5e3 1e4 5e4 1e5 5e5 1e6];
    nmax = 1e6;
    npoints = 10;
    nvals = round( exp( linspace(log(10), log(nmax), npoints ) ) );

    clear err
    err = zeros(length(nvals), 4); %RN, RS avg, RS worst, RS best

    % alpha = 0.9;
    alpha = 1;

    for i = 1:length(nvals)
        n = nvals(i);
        x = chop(ones(1,n) ./ (1:n).^alpha);
        sx = sum(x);
        for k = [1 5]

            options.round = k;
            % Initialize: subsequent calls chop(x) reuse options.
            chop([],options)

            if (k==5)
                rep = 10;
            else
                rep = 1;
            end
            errors = zeros(1,rep);
            for repeat = 1:rep
                s = 0;
                for j = 1:n
                    s = chop(s + x(j));
                end
                errors(repeat) = abs((sx-s)/sx);
            end
            eavg = sum(errors)/rep;
            eworst = max(errors);
            ebest = min(errors);
            err(i,round(k/5)+1) = eavg;
            err(i,3) = eworst;
            err(i,4) = ebest;
        end

    end
    nearmax = err(:,1);
    stochmax = err(:,2);
    stochworst = err(:,3);
    stochbest = err(:,4);

    fs=12; ms=7; lw=1;

    dim = nvals;
    lam = 1;
    u = 2^(-11); % unit roundoff for fp16
    pbound = exp(lam*sqrt(dim).*(2*u) + dim.*(2*u)^2/(1-2*u)) - 1;
    pbound(pbound < 0) = 1; pbound(pbound > 1) = 1;

    dbound = (dim.*u)./(1-dim.*u);
    dbound(dbound < 0) = 1; dbound(dbound > 1) = 1;

    filename = sprintf('test_sum_t%d.dat', t(global_index));
    fid = fopen(filename, 'w');
    for i=1:length(dim)
        fprintf(fid, "%d %f %f %f %f %f %f\n",...
                dim(i), nearmax(i), stochmax(i), stochworst(i), ...
                stochbest(i), dbound(i), pbound(i));
    end
    fclose(fid);
end