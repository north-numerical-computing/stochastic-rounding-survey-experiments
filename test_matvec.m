% test_matvec.m   Test accuracy of matrix-vector products.
%  This code uses the CPFloat library at https://github.com/mfasi/cpfloat.
%
% References:
%   [1] M. P. Connolly, N. J. Higham, T. Mary. Stochastic rounding and its
%       probabilistic backward error analysis. SIAM J. Sci. Comput., 43(1),
%       pp. 566–585. February 2021. http://dx.doi.org/10.1137/20m1334796
%
%   [2] M. Croci, M. Fasi, N. J. Higham, T. Mary, M. Mikaitis.
%       Stochastic Rounding: Implementation, Error Analysis, and
%       Applications. Tech. Report 2021.17, Manchester Institute for
%       Mathematical Sciences, The University of Manchester, UK.
%       October 2022. Revised January 2022.

clear all
rng(1)

fs = 14; ms = 7;
m = 100;
nlist = round(logspace(1,6,20));
rep = 10;

formats = ['b', 'h'];
precisions = [8, 11];
options.explim = 1;

for k = 1:2
    t = precisions(k);
    u = 2^-t;
    options.format = formats(k);
    i = 0;
    for n = nlist
        i = i + 1;
        fprintf('i = %2d, n = %7d\n',i,n);
        A = 1e-3*rand(m,n);
        x = rand(n,1);
        ye = A*x;
        absAx = abs(A)*abs(x);
        options.round = 1;
        y1 = matvec(A,x,options);
        berr1(i) = max(abs(ye-y1)./absAx);
        for j=1:rep
            options.round = 5;
            y2 = matvec(A,x,options);
            berr2(i,j) = max(abs(ye-y2)./absAx);
        end
    end
    berr2avg = sum(berr2,2)/rep;
    berr2max = max(berr2,[],2);
    berr2min = min(berr2,[],2);

    filename = sprintf('RTN_vs_SR_t%d.dat', t);
    fid = fopen(filename, 'w');
    for i=1:length(nlist)
        fprintf(fid, "%d %f %f %f %f %f %f\n", nlist(i), berr1(i), ...
                berr2avg(i), berr2max(i), berr2min(i), min(1,nlist(i)*u), ...
                min(1,sqrt(nlist(i))*u));
    end
    fclose(fid);
end

function y = matvec(A,x,options)
    [m,n] = size(A);
    y = zeros(m,1);
    for i=1:n
        y = cpfloat(y + cpfloat(A(:,i)*x(i), options), options);
    end
end
