% unit_circle_ODE.m   Tests for Euler integration of unit circle ODE.
%  This script integrates the unit circle ODE system
%      u'(t) = v(t),    u(0) = 1,
%      v'(t) = u(t),    v(0) = 0,
%  using the forward Euler method with various floating-point arithmetics.
%  This code uses the function chop at https://github.com/higham/chop.
%
% References:
%   [1] N. J. Higham (ed.). Princeton Companion to Applied Mathematics.
%       p. 51. 2015.
%
%   [2] M. Fasi, M. Mikaitis. Algorithms for stochastically rounded
%       elementary arithmetic operations in IEEE 754 floating-point
%       arithmetic. IEEE Trans. Emerg. Topics Comput., 9(3), pp. 1451–1466.
%       July 2021.
%       http://dx.doi.org/10.1109/TETC.2021.3069165
%
%   [3] M. Croci, M. Fasi, N. J. Higham, T. Mary, M. Mikaitis.
%       Stochastic Rounding: Implementation, Error Analysis, and
%       Applications. Tech. Report 2021.17, Manchester Institute for
%       Mathematical Sciences, The University of Manchester, UK.
%       October 2021. Revised January 2022.

% Clear chop options and reset the PRNG seed.
clear coordinates_exact coordinates coordinates_temp
rng(500)

% Set the number of times to repeat the SR experiments.
rep = 10;

% Steps of ODE integration.
if ~exist('N', 'var')
    N = 65536;
end

% Target format.
if ~exist('format', 'var')
    format = 'fp16';
end

% Set up number of coordinates to sample.
points = min(100, N);

% Arrays for coordinates.
coordinates_exact = zeros(2, points+2);
coordinates = zeros(4, points+2, 2);
coordinates_temp = zeros(2, points+2, 2, rep);

% Get the coordinates of the exact solution.
h_dp = 2*pi/N;
j = 1;
coordinates_exact(:, 1) = [1,0];
for i=1:N
    if (mod(i, floor(N/points)) == 0)
        j = j+1;
        coordinates_exact(:, j) = [cos((i-1)*h_dp), -sin((i-1)*h_dp)];
    end
end
coordinates_exact(:, end) = [1,0];

options.round = 1; % RN
% Run two cases: with RN and with SR.
for k = 1:2
    switch k
        case 1, options.format = format; options.subnormal = 1; sr = 0;
        case 2, options.format = format; options.subnormal = 1; sr = 1;
    end

    fprintf('k = %1.0f, prec = %s, subnormal = %1.0f\n',...
            k,options.format,options.subnormal)
    chop([],options)

    % Timestep size
    h_fp = chop(2*pi/N, options);
    closest_avg = 999;
    closest_avg_i = 0;
    furthest_avg_outside = 0;
    furthest_avg_outside_i = 0;
    furthest_avg_inside = 0;
    furthest_avg_inside_i = 0;
    avg_err = 0;
    avg_err_signed = 0;
    if sr
        repeat = rep;
    else
        repeat = 1;
    end
    for l=1:repeat
        j=1; avg_err = 0;
        % Initial values.
        u_fp = chop(1);
        v_fp = chop(0);
        coordinates_temp(k, j, :, l) = [u_fp, v_fp];
        for i=1:N
            [u_fp, v_fp] = Euler(0, h_fp, u_fp, v_fp, sr, options);
            if (mod(i, floor(N/points)) == 0 && j < points+1)
                j = j+1;
                coordinates_temp(k, j, :, l) = [u_fp, v_fp];
            end
        end
        coordinates_temp(k, end, :, l) = [u_fp, v_fp];
    end
    if (sr)
        % Compute closest and farthest points.
        u_pos = squeeze(coordinates_temp(k,:,1,:));
        v_pos = squeeze(coordinates_temp(k,:,2,:));
        norms = sqrt(u_pos.^2 + v_pos.^2);
        [~, minpos] = min(norms, [], 2);
        [~, maxpos] = max(norms, [], 2);
        for i = 1:points+2
            coordinates(2, i, 1) = u_pos(i, maxpos(i));
            coordinates(2, i, 2) = v_pos(i, maxpos(i));
            coordinates(3, i, 1) = u_pos(i, minpos(i));
            coordinates(3, i, 2) = v_pos(i, minpos(i));
        end
        coordinates(4, :, :) = squeeze( ...
            sum(coordinates_temp(k, :, :, :), 4)) / repeat;
    else
        coordinates(1, :, :) = squeeze( ...
            sum(coordinates_temp(1, :, :, :), 4)) / repeat;
    end
end

fileNamePrefix = sprintf('unit_circle_%s', options.format);
fileNameSuffix = sprintf('_%d.dat', N);

fileID = fopen('unit_circle_exact.dat', 'w');
fprintf(fileID, 'u v \n');
fprintf(fileID, '%e %e \n', [coordinates_exact(1,:); ...
                             coordinates_exact(2,:)]);
fclose(fileID);

filename = sprintf('%s_rn%s', fileNamePrefix, fileNameSuffix);
fileID = fopen(filename, 'w');
fprintf(fileID, 'u v \n');
fprintf(fileID, '%e %e \n', [coordinates(1,:,1); coordinates(1,:,2)]);
fclose(fileID);

filename = sprintf('%s_sr_worst_outside%s', fileNamePrefix, ...
                   fileNameSuffix);
fileID = fopen(filename, 'w');
fprintf(fileID, 'u v \n');
fprintf(fileID, '%e %e \n', [coordinates(2,:,1); coordinates(2,:,2)]);
fclose(fileID);

filename = sprintf('%s_sr_worst_inside%s', fileNamePrefix, fileNameSuffix);
fileID = fopen(filename, 'w');
fprintf(fileID, 'u v \n');
fprintf(fileID, '%e %e \n', [coordinates(3,:,1); coordinates(3,:,2)]);
fclose(fileID);

filename = sprintf('%s_sr_best%s', fileNamePrefix, fileNameSuffix);
fileID = fopen(filename, 'w');
fprintf(fileID, 'u v \n');
fprintf(fileID, '%e %e \n', [coordinates(4,:,1); coordinates(4,:,2)]);
fclose(fileID);

function [u, v] = Euler(accurate, h, u_prev, v_prev, SR, options)
    if (accurate)
        u = u_prev + h*v_prev;
        v = v_prev - h*u_prev;
    elseif (SR)
        temp = options.round;
        options.round = 5;
        u = chop(u_prev + chop(h * v_prev, options), options);
        v = chop(v_prev - chop(h * u_prev, options), options);
        options.round = temp;
    else
        u = chop(u_prev + chop(h * v_prev, options), options);
        v = chop(v_prev - chop(h * u_prev, options), options);
    end
end