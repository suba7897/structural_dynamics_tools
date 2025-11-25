function rayleigh_shape_function()
% ================================================================
%   RAYLEIGH / STODOLA FIRST-MODE SHAPE ESTIMATION
%   - User enters mass + stiffness
%   - Initial shape is specified
%   - Iterations stored and plotted on one plot
%   - Initial vs Final shape plotted as subplot
% ================================================================

clc; close all;

disp('-------------------------------------------------------');
disp('     RAYLEIGH (STODOLA) MODE SHAPE ESTIMATION TOOL');
disp('-------------------------------------------------------');

%% ===================== USER INPUT =======================

m = input('Enter story masses m (kg) as a row vector (include base as zero):\n e.g. [0 263 263 ...]\n\nm = ');
k = input('\nEnter story stiffness k (N/m) as a row vector (include base as zero):\n e.g. [0 9391 64065 ...]\n\nk = ');

m = m(:);
k = k(:);

N = length(m);

%% =============== INITIAL SHAPE FUNCTION ==================

% For your specific case you used [0 1 1 1 1]
shape_f = [0, 1, 1, 1, 1]';     % initial assumed shape
% Alternative (linear) would be:
% shape_f = linspace(0,1,N)';

%% =============== STORAGE FOR ITERATIONS ==================
shape_history = shape_f;         % each column = one iteration shape

%% =============== INITIAL ESTIMATE ========================

omega_old = angular_frequency(m, shape_f, k);

fprintf('\nInitial angular frequency ω₀ = %.4f rad/s\n', omega_old);

%% ===================== ITERATION LOOP =====================

iter = 1;
tol = 1e-6;

while true

    % Rayleigh internal force (with current ω)
    F_I = m .* (omega_old^2) .* shape_f;

    % Shear force from top-down
    SF = flip(cumsum(flip(F_I)));

    % Story displacements
    story_dis = SF ./ k;
    story_dis(isinf(story_dis)) = 0;

    % Floor deflections
    tot_dis = cumsum(story_dis);

    % Normalize shape
    shape_new = tot_dis ./ tot_dis(end);

    % New frequency from Rayleigh quotient
    omega_new = angular_frequency(m, shape_new, k);

    fprintf('Iteration %d: ω = %.6f rad/s\n', iter, omega_new);

    % Store iteration
    shape_history = [shape_history, shape_new];

    % Convergence check
    if abs(omega_new - omega_old) < tol
        disp('Convergence achieved.');
        break;
    end

    % Prepare next iteration
    shape_f = shape_new;
    omega_old = omega_new;
    iter = iter + 1;

end

%% ===================== FINAL RESULTS =====================

period_final = 2*pi/omega_new;
fprintf('\nFINAL RESULTS\n');
fprintf('  Final angular frequency ω = %.6f rad/s\n', omega_new);
fprintf('  Final period T = %.6f s\n', period_final);

%% ======================== PLOTS ===========================
figure('Name','Mode Shape Iterations','Position',[200 80 1200 800]);

% ------------- LEFT SUBPLOT: All Iterations ---------------
subplot(1,2,1);
hold on;

num_iter = size(shape_history,2);
floors = (0:N-1)';      % 0 = base, 1 = 1st floor, etc.

for i = 1:num_iter
    plot(shape_history(:,i), floors, 'LineWidth', 1.5);
end

set(gca, 'YDir', 'normal');  % bottom floor at bottom
xlabel('Mode Shape Value');
ylabel('Floor Level');
title('Rayleigh–Stodola Shape Convergence (All Iterations)');
grid on;
legend_strings = arrayfun(@(x) sprintf('Iter %d', x-1), 1:num_iter, 'UniformOutput', false);
legend(legend_strings, 'Location', 'bestoutside');

hold off;

% ------------- RIGHT SUBPLOT: Initial vs Final ---------------
subplot(1,2,2);
plot(shape_history(:,1), floors, 'k--o', 'LineWidth',1.4); hold on;
plot(shape_history(:,end), floors, 'r-o', 'LineWidth',1.8);

set(gca,'YDir','normal');
xlabel('Mode Shape Value');
ylabel('Floor Level');
legend('Initial Shape', 'Final Converged Shape', 'Location', 'best');
title('Initial vs Final Shape Function');
grid on;

sgtitle('Rayleigh Method Mode Shape Determination');

end % end main function

% ===============================================================
%   SUBFUNCTION: Compute angular frequency from Rayleigh quotient
% ===============================================================
function omega = angular_frequency(m, shape_f, k)
    % Ensure column vectors
    m       = m(:);
    shape_f = shape_f(:);
    k       = k(:);

    % Inter-story drift (difference of shape between levels)
    delta_sf = [0; diff(shape_f)];

    % Strain energy and kinetic energy
    U_max   = sum(k .* (delta_sf.^2));   % ∝ strain energy
    KE_max  = sum(m .* (shape_f.^2));    % ∝ kinetic energy

    % Rayleigh frequency
    omega = sqrt(U_max / KE_max);
end
