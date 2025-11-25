function newmark_beta_sdoff()

clc; close all;

disp('-----------------------------------------------------------');
disp('        SDOF TIME INTEGRATION — NEWMARK β METHOD           ');
disp('-----------------------------------------------------------');

%% ====================== MASS UNIT MENU ==========================

mass_unit = menu('Select MASS unit:', ...
    'kg (mass)', ...
    'kip (weight)', ...
    'kN (weight)' );

switch mass_unit
    case 1
        m = input('Enter mass (kg): ');
    case 2
        W_kip = input('Enter weight (kip): ');
        W_N = W_kip * 4448.22;      % kip → N
        m = W_N / 9.80665;          % N / g
    case 3
        W_kN = input('Enter weight (kN): ');
        W_N  = W_kN * 1000;         % kN → N
        m = W_N / 9.80665;
end

%% ====================== STIFFNESS UNIT MENU =====================

k_unit = menu('Select STIFFNESS unit:', ...
    'N/m', ...
    'kN/m', ...
    'kip/ft' );

switch k_unit
    case 1
        k = input('Enter stiffness (N/m): ');

    case 2
        k_kNm = input('Enter stiffness (kN/m): ');
        k = k_kNm * 1000;

    case 3
        k_kipft = input('Enter stiffness (kip/ft): ');
        k = k_kipft * 14593.9;    % kip/ft → N/m
end

%% ====================== DAMPING RATIO ===========================

zeta = input('Enter damping ratio ζ (e.g., 0.05): ');

% Damping coefficient
c = 2 * zeta * sqrt(k*m);

%% ====================== TIME STEP ===============================

dt = input('Enter time step Δt (seconds): ');

%% ====================== LOAD INPUT ===============================

disp('Select Excel file containing: time | load');
[file,path] = uigetfile({'*.xlsx;*.xls'});
raw = readmatrix(fullfile(path,file));

time = raw(:,1);
load_raw = raw(:,2);

%% ====================== LOAD UNIT MENU ==========================

load_unit = menu('Select LOAD type:', ...
    'Force: N', ...
    'Force: kN', ...
    'Force: kip', ...
    'Ground Acceleration: m/s^2', ...
    'Ground Acceleration: g');

% Ask if negative sign is needed (for relative motion)
neg_flag = menu('Apply negative sign for F = - m * a ?', ...
    'YES (-m a)', ...
    'NO (+m a)');

if neg_flag == 1
    sign_factor = -1;
else
    sign_factor = +1;
end

%% =============== CONVERT LOAD / ACCELERATION ====================

switch load_unit
    
    % ---------------- FORCE INPUTS ---------------------
    case 1   % N
        p = load_raw;      % already N

    case 2   % kN → N
        p = load_raw * 1000;

    case 3   % kip → N
        p = load_raw * 4448.22;

    % ---------------- ACCELERATION INPUTS --------------
    case 4   % m/s²
        a = load_raw;                     % acceleration
        p = sign_factor * m * a;          % F = ± m a

    case 5   % g → m/s² → force
        a = load_raw * 9.80665;           % g → m/s²
        p = sign_factor * m * a;          % F = ± m a
end

N = length(time);

%% =================== INITIAL CONDITIONS =======================

v  = zeros(N,1);
vd = zeros(N,1);
vdd = zeros(N,1);

vdd(1) = (p(1) - c*vd(1) - k*v(1)) / m;

%% =================== NEWMARK CONSTANTS =========================

beta  = 1/4;
gamma = 1/2;

a0 = 1/(beta*dt^2);
a1 = gamma/(beta*dt);
a2 = 1/(beta*dt);
a3 = 1/(2*beta) - 1;
a4 = gamma/beta - 1;
a5 = dt*(gamma/(2*beta) - 1);

k_eff = k + a1*c + a0*m;

%% ====================== TIME INTEGRATION =======================

for i = 1:N-1

    p_eff = p(i+1) ...
            + m*( a0*v(i) + a2*vd(i) + a3*vdd(i) ) ...
            + c*( a1*v(i) + a4*vd(i) + a5*vdd(i) );

    % Solve for displacement
    v(i+1) = p_eff / k_eff;

    % Update acceleration and velocity
    vdd(i+1) = a0*(v(i+1)-v(i)) - a2*vd(i) - a3*vdd(i);
    vd(i+1)  = vd(i) + dt*((1-gamma)*vdd(i) + gamma*vdd(i+1));

end

%% ======================= PLOTTING ===============================

figure;

% -------- DISPLACEMENT PLOT --------
subplot(3,1,1);
plot(time, v, 'LineWidth',1.5); hold on;
[ymax, idx_max] = max(abs(v));
plot(time(idx_max), v(idx_max), 'ro', 'MarkerSize',8, 'LineWidth',1.5);
text(time(idx_max), v(idx_max), sprintf('  Max = %.4g m', ymax), ...
    'FontSize',10, 'Color','r', 'FontWeight','bold');
xlabel('Time (s)'); ylabel('Displacement (m)');
title('SDOF Response using Newmark-\beta'); grid on;

% -------- VELOCITY PLOT --------
subplot(3,1,2);
plot(time, vd, 'LineWidth',1.5); hold on;
[ymax_v, idx_max_v] = max(abs(vd));
plot(time(idx_max_v), vd(idx_max_v), 'ro', 'MarkerSize',8, 'LineWidth',1.5);
text(time(idx_max_v), vd(idx_max_v), sprintf('  Max = %.4g m/s', ymax_v), ...
    'FontSize',10, 'Color','r', 'FontWeight','bold');
xlabel('Time (s)'); ylabel('Velocity (m/s)');
grid on;

% -------- ACCELERATION PLOT --------
subplot(3,1,3);
plot(time, vdd, 'LineWidth',1.5); hold on;
[ymax_a, idx_max_a] = max(abs(vdd));
plot(time(idx_max_a), vdd(idx_max_a), 'ro', 'MarkerSize',8, 'LineWidth',1.5);
text(time(idx_max_a), vdd(idx_max_a), sprintf('  Max = %.4g m/s^2', ymax_a), ...
    'FontSize',10, 'Color','r', 'FontWeight','bold');
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
grid on;

Period = 2*pi*sqrt(m/k);

sgtitle(sprintf('Newmark-\\beta Response (m = %.2f kg,  T = %.2f s)', m, Period), ...
    'FontSize',14, 'FontWeight','bold');

end
