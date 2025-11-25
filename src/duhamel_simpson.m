function duhamel_simpson_flex()

clc; close all;

disp('-----------------------------------------------------------');
disp('              D U H A M E L   I N T E G R A L');
disp('      Flexible Units: Mass | Stiffness | Load/Accel');
disp('-----------------------------------------------------------');

%% ======================= UNIT MENUS ==========================

mass_unit = menu('Select MASS unit (weight):', ...
    'kg', 'lb', 'kip', 'kN');

stiff_unit = menu('Select STIFFNESS unit:', ...
    'kN/m', 'kip/ft');

load_unit = menu('Select LOAD / ACCELERATION unit:', ...
    'kN  (force)', ...
    'kip (force)', ...
    'm/s^2 (acc)', ...
    'g (acc)', ...
    'ft/s^2 (acc)' );

dt      = input('Enter Δt (seconds): ');
T_final = input('Enter total simulation duration (seconds): ');
zeta    = input('Enter damping ratio (e.g., 0.05): ');

%% ================= MASS (WEIGHT) → kg =========================

switch mass_unit
    case 1      % kg (weight)
        W = input('Enter weight (kg): ');
        mass = W;  mass = mass;  % (kg is already mass)

    case 2      % lb → N → kg
        W_lb = input('Enter weight (lb): ');
        W_N  = W_lb * 4.44822;
        mass = W_N / 9.80665;

    case 3      % kip → N → kg
        W_kip = input('Enter weight (kip): ');
        W_N   = W_kip * 4448.22;
        mass  = W_N / 9.80665;

    case 4      % kN → N → kg
        W_kN = input('Enter weight (kN): ');
        W_N  = W_kN * 1000;
        mass = W_N / 9.80665;
end

%% ================== STIFFNESS → N/m ===========================

switch stiff_unit
    case 1      % kN/m → N/m
        k_input  = input('Enter stiffness (kN/m): ');
        stiffness = k_input * 1000;

    case 2      % kip/ft → N/m
        k_kip_ft  = input('Enter stiffness (kip/ft): ');
        stiffness = k_kip_ft * 14593.9;
end

%% ==================== LOAD FILE ===============================

disp('Select Excel file with: Time | Load/Accel');
[file, path] = uigetfile({'*.xlsx;*.xls'});
raw = readmatrix(fullfile(path, file));

time_raw = raw(:,1);
load_raw = raw(:,2);

%% =============== CONVERT LOAD/ACCELERATION =====================

switch load_unit

    case 1      % kN (force)
        load_conv = load_raw * 1000;         % → N

    case 2      % kip (force)
        load_conv = load_raw * 4448.22;      % → N

    case 3      % m/s² acceleration
        load_conv = -(mass * load_raw);      % → N  (negative sign)

    case 4      % g acceleration
        load_conv = -(mass * (load_raw * 9.80665));  % g → m/s² → force

    case 5      % ft/s² acceleration
        load_conv = -(mass * (load_raw * 0.3048));   % ft/s² → m/s² → force

end

%% ====================== INTERPOLATE LOAD =======================

time = (0:dt:T_final)';
force = interp1(time_raw, load_conv, time, 'linear', 0);

%% ===================== SYSTEM PROPERTIES =======================

omega_n = sqrt(stiffness / mass);
omega_d = omega_n * sqrt(1 - zeta^2);

%% ====================== DUHAMEL INTEGRAL =======================

Npts = length(time);
p = force ./ mass;

A = zeros(Npts,1);
B = zeros(Npts,1);

for N = 3:2:Npts
    tN = time(N);
    tau = time(1:N);
    p_tau = p(1:N);

    expTerm = exp( zeta*omega_n*(tau - tN) );

    w = ones(N,1);
    w(2:2:end-1) = 4;
    w(3:2:end-2) = 2;

    A(N) = (dt/(3*omega_d)) * sum( w .* p_tau .* expTerm .* cos(omega_d*tau) );
    B(N) = (dt/(3*omega_d)) * sum( w .* p_tau .* expTerm .* sin(omega_d*tau) );
end

A = fillmissing(A,'linear');
B = fillmissing(B,'linear');

response_m = A.*sin(omega_d*time) - B.*cos(omega_d*time);

%% ========== CONVERT RESPONSE BACK TO ORIGINAL DISPLACEMENT UNITS ========

switch stiff_unit
    case 1  % kN/m → displacement in meters
        response_plot = response_m;
        resp_label = 'Response (m)';

    case 2  % kip/ft → feet
        response_plot = response_m * 3.28084;   % m → ft
        resp_label = 'Response (ft)';
end

idx_start = find(response_plot ~= 0, 1, 'first');
if isempty(idx_start), idx_start = 1; end

%% ============================== PLOTS ===========================

figure;

subplot(2,1,1);
switch load_unit
    case 1
        plot(time_raw, load_raw, 'LineWidth',1.6); ylabel('Load (kN)');
    case 2
        plot(time_raw, load_raw, 'LineWidth',1.6); ylabel('Load (kip)');
    case 3
        plot(time_raw, load_raw, 'LineWidth',1.6); ylabel('Accel (m/s^2)');
    case 4
        plot(time_raw, load_raw, 'LineWidth',1.6); ylabel('Accel (g)');
    case 5
        plot(time_raw, load_raw, 'LineWidth',1.6); ylabel('Accel (ft/s^2)');
end
xlabel('Time (s)');
title('Input Load / Ground Acceleration');
grid on;

subplot(2,1,2);
plot(time(idx_start:end), response_plot(idx_start:end), 'r', 'LineWidth',1.6);
xlabel('Time (s)');
ylabel(resp_label);
title('SDOF Response (Duhamel Integral)');
grid on;

disp('Duhamel analysis complete.');

end
