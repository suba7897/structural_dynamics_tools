function duhamel_response_spectrum_mass()

clc; close all;

disp('-----------------------------------------------------------');
disp('              S D O F   R E S P O N S E   S P E C T R U M');
disp('      Flexible Units: Mass | Load/Accel | Period Array');
disp('-----------------------------------------------------------');

%% ======================= UNIT MENUS ==========================

mass_unit = menu('Select MASS unit (weight):', ...
    'kg', 'lb', 'kip', 'kN');

load_unit = menu('Select LOAD / ACCELERATION unit:', ...
    'kN  (force)', ...
    'kip (force)', ...
    'm/s^2 (acc)', ...
    'g (acc)', ...
    'ft/s^2 (acc)' );

dt      = input('Enter Δt (seconds): ');
T_final = input('Enter total simulation duration (seconds): ');
zeta    = input('Enter damping ratio (e.g., 0.05): ');
T_array = input('Enter periods array (e.g., [0.1 0.2 0.3 1]): ');

%% ================= MASS INPUT (WEIGHT) → kg ====================

switch mass_unit
    case 1
        W = input('Enter weight (kg): ');
        mass = W;
    case 2
        W_lb = input('Enter weight (lb): ');
        W_N  = W_lb * 4.44822;
        mass = W_N / 9.80665;
    case 3
        W_kip = input('Enter weight (kip): ');
        W_N   = W_kip * 4448.22;
        mass  = W_N / 9.80665;
    case 4
        W_kN = input('Enter weight (kN): ');
        W_N  = W_kN * 1000;
        mass = W_N / 9.80665;
end

%% ==================== LOAD FILE ===============================

disp('Select Excel file with: Time | Load/Accel');
[file, path] = uigetfile({'*.xlsx;*.xls'});
raw = readmatrix(fullfile(path, file));

time_raw = raw(:,1);
load_raw = raw(:,2);

%% =============== CONVERT LOAD/ACCELERATION =====================

switch load_unit
    case 1
        load_conv = load_raw * 1000;
    case 2
        load_conv = load_raw * 4448.22;
    case 3
        load_conv = -(mass * load_raw);
    case 4
        load_conv = -(mass * (load_raw * 9.80665));
    case 5
        load_conv = -(mass * (load_raw * 0.3048));
end

%% ====================== INTERPOLATE LOAD =======================

time = (0:dt:T_final)';
force = interp1(time_raw, load_conv, time, 'linear', 0);

%% =================== PRE-ALLOCATE STORAGE ======================

N_per = length(T_array);
max_disp = zeros(N_per,1);
max_velocity = zeros(N_per,1);
max_acceleration = zeros(N_per,1); 
response_history = zeros(length(time), N_per);

%% ===============================================================
%   LOOP THROUGH PERIODS – Compute Response Spectrum & Responses
% ===============================================================

for p_i = 1:N_per

    Tn = T_array(p_i);

    omega_n = 2*pi / Tn;
    omega_d = omega_n * sqrt(1 - zeta^2);

    stiffness = mass * omega_n^2;

    Npts = length(time);
    p = force ./ mass;

    A = zeros(Npts,1);
    B = zeros(Npts,1);

    for N = 3:2:Npts
        tN = time(N);
        tau = time(1:N);
        p_tau = p(1:N);

        expTerm = exp(zeta*omega_n*(tau - tN));

        w = ones(N,1);
        w(2:2:end-1) = 4;
        w(3:2:end-2) = 2;

        A(N) = (dt/(3*omega_d)) * sum( w .* (p_tau .* expTerm .* cos(omega_d*tau)) );
        B(N) = (dt/(3*omega_d)) * sum( w .* (p_tau .* expTerm .* sin(omega_d*tau)) );
    end

    A = fillmissing(A,'linear');
    B = fillmissing(B,'linear');

    response_m = A.*sin(omega_d*time) - B.*cos(omega_d*time);

    % Store max
    max_disp(p_i) = max(abs(response_m));
    max_velocity(p_i) = omega_n*max(abs(response_m));
    max_acceleration(p_i) = omega_n^2*max(abs(response_m));

    % Store full response
    response_history(:,p_i) = response_m;

end

%% ======================= SAVE ONCE (FIXED) =====================

output_full = [time, response_history];

[savefile, savepath] = uiputfile({'*.xlsx'; '*.csv'}, ...
    'Save FULL response time-history file as');

if savefile ~= 0
    fullpath = fullfile(savepath, savefile);

    if contains(savefile,'.xlsx')
        writematrix(output_full, fullpath, 'Sheet', 1);
    else
        writematrix(output_full, fullpath);
    end

    fprintf('\nSaved single output file:\n%s\n', fullpath);
end

%% ========================== PLOT ONCE ===============================

figure;
subplot(3,1,1)
plot(T_array, max_disp, 'r-o', 'LineWidth', 2);
xlabel('Period (s)');
ylabel('Spectral Displacement (SD) (m)');
title('SDOF Displacement Response Spectrum');
grid on;

subplot(3,1,2)
plot(T_array, max_velocity, 'r-o', 'LineWidth', 2);
xlabel('Period (s)');
ylabel('Spectral Velocity (SD) (m/s)');
title('SDOF Velocity Response Spectrum');
grid on;

subplot(3,1,3)
plot(T_array, max_acceleration, 'r-o', 'LineWidth', 2);
xlabel('Period (s)');
ylabel('Spectral Acceleration (SD) (m/s2)');
title('SDOF Acceleration Response Spectrum');
grid on;


end
