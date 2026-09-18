%% ME 30801 | Lab 1 - Control Volume Analysis | Fall 26
% Team #: 04
% 
% Team Members: Arnav Basu, Arnav Parag Shah, Matthew Ginart, Lucas
% Rodrigues
% 
% Date: 09/18/2026
% 
% Short Description: Control volume (momentum integral) analysis of drag 
% on a circular cylinder in cross-flow, computing drag force and 
% coefficient at six Reynolds numbers (5–102,700) from measured velocity 
% and pressure fields on a rectangular control surface.

%% 1. Initialize MATLAB

clear;
clc;
close all;

%% 2. Inputs and Flow Conditions

rho = 1.225;    % kg/m^3
mu = 1.81e-5;   % Pa-s
D = 0.1;        % m, cylinder diameter

% Add Re values here as more data comes in from the team
Re_list = [5, 50, 500, 3900, 34200, 102700];

% Preallocate storage
n = length(Re_list);
Cd_list = zeros(1,n);
Uinf_list = zeros(1,n);
CD_data = cell(1,n);   % stores the CD-face table for each Re, for plotting later

filepath = fileparts(matlab.desktop.editor.getActiveFilename);

%% 3. Loop over all Re cases: import, integrate, compute drag

for k = 1:n
    Re = Re_list(k);

    % ----- Import all 4 faces for this Re -----
    data_AB = readtable(fullfile(filepath, sprintf('Re%dAB.csv', Re)));
    data_BC = readtable(fullfile(filepath, sprintf('Re%dBC.csv', Re)));
    data_CD = readtable(fullfile(filepath, sprintf('Re%dCD.csv', Re)));
    data_AD = readtable(fullfile(filepath, sprintf('Re%dAD.csv', Re)));

    % ----- CS integral side AB -----
    y_AB = data_AB.Y_m_;
    u_AB = data_AB.XVelocityAvg_m_s_;
    momentum_flux_AB = trapz(y_AB, rho*u_AB.^2);

    % ----- CS integral side BC -----
    x_BC = data_BC.X_m_;
    u_BC = data_BC.XVelocityAvg_m_s_;
    v_BC = data_BC.YVelocityAvg_m_s_;
    momentum_flux_BC = trapz(x_BC, rho.*u_BC.*v_BC);

    % ----- CS integral side CD -----
    y_CD = data_CD.Y_m_;
    u_CD = data_CD.XVelocityAvg_m_s_;
    momentum_flux_CD = trapz(y_CD, rho*u_CD.^2);

    % ----- CS integral side AD -----
    x_AD = data_AD.X_m_;
    u_AD = data_AD.XVelocityAvg_m_s_;
    v_AD = data_AD.YVelocityAvg_m_s_;
    momentum_flux_AD = trapz(x_AD, rho.*u_AD.*v_AD);

    % ----- Surface forces (pressure), AB and CD only -----
    p_AB = data_AB.StaticPressureAvg_Pa_;
    pressure_force_AB = trapz(y_AB, p_AB);

    p_CD = data_CD.StaticPressureAvg_Pa_;
    pressure_force_CD = trapz(y_CD, p_CD);

    % ----- Drag force -----
    Fs = (pressure_force_AB + momentum_flux_AB) - (pressure_force_CD + momentum_flux_CD) ...
        + momentum_flux_AD - momentum_flux_BC;

    % ----- Freestream velocity: measured from AB inlet, not the Re formula -----
    % (formula-predicted U_inf did not match measured values consistently
    % across Re cases during testing, so we use the actual measured inlet
    % velocity for each Re instead)
    U_inf = mean(u_AB);

    % ----- Drag coefficient -----
    Cd = Fs / (0.5*rho*U_inf^2*D);

    % Store results
    Uinf_list(k) = U_inf;
    Cd_list(k) = Cd;
    CD_data{k} = data_CD;

    fprintf('Re=%d: U_inf=%.4e m/s, Fs=%.4e N/m, Cd=%.4f\n', Re, U_inf, Fs, Cd);
end

%% 4. Plotting

% Figure 1: y/D vs V/U_inf at the CD face, for all Re cases
figure;
hold on;
for k = 1:n
    data_CD = CD_data{k};
    y_over_D = data_CD.Y_m_ / D;
    u = data_CD.XVelocityAvg_m_s_;
    v = data_CD.YVelocityAvg_m_s_;
    V_mag = sqrt(u.^2 + v.^2);
    V_over_U = V_mag / Uinf_list(k);
    plot(V_over_U, y_over_D, '-o', 'DisplayName', sprintf('Re = %d', Re_list(k)));
end
hold off;
xlabel('V / U_\infty');
ylabel('y / D');
title('Velocity Profile at CD Face');
legend show;
grid on;

% Figure 2: Cd vs Re
figure;
semilogx(Re_list, Cd_list, 'o-', 'LineWidth', 1.5);
xlabel('Re');
ylabel('C_D');
title('Drag Coefficient vs Reynolds Number');
grid on;
