%% ME 30801 | Lab 1 - Control Volume Analysis | Fall 26
% Team #: 
% 
% Team Members: Member #1, Member #2, Member #3, Member #4
% 
% Date:
% 
% Short Description: 

%% 1. Initialize MATLAB

% Start from a clean workspace.
clear;
clc;
close all;

%% 2. Importing CSV files
% Use this section to import your CSV files into MATLAB

[filepath, filename, ext] = fileparts(matlab.desktop.editor.getActiveFilename);

% load AB side data
datafile = 'Re34200AB.csv';
inputfile_AB = fullfile(filepath, datafile);

% load BC side data
datafile = 'Re34200BC.csv';
inputfile_BC = fullfile(filepath, datafile);

% load CD side data
datafile = 'Re34200CD.csv';
inputfile_CD = fullfile(filepath, datafile);

% load AD side data
datafile = 'Re34200AD.csv';
inputfile_AD = fullfile(filepath, datafile);

%% 3. Inputs and Flow Conditions
% Define all physical quantities required for the analysis.

rho = 1.225;    % kg/m^3
mu = 1.81e-5;   % Pa-s

D = 0.1;        % m, cylinder diameter

Re = 34200;
U_int_Re34200 = Re*mu/rho/D

%% 4. Extracting data from CSV files

data_AB = readtable(inputfile_AB);
data_BC = readtable(inputfile_BC);
data_CD = readtable(inputfile_CD);
data_AD = readtable(inputfile_AD);

%% 4. Apply the COLM in the x-direction to calculate the drag force at each Re
% dP/dt*int_CV(rho*V)dV + int_CS(rho*V)(V.n)dA = F_B + F_S

% Write down the assumptions:
% 1. Steady state
% 2. Long cylinder (2D flow)
% 3. Air properties at standard atmospheric conditions
% 4. No body forces in x-direction

% ----- CS integral side AB at each Re -----
y_AB = data_AB.Y_m_;
u_AB = data_AB.XVelocityAvg_m_s_;

momentum_flux_AB = trapz(y_AB, rho*u_AB.^2); % (N)

% ----- CS integral side BC at each Re -----
x_BC = data_BC.X_m_;
u_BC = data_BC.XVelocityAvg_m_s_;
v_BC = data_BC.YVelocityAvg_m_s_;

momentum_flux_BC = trapz(x_BC, rho.*u_BC.*v_BC); % (N)

% ----- CS integral side CD at each Re -----
y_CD = data_CD.Y_m_;
u_CD = data_CD.XVelocityAvg_m_s_;

momentum_flux_CD = trapz(y_CD, rho*u_CD.^2); % (N)

% ----- CS integral side AD at each Re -----
x_AD = data_AD.X_m_;
u_AD = data_AD.XVelocityAvg_m_s_;
v_AD = data_AD.YVelocityAvg_m_s_;

momentum_flux_AD = trapz(x_AD, rho.*u_AD.*v_AD); % (N)

% ----- Surface forces side AB at each Re -----
p_AB = data_AB.StaticPressureAvg_Pa_;   % verify exact column name first
pressure_force_AB = trapz(y_AB, p_AB);

% ----- Surface forces side CD at each Re -----
p_CD = data_CD.StaticPressureAvg_Pa_;
pressure_force_CD = trapz(y_CD, p_CD);

% ----- Drag force at each Re -----
Fs = (pressure_force_AB + momentum_flux_AB) - (pressure_force_CD + momentum_flux_CD) ...
    + momentum_flux_AD - momentum_flux_BC

%% 5. Calculating drag coefficient at each Re

Cd_Re34200 = Fs/(0.5*rho*U_int_Re34200^2*D)

%% 6. Plotting
% 1) Present a plot of the right side of your control volume (point C-D), 
% showing the y-coordinates of the probe divided by the diameter of the 
% cylinder vs the dimensionless velocity, defined as the total velocity 
% magnitude V divided by the freestream velocity U, for all Reynolds number 
% cases.

% Figure 1


% 2) Present the Cd vs Re plot over different Re with meaningful descriptions

% Figure 2