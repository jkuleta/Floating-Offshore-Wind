clear all; close all; clc;

%% Input Parameters
% Flow
FLOW.V0 = 10;                                                              % Inflow velocity [m/s]   
FLOW.rho = 1.225;                                                          % Air density [kg/m2]
FLOW.omega = 9*2*pi/60;                                                    % Rotational speed [rad/s]
         
% Rotor
ROTOR.TurbineInput = 'NREL5MW.xlsx';                                       % Turbine input file [.xlsx]

[ROTOR.r,~,~] = xlsread(ROTOR.TurbineInput,'NREL5MW','A3:A21');            % Radial positions [-]
[ROTOR.beta,~,~] = xlsread(ROTOR.TurbineInput,'NREL5MW','B3:B21');         % Blade twist [deg]
[ROTOR.chord,~,~] = xlsread(ROTOR.TurbineInput,'NREL5MW','C3:C21');        % Blade chord [m]
[~,~,ROTOR.airfoil] = xlsread(ROTOR.TurbineInput,'NREL5MW','D3:D21');      % Blade Airfoils [-]
ROTOR.R = 63;                                                              % Diameter [m]
ROTOR.D = 2*ROTOR.R;                                                       % Radius [m]
ROTOR.H = 90;                                                              % Hub height [m]
ROTOR.B = 3;                                                               % Number of blades [-]
ROTOR.theta_pitch = 0;                                                     % Blade (collective) pitch angle [deg]
ROTOR.sigma = ROTOR.chord*ROTOR.B./(2*pi*ROTOR.r);                         % Rotor solidity [-]   

% Airfoil
AIRFOIL.Cylinder1 = xlsread(ROTOR.TurbineInput,'Cylinder1','A3:D5');       % Root airfoil: alpha, Cl, Cd, Cm
AIRFOIL.Cylinder2 = xlsread(ROTOR.TurbineInput,'Cylinder2','A3:D5');
AIRFOIL.DU40 = xlsread(ROTOR.TurbineInput,'DU40','A3:D138');
AIRFOIL.DU35 = xlsread(ROTOR.TurbineInput,'DU35','A3:D137');
AIRFOIL.DU30 = xlsread(ROTOR.TurbineInput,'DU30','A3:D145');
AIRFOIL.DU25 = xlsread(ROTOR.TurbineInput,'DU25','A3:D142');
AIRFOIL.DU21 = xlsread(ROTOR.TurbineInput,'DU21','A3:D144');
AIRFOIL.NACA64 = xlsread(ROTOR.TurbineInput,'NACA64','A3:D129');           % Tip airfoil: alpha, Cl, Cd, Cm

% Simulation options
SIMULATION.error = 0.01;                                                   % Convergence criteria BEM
SIMULATION.dt = 0.1;                                                       % Time step [s]
SIMULATION.time = 0:SIMULATION.dt:300;                                     % Time series [s]
SIMULATION.taustar_nw = 0.5;                                               % Constants for dynamic inflow model 
SIMULATION.taustar_fw = 2;                                                 % Constants for dynamic inflow model 


%% Solve (steady) BEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TSR_range = linspace(1, 15, 50); % Tip Speed Ratio range
CP_values = zeros(size(TSR_range));
CT_values = zeros(size(TSR_range));

for i = 1:length(TSR_range)
    FLOW.omega = TSR_range(i) * FLOW.V0 / ROTOR.R; % Update rotational speed for each TSR
    [P, T, CP, CT, ~, ~] = function_BEM(ROTOR, AIRFOIL, FLOW, SIMULATION, []);
    CP_values(i) = CP;
    CT_values(i) = CT;
end

% Plot TSR-CP and TSR-CT curves
figure;
subplot(2, 1, 1);
plot(TSR_range, CP_values, 'b-', 'LineWidth', 1.5);
xlabel('Tip Speed Ratio (TSR)');
ylabel('Power Coefficient (C_P)');
%title('TSR vs C_P');
grid on;

subplot(2, 1, 2);
plot(TSR_range, CT_values, 'r-', 'LineWidth', 1.5);
xlabel('Tip Speed Ratio (TSR)');
ylabel('Thrust Coefficient (C_T)');
%title('TSR vs C_T');
grid on;

%% Solve (unsteady) BEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSteps = length(SIMULATION.time);

% Floating motion parameters from literature
A_surge = 8;              f_surge = 1 / 107;    % ±4 m, 107 s period
A_pitch = 3 * pi/180;     f_pitch = 1 / 32.5;   % ±3 deg, 32.5 s period
A_yaw   = 3 * pi/180;     f_yaw   = 1 / 80.8;     % ±3 deg, 80.8 s period

% Convert pitch/yaw to linear velocity at rotor radius
V_surge = A_surge * 2*pi*f_surge * cos(2*pi*f_surge*SIMULATION.time);                    % [m/s]
V_pitch = ROTOR.R * A_pitch * 2*pi*f_pitch * cos(2*pi*f_pitch*SIMULATION.time);          % [m/s]
V_yaw   = ROTOR.R * A_yaw   * 2*pi*f_yaw   * cos(2*pi*f_yaw*SIMULATION.time);             % [m/s]

% Initialize storage
CP_time = zeros(1, nSteps);
CT_time = zeros(1, nSteps);
P_time  = zeros(1, nSteps);
T_time  = zeros(1, nSteps);
a_hist  = zeros(length(ROTOR.r), nSteps);
ap_hist = zeros(length(ROTOR.r), nSteps);

% Initial conditions for induction
PREVIOUSTIME.a_new  = 0.3 * ones(length(ROTOR.r),1);
PREVIOUSTIME.ap_new = 0.0 * ones(length(ROTOR.r),1);

% Time loop
for t = 1:nSteps
    FLOW.V_surge = V_surge(t);
    FLOW.V_pitch = V_pitch(t);
    FLOW.V_yaw   = V_yaw(t);

    [P, T, CP, CT, a_new, ap_new] = function_BEM(ROTOR, AIRFOIL, FLOW, SIMULATION, PREVIOUSTIME);

    % Store time history
    P_time(t)  = P;
    T_time(t)  = T;
    CP_time(t) = CP;
    CT_time(t) = CT;
    a_hist(:,t)  = a_new;
    ap_hist(:,t) = ap_new;

    % Update induction history
    PREVIOUSTIME.a_new  = a_new;
    PREVIOUSTIME.ap_new = ap_new;
end

% Plot aerodynamic response to floating motion
figure;
subplot(2,1,1);
plot(SIMULATION.time, P_time/1e6, 'b'); % MW
xlabel('Time [s]');
ylabel('Power [MW]');
%title('Aerodynamic Power Under Floating Motion');
grid on;

subplot(2,1,2);
plot(SIMULATION.time, T_time/1e3, 'r'); % kN
xlabel('Time [s]');
ylabel('Thrust [kN]');
%title('Aerodynamic Thrust Under Floating Motion');
grid on;


% Plot surge, pitch, and yaw movements
figure;
subplot(3, 1, 1);
plot(SIMULATION.time, V_surge, 'g', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Surge Velocity [m/s]');
grid on;

subplot(3, 1, 2);
plot(SIMULATION.time, V_pitch, 'm', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Pitch Velocity [m/s]');
grid on;

subplot(3, 1, 3);
plot(SIMULATION.time, V_yaw, 'c', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Yaw Velocity [m/s]');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


