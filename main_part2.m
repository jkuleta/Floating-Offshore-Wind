clear all; close all; clc;

%% Input Parameters
% Flow
FLOW.V0 = 10;                                                              % Inflow velocity [m/s]   
FLOW.rho = 1.225;                                                          % Air density [kg/m2]
FLOW.omega = 9*2*pi/60;                                                    % Rotational speed [rad/s]

% Waves - computed later baed on the sea-state
sea_state = input(['Enter the number of the sea state to use (1-8):\n' ...
    '1: H=0.09m, T=2.0s\n' ...
    '2: H=0.67m, T=4.8s\n' ...
    '3: H=1.40m, T=6.5s\n' ...
    '4: H=2.44m, T=8.1s\n' ...
    '5: H=3.66m, T=9.7s\n' ...
    '6: H=5.49m, T=11.3s\n' ...
    '7: H=9.14m, T=13.6s\n' ...
    '8: H=15.24m, T=17s\n' ...
    'Your choice: ']);
fprintf('Selected sea state: %d\n', sea_state);

if sea_state == 1
    load('OUTPUT_S1.mat');
    WAVES.H = 0.09;                                                            % Wave height [m]
    WAVES.T = 2.0;                                                             % Wave period [s]
elseif sea_state == 2
    load('OUTPUT_S2.mat');
    WAVES.H = 0.67;
    WAVES.T = 4.8;
elseif sea_state == 3
    load('OUTPUT_S3.mat');
    WAVES.H = 1.40;                                                            % Wave height [m]
    WAVES.T = 6.5;
elseif sea_state == 4
    load('OUTPUT_S4.mat');
    WAVES.H = 2.44;                                                            % Wave height [m]
    WAVES.T = 8.1;  
elseif sea_state == 5
    load('OUTPUT_S5.mat');
    WAVES.H = 3.66;                                                            % Wave height [m]
    WAVES.T = 9.7;      
elseif sea_state == 6
    load('OUTPUT_S6.mat');
    WAVES.H = 5.49;                                                            % Wave height [m]
    WAVES.T = 11.3;
elseif sea_state == 7
    load('OUTPUT_S7.mat');
    WAVES.H = 9.14;                                                            % Wave height [m]
    WAVES.T = 13.6;
elseif sea_state == 8
    load('OUTPUT_S8.mat');
    WAVES.H = 15.24;                                                            % Wave height [m]
    WAVES.T = 17;% Wave height [m]
else
    error('Invalid sea state selected. Please choose below 8.');
end

use_fast = input('Do you want to use hydrodynamic data from OpenFAST? (y/n): ', 's');
if strcmpi(use_fast, 'y')
    hydro_from_FAST = true;
    disp("Using hydrodynamic data from OpenFAST.");
else
    hydro_from_FAST = false;
    disp("Own hydrodynamic model will be used.");
end


WAVES.d = 200;                                                             % Water depth (positive) [m]
WAVES.WA = WAVES.H/2;                                                      % Wave amplitude [m]
WAVES.omega = 2*pi/WAVES.T;                                                % Wave frequency [rad/s]
WAVES.f = 1/WAVES.T;                                                       % Wave frequency [Hz]
WAVES.k = fsolve (@(ks) WAVES.omega.^2-9.81*ks*tanh(ks*WAVES.d),...        % Wave number [1/m]
    0.01,optimset('Display','off'));  
WAVES.L = 2*pi/WAVES.k;                                                    % Wave length [m]

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

% Floater
[FLOATER.M1,...                                                            % Normal mass matrix
    FLOATER.A1,...                                                         % Added-mass matrix (depends on the wave frequency)
    FLOATER.B,...                                                          % Damping matrix (depends on the wave frequency)
    FLOATER.C1,...                                                         % Hydrostatic matrix 
    FLOATER.K11,...                                                        % Mooring line matrix 
    FLOATER.Fhydrodynamic,FLOATER.Fhydrodynamicphase]...                   % Hydrodynamic forces 
    = function_floater(WAVES);

% Simulation options
SIMULATION.error = 0.01;                                                   % Convergence criteria BEM
SIMULATION.dt = 0.1;                                                       % Time step [s]
SIMULATION.time = 0:SIMULATION.dt:600;                                     % Time series [s]
SIMULATION.taustar_nw = 0.5;                                               % Constants for dynamic inflow model 
SIMULATION.taustar_fw = 2;                                                 % Constants for dynamic inflow model 

% Previous time for BEM
PREVIOUSTIME.a_new = 0; % Initialize previous time array
PREVIOUSTIME.ap_new = 0; % Initialize previous time array

%% State-space solution (TASK 8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The equation of motion of this floater-turbine system has the shape:
%%% (M1+A1)x'' + Bx' + (C1+K1)x = Fext
%%%      M1 is the mass matrix
%%%      A1 is the added-mass matrix, and depends on the wave frequency
%%%      B is the damping matrix, and depends on the wave frequency
%%%      C1 is the hydrostatic matrix
%%%      K1 is the mooring line matrix
%%%      Fext is the external forces: aerodynamic and the hydrodynamic forces
% Initial conditions
x0 = zeros(6, 1); % Initial states (3 translations and 3 rotations)
x_dot0 = zeros(6, 1); % Initial velocities

% Time integration
n = length(SIMULATION.time);
x = zeros(6, n); % States
x_dot = zeros(6, n); % Velocities
x_ddot = zeros(6, n); % Accelerations

% Assign initial conditions
x(:, 1) = x0;
x_dot(:, 1) = x_dot0;

% External forces (assumed zero for now)
Fext = zeros(6, 1);

% Effective matrices
M_eff = FLOATER.M1 + FLOATER.A1;
K_eff = FLOATER.C1 + FLOATER.K11;

% Time-stepping loop
for i = 1:n-1
    % Solve for acceleration
    x_ddot(:, i) = M_eff \ (Fext - FLOATER.B * x_dot(:, i) - K_eff * x(:, i));
    
    % Update velocity and position using explicit Euler method
    x_dot(:, i+1) = x_dot(:, i) + SIMULATION.dt * x_ddot(:, i);
    x(:, i+1) = x(:, i) + SIMULATION.dt * x_dot(:, i);
end

fprintf("Check for 0 displacement oscillation.\n");

% Plot results
figure;
subplot(3, 1, 1);
plot(SIMULATION.time, x(1, :), "LineWidth", 1.5);
xlabel('Time [s]');
ylabel('Surge [m]');
title('Surge Motion');
grid on;

subplot(3, 1, 2);
plot(SIMULATION.time, x(2, :), "LineWidth", 1.5);
xlabel('Time [s]');
ylabel('Sway [m]');
title('Sway Motion');

subplot(3, 1, 3);
plot(SIMULATION.time, x(3, :),  "LineWidth", 1.5);
xlabel('Time [s]');
ylabel('Heave [m]');
title('Heave Motion');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Free Decay Tests (TASK 10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions for free decay tests for all 6 degrees of freedom
initial_conditions = {
    'Surge [m]',    [1; 0; 0; 0; 0; 0];  % Initial displacement in surge
    'Sway [m]',     [0; 1; 0; 0; 0; 0];  % Initial displacement in sway
    'Heave [m]',    [0; 0; 1; 0; 0; 0];  % Initial displacement in heave
    'Roll [rad]',     [0; 0; 0; 1; 0; 0];  % Initial displacement in roll
    'Pitch [rad]',    [0; 0; 0; 0; 1; 0];  % Initial displacement in pitch
    'Yaw [rad]',      [0; 0; 0; 0; 0; 1];  % Initial displacement in yaw
};

% Initialize array to store periods
periods = zeros(length(initial_conditions), 1);
damping_linear = zeros(length(initial_conditions), 1);
damping_quadratic = zeros(length(initial_conditions), 1);

% Initialize state variables
x = zeros(6, n, length(initial_conditions));
x_dot = zeros(6, n, length(initial_conditions));
x_ddot = zeros(6, n, length(initial_conditions));

fprintf("Running free-decay tests...");


% Loop through each test to compute motion
for test = 1:length(initial_conditions)
    % Extract initial condition
    x0 = initial_conditions{test, 2};
    x_dot0 = zeros(6, 1); % Initial velocities

    % Assign initial conditions
    x(:, 1, test) = x0;
    x_dot(:, 1, test) = x_dot0;

    % Time-stepping loop
    for i = 1:n-1
        % Solve for acceleration
        x_ddot(:, i, test) = M_eff \ (-FLOATER.B * x_dot(:, i, test) - K_eff * x(:, i, test));

        % Update velocity and position using explicit Euler method
        x_dot(:, i+1, test) = x_dot(:, i, test) + SIMULATION.dt * x_ddot(:, i, test);
        x(:, i+1, test) = x(:, i, test) + SIMULATION.dt * x_dot(:, i, test);
    end

    % Find the period of oscillation for the current test
    [pks, locs] = findpeaks(x(find(x0, 1), :, test), SIMULATION.time);
    if length(locs) > 1
        periods(test) = mean(diff(locs));
    else
        % If no peaks are found, take the minimum (crest) and double the amount
        [~, minIdx] = min(x(find(x0, 1), :, test));
        periods(test) = 2 * SIMULATION.time(minIdx);
    end


    % Use only peaks for damping estimation, use only first 5 peaks or less
    % Add artificial peak at time 0 with amplitude 1
    eta_peaks = [1; pks(:)];
    t_peaks = [0; locs(:)];
    if length(eta_peaks) > 9
        eta_peaks = eta_peaks(1:9);
        t_peaks = t_peaks(1:9);
    end

    % Logarithmic decrement for linear damping estimation
    if length(eta_peaks) > 2
        deltas = log(eta_peaks(1:end-1) ./ eta_peaks(2:end));
        delta_mean = mean(deltas);
        damping_linear(test) = delta_mean / sqrt(4*pi^2 + delta_mean^2);
    else
        damping_linear(test) = NaN;
    end
end

% Plot results for all tests in a single figure
figure;
for test = 1:length(initial_conditions)
    subplot(2, 3, test);
    plot(SIMULATION.time, x(find(initial_conditions{test, 2}, 1), :, test),'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel(initial_conditions{test, 1});
    %title(initial_conditions{test, 1});
end

fprintf("Completed!\n");

% Display the periods and damping coefficients for each degree of freedom in one line
disp('Periods and damping coefficients for each degree of freedom:');
for test = 1:length(initial_conditions)
    fprintf('%s: Period = %.2f s, Linear Damping = %.4f, Quadratic Damping = %.4f\n', ...
        initial_conditions{test, 1}, periods(test), damping_linear(test), damping_quadratic(test));
end

%error("STOP HERE!"); % Stop execution here to check the results

%% TASK 10 - HYDRO + AERO LOADS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [5.1; 0; 0; 0.03; 0; 0]; % Initial states (3 translations and 3 rotations)
x_dot0 = zeros(6, 1); % Initial velocities

% Time integration
n = length(SIMULATION.time);
x = zeros(6, n); % States
x_dot = zeros(6, n); % Velocities
x_ddot = zeros(6, n); % Accelerations

Thrust = zeros(1, n);
Power = zeros(1, n);

% Assign initial conditions
x(:, 1) = x0;
x_dot(:, 1) = x_dot0;

% External forces (assumed zero for now)
Fext = zeros(6, 1);

fprintf("Running with aero/hydro loads...");

% Time-stepping loop
for i = 1:n-1
    % hydordynamic forces
    if hydro_from_FAST
            Fext_hydro = [OUTPUT.B1WvsFxi(i); OUTPUT.B1WvsFyi(i); OUTPUT.B1WvsFzi(i); ...
            OUTPUT.B1WvsMxi(i); OUTPUT.B1WvsMyi(i); OUTPUT.B1WvsMzi(i)];
    else
        Fext_hydro = FLOATER.Fhydrodynamic.* WAVES.WA .* cos(WAVES.k*x(1,i) - WAVES.omega .* SIMULATION.time(i)- FLOATER.Fhydrodynamicphase);
        Fext_hydro = Fext_hydro';
    end


    % aerodynamic forces
    [Power(i), Thrust(i), CP, CT, PREVIOUSTIME.a_new, PREVIOUSTIME.ap_new] = function_BEM(ROTOR, AIRFOIL, FLOW, SIMULATION, PREVIOUSTIME);
    
    % Calculate aerodynamic forces

    R_z = [cos(x(6,i)), -sin(x(6,i)), 0; ...
             sin(x(6,i)), cos(x(6,i)), 0; ...
             0, 0, 1]; % Rotation matrix for yaw

    
    R_y = [cos(x(5,i)), 0, sin(x(5,i));
            0, 1, 0;
            -sin(x(5,i)), 0, cos(x(5,i))]; % Rotation matrix for pitch
    
    R_x = [1 0 0;
            0 cos(x(4,i)) -sin(x(4,i));
            0 sin(x(4,i)) cos(x(4,i))]; % Rotation matrix for roll

    R = R_z * R_y * R_x; % Combined rotation matrix

    Forces_aero = R * [Thrust(i); 0; 0];

    arm = R*[0; 0; ROTOR.H]; % Position vector (surge, sway, 90)
    Moment_aero = cross(arm, Forces_aero); % Moment due to aerodynamic forces
    
    Fext_aero = [Forces_aero; Moment_aero];

    % combined external loads
    Fext = Fext_hydro + Fext_aero;

    % Solve for acceleration
    x_ddot(:, i) = M_eff \ (Fext - FLOATER.B * x_dot(:, i) - K_eff * x(:, i));
    
    % Update velocity and position using explicit Euler method
    x_dot(:, i+1) = x_dot(:, i) + SIMULATION.dt * x_ddot(:, i);
    x(:, i+1) = x(:, i) + SIMULATION.dt * x_dot(:, i);

    % Pass floater velocities
    FLOW.V_surge = x_dot(1, i+1);
    FLOW.V_pitch = ROTOR.H * x_dot(5, i+1); % Linear velocity at rotor radius due to pitch 
    FLOW.V_yaw   = ROTOR.H * sin(x(6, i+1)); % Linear velocity at rotor radius due to yaw
end

fprintf('Completed!\n');

figure;
dof_labels = {'Surge [m]', 'Sway [m]', 'Heave [m]', 'Roll [rad]', 'Pitch [rad]', 'Yaw [rad]'};
for j = 1:6
    subplot(3,2,j);
    plot(SIMULATION.time, x(j,:), 'LineWidth', 1.5);
    xlabel('Time [s]');
    ylabel(dof_labels{j});
    title(dof_labels{j});
    grid on;
end
sgtitle('Time Series Response of Each Degree of Freedom');

RESULTS.Time = SIMULATION.time;
RESULTS.Surge = x(1,:);
RESULTS.Sway = x(2,:);
RESULTS.Heave = x(3,:);
RESULTS.Roll = x(4,:);
RESULTS.Pitch = x(5,:);
RESULTS.Yaw = x(6,:);

figure;
subplot(3,2,1);
plot(OUTPUT.Time, OUTPUT.B1Surge, 'LineWidth', 1.5); hold on;
plot(RESULTS.Time, RESULTS.Surge, 'LineWidth', 1.5);
grid on;
legend('OpenFAST', 'Own code');
xlabel('Time [s]');
ylabel('Surge [m]');
title('Surge');

subplot(3,2,2);
plot(OUTPUT.Time, OUTPUT.B1Sway, 'LineWidth', 1.5); hold on;
plot(RESULTS.Time, RESULTS.Sway, 'LineWidth', 1.5);
grid on;
legend('OpenFAST', 'Own code');
xlabel('Time [s]');
ylabel('Sway [m]');
title('Sway');

subplot(3,2,3);
plot(OUTPUT.Time, OUTPUT.B1Heave, 'LineWidth', 1.5); hold on;
plot(RESULTS.Time, RESULTS.Heave, 'LineWidth', 1.5);
legend('OpenFAST', 'Own code');
grid on;
xlabel('Time [s]');
ylabel('Heave [m]');
title('Heave');

subplot(3,2,4);
plot(OUTPUT.Time, rad2deg(OUTPUT.B1Roll), 'LineWidth', 1.5); hold on;
plot(RESULTS.Time, rad2deg(RESULTS.Roll), 'LineWidth', 1.5);
legend('OpenFAST', 'Own code');
grid on;
xlabel('Time [s]');
ylabel('Roll [deg]');
title('Roll');

subplot(3,2,5);
plot(OUTPUT.Time, rad2deg(OUTPUT.B1Pitch), 'LineWidth', 1.5); hold on;
plot(RESULTS.Time, rad2deg(RESULTS.Pitch), 'LineWidth', 1.5);
legend('OpenFAST', 'Own code');
grid on;
xlabel('Time [s]');
ylabel('Pitch [deg]');
title('Pitch');

subplot(3,2,6);
plot(OUTPUT.Time, rad2deg(OUTPUT.B1Yaw), 'LineWidth', 1.5); hold on;
plot(RESULTS.Time, rad2deg(RESULTS.Yaw), 'LineWidth', 1.5);
legend('OpenFAST', 'Own code');
grid on;
xlabel('Time [s]');
ylabel('Yaw [deg]');
title('Yaw');

figure;
subplot(2,1,1);
plot(OUTPUT.Time, OUTPUT.RtAeroFxh, 'LineWidth', 1.5); hold on;
plot(SIMULATION.time, Thrust, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Thrust [N]');
title('Thrust');
legend('OpenFAST', 'Own code');

subplot(2,1,2);
plot(OUTPUT.Time, OUTPUT.RtAeroMxh, 'LineWidth', 1.5); hold on;
plot(SIMULATION.time, Power, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Power [W]');
title('Power');
legend('OpenFAST', 'Own code');
