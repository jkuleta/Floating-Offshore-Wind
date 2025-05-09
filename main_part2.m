clear all; close all; clc;

%% Input Parameters
% Flow
FLOW.V0 = 10;                                                              % Inflow velocity [m/s]   
FLOW.rho = 1.225;                                                          % Air density [kg/m2]
FLOW.omega = 9*2*pi/60;                                                    % Rotational speed [rad/s]

% Waves
WAVES.H = 2.44;                                                            % Wave height [m]
WAVES.T = 8.1;                                                             % Wave period [s]
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
SIMULATION.time = 0:SIMULATION.dt:100;                                     % Time series [s]
SIMULATION.taustar_nw = 0.5;                                               % Constants for dynamic inflow model 
SIMULATION.taustar_fw = 2;                                                 % Constants for dynamic inflow model 

%% State-space solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The equation of motion of this floater-turbine system has the shape:
%%% (M1+A1)x'' + Bx' + (C1+K1)x = Fext
%%%      M1 is the mass matrix
%%%      A1 is the added-mass matrix, and depends on the wave frequency
%%%      B is the damping matrix, and depends on the wave frequency
%%%      C1 is the hydrostatic matrix
%%%      K1 is the mooring line matrix
%%%      Fext is the external forces: aerodynamic and the hydrodynamic forces
%%%
%%%
%%% To be completed ...
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
