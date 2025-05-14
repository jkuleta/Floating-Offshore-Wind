function [P,T,CP,CT,a_new,ap_new] = function_BEM(ROTOR,AIRFOIL,FLOW,SIMULATION,PREVIOUSTIME)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function_BEM is a function solving the blade element momentum (BEM)
%%% theory. It can be used to determine the loads and performance of a wind
%%% turbine.
%%%
%%% Input parameters:
%%% - ROTOR: contains the rotor geometry   
%%%      radius sections, r(N) [m]
%%%      twist, beta(N) [deg]
%%%      chord, chord(N) [m]
%%%      airfoil shape, airfoil [string]
%%%      radius, R [m]
%%%      diameter, D [m]
%%%      tower height, H [m]
%%%      number of blades, B [-]
%%%      pitch angle, theta_pitch(N) [deg]
%%%      rotor solidity, sigmna(N), [-]
%%% - AIRFOIL: contains polars of each airfoil 
%%%      airfoil shape, [Alpha, Cl, Cd, Cm] [deg, -]
%%% - FLOW: contains all flow and operational properties
%%%      inflow velocity, V0 [m/s]
%%%      air density, rho [kg/m3]
%%%      rotor rotational speed [rad/s]
%%% - SIMULATION: contains all simulation properties
%%%      convergence criteria, error [-]
%%%      constants for dynamic inflow model (only used for unsteady), taustar_nw [-]
%%%      constants for dynamic inflow model (only used for unsteady), taustar_fw [-]
%%% - PREVIOUSTIME: contains induction parameters of the previous time step 
%%%   in case of unsteady simulation in time. For steady simulations or during the 
%%%   first timestep of a time-resolved solution, use PREVIOUSTIME = [].
%%%      axial induction of previous timestep, PREVIOUSTIME.a_new(N), [-]
%%%      tangential induction of previous timestep, PREVIOUSTIME.ap_new(N), [-]
%%%
%%% Output parameters:
%%%      power, P [W]
%%%      thrust, T [N]
%%%      power coefficient, CP [-]
%%%      thrust coefficeint, CT [-]
%%%      axial induction, a(N) [-]
%%%      tangential induction, ap(N) [-]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Pre-allocate variables
a_new = 0.3*ones(length(ROTOR.r),1);
ap_new = 0*ones(length(ROTOR.r),1);
a_old = 0.0*ones(length(ROTOR.r),1);
ap_old = 0*ones(length(ROTOR.r),1);

% Default values if not present
if ~isfield(FLOW, 'V_surge'); FLOW.V_surge = 0; end
if ~isfield(FLOW, 'V_pitch'); FLOW.V_pitch = 0; end
if ~isfield(FLOW, 'V_yaw');   FLOW.V_yaw   = 0; end

iteration = 0;

while sum(abs(a_new-a_old))>SIMULATION.error || ...                        % Condition for convergence
      sum(abs(ap_new-ap_old))>SIMULATION.error
    
    iteration = iteration+1;                                               % Number of iterations [-]
    
    % Update induction
    a_old = a_new;                                                         % Axial induction factor [-]
    ap_old = ap_new;                                                       % Tangential induction factor [-]
    
    % Velocity component
    Vn = (1 - a_new) * FLOW.V0 - FLOW.V_surge - FLOW.V_pitch - FLOW.V_yaw; % Normal velocity (to rotor plane) [m/s]
    %Vn = (1-a_new)*FLOW.V0;                                               % Normal velocity (to rotor plane) [m/s]
    Vt = (1+ap_new)*FLOW.omega.*ROTOR.r;                                   % Tangential velocity (to rotor plane) [m/s]
    Vrel = sqrt(Vn.^2+Vt.^2);                                              % Local relative velocity [m/s]
    
    % Inflow angles
    Phi = atand(Vn./Vt);                                                   % Local inflow angle [deg]
    Theta = ROTOR.theta_pitch+ROTOR.beta;                                  % Local pitch angle [deg]
    Alfa = Phi-Theta;                                                      % Local angle of attack [deg]
    
    % Lift and drag coefficient
    Cl = zeros(length(ROTOR.r),1);                                         % Lift coefficient [-]
    Cd = zeros(length(ROTOR.r),1);                                         % Drag coefficient [-]
    for i = 1:length(ROTOR.r)
        Cl(i) = interp1(AIRFOIL.(ROTOR.airfoil{i})(:,1), AIRFOIL.(ROTOR.airfoil{i})(:,2), Alfa(i)); 
        Cd(i) = interp1(AIRFOIL.(ROTOR.airfoil{i})(:,1), AIRFOIL.(ROTOR.airfoil{i})(:,3), Alfa(i)); 
    end
    
    % Normal and tangential coefficient (to rotor plane)
    cn = Cl.*cosd(Phi)+Cd.*sind(Phi);                                      % Normal force coefficient [-]
    ct = Cl.*sind(Phi)-Cd.*cosd(Phi);                                      % Tangential force coefficient [-]
    
    % Thrust and torque coefficient
    CT = Vrel.^2./FLOW.V0.^2.*ROTOR.sigma.*cn;                             % Thrust coefficient (normal to rotor plane) [-]
    CQ = Vrel.^2./FLOW.V0.^2.*ROTOR.sigma.*ct;                             % Torque coefficient (in rotor plane) [-]
    
    % Axial and tangential induction
    %%% Induction including prantl tip correction    
    f = ROTOR.B*(ROTOR.R-ROTOR.r)./(2*ROTOR.r.*sind(Phi));                
    F = 2/pi*acos(exp(-f));
    Ctb = CT./F;
    %%% Madsen's polynomial relation for heavily loaded rotors
    k = [0.00 ,0.251163 ,0.0544955 ,0.0892074];
    a_n = (k(4)*Ctb.^3+k(3)*Ctb.^2+k(2)*Ctb+k(1));                                   
    ap_n = CQ./(4*F.*(1-a_n).*FLOW.omega.*ROTOR.r/FLOW.V0); 
    %%% Relaxation
    a_new = (0.2*a_n+0.8*a_old);                                           % Steady axial induction factor [-]
    ap_new = 0.2*ap_n+0.8*ap_old;                                          % Steady tangential induction factor [-]  
    
    %  Dynamic inflow model
    if isstruct(PREVIOUSTIME)
        Vwake = (1-2*a_new)*FLOW.V0;
        a_nw = PREVIOUSTIME.a_new.*exp(-SIMULATION.dt./(SIMULATION.taustar_nw*ROTOR.R./Vwake)) + a_new.*(1-exp(-SIMULATION.dt./(SIMULATION.taustar_nw*ROTOR.R./Vwake)));
        a_fw = PREVIOUSTIME.a_new.*exp(-SIMULATION.dt./(SIMULATION.taustar_fw*ROTOR.R./Vwake)) + a_new.*(1-exp(-SIMULATION.dt./(SIMULATION.taustar_fw*ROTOR.R./Vwake)));
        ap_nw = PREVIOUSTIME.ap_new.*exp(-SIMULATION.dt./(SIMULATION.taustar_nw*ROTOR.R./Vwake)) + ap_new.*(1-exp(-SIMULATION.dt./(SIMULATION.taustar_nw*ROTOR.R./Vwake)));
        ap_fw = PREVIOUSTIME.ap_new.*exp(-SIMULATION.dt./(SIMULATION.taustar_fw*ROTOR.R./Vwake)) + ap_new.*(1-exp(-SIMULATION.dt./(SIMULATION.taustar_fw*ROTOR.R./Vwake)));
        a_new = 0.6*a_nw+0.4*a_fw;                                         % Unsteady axial induction factor after dynamic inflow correction [-]
        ap_new = 0.6*ap_nw+0.4*ap_fw;                                      % Unsteady angential induction factor after dynamic inflow correction [-]  
    end    
    
end
    
% Blade forces
pn = 1/2*FLOW.rho*Vrel.^2.*ROTOR.chord.*cn;                                % Local normal force [N]
pt = 1/2*FLOW.rho*Vrel.^2.*ROTOR.chord.*ct;                                % Local tangential force [N]
 
% Rotor performance
M = ROTOR.B* trapz(ROTOR.r,ROTOR.r.* pt(:));                               % Moment [Nm]
P = M*FLOW.omega;                                                          % Power [W]
CP = P/(1/2*FLOW.rho*FLOW.V0^3*pi*ROTOR.R^2);                              % Power coefficient [-]
T = ROTOR.B*trapz(ROTOR.r,pn(:));                                          % Thrust [N]
CT = T/(1/2*FLOW.rho*FLOW.V0^2*pi*ROTOR.R^2);                              % Thrust coefficient [-]   


