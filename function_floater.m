function [M1,A1,B,C1,K11,Fhydrodynamic,Fhydrodynamicphase] = function_floater(WAVES)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function_floater is a function used to extract the mass, damping and
%%% stiffness matrices for the OC5 floater. It also extracts the amplitude
%%% and phase of the hydrodynamic forces
%%%
%%% Input parameters
%%% - WAVES: contains the information of the waves   
%%%      wave period, T [s]
%%%
%%% Output parameters
%%%      mass matrix, M1 [-]
%%%      added-mass matrix, A1 [-]
%%%      damping matrix, B [-]
%%%      hydrostatic matrix, C1 [-]
%%%      mooring line matrix, K11[6x6] [-]
%%%      amplitude of hydrodynamic force[6x1] [N,Nm]
%%%      phase of hydrodynamic force[6x1] [deg]
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Read Wadam files
load('WADAM.mat')

ADDMASS = WADAM.ADDMASS;
DAMPING = WADAM.DAMPING;
MASS = WADAM.MASS;
WAVEDATA1 = WADAM.WAVEDATA1;
WAVEDATA2 = WADAM.WAVEDATA2;
FEXT = WADAM.FEXT;

%% Wave period for added mass and damping

% Find two closest wave periods
periodvector1 = WAVEDATA1(:,4)-WAVES.T;
periodvector_pos1 = periodvector1(periodvector1>=0);
period_above1 = find(periodvector1 == min(periodvector_pos1));
period_below1 = period_above1+1;

if period_below1>length(periodvector1)
    period_below1 = period_above1;
end

% Take multiplyfactor of period_above and 1-multiplyfactor of period_below
% for interpolation
multiplyfactor1 = (WAVES.T-WAVEDATA1(period_below1,4))/(WAVEDATA1(period_above1,4)-WAVEDATA1(period_below1,4)); 

if period_below1 == period_above1
    multiplyfactor1 = 1;
end

%% Wave period for force

% Find two closest wave periods
periodvector2 = WAVEDATA2(:,4)-WAVES.T;
periodvector_pos2 = periodvector2(periodvector2>=0);
period_above2 = find(periodvector2 == min(periodvector_pos2));
period_below2 = period_above2+1;

if period_below2>length(periodvector2)
    period_below2 = period_above2;
end

% Take multiplyfactor of period_above and 1-multiplyfactor of period_below
multiplyfactor2 = (WAVES.T-WAVEDATA2(period_below2,4))/(WAVEDATA2(period_above2,4)-WAVEDATA2(period_below2,4)); 

if period_below2 == period_above2
    multiplyfactor2 = 1;
end

%% Extract mass-matrix (M)
% Normal mass matrix (M1)
M1 = MASS;
 
% Added-mass matrix (A) ~ Depends on wave frequency
A1 = multiplyfactor1*ADDMASS(period_above1,:,:)+(1-multiplyfactor1)*ADDMASS(period_below1,:,:);
A1 = squeeze(A1);

%% Extract damping matrix (D)
B = multiplyfactor1*DAMPING(period_above1,:,:)+(1-multiplyfactor1)*DAMPING(period_below1,:,:);
B = squeeze(B);

B(1,:)=0.25*B(1,:);
B(5,:)=0.25*B(5,:);

%% Extract Stiffness matrix (K)
% Hydrostatic matrix (C)
C1 = zeros(6,6);
C1(3,3) = 3.836E6;
C1(3,5) = 7.73028E+01;
C1(4,4) = 7.72139E+08;
C1(5,5) = 8.05557E+08;
C1(4,6) = -1.52128E+07;
C1(5,6) = 3.83446E-02;

% % Mooring from OC4 > gives good decay period
K11 = zeros(6,6);
K11(1,1) = 7.08E4*1.10;
K11(1,5) = -1.08E5*1.25;
K11(2,2) = 7.08E4;
K11(2,4) = 1.08E5;
K11(3,3) = 1.91E4;
K11(4,2) = 1.07E5;
K11(4,4) = 8.73E7;
K11(5,1) = -1.07E5*3;
K11(5,5) = 8.73E7*3;
K11(6,6) = 1.17E8;


%% Forces F
Fhydrodynamic = multiplyfactor2*FEXT(period_above2,:,3)+(1-multiplyfactor2)*FEXT(period_below2,:,3); %per wave amplitude
Fhydrodynamicphase = multiplyfactor2*FEXT(period_above2,:,4)+(1-multiplyfactor2)*FEXT(period_below2,:,4); %per wave amplitude



