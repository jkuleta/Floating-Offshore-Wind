clc; clear all; close all

%% Input paramaters
casename = '5MW_OC4Semi_WSt_WavesWN';                                      % Name of the folder

%% Run openFAST
cd([casename])

% Run openfast to calculate loads and performance
system(['..\openfast\openfast_x64 ' '5MW_OC4Semi_WSt_WavesWN.fst'])
cd(['..'])

%% Read output file
% load file
file = fileread([casename '/5MW_OC4Semi_WSt_WavesWN.out']);            
lines = regexp(file, '\n', 'split');

Variable_list = strsplit(char(lines(7)),'\t');
Unit_list = strsplit(char(lines(8)),'\t');
Data_list = str2num(char(lines(10:4:length(lines))));
  
% Sort data in a table
OUTPUT = array2table(Data_list);
OUTPUT.Properties.VariableNames = Variable_list;

save("OUTPUT_S3.mat", "OUTPUT");

%% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%
% 
% load('RESULTS_S4.mat');
% figure;
% subplot(3,2,1);
% plot(OUTPUT.Time, OUTPUT.B1Surge, 'LineWidth', 1.5); hold on;
% plot(RESULTS.Time, RESULTS.Surge, 'LineWidth', 1.5);
% grid on;
% legend('OpenFAST', 'Own code');
% xlabel('Time [s]');
% ylabel('Surge [m]');
% title('Surge');
% 
% subplot(3,2,2);
% plot(OUTPUT.Time, OUTPUT.B1Sway, 'LineWidth', 1.5); hold on;
% plot(RESULTS.Time, RESULTS.Sway, 'LineWidth', 1.5);
% grid on;
% legend('OpenFAST', 'Own code');
% xlabel('Time [s]');
% ylabel('Sway [m]');
% title('Sway');
% 
% subplot(3,2,3);
% plot(OUTPUT.Time, OUTPUT.B1Heave, 'LineWidth', 1.5); hold on;
% plot(RESULTS.Time, RESULTS.Heave, 'LineWidth', 1.5);
% legend('OpenFAST', 'Own code');
% grid on;
% xlabel('Time [s]');
% ylabel('Heave [m]');
% title('Heave');
% 
% subplot(3,2,4);
% plot(OUTPUT.Time, OUTPUT.B1Roll, 'LineWidth', 1.5); hold on;
% plot(RESULTS.Time, RESULTS.Roll, 'LineWidth', 1.5);
% legend('OpenFAST', 'Own code');
% grid on;
% xlabel('Time [s]');
% ylabel('Roll [rad]');
% title('Roll');
% 
% subplot(3,2,5);
% plot(OUTPUT.Time, OUTPUT.B1Pitch, 'LineWidth', 1.5); hold on;
% plot(RESULTS.Time, RESULTS.Pitch, 'LineWidth', 1.5);
% legend('OpenFAST', 'Own code');
% grid on;
% xlabel('Time [s]');
% ylabel('Pitch [rad]');
% title('Pitch');
% 
% subplot(3,2,6);
% plot(OUTPUT.Time, OUTPUT.B1Yaw, 'LineWidth', 1.5); hold on;
% plot(RESULTS.Time, RESULTS.Yaw, 'LineWidth', 1.5);
% legend('OpenFAST', 'Own code');
% grid on;
% xlabel('Time [s]');
% ylabel('Yaw [rad]');
% title('Yaw');


%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%