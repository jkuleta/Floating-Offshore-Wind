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

%% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%% To be completed ...
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%