% PID_Plot 
% single PID file in single folder
clear all; close all; 

% Variables user can change: 
x_limits = [0000,16000];
%% Import data
PID_file = uigetfile('*.csv');
% importer interface --> importable range A2:B1001
% variables user can alter:
rowOffset = 1; 
colOffset = 0;
%
M = csvread(PID_file,rowOffset,colOffset);
%% Plot data
fig = figure;
Mx = M(:,1); 
My = M(:,2);
plot(Mx,My)
xlabel('Time')
ylabel('PID Reading')
set(gca, 'xlim', x_limits)
%% Save figure
figName = ['PID_',PID_file(1:end-4)];
saveas(fig, figName, 'fig'); 
saveas(fig, figName, 'png'); 