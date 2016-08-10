% PID_Plot_func
% can call from another script to process a PID file

function [fig, figName] = PId_Plot_func(PID_file);
% Data is taken from PID_file.
% importer interface --> importable range A2:B1001
% variables user can alter:
    rowOffset = 1;
    colOffset = 0;
%
    x_limits = [2000,7000];
    y_limits = [0,1]; % default y limits for plots
%
    M = csvread(PID_file,rowOffset,colOffset); % reads data
% Plot data
    fig = figure;
    Mx = M(:,1);
    My = M(:,2);
    plot(Mx,My)
    xlabel('Time')
    ylabel('PID Reading')
    set(gca, 'xlim', x_limits)
    set(gca, 'ylim', y_limits)
    % Save figure
    figName = ['PID_',PID_file(1:end-4)];
    saveas(fig, figName, 'fig');
    saveas(fig, figName, 'png');
end