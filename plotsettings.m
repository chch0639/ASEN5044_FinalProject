%% Assign defaults
function plotsettings(fontsize,linewidth)
set(0, 'defaultAxesXGrid', 'on')
set(0, 'defaultAxesYGrid', 'on')
set(0, 'defaultAxesZGrid', 'on')
set(0, 'RecursionLimit', 50);
set(0, 'DefaultFigurePaperType', 'A4');
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultTextFontSize', fontsize);
set(0, 'defaultlinelinewidth', linewidth);
set(0, 'DefaultAxesFontSize', fontsize);
set(0,'defaultfigurecolor',[1 1 1])
%set(0, 'DefaultUicontrolFontSize', 10);
set(groot,'defaultFigurePaperPositionMode','auto')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
