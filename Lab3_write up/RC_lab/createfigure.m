function createfigure(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 16-Sep-2015 16:50:33

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1);
set(plot1(1),'DisplayName','recorded voltages','Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
set(plot1(2),'DisplayName','theoretical values','LineWidth',2,...
    'Color',[1 0 0]);

% Create xlabel
xlabel('time (ms)','FontSize',11);

% Create ylabel
ylabel('v_c (V)','FontSize',11);

% Create title
title('voltage vs time','FontSize',11);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location','southeast');

