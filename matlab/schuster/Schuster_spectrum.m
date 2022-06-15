function [T_test,log_prob] = Schuster_spectrum(t_serie,T_min,T_max)


% SCHUSTER_SPECTRUM computes the Schuster spectrum of a given timeseries
% between 2 given periods.
%
% INPUT
%
% T_SERIE: The timeseries on which we want to test any periodicity
% T_MIN: The minimum period to be tested (same time unit as t_serie)
% T_MAX: The maximum period to be tested (same time unit as t_serie)
%
% OUTPUT
%
% T_TEST: Vector containing all the periods that have been tested  (same
% time unit as t_serie)
% LOG_PROB: Vector of the same length as 'T_test' giving for each period the
% log of the Schuster p-value. If D is the distance covered by the Schuster
% walk, and N the number of events in the catalog, log_prob = -D^2/N.
%
% COMMENTS
%
% The code uses epsilon = 1. This value has to be modified manually
% inside the code (first non-commented line) for other values.
%
% FUNCTIONS CALLED
%
% Schuster_test_log
%
%
% By Thomas Ader
% Copyright 2011-2012 Tectonics Observatory
% Created 05/01/2012
% Modified 09/01/2012 - 06/19/2013
% Add points at all tidal cycles -- Pond Sirorattanakul (29 MAR 2021)


% Generate the vector of periods to test
eps = 1;

t_serie = t_serie - min(t_serie); % So the time series starts at time 0.
t_span = max(t_serie); % Time spanned by the catalog
% Periods to test
T_test = 1./(1/T_max:eps/t_span:1/T_min);
% Also test periods at various tidal cycles
T_test2 = [0.4986/365.25 0.5/365.25 0.5175/365.25 0.5274/365.25 ...
            0.9973/365.25 1.0028/365.25 1.0758/365.25 13.661/365.25 ...
            14.765/365.25 27.55/365.25 31.812/365.25 182.621/365.25 ...
            365.26/365.25];

% Compute the Schuster p-values at the periods to test
log_prob=Schuster_test_log(t_serie,T_test);
log_prob2 = Schuster_test_log(t_serie,T_test2);

% Plot the Schuster spectrum
% figure;
% Spectrum
p_val = exp(log_prob);
p_val2 = exp(log_prob2);
eps_th = 1;
loglog(T_test,p_val,'o:', 'color', 0.8*[1 1 1], ...
    'markersize', 3, 'markerfacecolor', 0.8*[1 1 1], ...
    'MarkerEdgeColor', 0.7*[1 1 1]);
hold on
loglog(T_test2,p_val2,'o', 'color', 0.8*[0 0 1], ...
    'markersize', 3, 'markerfacecolor', 0.8*[0 0 1], ...
    'MarkerEdgeColor', 0.7*[0 0 1]);

% Expected value
loglog(T_test, eps_th*T_test/t_span, '--k', 'LineWidth', 1.5)
% Points above the expected value: change size and color
sig_pts = p_val < eps_th*T_test/t_span;
T_test_sig = T_test(sig_pts);
p_val_sig = p_val(sig_pts);
c_diff = log10(p_val_sig) - log10(eps_th*T_test_sig/t_span);
c_diff(c_diff<-2) = -2;
c_diff = 2*c_diff/10;

if ~isempty(T_test_sig)
    for i = 1:length(T_test_sig)
        loglog(T_test_sig(i),p_val_sig(i),'o', 'color', ...
            (0.8+c_diff(i))*[1 1 1], 'markersize', 5, ...
            'markerfacecolor', (0.8+c_diff(i))*[1 1 1], ...
            'MarkerEdgeColor', (0.7+1.5*c_diff(i))*[1 1 1]);
    end
end

% Points above the expected value: change size and color (for tidal
% periods)
sig_pts = p_val2 < eps_th*T_test2/t_span;
T_test_sig = T_test2(sig_pts);
p_val_sig = p_val2(sig_pts);
c_diff = log10(p_val_sig) - log10(eps_th*T_test_sig/t_span);
c_diff(c_diff<-2) = -2;
c_diff = 2*c_diff/10;

if ~isempty(T_test_sig)
    for i = 1:length(T_test_sig)
        loglog(T_test_sig(i),p_val_sig(i),'o', 'color', ...
            (0.8+c_diff(i))*[1 0 0], 'markersize', 5, ...
            'markerfacecolor', (0.8+c_diff(i))*[0 0 1], ...
            'MarkerEdgeColor', (0.7+1.5*c_diff(i))*[0 0 1]);
    end
end


% % 95% confidence level
% loglog(T_test, 0.05*eps_th*T_test/t_span, '--k', 'LineWidth', 1.5)

% 99% confidence level
loglog(T_test, 0.01*eps_th*T_test/t_span, '--k', 'LineWidth', 1.5)

% Points above the 99% confidence level
sig_pts_99 = p_val < 0.01*eps_th*T_test/t_span;

T_test_99 = T_test(sig_pts_99);
p_val_99 = p_val(sig_pts_99);
T_test_99(p_val_99 == 0) = [];
p_val_99(p_val_99 == 0) = [];
m_size = 5 - log10(p_val_99) + log10(0.01*eps_th*T_test_99/t_span);

if ~isempty(T_test_99)
    for i = 1:length(T_test_99)
        loglog(T_test_99(i),p_val_99(i),'o', 'color', 0.4*[1 1 1], ...
            'markersize', min([m_size(i), 10]), 'markerfacecolor', 0.4*[1 1 1], ...
            'MarkerEdgeColor', 0.1*[1 1 1]);
    end
end

% Points above the 99% confidence level (for tidal periods)
sig_pts_99 = p_val2 < 0.01*eps_th*T_test2/t_span;

T_test_99 = T_test2(sig_pts_99);
p_val_99 = p_val2(sig_pts_99);
T_test_99(p_val_99 == 0) = [];
p_val_99(p_val_99 == 0) = [];
m_size = 5 - log10(p_val_99) + log10(0.01*eps_th*T_test_99/t_span);

if ~isempty(T_test_99)
    for i = 1:length(T_test_99)
        loglog(T_test_99(i),p_val_99(i),'o', 'color', 0.4*[0 0 1], ...
            'markersize', min([m_size(i), 10]), 'markerfacecolor', 0.4*[0 0 1], ...
            'MarkerEdgeColor', 0.1*[0 0 1]);
    end
end

xlim([T_min T_max])
grid on; grid minor;
xlabel('Period tested','FontSize',14,...
    'FontWeight','bold','FontName','Times')
ylabel('Schuster p-value','FontSize',14,'FontWeight','bold',...
    'FontName','Times')
set(gca,'YDir','reverse')
set(gcf,'pos', [0 450 800 350])

% Write expected value and 99% confidence level
per_txt = 10^(0.75*log10(T_min) + 0.25*log10(T_max));
txt_pos = [per_txt eps_th*per_txt/t_span];
axes_size = get(gca,'pos').*get(gcf,'pos');
XL = axes_size(3);
YL = axes_size(4);
y_lim = get(gca,'Ylim');
angle_rot = -atand(YL/XL*...
    log10(T_max/T_min)/log10(y_lim(2)/y_lim(1)));

% Plot tides/year

% Tides periods in days
K2=[0.4986; 0.4986];
S2=[0.5; 0.5];
M2=[0.5175; 0.5175];
N2=[0.5274; 0.5274];
K1=[0.9973;0.9973];
P1=[1.0028;1.0028];
O1=[1.0758;1.0758];
T=[13.661;13.661];
S=[14.765;14.765];
A=[27.555;27.555];
E=[31.812;31.812];
SY=[182.621;182.621];
Y=[365.260;365.260];

% Tides periods in years
K2=K2/365.25;
S2=S2/365.25;
M2=M2/365.25;
N2=N2/365.25;
K1=K1/365.25;
P1=P1/365.25;
O1=O1/365.25;
T=T/365.25;
S=S/365.25;
A=A/365.25;
E=E/365.25;
SY=SY/365.25;
Y=Y/365.25;


Y1 = y_lim;
loglog(K2,Y1,'--b')
loglog(S2,Y1,'--b')
loglog(M2,Y1,'--b')
loglog(N2,Y1,'--b')
loglog(K1,Y1,'--b')
loglog(P1,Y1,'--b')
loglog(O1,Y1,'--b')
loglog(T,Y1,'--b')
loglog(S,Y1,'--b')
loglog(A,Y1,'--b')
loglog(E,Y1,'--b')
loglog(SY,Y1,'--b')
loglog(Y,Y1,'--b')

text(txt_pos(1),txt_pos(2), '\textbf{expected value}',...
    'FontSize',14,'FontWeight','bold','interpreter','latex',...
    'Color', 'k', 'VerticalAlignment', 'bottom', ...
    'Rotation', angle_rot);

% text(txt_pos(1),txt_pos(2)*0.05, '95% confidence level','FontSize',12,...
%     'FontWeight','bold','FontName','Times', 'Color', 'k',...
%     'VerticalAlignment', 'bottom','Rotation', angle_rot);

text(txt_pos(1),txt_pos(2)*0.01, '\textbf{99\% confidence level}',...
    'FontSize',14,'FontWeight','bold','interpreter','latex',...
    'Color', 'k', 'VerticalAlignment', 'bottom', ...
    'Rotation', angle_rot);

% % If we want to plot the 99.99% confidence level
% loglog(T_test, 0.0001*eps_th*T_test/t_span, '--k', 'LineWidth', 1.5)
% per_txt_2 = 10^(0.3*log10(T_min) + 0.7*log10(T_max));
% txt_pos_2 = [per_txt_2 eps_th*per_txt_2/t_span*1e-4];
% text(txt_pos_2(1),txt_pos_2(2), '99.99% confidence level','FontSize',12,...
%     'FontWeight','bold','FontName','Times', 'Color', 'k',...
%     'VerticalAlignment', 'bottom','Rotation', angle_rot);
% ylim(y_lim)

%%%% NEED FIXING
T_test = T_test2;
log_prob = log_prob2;

end








