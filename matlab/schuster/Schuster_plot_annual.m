function Schuster_plot_annual(t_serie)

% Schuster_plot_annual.m plots the Schuster walk associated with to a time
% series and a period of 1 year.
% -------------------------------------------------------------------------
% INPUT:
%   t_serie: a time series in datetime format
% -------------------------------------------------------------------------
% OUTPUT:
%   Plot of the Schuster walk associated to annual period.
% -------------------------------------------------------------------------
% AUTHOR:   Krittanon "Pond" Sirorattanakul
%           Seismological Laboratory
%           California Institute of Technology
%          
%           Adpated from Schuster_plot.m by Thomas Ader for the case of 
%           annual period.
% -------------------------------------------------------------------------
% REVISION: 1.0.0   30 MAR 2021     Initial creation
% -------------------------------------------------------------------------


if size(t_serie,2)>1 % To make sure we have column vectors
    t_serie = t_serie';
end

% Convert datetime to year and phase it based on Jan 1st of the first year
year_min = floor(decyear(min(t_serie)));
t = years(t_serie - datetime(year_min,1,1));

% Phase all the times from the timeseries with respect to the period T.
T = 1;
phase=mod(t,T)*2*pi/T; 

% Perform Schuster walk
end_walk = [cumsum(cos(phase)),cumsum(sin(phase))];

% Make figure
% figure; hold on;

% Plot the walk itself
plot(end_walk(:,1),end_walk(:,2))

% Plot big dots at the beginning of each year
% Include last year in the interpolation only if we have data at least 85%
% of the last year.
if t(end)-floor(t(end)) > 0.9
    x_year = interp1(t,end_walk(:,1),[0:1:floor(t(end))+1],'linear','extrap');
    y_year = interp1(t,end_walk(:,2),[0:1:floor(t(end))+1],'linear','extrap');
else
    x_year = interp1(t,end_walk(:,1),[0:1:floor(t(end))],'linear','extrap');
    y_year = interp1(t,end_walk(:,2),[0:1:floor(t(end))],'linear','extrap');
end
plot(x_year,y_year,'.','Marker','o',...
    'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5)

% Plot circles at probabilities 10^(-2) and 10^(-5)
p=[0.01,0.00001]; % probability to get there by random walk
% <=> 1-p: confidence that being that far isn't due to random walk
theta = 0:0.01:(2*pi);

for j=1:length(p)
    D = sqrt(length(t)*log(1/p(j)));
    
    x_circ = D*cos(theta);
    y_circ = D*sin(theta);
    
    plot(x_circ,y_circ,'--k')
    conf = [num2str(p(j)*100) '%'];
    text(x_circ(60)+5,y_circ(60)+5,conf,'HorizontalAlignment','left','FontSize',10)
end

% Plot labels for months
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
theta = 0:(2*pi)/12:(11/12)*2*pi;
D = sqrt(length(t)*log(1/p(2)));
for ii=1:length(months)
    plot([0 D*cos(ii*2*pi/12)],[0 D*sin(ii*2*pi/12)],'k:','LineWidth',0.5)
    text(D*cos((ii-0.5)*2*pi/12)*1.1-0.05*D,D*sin((ii-0.5)*2*pi/12)*1.1,months(ii),'FontSize',10)
end

% Adjust the axes
axis equal
xlim(1.1*get(gca,'XLim'))
ylim(1.1*get(gca,'YLim'))

% Add text labels to beginning of each year
dx = 0.01*diff(get(gca,'XLim'));
dy = 0.01*diff(get(gca,'YLim'));
for ii=1:length(x_year)
    text(x_year(ii)+dx,y_year(ii)-dy,num2str(year_min+ii-1),'FontSize',5)
end
   

end