function Schuster_plot(t_serie,periods)


% SCHUSTER_PLOT plots the Schuster walk associated to a timeseries for a
% given set of periods.
%
% INPUT
%
% T_SERIE: The timeseries on which we want to test any periodicity
% PERIODS: Vector containing all the periods to be tested
%
% OUTPUT
%
% Plot of the Schuster walk.
%
%
%
% Thomas Ader
% June 2010
%

if size(t_serie,2)>1 % To make sure we have column vectors
    t_serie = t_serie';
end

n_per = length(periods); % Number of periods tested
t_serie = t_serie - t_serie(1); % So the timeseries starts at time 0.

for i = 1:n_per

    T = periods(i);
    tlim=t_serie(end)-mod(t_serie(end),T);%To have a round number of cycles.
    % If not it induces artefacts to Schuster.
    t=t_serie(t_serie<=tlim); % selects a round number of cycles from the
    % timeseries.

    phase=mod(t,T)*2*pi/T; % Phase all the times from the timeseries
    % with respect to the period T.

    end_walk = [cumsum(cos(phase)),cumsum(sin(phase))]; % Point where we end the
    % Schuster walk

    figure; hold on
    %     Plot the walk itself
    plot(end_walk(:,1),end_walk(:,2))

    %     Plot a big dot at the starting and end points
    plot([0 end_walk(end,1)],[0 end_walk(end,2)],'.','Marker','o',...
        'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5)

    %     Plot circles at probabilities 10^(-2) and 10^(-5)
    p=[0.01,0.00001]; % probability to get there by random walk
    % <=> 1-p: confidence that being that far isn't due to random walk
    theta = 0:0.01:(2*pi);

    for j=1:length(p)
        D = sqrt(length(t)*log(1/p(j)));

        x_circ = D*cos(theta);
        y_circ = D*sin(theta);

        plot(x_circ,y_circ,'--k')
        conf = [num2str(p(j)*100) '%'];
        text(x_circ(120)+5,y_circ(120)+5,conf,'HorizontalAlignment','left')
    end

    axis equal

end

end