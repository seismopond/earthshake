function log_prob=Schuster_test_log(t_serie,periods)


% SCHUSTER_TEST_LOG performs the Schuster test on a given timeseries
%
% INPUT
%
% T_SERIE: The timeseries on which we want to test any periodicity
% PERIODS: Vector containing all the periods to be tested
%
% OUTPUT
%
% LOG_PROB: Vector of the same length as 'period' giving for each period the
% log of probability that the distance covered by the Schuster walk is due to a
% random walk. The smaller the probability, the more likely it is that
% there is a periodicity in the data at the given period.
%
%
% By Thomas Ader
% Copyright 2010-2011 Tectonics Observatory
% Created 06/10/2012
% Modified 05/01/2012


log_prob = zeros(size(periods));
n_per = length(periods); % Number of periods tested
t_serie = t_serie - min(t_serie); % So the time series starts at time 0.

for i = 1:n_per

    T = periods(i);
    tlim=max(t_serie)-mod(max(t_serie),T);%To have a round number of cycles.
    % If not it induces artefacts to Schuster.
    t=t_serie(t_serie<=tlim); % selects a round number of cycles from the
    % timeseries.

    phase=mod(t,T)*2*pi/T; % Phase all the times from the timeseries
    % with respect to the period T.

    % Point where the Schuster walk ends
    end_walk = [sum(cos(phase)),sum(sin(phase))]; 
    % Distance from the origin
    D = norm(end_walk); 
    % Probability to reach the same point by random walk
    log_prob(i) = -D^2/length(t);
    % length(t) is the number of steps, i.e. the number of events in the
    % catalog.
end

end