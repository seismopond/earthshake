function calc_Poisson_ts(time,params,options)

%%% NEED TO ADD COMPATIBILITY FOR NON-DATETIME
%%% CURRENTLY ONLY SUPPORT MAXC AND MLE METHOD
% Last update May 10, 2021

%% Setup
% Convert both time and mag to column vector if they aren't already a
% column vector
time = time(:);

% Check if the format of time is datetime
if ~isa(time,'datetime')
    error(['"time" variable is not in a datetime format. Only datetime format is supported for calc_bval_ts.m.'])
end

% Initialize variables if not specified by the user
% Return errors for the following
if ~exist('params','var')
    params.start_date = time(1);
end

% Pre-allocate default values for the following
if ~isfield(params,'start_date')
    params.start_date = time(1);
end
if ~isfield(params,'end_date')
    params.end_date = time(end);
end

if ~isfield(params,'num_points')
    params.num_points = 100;
end
if ~isfield(params,'sliding_windows')
    params.sliding_windows = 20;
end

% if ~exist('options','var')
%     options.plot_output = 1;
% end
% if ~isfield(options,'plot_output')
%     options.plot_output = 1;
% end
% if ~isfield(options,'plot_mc')
%     options.plot_mc = 1;
% end
% if ~isfield(options,'plot_mmax')
%     options.plot_mmax = 0;
% end


%% Pick data range of interest
index = ((time > params.start_date) & (time < params.end_date));
time = time(index);



%% Calculate b values, Mc, and Mmax with sliding window
ii = 1;

while(params.sliding_windows*ii + params.num_points < length(time))      
    cov_time(ii) = time(params.sliding_windows*ii) + ...
        0.5 * (time(params.sliding_windows*ii+params.num_points) - time(params.sliding_windows*ii));
    [cov(ii),cov_min(ii),cov_max(ii)] = calc_cov(...
        time(params.sliding_windows*ii:params.sliding_windows*ii+params.num_points));
    [bg(ii),bg_rate(ii)] = calc_branching_ratio(...
        time(params.sliding_windows*ii:params.sliding_windows*ii+params.num_points));

    ii = ii + 1;
end

cov_time = cov_time(:);
cov = cov(:);
cov_min = cov_min(:);
cov_max = cov_max(:);
bg = bg(:);
bg_rate = bg_rate(:);


% Make plots
figure;
set(gcf,'position',[0,-200,600,800])

% Subplot 1 - Cumulative number of events vs. time
subplot(4,1,1)
plot(time,[1:length(time)]','k','LineWidth',1)
grid on; grid minor; hold on; box on;
ylabel({'Cumulative number'; 'of events'})
set(gca,'xticklabel',{[]},'FontSize',10)
set(gca,'FontSize',10)
xlim([time(1) time(end)])

% Subplot 2 - Coefficient of variations (CoV) vs. time
subplot(4,1,2)
plot(cov_time,cov,'k','LineWidth',1)
grid on; grid minor; hold on; box on;
plot(cov_time,cov_min,'k--','LineWidth',1)
plot(cov_time,cov_max,'k--','LineWidth',1)
ylabel({'Coefficient'; 'of Variation'})
set(gca,'xticklabel',{[]},'FontSize',10)
set(gca,'FontSize',10)
xlim([time(1) time(end)])

% Subplot 3 - Branching ratio (fraction of mainshock)
subplot(4,1,3)
grid on; grid minor; hold on; box on;
plot(cov_time,bg,'k','LineWidth',1)
ylabel({'Fraction of'; 'Mainshocks'})
set(gca,'xticklabel',{[]},'FontSize',10)
set(gca,'FontSize',10)
xlim([time(1) time(end)])

% Subplot 4 - Background rate (gamma distribution fitting) 
subplot(4,1,4)
grid on; grid minor; hold on; box on;
plot(cov_time,bg_rate,'k','LineWidth',1)
ylabel({'Background rate'; '(event per day)'})
set(gca,'FontSize',10)
xlim([time(1) time(end)])



end