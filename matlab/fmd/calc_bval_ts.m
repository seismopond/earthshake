function [b_time,b_val,b_err,mc,Mmax_fmd] = calc_bval_ts(time,mag,params,options)

%%% NEED TO ADD COMPATIBILITY FOR NON-DATETIME
%%% CURRENTLY ONLY SUPPORT MAXC AND MLE METHOD


%% Setup
% Convert both time and mag to column vector if they aren't already a
% column vector
time = time(:);
mag = mag(:);

% Ensure that the values of time is properly sorted
% Need sorting? Check if these two quantities are the same
if sum(time == sort(time)) ~= length(time)
    [time,index] = sort(time);
    mag = mag(index);
end

% Check if time and mag have equal length
if length(time) ~= length(mag)
    error('"time" and "mag" do not have the same length.')
end

% Check if the format of time is datetime
if ~isa(time,'datetime')
    error(['"time" variable is not in a datetime format. Only datetime format is supported for calc_bval_ts.m.'])
end



% Initialize variables if not specified by the user
% Return errors for the following
if ~exist('params','var')
    error(['Please input related parameters for b-value calculations.'])
end
if ~isfield(params,'bin_catalog')
    error('Please input binning value of catalog magnitude to params.bin_catalog')
end

% Pre-allocate default values for the following
if ~isfield(params,'bin_hist')
    params.bin_hist = 0.2;
    disp('Using default bin_hist = 0.2 ...')
end
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

if ~exist('options','var')
    options.plot_output = 1;
end
if ~isfield(options,'plot_output')
    options.plot_output = 1;
end
if ~isfield(options,'plot_mc')
    options.plot_mc = 1;
end
if ~isfield(options,'plot_mmax')
    options.plot_mmax = 0;
end


%% Pick data range of interest
index = ((time > params.start_date) & (time < params.end_date));
mag = mag(index);
time = time(index);



%% Calculate b values, Mc, and Mmax with sliding window
ii = 1;

while(params.sliding_windows*ii + params.num_points < length(mag))
    b_time(ii) = time(params.sliding_windows*ii) + ...
        0.5 * (time(params.sliding_windows*ii+params.num_points) - time(params.sliding_windows*ii));
    [b_val(ii),b_err(ii),mc(ii)] = calc_bval_maxc(...
        mag(params.sliding_windows*ii:params.sliding_windows*ii+params.num_points), ...
        params.bin_catalog, params.bin_hist, 0);
    [Mmax_fmd(ii),Mmax(ii)] = estimate_Mmax_GR(...
        mag(params.sliding_windows*ii:params.sliding_windows*ii+params.num_points),...
        b_val(ii),mc(ii));
    ii = ii + 1;
end

b_time = b_time(:);
b_val = b_val(:);
b_err = b_err(:);
mc = mc(:);
Mmax_fmd = Mmax_fmd(:);
Mmax = Mmax(:);


% Make plots
figure;
set(gcf,'position',[100,0,800,600])

% Subplot 1 - Cumulative number of events vs. time
subplot(3,1,1)
plot(time,[1:length(time)]','k','LineWidth',1)
grid on; grid minor; hold on; box on;
ylabel({'Cumulative number'; 'of events'})
set(gca,'xticklabel',{[]})
xlim([b_time(1) b_time(end)])
set(gca,'FontSize',12)

% Subplot 2 - Mc and b-values vs. time
subplot(3,1,2)
grid on; grid minor; hold on; box on;
yyaxis left
plot(b_time,b_val,'k','LineWidth',1)
plot(b_time,b_val-b_err,'k--','LineWidth',1)
plot(b_time,b_val+b_err,'k--','LineWidth',1)
ylabel('b-value')
yyaxis right
plot(b_time,mc,'r','LineWidth',1)
ylabel('M_c')
set(gca,'xticklabel',{[]})
xlim([b_time(1) b_time(end)])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'r';
set(gca,'FontSize',12)

% Subplot 3 - Mmax observed and predicted from GR law
subplot(3,1,3)
grid on; grid minor; hold on; box on;
plot(b_time,Mmax,'k','LineWidth',1)
plot(b_time,Mmax_fmd,'k--','LineWidth',1)
legend('Observed','G-R law')
ylabel('M_{max}')
xlim([b_time(1) b_time(end)])
set(gca,'FontSize',12)



end