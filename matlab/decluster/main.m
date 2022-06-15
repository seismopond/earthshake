close all
clear all

%close all; clear all;
addpath(genpath('Code'));
fig_dir = 'Results/';
mkdir(fig_dir)

%% Load seismicity data
ncsn_file = 'test_data';
options.lat_max = 38;
seis = load_ncsn(ncsn_file, options);

figure;
%%%%%%%%% Plot seismicity %%%%%%%%%
scatter(seis.lon,seis.lat,5.^(seis.mag)/100,'r')

%%%%%%%%% Axes labels and tick marks %%%%%%%%%
ax = gca;
ax.FontSize = 12;
xlabel('Longitude')
ylabel('Latitude')
saveas(gca,[fig_dir,'map.png'])


%% Frequency-Magnitude Distribution
% Choose completeness magnitude
mc = 1.4;
% Calculate b-value using maximum likelihood method
b = log10(exp(1))/(mean(seis.mag(seis.mag>mc))-mc);
% Generate counts for plotting histogram
nedges = [0:0.1:8];
counts = histcounts(seis.mag,nedges);
counts_cum = zeros(1,length(counts));
dummy = 0;
for ii=1:length(counts)
    dummy = dummy + counts(end-ii+1);
    counts_cum(end-ii+1) = dummy;
end

% Figure
figure; grid on; grid minor; hold on; box on;
set(gcf,'position',[100,100,600,400])
scatter(nedges(1:end-1)+0.05,counts,50,'or','filled')
scatter(nedges(1:end-1)+0.05,counts_cum,50,'^k','filled')
xlim([0 6.2])
set(gca,'YScale','log')
plot([mc mc],get(gca,'YLim'),'k--')
plot([mc mc],get(gca,'YLim'),'k--')
plot([mc max(seis.mag)],[length(seis.mag(seis.mag>=mc)) 10^(log10(length(seis.mag(seis.mag>mc)))-b*(max(seis.mag)-mc))],'k','LineWidth',1)
text(min(get(gca,'XLim')) + 0.8*(max(get(gca,'XLim'))-min(get(gca,'XLim'))),...
        10^(log10(min(get(gca,'YLim'))) + 0.93*(log10(max(get(gca,'YLim')))-log10(min(get(gca,'YLim'))))),...
        ['M_c = ', num2str(mc,'%.2f')],'FontSize',12)
text(min(get(gca,'XLim')) + 0.8*(max(get(gca,'XLim'))-min(get(gca,'XLim'))),...
        10^(log10(min(get(gca,'YLim'))) + 0.87*(log10(max(get(gca,'YLim')))-log10(min(get(gca,'YLim'))))),...
        ['b = ', num2str(b,'%.2f')],'FontSize',12)
xlabel('Magnitude')
ylabel('Counts')
set(gca,'FontSize',14)
saveas(gca,[fig_dir,'fmd.png'])


%% Magnitude-time and Cumulative # of events
% Make a figure
figure; box on; hold on;
set(gcf,'position',[100,100,1000,300])
scatter(seis.time,seis.mag,10.^(seis.mag)/1000,'k')
ylabel('Magnitude')
ylim([0 6.2])
yyaxis right
plot(seis.time,[1:length(seis.time)],'LineWidth',2)
ylabel('Cumulative # of events')
ax = gca;
ax.FontSize = 14;
saveas(gca,[fig_dir,'mag_time.png'])


%% Latitude-time
% Make a figure
figure; box on; hold on;
set(gcf,'position',[100,100,1000,300])
scatter(seis.time,seis.lat,0.25,'k')
ylabel('Latitude')
ax = gca;
ax.FontSize = 14;
saveas(gca,[fig_dir,'lat_time.png'])


%% Calculating Nearest-Neighbor Distances (NND)
% Declustering parameters (mc and b are taken from the fit to the
% frequency-magnitude distribution)
D = 1.6;

% Remove seismicity with magnitude smaller than mc
index = (seis.mag>mc);
seis_crop.time = seis.time(index);
seis_crop.lat = seis.lat(index);
seis_crop.lon = seis.lon(index);
seis_crop.mag = seis.mag(index);
idx_Geysers = seis_crop.lat > 38.7 & seis_crop.lat < 38.9 & seis_crop.lon > -122.95 & seis_crop.lon < -122.65;

% Convert lat/lon to UTM x/y
[seis_x,seis_y,~] = deg2utm(seis_crop.lat,seis_crop.lon);
seis_x = 0.001*(seis_x-mean(seis_x));
seis_y = 0.001*(seis_y-mean(seis_y));
seis_time = years(seis_crop.time-seis_crop.time(1));

% Translated from python code used in Ge193/271 (Dynamics of Seismicity)
nnd = [];
R = [];
T = [];
for ii=1:length(seis_time)
    if rem(ii,100)==0
        disp([num2str(ii)])
    end
    tau_ij = [];
    r_ij = [];
    eta_ij = [];
    mag_ij = [];
    for kk=1:length(seis_time)
        if ii~=kk
            dt = seis_time(ii)-seis_time(kk);
            if dt>0
                tau_ij(end+1) = dt;
            else
                dt = inf;
                tau_ij(end+1) = dt;
            end
            d = sqrt((seis_x(ii)-seis_x(kk))^2 + (seis_y(ii)-seis_y(kk))^2);
            r_ij(end+1) = d;
            eta_ij(end+1) = dt*(d^D)*10^(-b*(seis_crop.mag(kk)-mc));
            mag_ij(end+1) = seis_crop.mag(kk);
        end
    end
    [eta_min,index] = min(eta_ij);
    nnd(end+1) = eta_min;
    R(end+1) = r_ij(index)*10^(-b*0.5*(mag_ij(index)-mc));
    T(end+1) = tau_ij(index)*10^(-b*0.5*(mag_ij(index)-mc));  
end
rescaled_T = log10(T);
rescaled_R = log10(R);

% Remove inf values
index = (isinf(rescaled_T) | isinf(rescaled_R));
rescaled_T = rescaled_T(~index);
rescaled_R = rescaled_R(~index);
select = ~idx_Geysers(~index)';

% Figure 1: rescaled distance-time
figure;
[values,centers] = hist3([rescaled_T(select)' rescaled_R(select)'],[41 41]);
s = pcolor(centers{:}, values.'); hold on;
s.FaceColor = 'interp';
h = colorbar;
ylabel(h,'Counts','FontSize',14)
colormap(jet)
%caxis([0 500])
limit_x = get(gca,'XLim');
limit_y = get(gca,'YLim');
plot(get(gca,'XLim'), -1*get(gca,'XLim')-2.6,'w')
xlim(limit_x)
ylim(limit_y)
xlabel('log(Rescaled time, T [yr])');
ylabel('log(Rescaled distance, R [km])');
set(gca,'FontSize',14)
saveas(gca,[fig_dir,'rescaled_distance_time.png'])

% Figure 2: Histogram of nearest-neighbor distance
figure;
h = histogram(log10(nnd(select)),50,'Normalization','pdf');
set(gcf,'position',[100,100,600,400])
hold on; box on;
plot([-2.6 -2.6],get(gca,'YLim'),'k--')
plot([-2.6 -2.6],get(gca,'YLim'),'k--')
xlabel('Nearest-neighbor distance, log \eta')
ylabel('Probability Density')
set(gca,'FontSize',14)
saveas(gca,[fig_dir,'nnd_histogram.png'])


%% Declustering
% Declustering ---- MUST PICK A CUTOFF VALUES HERE
declust = (log10(nnd)>-4 & nnd~=inf & ~idx_Geysers');
seis_declust.time = seis_crop.time(declust);
seis_declust.lon = seis_crop.lon(declust);
seis_declust.lat = seis_crop.lat(declust);
seis_declust.mag = seis_crop.mag(declust);


%% Latitude-time (Declustered)
% Make a figure
figure; box on; hold on;
set(gcf,'position',[100,100,1000,300])
scatter(seis_declust.time,seis_declust.lat,0.25,'k')
ylabel('Latitude')
ax = gca;
ax.FontSize = 14;
saveas(gca,[fig_dir,'lat_time_declustered.png'])

