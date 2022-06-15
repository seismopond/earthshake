function [b_val,b_err,mc] = calc_bval_maxc(mag,bin_catalog,bin_hist,plot_output)

% calc_bval calculates b-values in the Gutenberg-Richter Frequency
% Mangnitude Distribution (FMD). 
% 
% The method first estimates the magnitude of completeness using the 
% maximum curvature method (MAXC) with a simple corretion value of 0.2 
% (Woessner and Wiember, 2005, BSSA), i.e.
%           Mc = Mc (MAXC) + 0.2
%
% On top of this correction, we also apply the correction due to effects of
% binning in catalog magnitude. This correction is very small for modern
% catalogs, but can be quite large for historical data. This correction is
% described in Marzocchi and Sandri (2003), Annals Geophysics:
%           M_min = min(M) - bin_catalog/2
%
% After determining M_min, we calculate the b-value by using maximum
% likelihood estimate (MLE) (see Aki, 1965), i.e.
%           b = log10(e) / (mean(M) - M_min)
%
% The uncertainty of b-value is estimated with the formula derived by Shi
% and Bolt (1982), BSSA, i.e.
%           sigma_b = 2.30 * b^2 * sqrt( sum((M-mean(M))^2) / N(N-1) )
% -------------------------------------------------------------------------
% INPUT:
%   mag: A magnitude time series from the catalog that we want to compute 
%        mc and b-values
%   bin_catalog: Binning value of catalog magnitude. In other words, this
%                represents the rounding errors of catalog magnitudes.
%   bin_hist [OPTIONAL]: A user-specified bin size to be used to calculate
%                        histogram. If there is no input, this value would 
%                        default to 0.1.
%   plot_output [OPTIONAL]: A user-specified input whether the user wants
%                           the plot to display the outputs (0 = no plot, 
%                           1 = plot). Default to 1 (plot).
% -------------------------------------------------------------------------
% OUTPUT:
%   b_val: b-value of the frequency magnitude distribution
%   b_err: error estimate of the b-value calculated
%   mc: magnitude of completeness after applying the correction value of
%       0.2
%   Plot of the frequency magnitude distribution if plot_output = 1
% -------------------------------------------------------------------------
% AUTHOR:   Krittanon "Pond" Sirorattanakul
%           Seismological Laboratory
%           California Institute of Technology
% -------------------------------------------------------------------------
% REVISION: 1.0.0   19 APR 2021     Initial creation
% -------------------------------------------------------------------------


%% Initialization
% Initialize parameters if there is no input of specific values
if ~exist('bin_hist')
    bin_hist = 0.1;
end
if ~exist('plot_output')
    plot_output = 1;
end


%% Generate histogram
% Ensure that magnitude time series is a column vector
mag = mag(:);

% Generate counts for histogram
nedges = [floor(min(mag))-bin_hist/2 : bin_hist : max(mag)+bin_hist/2];
bin_centers = nedges(1:end-1) + bin_hist/2;
counts = histcounts(mag,nedges);

% Calculate cumulative number of events with magnitude larger than M_bin
counts_cum = zeros(1,length(counts));
for ii=1:length(counts)
    counts_cum(end-ii+1) = length(mag(mag>bin_centers(end-ii+1)));
end


%% Calculate magnitude of completeness
% Maximum curvature in cumulative number of events is at the bins with the
% largest counts
[~,idx] = max(counts);
mc = bin_centers(idx);

% Add a correction value of 0.2
mc = mc + 0.2;


%% Calculate b-value and associated estimated error
% Extract magnitude time serie with magnitude larger than mc to be used for
% b-value estimation
mag_b = mag(mag>=mc);
mag_mean = mean(mag_b);
N = length(mag_b);

b_val = log10(exp(1))/(mag_mean - (mc-bin_catalog/2));

b_err = 2.3 * (b_val^2) * sqrt(sum((mag_b - mag_mean).^2) / (N*(N-1)));


%% Plot frequency magnitude distribution

if plot_output == 1
    % Figure
    figure; grid on; grid minor; hold on; box on;
    set(gcf,'position',[100,100,600,400])
    scatter(bin_centers,counts,50,'or','filled')
    scatter(bin_centers,counts_cum,50,'^k','filled')
    set(gca,'YScale','log')
    ylim([0.8 1.5*max(counts_cum)])
    plot([mc max(mag)],[length(mag(mag>=mc)) 10^(log10(length(mag(mag>=mc)))-b_val*(max(mag)-mc))],'k','LineWidth',1)
    text(min(get(gca,'XLim')) + 0.65*(max(get(gca,'XLim'))-min(get(gca,'XLim'))),...
            10^(log10(min(get(gca,'YLim'))) + 0.72*(log10(max(get(gca,'YLim')))-log10(min(get(gca,'YLim'))))),...
            ['M_c = ', num2str(mc,'%.2f')],'FontSize',12)
    text(min(get(gca,'XLim')) + 0.65*(max(get(gca,'XLim'))-min(get(gca,'XLim'))),...
            10^(log10(min(get(gca,'YLim'))) + 0.66*(log10(max(get(gca,'YLim')))-log10(min(get(gca,'YLim'))))),...
            ['b = ', num2str(b_val,'%.2f'), ' \pm ', num2str(b_err,'%.2f')],'FontSize',12)
    xlabel('Magnitude')
    ylabel('Number of events')
    plot([mc mc],get(gca,'YLim'),'k--')
    legend('Non-cumulative','Cumulative','G-R law')
    set(gca,'FontSize',12)
end




end