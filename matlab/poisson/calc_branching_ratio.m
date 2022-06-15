function [bg_frac, bg_rate] = calc_branching_ratio(time)

% calc_branching_ratio.m estimates the fraction of seismicity that is
% considered non-interacting background seismicity using Gamma distribution
% (see Hainzl et al., 2006, BSSA for details). The script also estimates 
% the background seismicity rate.
% ------------------------------------------------------------------------
% INPUT:
%   time: A vector containings times in either datetime or double format.
% ------------------------------------------------------------------------
% OUTPUT:
%   bg_frac: Fraction of seismicity that is considered non-interacting
%            background seismicity.
%   bg_rate: Estimated background seismicity rate.
% ------------------------------------------------------------------------
% DEPENDENCIES:
%   N/A
% -------------------------------------------------------------------------
% AUTHOR:   Krittanon "Pond" Sirorattanakul
%           Seismological Laboratory
%           California Institute of Technology
% -------------------------------------------------------------------------
% REVISION: 1.0.0   20 APR 2021     Initial creation
%           1.0.1   10 MAY 2021     Calculate background rate
%           1.0.2   14 JUN 2022     Check input type and clean up
% -------------------------------------------------------------------------


% Evaluate the type of the input and calculate the interevent times
if class(time) == "datetime"
    dt = seconds(diff(time));
elseif class(time) == "double"
    dt = diff(time);
else
    error("Please input time in either datetime or double format...")
end

% Calculate mean seismicity rate
rate = (length(time)-1)/seconds(max(time)-min(time));

% Normalized the interevent times so that mean is one
dt_norm = rate*dt;

% Calculate beta parameter from Gamma distribution and the fraction of
% background seismicity (see Hainzl et al., 2006, BSSA for details)
beta = (std(dt_norm))^2;
bg_frac = 1/beta;
bg_rate = bg_frac*rate;


end