function [cov,cov_min,cov_max] = calc_cov(time)

% calc_cov.m calculates Coefficient of Variations (CoV) and estimate 95
% percent confidence interval range of CoV in which the process would be
% homogeneous Poisson.
% ------------------------------------------------------------------------
% INPUT:
%   time: A vector containings times in either datetime or double format.
% ------------------------------------------------------------------------
% OUTPUT:
%   cov: Coefficient of Variations (CoV) calculated for the input dataset.
%   cov_min: Minimum Coefficient of Variations (CoV) for the process to be
%            considered Poissonian at 95 percent confidence interval.
%   cov_max: Maximum Coefficient of Variations (CoV) for the process to be
%            considered Poissonian at 95 percent confidence interval.
% ------------------------------------------------------------------------
% DEPENDENCIES:
%   N/A
% -------------------------------------------------------------------------
% AUTHOR:   Krittanon "Pond" Sirorattanakul
%           Seismological Laboratory
%           California Institute of Technology
% -------------------------------------------------------------------------
% REVISION: 1.0.0   20 APR 2021     Initial creation
%           1.0.1   14 JUN 2022     Check input type and clean up
% -------------------------------------------------------------------------

% Evaluate the type of the input and calculate the interevent times
if class(time) == "datetime"
    dt = seconds(diff(time));
elseif class(time) == "double"
    dt = diff(time);
else
    error("Please input time in either datetime or double format...")
end

% Calculate the Coefficient of Variations (CoV)
cov = std(dt) / mean(dt);

% Generate 200 realizations of simulated Poisson process and estimate 95
% percent confidence interval arising from given Poisson rate
rate = 1 / mean(dt);
data_sim = exprnd(1/rate,length(time),200);
cov_sim = std(data_sim)./mean(data_sim);
cov_min = prctile(cov_sim,2.5);
cov_max = prctile(cov_sim,97.5); 


end