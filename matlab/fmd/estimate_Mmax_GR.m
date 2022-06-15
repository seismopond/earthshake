function [Mmax_fmd,Mmax] = estimate_Mmax_GR(mag,b_val,mc)

% estimate_Mmax_GR.m calculates maximum earthquake magitude expected from
% Gutenberg-Richter law (Poisson point process).
% -------------------------------------------------------------------------
% INPUT:
%   mag: A magnitude time series from the catalog that we want to compute 
%        mc and b-values
%   b_val: b-value of the frequency magnitude distribution
%   mc: magnitude of completeness after applying the correction value of
%       0.2
% -------------------------------------------------------------------------
% OUTPUT:
%   Mmax_fmd: Maximum earthquake magitude expected from Gutenberg-Richter 
%             law
%   Mmax: Magnitude of the largest earthquake in the catalog
% -------------------------------------------------------------------------
% AUTHOR:   Krittanon "Pond" Sirorattanakul
%           Seismological Laboratory
%           California Institute of Technology
% -------------------------------------------------------------------------
% REVISION: 1.0.0   19 APR 2021     Initial creation
% -------------------------------------------------------------------------

if max(mag) <= mc
    error('Largest magnitude is smaller than mc')
end
a = log10(length(mag(mag>=mc))) + b_val*mc;
Mmax_fmd = a / b_val;
Mmax = max(mag);

end