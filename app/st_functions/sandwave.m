function [delta,lambda] = sandwave(dep,tau,taucr,d50)
%
%-------header-------------------------------------------------------------
% NAME
%   sandwave.m
% PURPOSE
%   compute the height and wave length of bed ripples
% USAGE
%    [delta,lambda] = sandwave(d,tau,taucr,d50)
% INPUTS
%   dep - water depth (m)
%   tau - wave-current induced bed shear stress (N/m2)
%   taucr - critical bed shear stress (N/m2)
%   d50 - median grain size (m)
% RESULTS
%   delta - sandwave height (m)
%   lambda - sandwave wavelength (m)
% NOTES
%   uses van Rijn (1984) as in Soulsby, Dynamics of Marine Sands, Eq(83)
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2023
%--------------------------------------------------------------------------
%
lambda = 7.3*dep;   %ripple length

idx = tau>taucr & tau<26*taucr;
delta = zeros(size(dep)); Ts = delta;
%ripple height
Ts(idx) = (tau(idx)-taucr)/taucr;
delta(idx) = 0.11*dep(idx).*(d50./dep(idx)).^0.3.*(1-exp(-0.5*Ts(idx))).*(25-Ts(idx));