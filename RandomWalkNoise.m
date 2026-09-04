%% ========================================================================
%  GNSS-TS-NoiseVCE : Noise Analysis of GNSS Coordinate Time Series
%  ------------------------------------------------------------------------
%  Function: RandomWalkNoise.m
%  Purpose : Constructs the random walk noise cofactor matrix Qrw for the
%            LS-VCE stochastic model, scaled by the sampling frequency.
%  Inputs  : t   - time vector (decimal years)
%            m   - number of observations
%  Outputs : Qrw - m x m random walk noise cofactor matrix
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [Qrw] = RandomWalkNoise(t, m)
T = t(end) - t(1) ;
fs = (m-1)/T ;
h = waitbar(0, 'Random walk noise cofactor matrix construction...') ;
for k = 1:m
    Qrw(k:m, k) = k ;
    Qrw(k, k:m) = k ;
    waitbar(k/m)
end
close(h)
Qrw = fs^(-1)*Qrw ;
end