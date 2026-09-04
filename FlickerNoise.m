%% ========================================================================
%  GNSS-TS-NoiseVCE : Noise Analysis of GNSS Coordinate Time Series
%  ------------------------------------------------------------------------
%  Function: FlickerNoise.m
%  Purpose : Constructs the flicker (1/f) noise cofactor matrix Qf for the
%            LS-VCE stochastic model, following the approximation of
%            Zhang et al. (1997) / Amiri-Simkooei (2007).
%  Inputs  : t  - time vector (decimal years)
%            m  - number of observations
%  Outputs : Qf - m x m flicker noise cofactor matrix
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [Qf] = FlickerNoise(t, m)
Qf = (9/8)*eye(m) ;
h = waitbar(0, 'Flicker noise cofactor matrix construction...') ;
for k = 1:m
    tau = abs(t(k+1:end) - t(k)) ;
    tau = tau' ;
    Qf(k, k+1:end) = (9/8)*(1 - (log10(tau)/log10(2)+2)/24) ;
    tau = tau' ;
    Qf(k+1:end, k) = (9/8)*(1 - (log10(tau)/log10(2)+2)/24) ;
    waitbar(k/m)
end
close(h)
end