%% ========================================================================
%  GNSS-TS-NoiseVCE : Noise Analysis of GNSS Coordinate Time Series
%  ------------------------------------------------------------------------
%  Function: WhiteNoise.m
%  Purpose : Returns the white noise cofactor matrix (identity) for the
%            LS-VCE stochastic model.
%  Inputs  : m  - number of observations
%  Outputs : Qw - m x m identity matrix
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [Qw] = WhiteNoise(m)
Qw = eye(m) ;
end