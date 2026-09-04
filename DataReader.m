%% ========================================================================
%  GNSS-TS-NoiseVCE : Noise Analysis of GNSS Coordinate Time Series
%  ------------------------------------------------------------------------
%  Function: DataReader.m
%  Purpose : Reads a GNSS coordinate time series text file (fixed-width
%            format) and returns time, observation, and sigma vectors.
%  Inputs  : fid - file identifier from fopen()
%  Outputs : t   - time vector (decimal years)
%            y   - coordinate time series
%            s   - formal standard deviations
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [t, y, s] = DataReader(fid)

% clear
% clc
% close all
% format long g

%fid = fopen('AJAC.lat.txt') ;

c = 1 ;
while ~feof(fid)
    
    txt = fgetl(fid) ;
    txt = txt(1:57) ;
    txt = str2num(txt) ;
    
    t(c, :) = txt(1) ;
    y(c, :) = txt(2) ;
    s(c, :) = txt(3) ;
    
    c = c + 1 ;
      
end