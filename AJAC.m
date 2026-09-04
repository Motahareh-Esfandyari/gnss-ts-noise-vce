%% ========================================================================
%  GNSS-TS-NoiseVCE : Noise Analysis of GNSS Coordinate Time Series
%  ------------------------------------------------------------------------
%  Script  : AJAC.m
%  Purpose : Main driver script. Runs LS-VCE on the AJAC station (Corsica,
%            France) latitude and longitude time series using four noise
%            model combinations:
%              I)   White noise only
%              II)  White + Flicker noise
%              III) White + Random Walk noise
%              IV)  White + Flicker + Random Walk noise
%            The functional model includes a linear trend plus the
%            harmonics detected by LS-HE (see LS_HE.m).
%  Inputs  : AJAC.lat.txt, AJAC.lon.txt, AJAC.rad.txt
%            AJAC_Lat_Harmonics.mat, AJAC_Lon_Harmonics.mat
%  Outputs : AJAC_Lat_I..IV.mat, AJAC_Lon_I..IV.mat
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

clear
clc
close all
format long g

fid_lat = fopen('AJAC.lat.txt') ;
fid_lon = fopen('AJAC.lon.txt') ;
fid_rad = fopen('AJAC.rad.txt') ;

[t, lat, s_lat] = DataReader(fid_lat) ;
[~, lon, s_lon] = DataReader(fid_lon) ;

sample = 0.0027 ; Alpha = 0.01 ;

%[DW, DT, T, P, xcap, ycap] = LS_HE(t, lat, sample, Alpha) ;

%save AJAC_Lat_Harmonics DW DT T P xcap ycap

%[DW, DT, T, P, xcap, ycap] = LS_HE(t, lon, sample, Alpha) ;

%save AJAC_Lon_Harmonics DW DT T P xcap ycap

%% Lat:

load AJAC_Lat_Harmonics

y = lat ;
m = length(y) ;

A = [ones(m, 1), t] ;
for k = 1:length(DW)
    A = [A, cos(DW(k)*t), sin(DW(k)*t)] ;
end

Epsilon = 1e-3 ;

Qw  = WhiteNoise(m)         ;
Qf  = FlickerNoise(t, m)    ;
Qrw = RandomWalkNoise(t, m) ; 

%I)White:

s0 = [1]' ;
Q(:, :, 1) = Qw ;

[scap1, c1, Error1, N1] = Linear_LS_VCE(s0, Q, Epsilon, y, A) ;

save AJAC_Lat_I scap1 c1 Error1 N1

%II)White+Flicker:

s0 = [1, 1]' ;
Q(:, :, 1) = Qw ;
Q(:, :, 2) = Qf ;

[scap2, c2, Error2, N2] = Linear_LS_VCE(s0, Q, Epsilon, y, A) ;

save AJAC_Lat_II scap2 c2 Error2 N2

%III)White+Random Walk:

s0 = [1, 1]' ;
Q(:, :, 1) = Qw  ;
Q(:, :, 2) = Qrw ;

[scap3, c3, Error3, N3] = Linear_LS_VCE(s0, Q, Epsilon, y, A) ;

save AJAC_Lat_III scap3 c3 Error3 N3

%White+Flicker+Random Walk:

s0 = [1, 1, 1]' ;
Q(:, :, 1) = Qw  ;
Q(:, :, 2) = Qf  ;
Q(:, :, 3) = Qrw ;

[scap4, c4, Error4, N4] = Linear_LS_VCE(s0, Q, Epsilon, y, A) ;

save AJAC_Lat_IV scap4 c4 Error4 N4

%% Lon:

load AJAC_Lon_Harmonics

A = [] ;

y = lon ;
m = length(y) ;

A = [ones(m, 1), t] ;
for k = 1:length(DW)
    A = [A, cos(DW(k)*t), sin(DW(k)*t)] ;
end

Epsilon = 1e-3 ;

Qw  = WhiteNoise(m)         ;
Qf  = FlickerNoise(t, m)    ;
Qrw = RandomWalkNoise(t, m) ; 

%I)White:

s0 = [1]' ;
Q(:, :, 1) = Qw ;

[scap1, c1, Error1, N1] = Linear_LS_VCE(s0, Q, Epsilon, y, A) ;

save AJAC_Lon_I scap1 c1 Error1 N1

%II)White+Flicker:

s0 = [1, 1]' ;
Q(:, :, 1) = Qw ;
Q(:, :, 2) = Qf ;

[scap2, c2, Error2, N2] = Linear_LS_VCE(s0, Q, Epsilon, y, A) ;

save AJAC_Lon_II scap2 c2 Error2 N2

%III)White+Random Walk:

s0 = [1, 1]' ;
Q(:, :, 1) = Qw  ;
Q(:, :, 2) = Qrw ;

[scap3, c3, Error3, N3] = Linear_LS_VCE(s0, Q, Epsilon, y, A) ;

save AJAC_Lon_III scap3 c3 Error3 N3

%White+Flicker+Random Walk:

s0 = [1, 1, 1]' ;
Q(:, :, 1) = Qw  ;
Q(:, :, 2) = Qf  ;
Q(:, :, 3) = Qrw ;

[scap4, c4, Error4, N4] = Linear_LS_VCE(s0, Q, Epsilon, y, A) ;

save AJAC_Lon_IV scap4 c4 Error4 N4




