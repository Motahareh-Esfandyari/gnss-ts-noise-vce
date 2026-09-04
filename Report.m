%% ========================================================================
%  GNSS-TS-NoiseVCE : Noise Analysis of GNSS Coordinate Time Series
%  ------------------------------------------------------------------------
%  Script  : Report.m
%  Purpose : Collects the LS-VCE results of all four noise model scenarios
%            (I-IV) for the AJAC latitude and longitude series, computes
%            the standard deviations of the estimated variance components
%            from inv(N), and exports summary tables to Excel.
%  Inputs  : AJAC_Lat_I..IV.mat, AJAC_Lon_I..IV.mat
%  Outputs : AJAC_Lat_Noise.xls, AJAC_Lon_Noise.xls
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

clear
clc
close all
format long g

Station = {'AJAC'} ;

%% Lat:

load AJAC_Lat_I
load AJAC_Lat_II
load AJAC_Lat_III
load AJAC_Lat_IV

StD1 = diag(sqrt(inv(N1))) ;
StD2 = diag(sqrt(inv(N2))) ;
StD3 = diag(sqrt(inv(N3))) ;
StD4 = diag(sqrt(inv(N4))) ;

Report1 = num2cell(zeros(5, 7)) ;
Report1(1, :) = [{'Station'},...
    {'Sw'}, {'Sf'}, {'Srw'},...
    {'StD of Sw'}, {'StD of Sf'}, {'StD of Srw'}] ;

Report1(2, :) = [Station, num2cell(scap1), {'-'}, {'-'},...
    num2cell(StD1), {'-'}, {'-'}] ;

Report1(3, :) = [Station, num2cell(scap2(1)), num2cell(scap2(2)), {'-'},...
    num2cell(StD2(1)), num2cell(StD2(2)), {'-'}] ;

Report1(4, :) = [Station, num2cell(scap3(1)), {'-'}, num2cell(scap3(2)),...
    num2cell(StD3(1)), {'-'}, num2cell(StD3(2))] ;

Report1(5, :) = [Station, num2cell(scap4(1)), num2cell(scap4(2)), num2cell(scap4(3)),...
    num2cell(StD4(1)), num2cell(StD4(2)), num2cell(StD4(3))] ;

%% Lon:

load AJAC_Lon_I
load AJAC_Lon_II
load AJAC_Lon_III
load AJAC_Lon_IV

StD1 = diag(sqrt(inv(N1))) ;
StD2 = diag(sqrt(inv(N2))) ;
StD3 = diag(sqrt(inv(N3))) ;
StD4 = diag(sqrt(inv(N4))) ;

Report2 = num2cell(zeros(5, 7)) ;
Report2(1, :) = [{'Station'},...
    {'Sw'}, {'Sf'}, {'Srw'},...
    {'StD of Sw'}, {'StD of Sf'}, {'StD of Srw'}] ;

Report2(2, :) = [Station, num2cell(scap1), {'-'}, {'-'},...
    num2cell(StD1), {'-'}, {'-'}] ;

Report2(3, :) = [Station, num2cell(scap2(1)), num2cell(scap2(2)), {'-'},...
    num2cell(StD2(1)), num2cell(StD2(2)), {'-'}] ;

Report2(4, :) = [Station, num2cell(scap3(1)), {'-'}, num2cell(scap3(2)),...
    num2cell(StD3(1)), {'-'}, num2cell(StD3(2))] ;

Report2(5, :) = [Station, num2cell(scap4(1)), num2cell(scap4(2)), num2cell(scap4(3)),...
    num2cell(StD4(1)), num2cell(StD4(2)), num2cell(StD4(3))] ;

%% Write to Excel:

xlswrite('AJAC_Lat_Noise', Report1)
xlswrite('AJAC_Lon_Noise', Report2)