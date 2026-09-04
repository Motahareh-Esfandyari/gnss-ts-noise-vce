%% ========================================================================
%  GNSS-TS-NoiseVCE : Noise Analysis of GNSS Coordinate Time Series
%  ------------------------------------------------------------------------
%  Script  : LSHE.m
%  Purpose : Standalone script version of the LS-HE algorithm (spectral
%            analysis + Fisher-test signal detection) applied to the AJAC
%            latitude series. See LS_HE.m for the function version.
%  Inputs  : AJAC.lat.txt
%  Outputs : DW, DT (detected frequencies/periods), xcap, ycap
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

fid = fopen('AJAC.lat.txt') ;
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

T(1, :) = 2 * sample ; m = length(y) ;
TT = t(end) - t(1) ; j = 1 ; alpha = 0.5 ;
while T(j, :) < TT
    T(j + 1, :) = T(j, :) * (1 + alpha*T(j, :)/TT) ; j = j + 1 ;
end

T(end, :) = [] ;

w = 2*pi./T ;
a = [ones(size(y, 1), 1), t] ; A = a ;
h = waitbar(0, 'Spectral values are being computed, please wait...') ;
for i = 1:size(w, 1)
    if rcond(A'*A) < 1e-20
        Aj = [cos(w(i, :)*t), sin(w(i, :)*t)] ;
        P(i, :) = NaN ; 
    else
        PAo = eye(m) - A*(A'*A)^(-1)*A' ; ecap0 = PAo*y ;
        Aj = [cos(w(i, :)*t), sin(w(i, :)*t)] ;
        P(i, :) = ecap0'*Aj*(Aj'*PAo*Aj)^(-1)*Aj'*ecap0 ;
    end
    waitbar(i/size(w, 1)) ;
end
close(h)

[idx, ~] = find(isnan(P)==1) ; P(idx, :) = [] ; T(idx, :) = [] ;

P_sort = sort(P) ; [~, IDX] = ismember(P_sort, P) ;
w2 = w(IDX) ; T_sort = T(IDX) ;

AA = A ; DW = [] ; DT = [] ; i = 1 ;
h = waitbar(0, 'Signal detection is in proccesing, please wait...') ;
for kk = 1:size(w, 1)
    k = size(w, 1) - kk + 1 ;
    Ak = [cos(w2(k, :)*t), sin(w2(k, :)*t)] ;
    if rcond([A, Ak]'*[A, Ak]) < 1e-15
        A = A ;
    else
        [m, n] = size(A) ;
        PAo  = eye(m) - A*(A'*A)^(-1)*A' ; ecap_0 = PAo*y ;
        PAko = eye(m) - [A, Ak]*([A, Ak]'*[A, Ak])^(-1)*[A, Ak]' ;
        ecap_a = PAko*y ;
        sigcap_a = ecap_a'*ecap_a/(m - n - 2) ;
        T2(i, :) = ecap_0'*Ak*(Ak'*PAo*Ak)^(-1)*Ak'*ecap_0/...
            (2*sigcap_a);
        F(i, :) = finv(1-Alpha/2, 2, m - n - 2) ;
        if T2(i) > F(i)
            A = [A, Ak] ;
            DW = [DW; w2(k, :)] ; DT = [DT; T_sort(k, :)] ;%break
        end
         i = i + 1 ;
    end
    waitbar(kk/size(w, 1))
end
close(h)
xcap = inv(A'*A)*A'*y ;
ycap = A*xcap ;