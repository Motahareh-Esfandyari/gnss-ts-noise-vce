%% ========================================================================
%  GNSS-TS-NoiseVCE : Noise Analysis of GNSS Coordinate Time Series
%  ------------------------------------------------------------------------
%  Function: LS_HE.m
%  Purpose : Least-Squares Harmonic Estimation (LS-HE). Computes the
%            power spectrum of a time series and detects statistically
%            significant frequencies via a Fisher test (Type I error
%            Alpha), following Amiri-Simkooei (2007).
%  ------------------------------------------------------------------------
%  Adapted and applied by : Motahareh Esfandyari-Kaloukan
%  Original implementation: Hamed Karimi (see credit below)
%  ========================================================================

 function [DW, DT, T, P, xcap, ycap] = LS_HE(t, y, sample, Alpha)
%% LS_HE.M
% Description: LS_HE performs Least Squares - Harmonic Estimation on a time
% series. This is uni-variate one of LS-HE(we have multi-variate too).
% Input:
%   t : Time argument related to the time series in y
%   y : The vector of time series sort by time in the vector t
%   Qy : (Co)variance matrix of time series, it might consist of just white
%   noise or other noises such as colored noise or flicker noise or random
%   walk noise etc.
%   sample : The sample of time(or time variation rate)
%   Alpha : Type I error, use in statistical testing for detection the
%   frequencies exist in the time series(use in finv())
% Output:
%   DW : The frequencies that are significant in signal detection test by
%   fisher test, under the type I error Alpha
%   DT : The detected periods related to the detected frequencies DW
%   T : This is all entire periods in LSHE algorithm, this is related to
%   the spectral values vector P generates in LSHE
%   P : Spectral values related to the periods in the vector T
%   xcap : Least squares estimation of harmonic coefficients related to
%   each frequency in the vector DW, of course the length of xcap equals to
%   2 + 2*length(DW), because of 2 number of trend coefficients and
%   2*length(DW) related to [cos(DW(i)*t), sin(DW(i)*t)]
%   ycap : Reconstructed time series(reconstructed y)
% Note: This methodology is according to:
% Amiri-Simkooei, A. (2007). Least-squares variance component estimation:
% theory and GPS applications. 
%==========================================================================
% Written by: Hamed Karimi, Geodesy student at University of Isfahan
% PZSD-DANACO.
%==========================================================================

%% START
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
        P(i, :) = ecap0'*Aj*(Aj'*PAo*Aj)^(-1)*Aj'*ecap0 ;    %power omega
    end
    waitbar(i/size(w, 1)) ;
end
close(h)

[idx, ~] = find(isnan(P)==1) ; P(idx, :) = [] ; T(idx, :) = [] ;

P_sort = sort(P) ; [~, IDX] = ismember(P_sort, P) ;
wn2 = w(IDX) ; T_sort = T(IDX) ;

AA = A ; DW = [] ; DT = [] ; i = 1 ;
h = waitbar(0, 'Signal detection is in proccesing, please wait...') ;
for kk = 1:size(w, 1)
    k = size(w, 1) - kk + 1 ;
    Ak = [cos(wn2(k, :)*t), sin(wn2(k, :)*t)] ;
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
            DW = [DW; wn2(k, :)] ; DT = [DT; T_sort(k, :)] ;%break
        end
         i = i + 1 ;
    end
    waitbar(kk/size(w, 1))
end
close(h)
xcap = inv(A'*A)*A'*y ;
ycap = A*xcap ;

end