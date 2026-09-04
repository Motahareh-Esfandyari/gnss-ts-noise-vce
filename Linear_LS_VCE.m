%% ========================================================================
%  GNSS-TS-NoiseVCE : Noise Analysis of GNSS Coordinate Time Series
%  ------------------------------------------------------------------------
%  Function: Linear_LS_VCE.m
%  Purpose : Iterative Least-Squares Variance Component Estimation (LS-VCE)
%            for a linear model y = Ax + e with stochastic model
%            Qy = sum_k( s_k * Q_k ). Estimates the variance components
%            s_k and their precision, following Teunissen & Amiri-Simkooei
%            (2008).
%  Inputs  : s0      - initial variance components (p x 1)
%            Q       - m x m x p stack of cofactor matrices
%            Epsilon - convergence threshold on ||s_new - s_old||
%            y       - observation vector (m x 1)
%            A       - design matrix (m x n)
%  Outputs : scap    - estimated variance components
%            c       - number of iterations
%            Error0  - convergence history
%            N       - normal matrix (inv(N) gives covariance of scap)
%  ------------------------------------------------------------------------
%  Author  : Motahareh Esfandyari-Kaloukan
%  ========================================================================

function [scap, c, Error0, N] = Linear_LS_VCE(s0, Q, Epsilon, y, A)
m = length(y) ; p = length(s0) ; s = s0 ;
Error = 1000 ; Error0 = [] ;
c = 0
while Error>Epsilon
    Qy = zeros(m) ;
    for k = 1:p
        Qy = Qy + s(k)*Q(:, :, k) ;
    end
    %invQy = Qy^(-1) ;
    invQy = inv(Qy) ;
    %PAo = eye(m) - A*(A'*invQy*A)^(-1)*A'*invQy ;
    PAo = eye(m) - A*inv(A'*invQy*A)*A'*invQy ;
    ecap = PAo*y ;
    for k = 1:p
        L0 = ecap'*invQy*Q(:, :, k)*invQy*ecap ;
        L(k, :) = (1/2)*trace(L0) ;
        for l = 1:p
            N0 = Q(:, :, k)*invQy ;
            N0 = N0*PAo ;
            N0 = N0*Q(:, :, l) ;
            N0 = N0*invQy ;
            N0 = N0*PAo ;
            %N0 = Q(:, :, k)*invQy*PAo*Q(:, :, l)*invQy*PAo ;
            N(k, l) = (1/2)*trace(N0) ;
        end
    end
    scap = inv(N)*L ;
    Error = norm(scap-s0)
    s0 = scap ;
    s = scap ;
    c = c + 1
    Error0(c, :) = Error ;
end
end