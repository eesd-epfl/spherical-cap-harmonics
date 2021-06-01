%*************************************************************************%
% All rights reserved (C) to the authors: Mahmoud SHAQFA, Gary CHOI       %
%		 and Katrin BEYER                                                 %
%                                                                         %
% M. Shaqfa Contact:                                                      %
% Earthquake Engineering and Structural Dynamics Laboratory (EESD),       %
% School of Architecture, Civil and Environmental Engineering (ENAC),     %
% Ecole polytechnique federale de Lausanne (EPFL),                        %
% CH-1015 Lausanne, Switzerland.                                          %
%               Tel.: +41 21 69 33297                                     %
%               Email: mahmoud.shaqfa@epfl.ch                             %
%                                                                         %
% G. Choi Contact:                                                        %
% Department of Mathematics, Massachusetts Institute of Technology (MIT)  %
% Cambridge, MA, USA                                                      %
%               Email: ptchoi@mit.edu                                     %
%                                                                         %
% K. Beyer Contact:                                                       %
%               Email: katrin.beyer@epfl.ch                               %
%*************************************************************************%
% This code includes implementations for:                                 %
%				- Spherical Cap Harmonics (SCH)                           %
%				- Spherical Harmonics (SH)                                %
%				- HemiSpherical Harmonics (HSH)                           %
% This code is part of the paper: "Spherical Cap Harmonic Analysis(SCHA)..%
%	 for Characterising the Morphology of Rough Surface Patches"          %
%                                                                         %
%*************************************************************************%
% This library is free software; you can redistribute it and/or modify	  %
% it under the terms of the GNU Lesser General Public License as published%
% by the Free Software Foundation; either version 2.1 of the License, or  %
% (at your option) any later version.                                     %
%                                                                         %
% This library is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %
% See the GNU Lesser General Public License for more details.        	  %
% You should have received a copy of the GNU Lesser General Public License%
% along with this library; if not, write to the Free Software Foundation, %
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA       %
%*************************************************************************%
% Author of this file: Mahmoud S. Shaqfa
function legendre_basis = fractional_associated_Legendre_functions(theta_c, K, BC, eigen_table, thetas, normalization_method, N_eps)
% This function computes the fully normalized fractional associated
% Legendre fuctions at non-integral degrees and integer orders and
% indicies.

% theta_o: in rads
% thetas: one column vector (\theta) resulted from the parametrized verticies (rads).
% normalization_method: "Hainse" or "Hwang" normalization
if nargin == 0
    % Testing and benchmarking area
%     close all; clear; clc;
    tic;
    K = 39;
    theta_c = 10;
    thetas = deg2rad(theta_c:-0.1:0);
    theta_c = deg2rad(10);
    BC = "even";
    calculate_eigens = false;
    normalization_method = "Hainse";
%     N_eps = 400000; % The stable number for RK
    N_eps = 500000;
    if calculate_eigens
        eigen_table = Sturm_Liouville_eigenvalues(x_c, K, BC, 2, 30);
    else
        if BC == "even"
            % For testing purposes only
            eigen_table = load("eigenvalues_k_40_even_theta_10.mat");
        elseif BC == "odd"
            eigen_table = load("eigenvalues_k_8_odd.mat"); % For testing purposes only (not used in this paper)
        end
        eigen_table = eigen_table.eigen_table;
    end
end

% Calculate the cos of thetas involved
x_c = cos(theta_c);
x = cos(thetas);

% Check inputs' size (1D and 2D matrix)
sz = length(size(thetas)); % For reconstruction matricies only
if isempty(find(size(thetas) == 1, 1))
    sz = sz + 1;
end

if sz == 2
    legendre_basis = zeros(length(x), (K+1));
elseif sz == 3
    temp_sz = size(thetas); temp_sz(3) = (K+1);
    legendre_basis = zeros(temp_sz);
end
%% Calculate the associated non-normalized basis
for m = 0:K
    l = eigen_table(K+1, m+1);
    if normalization_method == "Hwang" %(Not used in this paper)
        normalization = fully_normalized_factor(m, l, x_c, BC);
        if sz == 2
            legendre_basis(:, m+1) = normalization .* associated_basis(m, l, x, x_c);
        elseif sz == 3
            legendre_basis(:, :, m+1) = normalization .* associated_basis(m, l, x, x_c);
        end
    elseif normalization_method == "Hainse" % The default basis used in this paper
        normalization = fully_normalized_factor2(m, l);
        if sz == 2
            legendre_basis(:, m+1) = normalization .* associated_basis2(K, m, l, x, N_eps);
        elseif sz == 3
            legendre_basis(:, :, m+1) = normalization .* associated_basis2(K, m, l, x, N_eps);
        end
    end
end

if nargin == 0
    fprintf("\nThe solver time: %1.5f min\n", toc/60)
    figure
    for m = 1:K+1
        plot(x, (legendre_basis(:, m))); xlim([min(x) 1]);
        hold on
    end
end
end
%% Other functions
function normalization = fully_normalized_factor(m, l, x_c, BC)
% Check the imag and real normalization factor, please!
if m == 0
    delta_m = 1;
else
    delta_m = 0;
end

if BC == "odd"
    G1 = G(m, l, x_c);
    H = (l-m) * hypergeom([m-l+1, m+l], m+1, x_c) - l * hypergeom([m-l, m+l+1], m+1, x_c);
    normalization = sqrt((2-delta_m)*(1-x_c)*(2*l+1)/abs(G1*H)); % with abs(GG*H) sign correction!
elseif BC == "even"
    F = hypergeom([m-l, m+l+1], m+1, x_c);
    G1 = G(m, l-1, x_c);
    G2 = G(m, l, x_c);
    J = (l-m) * G1 - l * x_c * G2 + ...
        ((l-m)/(l+m)) * hypergeom([m-l+1, m+l], m+1, x_c) - ...
        x_c * hypergeom([m-l, m+l+1], m+1, x_c);
    normalization = sqrt(((2-delta_m)*(1-x_c)*(2*l+1))/abs(J*F)); % with abs(J*F) sign correction!
end
% normalization = normalization * sqrt(1/(2*pi*(1-x_c))); % Pi normalization
end

function p_l_m = associated_basis(m, l, x, x_c)
% Hwang (1997) (not used in this paper)
x_ = (1-x)./2;
p_l_m = ((1-x.^2)./(1-x_c^2)).^(m/2) .* hypergeom([m-l, m+l+1], m+1, x_);
end

function p_l_m = associated_basis2(k, m, l, x, N_eps)
% Haines (1985)
x_ = (1-x)./2;
% p_l_m = (1-x.^2).^(m/2) .* hypergeom([m-l, m+l+1], m+1, x_);
p_l_m = (1-x.^2).^(m/2) .* hypergeom_2F1(k, m, l, x_, N_eps);
end

function gk_1 = G(m, l, x)
% This is a simple implementation for the hypergeometric series
% This function should not be used for integr values of l.
% This series diverges at x = -1 (Ordinary Spherical Harmonics)
% The values here only valid for theta != 90 or 180 degs.
% The values of m must be only postive integers
% Otherwise, this function (normalization) should be replaced with a proper
% one.

k = 1;
betak = 0;
ak = 1;
eps = inf;
bl_m = 0;

if m == 0
    bl_m = 0;
else
    for ii = 0:2*m-1
        bl_m = bl_m + (1/(l - m + 1 + ii));
    end
end
x_o = 0.5 * (1 - x);
gk_1 = 0;
gk = 0;
while eps > 10^-5 && k < 10^5
    ak = ak * (m - l + k - 1) * (m + l + k) / (k * (m + k));
    betak = betak + (1 / ((m + l + 1 + k) * (m - l + k)));
    temp_betak = betak * -1 * (2 * l + 1);
    if temp_betak == -inf || temp_betak == inf || isnan(temp_betak) || betak == inf
        break;
    else
        gk = ak * (bl_m + temp_betak) * (x_o)^k;
    end
    eps = abs(gk - gk_1);
    gk_1 = gk;
    k = k + 1;
end
end

function Km_n = fully_normalized_factor2(m, l)
% This is a simple implementation for the approximate normalization factor
% derived by Haines (1985a)

if m == 0
    Km_n =  1;
else
    p = (l/m)^2 - 1;
    e1 = (-1/(12 * m)) * (1 + (1/p));
    e2 = (1/(360 * m^3)) * (1 + (3/(p^2)) + (4/(p^3)));
    Km_n = (2^(-m) / (sqrt(m * pi))) * ((l+m)/(l-m))^(0.5*l + 0.25) *...
        p^(0.5*m) *  exp(e1 + e2);
end
end

function Km_n = fully_normalized_factor3(m, l, x_c)
% This is a simple implementation for the approximate normalization factor
% derived by Haines (1985a)

if m == 0
    Km_n =  1;
else
    Km_n = (sqrt(2)/(2^m * factorial(m))) * sqrt(gamma(m+l+1)/gamma(m+l-1));
    Km_n = Km_n * sqrt(1/(2*pi*(1-x_c))); % Pi normalization
end
end