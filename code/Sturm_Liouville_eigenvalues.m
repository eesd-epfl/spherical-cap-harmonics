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
function eigen_table = Sturm_Liouville_eigenvalues(x_o, K, BC, solver_increment_diag, solver_increment_tri, iterations_limit, N_eps1, N_eps2)
% This is a solver to return the first n-roots of the fraction associated
% Legendre functions by applying either Dirichlet or Neumann Boundary
% conditions
% x_o: the value of cos(x_0) that corresponds to the cap latitudinal
% half-angle.
% K: (Postive integer) the index of the m orders of the fraction associated Legendre
% function
% BC: the boundary conditions type. For Dirichlet we apply constant
% boundary of the expanded function (odd basis). For Neumann we apply
% boundaries for the first derivative and equate it by zero (even basis).
% BC is "odd" or "even".
% The function returns KxM table of roots.
% example:
%         eigen_table = Sturm_Liouville_eigenvalues(cosd(30), 15, 'odd')
% solver_increment_diag: the diagonal solver for all m = k; this should be
% smaller than the solver_increment_tri for the rest of the m > 0 roots.
% WARNING: Be careful with the roots extracted from this solver (check they make sense).
if nargin == 0
% This solves the example IN Table 1 Paper:
%     "Fully normalized spherical cap harmonics: application to the analysis
%     of sea-level data from TOPEX/POSEIDON and ERS-1", 
% Authors: Cheinway Hwang, Shin-Kuen Chen
    close all; clc; tic
%     N_eps = 10000;
    N_eps1 = 30000; % For Taylor (not used in this version)
    N_eps2 = 50000; % For RK     (not used in this version)
    K = 15 + 1;     % Number of degrees to be found
    x_o = cosd(10);
    BC = "even";
    solver_increment_diag = 1;
    solver_increment_tri = 5;
    iterations_limit = 50;
elseif nargin == 3
    solver_increment_diag = 5;
    solver_increment_tri = 15;
    iterations_limit = 20;
    N_eps1 = 30000; N_eps2 = 50000;
elseif nargin == 6
    N_eps1 = 30000; N_eps2 = 50000;
end
eigen_table = zeros(length(0:1:K));

if BC == "odd"
    f_BC = @(k, l, m, x) Dirichlet_BC(k, l, m, x, N_eps1, N_eps2);
elseif BC == "even"
    f_BC = @(k, l, m, x) Neumann_BC(k, l, m, x, N_eps1, N_eps2);
end
%% Start the solver hereJidong Zhao
% Find the first K-roots (eigenvalues) at \theta_o

% Solve for diagonals
for m = 1:K+1
    func = @(l) f_BC(m-1, l, m-1, x_o); % To find the final roots
    init_l = m; % Solve for the diagonal elements
    fun_val = fzero(func, init_l);
    idx = 0;
    while fun_val < 0 && idx <= iterations_limit
        init_l = init_l + solver_increment_diag;
        fun_val = fzero(func, init_l);
        idx = idx + 1;
    end
    eigen_table(m, m) = fun_val;
    fprintf("\nSolved the root l(m): %2.0f x%2.0f = %f\n", m-1, m-1, fun_val)
end

% Solve for the lower triangle
for m = 1:K+1
    for k = 1:K+1
        func = @(l) f_BC(k-1, l, m-1, x_o); % To find the final roots
        if m > k || m == k
            continue;
        else
            init_l = eigen_table(k-1, m) + solver_increment_tri;
            fun_val = fzero(func, init_l); idx = 0;
            while idx <= iterations_limit && round(fun_val, 4) <= round(init_l, 4)
                init_l = init_l + solver_increment_tri;
                fun_val = fzero(func, init_l);
                idx = idx + 1;
            end
            eigen_table(k, m) = fun_val;
            fprintf("\nSolved the root l(m): %2.0f x%2.0f = %f\n", k-1, m-1, fun_val)
        end
    end
end

if nargin == 0
    toc
%     save(["eigenvalues_k_", num2str(K), "_", BC, "_theta_",...
%         num2str(round(acosd(x_o), 0)), ".mat"], 'eigen_table.eigen_table')
    fprintf("\nThe solver time: %1.5f min\n", toc/60)
    % Define a color map similar to Paraview's shades
    red_color = zeros(1, 255); green_color = zeros(1, 255); blue_color = zeros(1, 255);
    red_color(1:127) = linspace(42,220,127); red_color(128:255) = linspace(220,174,127+1);
    green_color(1:127) = linspace(63,220,127); green_color(128:255) = linspace(220,0,127+1);
    blue_color(1:127) = linspace(181,220,127); blue_color(128:255) = linspace(220,22,127+1);
    ParaviewMap = [red_color', green_color', blue_color']./255;
    
    font = 'Amiri';
    figure('DefaultTextFontName', font, 'DefaultAxesFontName', font,'renderer','painters');
    b = bar3(eigen_table'); colorbar; xlabel("k"); ylabel("m"); zlabel("l(m)")
    for k = 1:length(b)
        zdata = b(k).ZData;
        b(k).CData = zdata;
        b(k).FaceColor = 'interp';
    end
    colormap(ParaviewMap);
    
    view(0,90)
    axis equal
    set(gca,'XTick', 0:1:K);
    set(gca,'YTick', 0:1:K);
end
end
%% Boundary conditions
function y = Dirichlet_BC(k, l, m, x, N_eps1, N_eps2)
    y = hypergeom([m-l, m+l+1], m+1, (1-x)./2);    % MATLAB's
%     y = hypergeom_2F1(k, m, l, (1-x)./2, N_eps1, N_eps2); % Mahmoud's
end

function y = Neumann_BC(k, l, m, x, N_eps1, N_eps2)
    y = hypergeom([m-l, m+l+1], m+1, (1-x)./2).*x.*l - (l-m).*hypergeom([m-l+1, m+l], m+1, (1-x)./2); % MATLAB's
%     y = hypergeom_2F1(k, m, l, (1-x)./2, N_eps1, N_eps2).*x.*l - (l-m).* hypergeom_2F1(k, m, l-1, (1-x)./2, N_eps1, N_eps2); % Mahmoud's
end