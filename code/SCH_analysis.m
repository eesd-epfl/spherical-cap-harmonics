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
function [qm_k, C_mat] = SCH_analysis(K_max, theta_c, eigen_table, thetas, phis, vertices, N_eps)
% This function to expand the coefficients from the provided surface
% patches.
% K_max: the maximum of index degree for the expnasion.
% All angles' inputs are in radians only.
C_mat = zeros(size(thetas,1), (K_max+1)^2);
fprintf("\n\n Compiling the basis functions...\n")
% Combine basis fucntions
for n = 0:K_max
    C_mat(:, n^2 + 1:n^2 + 2*n + 1) = spherical_cap_harmonic_basis(n, theta_c, eigen_table, thetas, phis, N_eps);
    fprintf("\nCompleted: %3.2f %% for n = %3.0f", ((n+1)^2/(K_max+1)^2*100), n)
end

% Least Square Fitting (LSF) Method
fprintf("\n\n Computing the Least Square Fit...\n")
qm_k = (C_mat' * C_mat) \ C_mat' * vertices;
end