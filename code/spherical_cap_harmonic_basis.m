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
function SCH_basis = spherical_cap_harmonic_basis(K, theta_c, eigen_table, thetas, phis, N_eps)
% This function returns the complete basis that includes Legendre and
% Fourier for \theta and \phi.
% It returns a table for positive and negative orders for a certain index K
% The number of the columns in the matrix is 2K+1 for (-m:1:m) orders
% Note: thetas and phis matrix msut be identical in size
% All the angles here must be in radians.
% For our applications the following variables shall be fixed:
% thetas and phis must be square matrix or column matrix (row matrix will give you an error!)
BC = "even";
normalization_method = "Hainse";

% Initializing and filling up
ms = 0:K;

% Check inputs' size (1D and 2D matrix)
sz = length(size(thetas)); % For reconstruction matricies only
if isempty(find(size(thetas) == 1, 1))
    sz = sz + 1;
end
if sz == 2 % For normal analysis or reconstruction using geodesic domes
    SCH_basis = zeros(length(thetas), (2*K+1));
    positive_ms = repmat(ms, length(thetas), 1);
    % Compute and fill up the postive orders
    SCH_basis(:, K+1:end) = fractional_associated_Legendre_functions(theta_c,...
        K, BC, eigen_table, thetas, normalization_method, N_eps) .* exp(1i .* positive_ms .* repmat(phis, 1, K+1));
    % Compute and fill up the negative orders
    SCH_basis(:, 1:K) = flip(conj(SCH_basis(:, K+2:end)) .* (-1).^ positive_ms(:, 2:end), 2);
elseif sz == 3 % For the meshgrid reconstruction case only or 3D render of the basis
    temp_sz = size(thetas); temp_sz(3) = (2*K+1);
    SCH_basis = zeros(temp_sz);
    % Prepare the ms
    positive_ms = zeros(temp_sz(1), temp_sz(2), K+1);
    for ii = ms
        positive_ms(:, :, ii+1) = ones(temp_sz(1), temp_sz(2)) .* ii;
    end
    % Compute and fill up the postive orders
    SCH_basis(:, :, K+1:end) = fractional_associated_Legendre_functions(deg2rad(theta_c),...
        K, BC, eigen_table, thetas, normalization_method, N_eps) .* exp(1i .* positive_ms .* phis);
    % Compute and fill up the negative orders
    SCH_basis(:, :, 1:K) = flip(conj(SCH_basis(:, :, K+2:end)) .* (-1).^ positive_ms(:, :, 2:end), 3);
end
end