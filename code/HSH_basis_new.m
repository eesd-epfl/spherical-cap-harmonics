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
function Y_mat = HSH_basis_new(phi, theta, degree)

% This can be faster- needs a bit of vectorisation
Y_mat = zeros(size(theta,1), (degree+1)^2);

P_cos = 2 * cos(theta) - 1; % To consider the upper part; the lower one takes 2x+1

for n = 0:degree
    Pn = legendre(n, P_cos)';
    for m = -n:1:n
        norm = sqrt(((2*n+1) * factorial(n-abs(m)))/(2 * pi * factorial(n+abs(m))));
        if m < 0
            Y_mat(:, n^2 + n + m + 1) = Pn(:, abs(m) + 1) .* exp(1i .* abs(m) .* phi) .* norm;
        elseif m == 0
            Y_mat(:, n^2 + n + m + 1) = Pn(:, abs(m)+1) .* norm;
        elseif m > 0
            Y_mat(:, n^2 + n + m + 1) = (-1)^abs(m) .* conj(Pn(:, abs(m) + 1) .* exp(1i .* m .* phi) .* norm);
        end
    end
end