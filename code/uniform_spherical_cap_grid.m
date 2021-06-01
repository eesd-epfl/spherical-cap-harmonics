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
% Authors of this file: Mahmoud S. Shaqfa and Gary Choi

function [v,f] = uniform_spherical_cap_grid(theta_c, edge_length, resolution, reconstruction_mesh_option)
% Authors Mahmoud S. Shaqfa, Gary Choi, Katrin Beyer
% This function to generate a uniform grid on a spherical cap or a
% We used the DISTMESH Library: http://persson.berkeley.edu/distmesh/
% for generating a uniform mesh on a unit disk.
% prescribed half-angle \theta_c.
% option 1 gives you reconstruction by edge length and 2 gives you
% reconstruction by number of even elements used.
% edge_length: the desired length of each edge on a unit disk
% resolution
%% Create a uniform spherical cap grid
if reconstruction_mesh_option == 1
    % Planar uniform mesh generation on a unit disk
    fd=@(p) sqrt(sum(p.^2,2))-1;
    [p,f]=distmesh2d(fd,@huniform,edge_length,[-1,-1;1,1],[]);
else
    samples = linspace(-1, 1, resolution);
    [sample_phi_1, sample_phi_2] = meshgrid(samples, samples);
    a = sample_phi_1; b = sample_phi_2;
    A = zeros(size(a)); B = zeros(size(a));
    
    indx = find(and(abs(a) >= abs(b), and(b ~= 0, a~=0)));
    A(indx) = 2.* a(indx) ./sqrt(pi) .* cos(b(indx).*pi./(4.*a(indx)));
    B(indx) = 2.* a(indx) ./sqrt(pi) .* sin(b(indx).*pi./(4.*a(indx)));
    
    indx2 = find(and(abs(a) < abs(b), and(b ~= 0, a~=0)));
    A(indx2) = 2.* b(indx2) ./sqrt(pi) .* sin(a(indx2).*pi./(4.*b(indx2)));
    B(indx2) = 2.* b(indx2) ./sqrt(pi) .* cos(a(indx2).*pi./(4.*b(indx2)));
    
    stlWrite('tmp_dome.stl', A, B, zeros(size(A)))
    [p, f] = stlRead('tmp_dome.stl');
    delete('tmp_dome.stl');
    
    p = p./max(abs(p(:, 1))); % Unit disk
end
% Rescale the disk based on the area of the desired spherical cap domain
% Note that surface area of spherical cap with height h = 2*pi*h, and hence
% the desired area is 2*pi*(1-z_lim)
p = p*sqrt(2*(1-cos(theta_c))); % Scale it for a spherical cap
% Lambert equi-area projection
v = -[sqrt(1-(p(:,1).^2+p(:,2).^2)/4).*p(:,1), ...
     sqrt(1-(p(:,1).^2+p(:,2).^2)/4).*p(:,2), ...
     -1+(p(:,1).^2+p(:,2).^2)/2];
end