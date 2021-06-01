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
% This is the main script interface for the Hemispherical Harmonic Analysis
%% Paths
% function Main_HSH
% close all;
clear; clc; tic
addpath('code') % for the hemispherical parmaeterization
addpath('stlTools') % for reading the STL files
addpath('input_geom')
addpath('distmesh')
%% Inputs
fprintf('The analysis started at: %s\n', datestr(now,'HH:MM:SS.FFF'))

% To appply spherical cap parametrization step (true)
parametrization = true;

surface_name = 'half_rock_mesh_manifold.stl';

para_surface_path = ''; % if the parametrization already computed
solve_roots = false; % To solve for the roots if not pre-calculated!

max_degree = 45;
max_reconstruction_degree = 45;
max_reconstruction_degree = min(max_reconstruction_degree, max_degree);

% Math grid reconstruction
reconstruction_resolution = 70;
math_grid_reconstruction = false;

% option 1 gives you reconstruction by edge length and 2 gives you
% reconstruction by number of even elements used.
reconstruction_mesh_option = 1;
% The desired length of each edge on a unit disk from (0.035 < edge_length < 1.0) 
edge_length = 0.035;
% The resolution or the number of vertices must be an even number!
resolution = 30*2;

plot_figures = true;
%% Spherical Cap Parametrization

% Load the STL file
[v,f] = stlRead(surface_name);

% Normalize the vertices
v = [v(:,1) - mean(v(:,1)), v(:,2) - mean(v(:,2)), v(:,3) - mean(v(:,3))];

if parametrization
    % ensure there is no triangle with all three vertices lying on the boundary
    % (such a triangle will be highly squeezed under the parameterization as
    % all three vertices will be constrained)
    v_old = v;
    f_old = f;
    [v,f] = clean_mesh(v_old,f_old,1);
    while length(v) ~= length(v_old)
        v_old = v;
        f_old = f;
        [v,f] = clean_mesh(v_old,f_old,1);
    end
    height_map = v(:, 3);
    map = spherical_cap_conformal_map(v, f, cos(pi/2));
    
    if plot_figures
        paraview_patch(v, f, height_map)
        view([50 20])
        title('Input surface')
    
        paraview_patch(map, f, height_map)
        view([50 20])
        title('Spherical cap conformal parameterization')
        axis equal tight on
    end
    % Save parametrised mesh
    stlWrite(['rec_output/para_output_HSH_', surface_name], f, map,'mode','ascii')
else
    [map,f] = stlRead(para_surface_path);
end
fprintf("\n\n Verts number: %3.0f", length(v(:, 1)))
%% The analysis step (HSH)
thetas = acos(map(:,3)); % [0, theta_c]
phis = atan2(map(:,2),map(:,1)); % [-pi, pi]
C_mat = HSH_basis_new(phis, thetas, max_degree); % hemispherical harmonics basis
qm_k = (C_mat'*C_mat)\C_mat'*v;
%% Shape descriptors and the fractal dimension
% Shape descriptors (2-norm) for frequency accumulates at a certain
% frequency degree..
Dl = zeros([3, max_degree+1]);
for  k = 1:3
    for l = 1:max_degree
        for m = -l:1:l
            Dl(k, l) = Dl(k, l) + (real(qm_k(l^2 + l + m + 1, k)))^2 + (imag(qm_k(l^2 + l + m + 1, k)))^2;
        end
        Dl(k, l) = sqrt(Dl(k, l)); % Export only the last object stone.stl
    end
    Dl(k, :) = Dl(k, :)/Dl(k, 1);
end
if plot_figures
   figure
   subplot(2,3,1)
   x_temp = 1:length(Dl(1, :));
   loglog(x_temp(2:end), Dl(1, 2:end), 'LineWidth', 2)
   title('The shape descriptors')
   xlabel('Freqency index (k)')
   ylabel('Normalized amplitude (Dx)')
   grid on
   
   subplot(2,3,2)
   loglog(x_temp(2:end), Dl(2, 2:end), 'LineWidth', 2)
   title('The shape descriptors')
   xlabel('Freqency index (k)')
   ylabel('Normalized amplitude (Dy)')
   grid on
   
   subplot(2,3,3)
   loglog(x_temp(2:end), Dl(3, 2:end), 'LineWidth', 2)
   title('The shape descriptors')
   xlabel('Freqency index (k)')
   ylabel('Normalized amplitude (Dz)')
   grid on
   
   subplot(2,3,4:6)
   loglog(x_temp(2:end), sqrt(Dl(1, 2:end).^2 + ...
       Dl(2, 2:end).^2 + Dl(3, 2:end).^2), 'LineWidth', 2)
   title('The shape descriptors')
   xlabel('Freqency index (k)')
   ylabel('Normalized amplitude (Dr)')
   grid on
end
%% The reconstruction step
fprintf("\n\n The reconstruction part has started...")
% Define colormap
red_color = zeros(1, 255); green_color = zeros(1, 255); blue_color = zeros(1, 255);
red_color(1:127) = linspace(42,220,127); red_color(128:255) = linspace(220,174,127+1);
green_color(1:127) = linspace(63,220,127); green_color(128:255) = linspace(220,0,127+1);
blue_color(1:127) = linspace(181,220,127); blue_color(128:255) = linspace(220,22,127+1);
ParaviewMap = [red_color', green_color', blue_color']./255;
if math_grid_reconstruction
    sample_phi_range = linspace(-pi,pi,reconstruction_resolution);
    sample_theta_range = linspace(0,pi/2-0.1, reconstruction_resolution);
    [sample_theta, sample_phi] = meshgrid(sample_theta_range, sample_phi_range);
    Z_reconstructed = HSH_basis_new(sample_phi(:), sample_theta(:), max_reconstruction_degree); % hemispherical harmonics
    reconstructed = real(Z_reconstructed * qm_k(1:(max_reconstruction_degree+1)^2, :));
    
    x_reconstructed = reshape(reconstructed(:,1), length(sample_phi_range), length(sample_theta_range));
    y_reconstructed = reshape(reconstructed(:,2), length(sample_phi_range), length(sample_theta_range));
    z_reconstructed = reshape(reconstructed(:,3), length(sample_phi_range), length(sample_theta_range));
    figure;
    surf(x_reconstructed, y_reconstructed, z_reconstructed);
    axis equal tight off
    title(['Reconstructed surface', ' with l = ', num2str(max_reconstruction_degree)]);
    view([50 20])
    colormap(ParaviewMap);
    s.CData = z_reconstructed;
    colorbar
else
    [v_rec, f_rec] = uniform_spherical_cap_grid(pi/2-0.02, edge_length, resolution,...
        reconstruction_mesh_option);
    thetas_rec = acos(v_rec(:,3));
    phis_rec = atan2(v_rec(:,2),v_rec(:,1));
    Z_reconstructed = HSH_basis_new(phis_rec, thetas_rec, max_reconstruction_degree);
    v_rec = real(Z_reconstructed * qm_k(1:(max_reconstruction_degree+1)^2, :));
    % This to export the final reconstruction degrees in one files
    stlWrite(['rec_output/rec_k_', num2str(max_reconstruction_degree),...
        '_HSH_', surface_name]...
        , f_rec, v_rec,'mode','ascii')
    if plot_figures
        paraview_patch(v_rec, f_rec, v_rec(:, 3))
        colormap(ParaviewMap);
    end
end
surface_name_wo_format = split(surface_name,'.');
surface_name_wo_format = surface_name_wo_format{1};
save(['rec_output/rec_k_', num2str(max_degree), '_HSH_',...
    surface_name_wo_format, '.mat'])
fprintf("\n\nFinished in %f min.\n", toc/60)
fprintf("\n\n Verts number: %3.0f", length(v(:, 1)))
fprintf("\n\n Equations number: %3.0f", length(v(:, 1)) * (max_degree+1)^2)
fprintf("\n\n (K+1)^2: %3.0f\n", (max_degree+1)^2)
fprintf("\n\n Total no. of coefficients 3(K+1)^2: %3.0f\n", 3*(max_degree+1)^2)
fprintf('The analysis ended at: %s\n', datestr(now,'HH:MM:SS.FFF'))
% end