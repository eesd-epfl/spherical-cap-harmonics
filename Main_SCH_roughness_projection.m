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
% This file is for the roughness projection on a spherical cap (including a
% hemisphere).
%% Paths
% function Main_SCH_roughness_projection
close all;
clear; clc; tic
addpath('code')
addpath('stlTools')
addpath('geodesic_domes')
addpath('input_geom/donor_meshes/')
addpath('analysis_results')
%% Reconstruction range for k \in [start, end]
starting_degree = 15;
ending_degree = 40;
N_eps = 1000;
scale_factor = 0.1;
touch_one_direction = true; % project vertices along one direction.
touched_direction = 3; % (1) -> x; (2) -> y; (3) -> z

plot_figures = true;
%% Load the analysis reults from HSH or SCH only
analysis_mat_file = 'rec_k_40_theta_10_rocks_patch_rough.mat';
load(analysis_mat_file);
%% The donor mesh is the mesh that provides us with the shape and we project
% roughness on it extracted from the SCH or HSH analyses.

donor_mesh = true; % if false it will load a unit sphere as a donor mesh
icosahedron_dome_refinement = 4;

% donor_mesh_name = 'grad_ref_cube.stl';
% donor_mesh_name = 'refined_sphere.stl';
% donor_mesh_name = 'cube_refined_sample_2.stl';
donor_mesh_name = 'cube_refined_sample_3.stl';

if donor_mesh
    [v, f] = stlRead(donor_mesh_name);
else
    [v, f] = icosphere(icosahedron_dome_refinement);
end
% v = [v(:,1) - mean(v(:,1)), v(:,2) - mean(v(:,2)), v(:,3) - mean(v(:,3))];
%% Specify (select) the group of vertices to be parametrized
% filter_query = strcat('selected_verts = find(v(v(:,1)>(0.8-0.0001)));'); % for a sphere
filter_query = strcat('selected_verts = find(v(:,3)>0.577-0.001);'); % for the donor mesh (cube)
% filter_query = strcat('selected_verts = find(v(:,3)>0.63-0.001);'); % for the donor mesh (refined sphere)

eval(filter_query);

selected_groups = zeros(length(v(:, 1)), 1);
selected_groups(selected_verts) = 1;
fprintf('\nNumber of vertices chosen: %f \n', length(selected_groups(selected_verts)))
%% Project roughness and render results
if plot_figures
    paraview_patch(v, f, selected_groups)
    title('Input donor mesh')
%     axis on
%     xlabel('x'); ylabel('y'); zlabel('z')
    view([50 20])
    colorbar off
end

if isempty(selected_verts)
    error('The selected vetices list must be specified.');
end

% Solve for \phi and \theta from the donor mesh patch
map = spherical_map_with_prescribed_cap(v, f, selected_verts, cos(theta_c));

if plot_figures
    paraview_patch(map, f, selected_groups)
    title('Input donor mesh parametrisation on a prescribed cap')
end

selected_map = map(selected_verts, :);
thetas = acos(selected_map(:,3)); % [0, theta_c]
phis = atan2(selected_map(:,2),selected_map(:,1)); % [-pi, pi]

k_range = [starting_degree, ending_degree];
recursive_reconstruction = scale_factor .* real(SCH_range_reconstruction_icosahedron_patch(...
    k_range, qm_k, theta_c, eigen_table, thetas, phis, N_eps));

if touch_one_direction
    v(selected_verts, touched_direction) = v(selected_verts, touched_direction)...
        + recursive_reconstruction(:, touched_direction, end);
else
    v(selected_verts, :) = v(selected_verts, :) + recursive_reconstruction(:, end);
end

if plot_figures
    paraview_patch(v, f, selected_groups)
    title('Output of the roughened mesh')
    colorbar off
end
%% Save the resulted STL file
stlWrite(['roughness_output/rec_l_', num2str(starting_degree), '_to_',...
        num2str(ending_degree), '_SH_', surface_name]...
        , f, v,'mode','ascii')

% end