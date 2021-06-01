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
% function Main_optimum_half_angle_finder
% This function to solve for the optimum half-angle (\theta_c) by
% minimizing area and angluar distortion
%% Load paths
format long
clc; close all; clear; tic
addpath('code')
addpath('stlTools')
addpath('input_geom')
%% The PSS inputs
silent_mode = false;
PSS_plot_results = true;
LB = deg2rad(7); % Lower bound in radians for \theta_c
UB = deg2rad(120);% Upper bound in radians for \theta_c
pop_size = 20;
iterations = 5;
acceptance_rate = 0.97;
%% Input surface
surface_name = 'half_rock_mesh_manifold.stl';  % Best answer 100 degrees...
surface_name = "test_rough_surface.stl"; % Best answer 8.6 degrees...

plot_figures = true;

% Load the STL file
[v,f] = stlRead(surface_name);
% Clean the STL mesh (if needed)
v_old = v;
f_old = f;
[v,f] = clean_mesh(v_old,f_old,1);
while length(v) ~= length(v_old)
    v_old = v;
    f_old = f;
    [v,f] = clean_mesh(v_old,f_old,1);
end
clear v_old f_old
% Normalize the vertices
v = [v(:,1) - mean(v(:,1)), v(:,2) - mean(v(:,2)), v(:,3) - mean(v(:,3))];
if plot_figures
    height_map = v(:, 3);
    height_map = sqrt(v(:, 3).^2 + v(:, 2).^2 + v(:, 1).^2);
    paraview_patch(v, f, height_map)
    view([50 20])
    ff = gcf;
    title('Input surface')
end
%% The core of the solver
objective_fun = @(theta_c) optimum_half_angle_objective(v, f, theta_c);

results_struct = PSS_optimizer(objective_fun, LB, UB, pop_size, iterations, ...
    acceptance_rate, PSS_plot_results, silent_mode);

% Save parametrised mesh
best_theta_c = results_struct.best_solution;
optimum_eval = results_struct.evaluations_history(end);
map = spherical_cap_conformal_map(v, f, cos(best_theta_c), false);
stlWrite(strcat('optimum_half_angle/para_th_', num2str(rad2deg(best_theta_c)), '_eval_',...
    num2str(optimum_eval), '_', surface_name), f, map,'mode','ascii')
if plot_figures
    height_map = sqrt(v(:, 3).^2 + v(:, 2).^2 + v(:, 1).^2);
    paraview_patch(map, f, height_map)
    view([50 20])
    title(['Output surface, \theta_{c} = ', num2str(rad2deg(best_theta_c))])
    ff = gcf;
end
disp(['Best half-angle obtained is: ', num2str(rad2deg(best_theta_c))])
toc
% end
%% The objective function definition
function optimum_eval = optimum_half_angle_objective(v, f, theta_c)
[~, optimum_eval] = spherical_cap_conformal_map(v, f, cos(theta_c), true);
end