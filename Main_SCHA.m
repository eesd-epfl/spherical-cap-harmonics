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
% This is the main script interface for the Spherical Cap Harmonic Analysis
%% Paths
% function Main_SCHA
% close all;
clear; clc; tic
addpath('code')
addpath('stlTools')
addpath('input_geom')
addpath('distmesh')
%% Inputs
fprintf('The analysis started at: %s\n', datestr(now,'HH:MM:SS.FFF'))

% To appply spherical cap parametrization step (true)
parametrization = true;

% surface_name = 'rocks_patch_rough.stl'; % This one should be downloaded
% from Zenodo
surface_name = 'test_rough_surface.stl'; % This is a simple and dummy surface (a fast example)

para_surface_path = ''; % if the parametrization already computed
solve_roots = false; % To solve for the roots if not pre-calculated!

 % The table is included already no need to calculate the roots
 % If the user wants to use different \theta_c or roots for K>39+1
 % uncomment the next line and comment the table with path so it could
 % recompute the S-L roots AGAIN. Be careful with the computed roots!(check them out befor usage)
% eigen_table_path = '';
eigen_table_path = "eigenvalues_k_40_even_theta_10.mat";
% eigen_table_path = "eigenvalues_k_35_even_theta_50.mat"; %(just an example)

theta_c = deg2rad(10); % The half-angle of the spherical cap (must be in radians)

% To try this method real quick, firts use max_degree = 12;
% For k <= 12 the hypergeometric function is computed really fast.
max_degree = 39 + 1; % +1 to account for the 0-degree.
max_reconstruction_degree = 39 + 1;
truncation_degree_fit = 40; % For fitting the fractal dimension
truncation_degree_fit = min(truncation_degree_fit, max_degree);
max_reconstruction_degree = min(max_reconstruction_degree, max_degree);

% In case you wanted to compute the roots again from here.
max_iterations_steps = 50;
solver_steps = 0.5;
solver_increment_diag = 0.1;
solver_increment_tri = 0.5;

% Math grid reconstruction
reconstruction_resolution = 50;
math_grid_reconstruction = false;

% Geodesic dome STL reconstruction
recursive_reconstruction = true;

% option 1 gives you reconstruction by edge length and 2 gives you
% reconstruction by number of even elements used.
reconstruction_mesh_option = 1;
% The desired length of each edge on a unit disk from (0.035 < edge_length < 1.0) 
edge_length = 0.032;
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
    if plot_figures
        height_map = v(:, 3);
        paraview_patch(v, f, height_map)
        view([50 20])
        title('Input surface')
    end
    map = spherical_cap_conformal_map(v, f, cos(theta_c));
    if plot_figures
        paraview_patch(map, f, height_map)
        view([50 20])
        title('Spherical cap conformal parameterization')
        axis equal tight on
    end
    % Save parametrised mesh
    stlWrite(['rec_output/para_output_', surface_name], f, map,'mode','ascii')
else
    [map,f] = stlRead(para_surface_path);
end
fprintf("\n\n Verts number: %3.0f", length(v(:, 1)))
%% Finding the Sturm–Liouville eigenvalues for the spherical cap basis.
% This step should be ideally run from the "Sturm_Liouville_eigenvalues.m"
% file and not here so it can be adjusted and calibrated.
if solve_roots 
    eigen_table = Sturm_Liouville_eigenvalues(cos(theta_c), max_degree...
        , "even", solver_increment_diag, solver_increment_tri, max_iterations_steps)
    % Add saving method for the table here (Mahmoud)
else
    % Load the table
    eigen_table = load(eigen_table_path);
    eigen_table = eigen_table.eigen_table;
end
%% The analysis step (SCHA)
thetas = acos(map(:,3)); % [0, theta_c]
phis = atan2(map(:,2),map(:,1)); % [-pi, pi]
% profile on
N_eps = 5000;
[qm_k, C_mat] = SCH_analysis(max_degree, theta_c, eigen_table, thetas, phis, v, N_eps);
% prof = profile('info');
% save prof prof
%% Shape descriptors and the fractal dimension
% Shape descriptors (2-norm) for frequency accumulates at a certain
% frequency degree..
Dl = zeros([3, max_degree+1]);
for  k = 1:3
    for l = 1:max_degree
        for m = -l:1:l
            Dl(k, l) = Dl(k, l) + (real(qm_k(l^2 + l + m + 1, k)))^2 + (imag(qm_k(l^2 + l + m + 1, k)))^2;
        end
        Dl(k, l) = sqrt(Dl(k, l));
    end
    Dl(k, :) = Dl(k, :)/Dl(k, 1);
end
if plot_figures
   figure
   subplot(3,3,1)
   x_temp = 1:length(Dl(1, :));
   loglog(x_temp(2:end), Dl(1, 2:end), 'LineWidth', 2)
   title('The shape descriptors')
   xlabel('Freqency index (k)')
   ylabel('Normalized amplitude (Dx)')
   grid on
   
   subplot(3,3,2)
   loglog(x_temp(2:end), Dl(2, 2:end), 'LineWidth', 2)
   title('The shape descriptors')
   xlabel('Freqency index (k)')
   ylabel('Normalized amplitude (Dy)')
   grid on
   
   subplot(3,3,3)
   loglog(x_temp(2:end), Dl(3, 2:end), 'LineWidth', 2)
   title('The shape descriptors')
   xlabel('Freqency index (k)')
   ylabel('Normalized amplitude (Dz)')
   grid on
   
   subplot(3,3,4:6)
   loglog(x_temp(2:end), sqrt(Dl(1, 2:end).^2 + ...
       Dl(2, 2:end).^2 + Dl(3, 2:end).^2), 'LineWidth', 2)
   title(['The shape descriptors, \theta_{c} = ', num2str(rad2deg(theta_c))])
   xlabel('Freqency index (k)')
   ylabel('Normalized amplitude (Dr)')
   grid on
   
   subplot(3,3,7:9)
   loglog(x_temp(2:truncation_degree_fit), (Dl(1, 2:truncation_degree_fit).^2 + ...
       Dl(2, 2:truncation_degree_fit).^2 + Dl(3, 2:truncation_degree_fit).^2), 'LineWidth', 2)
   title(['Surface Energy, \theta_{c} = ', num2str(rad2deg(theta_c))])
   xlabel('Freqency index (k)')
   ylabel('Energy of amplitude (Dr^{2})')
   grid on
   hold on
   % This fitting is not accurate use Excel instead!!!
   y_fit = (Dl(1, 2:truncation_degree_fit).^2 + ...
       Dl(2, 2:truncation_degree_fit).^2 + Dl(3, 2:truncation_degree_fit).^2)';
   fitting = fit(x_temp(2:truncation_degree_fit)',y_fit,'b*x^m', ...
       'StartPoint', [1, -1]);
   plot(fitting,'k')
   hold off
end
surface_name_wo_format = split(surface_name,'.');
surface_name_wo_format = surface_name_wo_format{1};
save(['rec_output/analysis_k_', num2str(max_degree), '_theta_', ...
    num2str(round(rad2deg(theta_c),0)), '_', surface_name_wo_format,...
    '.mat'], 'Dl')

%% Estimate the wavelengths by FDEC
FDEC = qm_k(2:4, :);
A = [(FDEC(:, 1) - FDEC(:, 3)), 1i .* (FDEC(:, 1) + FDEC(:, 3)), sqrt(2) .* FDEC(:, 2)];
[vv,dd] = eig(A*A');
d = sqrt(diag(abs(dd)));
omega = (2*pi/k) * sqrt( ((0.5*d(2))^2 + (0.5*d(3))^2)/2);
fprintf("\n\n The dimensions of FDEC: a = %3.3f, b = %3.3f, c = %3.3f.\n",...
    d(3), d(2), d(1))
if plot_figures
   figure
   x_omega = 1:(max_degree-2);
   y_omega = (2*pi) * sqrt(((0.5*d(2))^2 + (0.5*d(3))^2)/2) ./ x_omega;
   loglog(x_omega, y_omega, 'LineWidth', 2)
   xlabel('Freqency index (k)')
   ylabel('Wavelength (\omega_{k})')
   title('The corresponding wavelengths')
end
fprintf("\n\n The largest wavelength = %3.4f, and the smallest = %3.4f.\n",...
    max(y_omega), min(y_omega))
%% The reconstruction step
fprintf("\n\n The reconstruction part has started...")
if math_grid_reconstruction
    reconstructed = real(SCH_reconstruction_math(max_reconstruction_degree, qm_k, ...
        reconstruction_resolution, theta_c, eigen_table, N_eps));
    figure;
    surf(reconstructed(:, :, 1), reconstructed(:, :, 2), reconstructed(:, :, 3));
    axis equal tight off
    title(['Reconstructed surface', ' with K = ', num2str(max_reconstruction_degree)]);
    view([50 20])
    % Define colormap
    red_color = zeros(1, 255); green_color = zeros(1, 255); blue_color = zeros(1, 255);
    red_color(1:127) = linspace(42,220,127); red_color(128:255) = linspace(220,174,127+1);
    green_color(1:127) = linspace(63,220,127); green_color(128:255) = linspace(220,0,127+1);
    blue_color(1:127) = linspace(181,220,127); blue_color(128:255) = linspace(220,22,127+1);
    ParaviewMap = [red_color', green_color', blue_color']./255;
    colormap(ParaviewMap);
    s.CData = reconstructed(:, :, 3);
    colorbar
else
    [v_rec, f_rec] = uniform_spherical_cap_grid(theta_c, edge_length, resolution,...
        reconstruction_mesh_option);
    thetas_rec = acos(v_rec(:,3));
    phis_rec = atan2(v_rec(:,2),v_rec(:,1));
    v_rec = real(SCH_reconstruction_icosahedron_dome(max_reconstruction_degree, qm_k,...
        theta_c, eigen_table, thetas_rec, phis_rec, N_eps));
    if recursive_reconstruction
        % This to export all the reconstruction degrees in different files
        for kk = 0:max_reconstruction_degree-1
            stlWrite(['rec_output/rec_k_', num2str(kk),...
                '_theta_', num2str(round(rad2deg(theta_c),0)), '_', surface_name]...
                , f_rec, v_rec(:, :, kk+1),'mode','ascii')
            fprintf("\nFile: (%s) has been saved.", ['rec_k_', num2str(kk),...
                '_theta_', num2str(round(rad2deg(theta_c),0)), '_', surface_name]);
        end
    else
        % This to export the final reconstruction degrees in one files
        stlWrite(['rec_output/rec_k_', num2str(max_reconstruction_degree),...
            '_theta_', num2str(round(rad2deg(theta_c),0)), '_', surface_name]...
            , f_rec, v_rec(:, :, end),'mode','ascii')
    end
    rec_height_map = v_rec(:, 3, end);
    if plot_figures
        paraview_patch(v_rec(:, :, end), f_rec, rec_height_map)
        view([50 20])
    end
end
save(['rec_output/rec_k_', num2str(max_degree), '_theta_', ...
    num2str(round(rad2deg(theta_c),0)), '_', surface_name_wo_format,...
    '.mat'])
fprintf("\n\nFinished in %f min.\n", toc/60)
fprintf("\n\n Verts number: %3.0f", length(v(:, 1)))
fprintf("\n\n Equations number: %3.0f", length(v(:, 1)) * (max_degree+1)^2)
fprintf("\n\n (K+1)^2: %3.0f\n", (max_degree+1)^2)
fprintf("\n\n Total no. of coefficients 3(K+1)^2: %3.0f\n", 3*(max_degree+1)^2)
fprintf('The analysis ended at: %s\n', datestr(now,'HH:MM:SS.FFF'))
fprintf("\n The Hurst exponent: %3.3f, and the fractal dimesnion is: %3.3f (not accurate fit by MATLAB!)\n",...
    -0.5*fitting.m, (6 + fitting.m)*0.5)
% end