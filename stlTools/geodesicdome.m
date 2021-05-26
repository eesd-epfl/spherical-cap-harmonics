function [v, f] = geodesicdome(N, theta_c)
% Authors: Mahmoud Shaqfa, Gary Choi, Katrin Beyer
% Generate a Geodesic Dome for the spherical caps
% theta_c in radians
% N: the number of refinements

if nargin == 0
    close all
    % Test codes
    N = 3;
    theta_c = deg2rad(90);
end

[v, f] = icosphere(N);
v = [v(:,1) - mean(v(:,1)), v(:,2) - mean(v(:,2)), v(:,3) - mean(v(:,3))];

% Find unneeded vertices
x_c = cos(theta_c);
list = find(v(v(:, 3) < x_c));

% Delete unneeded vertices
[v, f] = stlDelVerts2(v, f, list');

% Rotate about the z-axis
theta=pi/2;
R = [cos(theta), 0, sin(theta)
    0 1 0
    -sin(theta), 0, cos(theta)];
v = v*R;

if nargin == 0
    % Test codes
    cd ..
    addpath('/Spherical Cap Harmonic Analysis')
    addpath('/stlTools')
    paraview_patch(v, f)
    axis on
    xlabel('x'); ylabel('y'); zlabel('z')
    cd stlTools/
end
end