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
% Author of this file: Gary Choi
function [map, optimum_eval] = spherical_cap_conformal_map(v,f,z_lim, silent_mode)

% Compute a spherical cap conformal mapping using the disk conformal 
% parameterization method in [Choi and Lui, J. Sci. Comput. 2015] together 
% with an optimal Mobius transformation for further reducing the area 
% distortion.
%
% Input:
% v: nv x 3 vertex coordinates of a simply-connected open triangle mesh
% f: nf x 3 triangulations of a simply-connected open triangle mesh
% z_lim: the z-limit of the target spherical cap (all points will be with z >= z_lim), must be in [-1,1)
% 
% Output:
% map: nv x 3 vertex coordinates of the spherical cap conformal parameterization
% 
% Remark:
% 1. Please make sure that the input mesh does not contain any 
%    unreferenced vertices/non-manifold vertices/non-manifold edges.
% 2. Please remove all valence 1 boundary vertices (i.e. vertices with 
%    only 1 face attached to them) before running the program.
%% Initialization
if nargin == 3
   silent_mode = false; 
end

% Parameter
parameter.north = 5;
parameter.south = 100;
parameter.threshold = 0.00001;

nv = length(v);
TR = TriRep(f,v);
B = freeBoundary(TR);
bdy_index = B(:,1);
% Remark: The above approach for getting the surface boundary may not work
% well in case the boundary contains vertices with valence 1
% In that case, use some other method to obtain the boundary first

bdy_length = sqrt((v(bdy_index,1) - v(bdy_index([2:end,1]),1)).^2 + ...
            (v(bdy_index,2) - v(bdy_index([2:end,1]),2)).^2 + ...
            (v(bdy_index,3) - v(bdy_index([2:end,1]),3)).^2);
partial_edge_sum = zeros(length(bdy_length),1);

% Arclength parameterization boundary constraint for the initial map
for i = 2:length(bdy_length)
    for j = 1:i-1
    partial_edge_sum(i) = partial_edge_sum(i) + bdy_length(j);
    end
end
theta = 2*pi.*partial_edge_sum/sum(bdy_length)+2*pi*0.01;
bdy = exp(theta*1i)';

% Disk harmonic map
M = cotangent_laplacian(v,f);

[mrow,mcol,mval] = find(M(bdy_index,:));
M = M - sparse(bdy_index(mrow),mcol,mval,nv, nv) + ...
        sparse(bdy_index,bdy_index,ones(length(bdy_index),1),nv, nv);
c = zeros(nv,1); c(bdy_index) = bdy;
z = M \ c;
map_disk = [real(z),imag(z)]; map_prev = map_disk;

if sum(sum(isnan(map_disk))) ~= 0
    % use tutte embedding instead
    map_disk = tutte_map(v,f,bdy_index,bdy); map_prev = map_disk;
    z = map_disk(:,1) + 1i*map_disk(:,2);
end

if ~silent_mode
    fprintf('Initialization completed.\n')
end
%% North Pole iteration
% Use the Cayley transform to map the disk to the upper half plane
% All boundary points will be mapped to the real line

mu = beltrami_coefficient(map_disk, f, v); 
mu_v = f2v(v,f)*mu;
bdy_index_temp = [bdy_index(2:end);bdy_index(1)];
[~, least] = min(abs(mu_v(bdy_index))+abs(mu_v(bdy_index_temp)));
z = z* exp(-1i*(angle(z(bdy_index(least)))+angle(z(bdy_index(mod(least,length(bdy_index))+1)))/2));
g = 1i*(1 + z)./(1 - z);

% fix the points near the puncture, i.e. near z = 1
[~, ind] = sort(-real(z));
fixed = setdiff(ind(1:max(round(length(v)/parameter.north),min(100,length(z)))), bdy_index);
fixed = [fixed; find(real(g) == max(real(g))); find(real(g) == min(real(g)))];

P = [real(g),imag(g),ones(length(g),1)];
mu = beltrami_coefficient(P, f, v); 

% compute the updated x coordinates
target = P(fixed,1);
A = generalized_laplacian(P,f,mu); Ax = A; Ay = A;
b = -Ax(:,fixed)*target;
b(fixed) = target;
Ax(fixed,:) = 0; Ax(:,fixed) = 0;
Ax = Ax + sparse(fixed,fixed,ones(length(fixed),1), size(A,1), size(A,2));
x = Ax\b;

% compute the updated y coordinates
target = P(fixed,2); 
fixed = [fixed;bdy_index];
target = [target;zeros(length(bdy_index),1)];
b = -Ay(:,fixed)*target;
b(fixed) = target;
Ay(fixed,:) = 0; Ay(:,fixed) = 0;
Ay = Ay + sparse(fixed,fixed,ones(length(fixed),1), size(Ay,1), size(A,2));
y = Ay\b;

g_new = complex(x,y);
z_new = (g_new - 1i)./(g_new + 1i);
map_disk = [real(z_new), imag(z_new)];

if sum(sum(isnan(map_disk))) ~= 0 
    % use the old result in case of getting NaN entries
    map_disk = map_prev;
    z_new = map_disk(:,1) + 1i*map_disk(:,2);
end
if ~silent_mode
    fprintf('North pole step completed.\n')
end
%% Reflection along the unit circle

f_temp = f + length(v);
a = sort(bdy_index + length(v));
for i = length(a):-1:1
    f_temp(f_temp == a(i)) = a(i) - length(v);
    f_temp = f_temp - (f_temp >a(i));
end
f_filled = [f;fliplr(f_temp)];

z_filled = [z_new;1./conj(z_new)];
z_filled(bdy_index + length(v)) = [];

energy_old = 0;
energy = mean(abs(beltrami_coefficient([real(z_new),imag(z_new),0.*z_new], f, v)));

iteration_count = 1; 

map_disk_opt = map_disk;

if ~silent_mode
    fprintf('Reflection completed.\n')
end
%% South pole iteration
% Iteratively compose the map with a quasi-conformal map,
% then normalize the boundary
while abs(energy_old-energy) > parameter.threshold

    energy_old = energy;
    
    mu = beltrami_coefficient([real(z_new),imag(z_new),0.*z_new], f, v);
    mu_filled = [mu;1/3*((z_new(f(:,1))./(conj(z_new(f(:,1))))).^2 + ...
        (z_new(f(:,2))./(conj(z_new(f(:,2))))).^2 + ...
        (z_new(f(:,3))./(conj(z_new(f(:,3))))).^2).*conj(mu)./...
        abs(((z_new(f(:,1))./(conj(z_new(f(:,1))))).^2 + ...
        (z_new(f(:,2))./(conj(z_new(f(:,2))))).^2 + ...
        (z_new(f(:,3))./(conj(z_new(f(:,3))))).^2))];

    % fix the points near infinity
    [~, ind] = sort(-abs(z_filled));
    fixed2 = ind(1:max(round(length(z_filled)/parameter.south),...
        min(100,length(z))));
    map_filled = linear_beltrami_solver(...
        [real(z_filled),imag(z_filled),0.*z_filled],f_filled,mu_filled,...
        fixed2,[real(z_filled(fixed2)), imag(z_filled(fixed2))]);

    z_big = complex(map_filled(:,1),map_filled(:,2));
    z_final = z_big(1:length(v));
    
    % normalization
    z_final = z_final - mean(z_final); % move centroid to zero
    if max(abs(z_final))>1
        z_final = z_final/(max(abs(z_final))); % map it into unit circle
    end
    mu_temp = beltrami_coefficient([real(z_final),imag(z_final),0.*z_final],f,v);
    map_temp = linear_beltrami_solver(...
        [real(z_final),imag(z_final),0.*z_final],f,mu_temp,...
        bdy_index,[real(z_final(bdy_index)./abs(z_final(bdy_index))), ...
        imag(z_final(bdy_index)./abs(z_final(bdy_index)))]);

    z_new = map_temp(:,1) + 1i*map_temp(:,2);

    z_filled = [z_new;1./conj(z_new)];
    z_filled(bdy_index + length(v)) = [];

    map_disk = [real(z_new), imag(z_new)];
    
    if sum(sum(isnan(map_disk))) ~= 0 
        % use the previous result in case of getting NaN entries
        break;
    end

    energy = mean(abs(beltrami_coefficient(map_disk, f, v)));
    map_disk_opt = map_disk;
    if ~silent_mode
        fprintf('Iteration %d: mean(|mu|) = %.4f\n',[iteration_count,energy]);
    end
    iteration_count = iteration_count+1;
    
    if iteration_count > 5
        % it usually converges within 5 iterations so we set 5 here
        break;
    end
end
map_disk = map_disk_opt;
map_disk(:,1) = -map_disk(:,1);

if ~silent_mode
    fprintf('South pole step completed.\n')
end
%% Compute the scaling factor for matching the target spherical cap shape
r = sqrt((1-z_lim)/(1+z_lim));

%% Search for an optimal Mobius transformation for reducing the area distortion while preserving conformality
% The optimal Mobius transformation is in the form:
%    f(z) = \frac{z-a}{1-\bar{a} z}
%    x(1): |a| (0 ~ 1)        magnitude of a
%    x(2): arg(a) (-pi ~ pi)   argument of a

% Compute the area with normalization
area_v = face_area(f,v); 
area_v = area_v/sum(area_v);

z = complex(map_disk(:,1),map_disk(:,2));

% Function for calculating the area after the Mobius transformation and the
% inverse stereographic projection
area_map = @(x) face_area(f,-stereographic(r*(z-x(1)*exp(1i*x(2)))./(1-x(1)*exp(-1i*x(2))*z)))/...
    sum(face_area(f,-stereographic(r*(z-x(1)*exp(1i*x(2)))./(1-x(1)*exp(-1i*x(2))*z))));

% objective function: mean(abs(log(area_map./area_v)))
d_area = @(x) finitemean(abs(log(area_map(x)./area_v)));

% Optimization setup
x0 = [0,0]; % initial guess, try something diferent if the result is not good
lb = [0,-pi]; % lower bound for the parameters
ub = [1,pi]; % upper bound for the parameters

% Use the Pareto-like optimisation algorithm (PSS algorithm)
if ~silent_mode
    fprintf('Finding the optimum Mobius transformation.\n')
end
PSS_algo = true;
if PSS_algo
    % Mahmoud S. Shaqfa
    plot_results =  false;
    acceptance_rate = 0.97;
    iterations = 20;
    pop_size = 30;
    PSS_results = PSS_optimizer(d_area, lb, ub, pop_size, iterations,...
        acceptance_rate, plot_results, silent_mode);
    x = PSS_results.best_solution;
%     best_runs(i, :) = x;
%     best_vals(i) = d_area(x);
else
    % Optimization (may further supply gradients for better result, not yet implemented)
    options = optimoptions('fmincon','Display','off');
    x = fmincon(d_area,x0,[],[],[],[],lb,ub,[],options);
end

% obtain the conformal parameterization with area distortion corrected
fz = r*(z-x(1)*exp(1i*x(2)))./(1-conj(x(1)*exp(1i*x(2)))*z);

% ensure everything is inside the target disk with radius r
fz = fz/max(abs(fz))*r*(1-eps); % (1-eps) avoids numerical error that maps boundary vertices to be below z = z_lim 

if ~silent_mode
    fprintf('Mobius transformation completed.\n')
end

% apply the inverse stereographic projection
map = -stereographic(fz); % add the negative sign so that the unit disk is mapped to the upper hemisphere

if ~silent_mode
    fprintf('Inverse stereographic projection completed.\n')
end
%% Objective function to be minimized it includes the area and angular
% distortions
% distortion = angle_distortion(v, f, map, plot_results);
% average_angular_distortion = mean(abs(distortion).^2);
% optimum_eval = d_area(x) + average_angular_distortion;
optimum_eval = d_area(x);
end

function m = finitemean(A)
% for avoiding the Inf values caused by division by a very small area
    m = mean(A(isfinite(A)));
end
