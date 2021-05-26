function [z , PixelWidth , PSD] = artificial_surf(sigma, H, Lx, m , n, qr)
% Modified to control the seed number and export stl files by Mahmoud S.
% Shaqfa.

% https://au.mathworks.com/matlabcentral/fileexchange/60817-surface-generator-artificial-randomly-rough-surfaces
% Generates artifial randomly rough surfaces with given parameters.
% In other words, generates fractal topographies with different fractal
% dimensions.
% ======================= inputs
% parameters (in SI units)
% sigma: standard deviation , i.e. root-mean-square roughness Rq(m)
% H: Hurst exponent (roughness exponent), 0<= H <= 1
% It relates to the fractal dimension of a surface by D = 3-H.
% Lx: length of topography in x direction.
% m: number of pixels in x
% n: number of pixels in y
% !!! please note that setting x and y to be a power of 2 will reduce
% computing time, like using the numbers (256,512,1024,2048,...)
% !!! a normal desktop computer cannot handle pixels above 2048
% !!! you can increase resolution (i.e. pixelwidth) by increasing n and m,
% or reducing size of final topography (reducing Lx)
% qr (((OPTIONAL parameter)))):
% cut-off wavevector or roll-off wavevector (1/m), it relates to
% roll-off wavelength by qr = (2*pi)/lambda_r, where lambda_r is the
% roll-off wavelength. Check the image that I uploaded to Mathworks for its
% meaning.
% !!! note that qr cannot be smaller than image size, i.e. qr > (2*pi/L)
% where L is the length of image in x and/or y direction ( Lx and Ly ).
% Also, qr cannot be bigger than Nyquist frequency qr < (pi/PixelWidth)
% ======================== outputs
% z : the surface height profile of a randomly rough surface (m)
% PixelWidth : spatial resolution, i.e. the distance between each two points
% in the generated height topography z (m). This directly relates to the size
% of generated surface by the relation Lx = (m-1) * PixelWidth
% PSD: this is a radially averaged power spectrum  of the generated surface.
% To check your results,when you have generated the surface z, again apply
% power spectrum technique to z topography by code:
% http://se.mathworks.com/matlabcentral/fileexchange/54297-radially-averaged-surface-roughness-power-spectrum--psd-
% then compare it the PSD output here. They must be the same!
% ======================== plot the topography
% [n,m] = size(z);
% x = linspace(0,(m-1) * PixelWidth , m);
% y = linspace(0,(n-1) * PixelWidth , n);
% [X,Y] = meshgrid(x,y);
% mesh(X,Y,z)
% colormap jet
% axis equal
% =========================== example inputs
% [z , PixelWidth, PSD] = artificial_surf(0.5e-3, 0.8, 0.1, 512 , 512);
% generates surface z with sigma = 0.5 mm, H = 0.8 Hurst exponent
% [i.e. a surface with fractal dimension D = 2.2], Lx = 10 cm [the size of
% topography in x direction], m = n = 512 [results in a square surface,
% i.e. Lx = Ly = 10 cm]. With these parameters simply the Pixel Width of
% the final topography is Lx/(m-1) = 1.96e-4.
% [z , PixelWidth] = artificial_surf(1e-3, 0.7, 0.1, 1024 , 512, 1000);
% generated surface z with sigma = 1 mm, H = 0.7 (i.e. D = 2.3) with Lx =
% 10 cm. Pixelwidth = Lx / (m-1) = 9.8e-5 m.
% Since m and n are not equal, the surface is rectangular (not square) and
% Ly = n * Pixelwidth = 5 cm. The surface has a roll-off wavevector
% ar qr = 1000 (1/m) which is equal to lambda_r = (2*pi/qr) = 6.3 mm.
% !!! Play with qr to realize its physical meaning in a topography
% =========================== general notes
% !!! note that the code assumes same pixel width in x and y direction,
% which is typical of measurement instruments
% !!! The generated surface could be rectangular (if n =~ m)
% !!! Lx = m * PixelWidth; % image length in x direction
% !!! Ly = n * PixelWidth; % image length in y direction
% =========================================================================
% Check number of inputs (Mahmoud S. Shaqfa)
if nargin < 5 || nargin > 6
    myFolders = split(pwd,"/");
    myFolder = myFolders{length(myFolders)};
    if myFolder == "generating_artificial_surfaces"
        cd ..
    end
    addpath('code')
    addpath('stlTools')
    addpath('input_geom')
    addpath('distmesh')
    
    close all; clc;
    %     error(['The code requires at least 5, and maximum 6, inputs.',...
    %         'Check the code description']);
    % [z , PixelWidth, PSD] = artificial_surf(0.5e-3, 0.8, 0.1, 512 , 512);
    addpath('stlTools')
    sigma = 0.7e-2; H = 0.40; Lx = 0.1; m = 128;
    n = m;
    seed = 10; % The random seed number used for reproducability
    rng(seed)
end
% =========================================================================
% make surface size an even number
if mod(n,2)
    n = n -1;
end
if mod(m,2)
    m = m -1;
end
% =========================================================================
PixelWidth = Lx/m;
Lx = m * PixelWidth; % image length in x direction
Ly = n * PixelWidth; % image length in y direction
% =========================================================================
% Wavevectors (note that q = (2*pi) / lambda; where lambda is wavelength.
% qx
qx = zeros(m,1);
for k = 0:m-1
    qx(k+1)=(2*pi/m)*(k);
end
qx = (unwrap((fftshift(qx))-2*pi))/PixelWidth;
% qy
qy = zeros(n,1);
for k = 0:n-1
    qy(k+1)=(2*pi/n)*(k);
end
qy = (unwrap((fftshift(qy))-2*pi))/PixelWidth;
% 2D matrix of q values
[qxx,qyy] = meshgrid(qx , qy);
[~,rho] = cart2pol(qxx,qyy);
% =========================================================================
% handle qr case
 if ~exist('qr','var')
     % default qr
      qr = 0; % no roll-off
 end
% 2D matrix of Cq values
Cq = zeros(n,m);
for i = 1:m
    for j = 1:n
        if rho(j,i) < qr
            Cq(j,i) = qr^(-2*(H+1));
        else
            Cq(j,i) = rho(j,i).^(-2*(H+1));
        end
    end
end
% =====================paraview_patch(X,Y, z)====================================================
% applying rms
Cq(n/2+1,m/2+1) = 0; % remove mean
RMS_F2D = sqrt((sum(sum(Cq)))*(((2*pi)^2)/(Lx*Ly)));
alfa = sigma/RMS_F2D;
Cq = Cq.*alfa^2;
% =========================================================================
% Radial Averaging : DATA WILL BE USEFUL FOR CHECKING RESULTS
rhof = floor(rho);
J = 200; % resolution in q space (increase if you want)
qrmin = log10(sqrt((((2*pi)/Lx)^2+((2*pi)/Ly)^2)));
qrmax = log10(sqrt(qx(end).^2 + qy(end).^2)); % Nyquist
q = floor(10.^linspace(qrmin,qrmax,J));
C_AVE = zeros(1,length(q));
ind = cell(length(q)-1 , 1);
for j = 1:length(q)-1
    ind{j} = find( rhof > q(j) & rhof <=(q(j+1)));
    C_AVE(j) = nanmean(Cq(ind{j}));
end
ind = ~isnan(C_AVE);
C = C_AVE(ind);
q = q(ind);
% =========================================================================
% reversing opertation: PSD to fft
Bq = sqrt(Cq./(PixelWidth^2/((n*m)*((2*pi)^2))));
% =========================================================================
% apply conjugate symmetry to magnitude
Bq(1,1) = 0;
Bq(1,m/2+1) = 0;
Bq(n/2+1,m/2+1) = 0;
Bq(n/2+1,1) = 0;
Bq(2:end,2:m/2) = rot90(Bq(2:end,m/2+2:end),2);
Bq(1,2:m/2) = rot90(Bq(1,m/2+2:end),2);
Bq(n/2+2:end,1) = rot90(Bq(2:n/2,1),2);
Bq(n/2+2:end,m/2+1) = rot90(Bq(2:n/2,m/2+1),2);
% =========================================================================
% defining a random phase between -pi and phi (due to fftshift,
% otherwise between 0 and 2pi)
phi =  -pi + (pi+pi)*rand(n,m);
% =========================================================================
% apply conjugate symmetry to phase
phi(1,1) = 0;
phi(1,m/2+1) = 0;
phi(n/2+1,m/2+1) = 0;
phi(n/2+1,1) = 0;
phi(2:end,2:m/2) = -rot90(phi(2:end,m/2+2:end),2);
phi(1,2:m/2) = -rot90(phi(1,m/2+2:end),2);
phi(n/2+2:end,1) = -rot90(phi(2:n/2,1),2);
phi(n/2+2:end,m/2+1) = -rot90(phi(2:n/2,m/2+1),2);
% =========================================================================
% Generate topography
[a,b] = pol2cart(phi,Bq);
Hm = complex(a,b); % Complex Hm with 'Bq' as abosulte values and 'phi' as
% phase components.
z = ifft2(ifftshift((Hm))); % generate surface
PSD.qx = qx;
PSD.qy = qy;
PSD.Cq = Cq;
PSD.q = q;
PSD.C = C;

% Plot and save meshes
[n,m] = size(z);
x = linspace(0,(m-1) * PixelWidth , m);
y = linspace(0,(n-1) * PixelWidth , n);
[X,Y] = meshgrid(x,y);
stl_name = ['input_geom/fractal_surface_H_', num2str(H), '.stl'];
stlWrite(stl_name, X, Y, z)
[v, f] = stlRead(stl_name);
paraview_patch(v,f, v(:,3))