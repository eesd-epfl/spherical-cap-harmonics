function [field1,field2,tx,ty]=Brownian_field(H,n)
%% simulate Fractional Brownian field on unit disk, with Hurst parameter 'H';
%  Note that the covariance function is isotropic, see reference below.
% INPUTS:
%        - 'H' is the Hurst parameter of the Gaussian process
%        - 'n' is the number of grid points, where 'n' is a power of 2;
%            if the 'n' supplied is not a power of two,
%            then we set n=2^ceil(log2(n)); default is n=2^8;
% OUTPUT:
%          - two statistically independent fields 'field1' and 'field2'
%            over unit disk; if not output requested, then function
%            outputs a figure of one of the fields
%          - vectors 'tx' and 'ty' so that the field is plotted via
%            surf(tx,ty,field1,'EdgeColor','none')
%
% Example:
%  [field1,field2,tx,ty]=Brownian_field(.9,2^10);
%   surf(tx,ty,field2,'EdgeColor','none'),colormap bone
%% Reference:
% Kroese, D. P., & Botev, Z. I. (2015). Spatial Process Simulation.
% In Stochastic Geometry, Spatial Statistics and Random Fields(pp. 369-404)
% Springer International Publishing, DOI: 10.1007/978-3-319-10064-7_12
if nargin==0
    H = 0.5;
    n=2^8; % default value of points
else
    n=2^ceil(log2(n));
end

if (H>1)|(H<0) % Hurst parameter error check
    error('Hurst parameter must be between 0 and 1')
end

R=2; % [0,R]^2 grid, may have to extract only [0,R/2]^2
m=n; % size of grid is m*n; covariance matrix is m^2*n^2
tx=[1:n]/n*R; ty=[1:m]/m*R; % create grid for field
Rows=zeros(m,n);
for i=1:n
    for j=1:m % rows of blocks of cov matrix
        Rows(j,i)=rho([tx(i),ty(j)],[tx(1),ty(1)],R,2*H);
    end
end
BlkCirc_row=[Rows, Rows(:,end-1:-1:2);
    Rows(end-1:-1:2,:), Rows(end-1:-1:2,end-1:-1:2)];
% compute eigen-values
lam=real(fft2(BlkCirc_row))/(4*(m-1)*(n-1));
lam=sqrt(lam);
% generate field with covariance given by block circular matrix
Z=complex(randn(2*(m-1),2*(n-1)),randn(2*(m-1),2*(n-1)));
F=fft2(lam.*Z);
F=F(1:m,1:n); % extract sub-block with desired covariance
[out,c0,c2]=rho([0,0],[0,0],R,2*H);
field1=real(F); field2=imag(F); % two independent fields
field1=field1-field1(1,1); % set field zero at origin
field2=field2-field2(1,1); % set field zero at origin
% make correction for embedding with a term c2*r^2
field1=field1 + kron(ty'*randn,tx*randn)*sqrt(2*c2);
field2=field2 + kron(ty'*randn,tx*randn)*sqrt(2*c2);
[X,Y]=meshgrid(tx,ty);
field1((X.^2+Y.^2)>1)=nan;field2((X.^2+Y.^2)>1)=nan;
field1=field1(1:n/2,1:m/2);field2=field2(1:n/2,1:m/2);
tx=tx(1:n/2);ty=ty(1:m/2);
if nargout==0
    surf(tx,ty,field1,'EdgeColor','none')
    title('Fractional Gaussian Field on Unit Disk')
    colormap bone
    stl_name = ['input_geom/fractal_surface_H_', num2str(H), '.stl'];
    stlWrite(stl_name, tx,ty, field1)
    [v, f] = stlRead(stl_name);
    paraview_patch(v,f, v(:,3))
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out,c0,c2]=rho(x,y,R,alpha)
% embedding of covariance function on a [0,R]^2 grid
if alpha<=1.5 % alpha=2*H, where H is the Hurst parameter
    beta=0;c2=alpha/2;c0=1-alpha/2;
else % parameters ensure piecewise function twice differentiable
    beta=alpha*(2-alpha)/(3*R*(R^2-1)); c2=(alpha-beta*(R-1)^2*(R+2))/2;
    c0=beta*(R-1)^3+1-c2;
end
% create continuous isotropic function
r=sqrt((x(1)-y(1))^2+(x(2)-y(2))^2);
if r<=1
    out=c0-r^alpha+c2*r^2;
elseif r<=R
    out=beta*(R-r)^3/r;
else
    out=0;
end
end