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
%	 for Characterising the Morphology of Rough Surface Patche            %
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
% This file was re-written and modified by: Mahmoud S. Shaqfa
function y1 = deivpdp2f1(a,b,c,xf,n)
% This is a modified version of the hypergeometric function 2F1 by Mahmoud
% S. Shaqfa. The license belgongs to the original authors John Pearson.
%-------------------------------------------------------------------------%
%                             By John Pearson                             %
%                           University of Oxford                          %
%    Part of MSc dissertation 'Computation of Hypergeometric Functions'   %
%-------------------------------------------------------------------------%

% Computes the Gauss hypergeometric function 2F1(a,b;c;z) using the       %
% Dormand-Prince method to solve the initial value problem (3.2), as in   %
% Appendix H.2.                                                        %

%-------------------------------------------------------------------------%
% Input:  a =Re(a)                                                        %
%         b =Re(b)                                                        %
%         c =Re(c)                                                        %
%         xf=Point where computation is desired                           %
%         n=Number of mesh points used                                    %
% Output: y1 =Computed value of 2F1(a,b;c;z)                              %
%-------------------------------------------------------------------------%
% The vectorized version of the code:

% Define coefficients of differential operators
A1 = @(x,y) x.*(1-x);
A2 = @(x,y) c-(a+b+1).*x ;
A3 = @(x,y) -a.*b;

f1 = @(x,y,z) z;
f2 = @(x,y,z) -1./A1(x,y).*(A2(x,y).*z+A3(x,y).*y);

h = xf./n; x0 = 0;
% Define ICs using another method or a Taylor series expansion
y1=taylor_2f1(a,b,c,h,200);
z1=a.*b./c.*taylor_2f1(a+1,b+1,c+1,h,200);

% RK4 method
for i = 2:n
    x=x0+(i-1).*h;
    k1y=f1(x,y1,z1);
    k1z=f2(x,y1,z1);
    k2y=f1(x+h./5,y1+k1y.*h./5,z1+k1z.*h./5);
    k2z=f2(x+h./5,y1+k1y.*h./5,z1+k1z.*h./5);
    k3y=f1(x+3./10.*h,y1+3./40.*k2y.*h,z1+3./40.*k2z.*h);
    k3z=f2(x+3./10.*h,y1+3./40.*k2y.*h,z1+3./40.*k2z.*h);
    k4y=f1(x+4./5.*h,y1+44./45.*k1y.*h-56./15.*k2y.*h+32./9.*k3y.*h,...
        z1+44./45.*k1z.*h-56./15.*k2z.*h+32./9.*k3z.*h);
    k4z=f2(x+4./5.*h,y1+44./45.*k1y.*h-56./15.*k2y.*h+32./9.*k3y.*h,...
        z1+44./45.*k1z.*h-56./15.*k2z.*h+32./9.*k3z.*h);
    k5y=f1(x+8./9.*h,y1+19372./6561.*k1y.*h-25360./2187.*k2y.*h...
        +64448./6561.*k3y.*h-212./729.*k4y.*h,z1+19372./6561.*k1z.*h...
        -25360./2187.*k2z.*h+64448./6561.*k3z.*h-212./729.*k4z.*h);
    k5z=f2(x+8./9.*h,y1+19372./6561.*k1y.*h-25360./2187.*k2y.*h...
        +64448./6561.*k3y.*h-212./729.*k4y.*h,z1+19372./6561.*k1z.*h...
        -25360./2187.*k2z.*h+64448./6561.*k3z.*h-212./729.*k4z.*h);
    k6y=f1(x+h,y1+9017./3168.*k1y.*h-355./33.*k2y.*h+46732./5247.*k3y.*h...
        +49./176.*k4y.*h-5103./18656.*k5y.*h,z1+9017./3168.*k1z.*h...
        -355./33.*k2z.*h+46732./5247.*k3z.*h+49./176.*k4z.*h-5103./18656.*k5z.*h);
    k6z=f2(x+h,y1+9017./3168.*k1y.*h-355./33.*k2y.*h+46732./5247.*k3y.*h...
        +49./176.*k4y.*h-5103./18656.*k5y.*h,z1+9017./3168.*k1z.*h...
        -355./33.*k2z.*h+46732./5247.*k3z.*h+49./176.*k4z.*h-5103./18656.*k5z.*h);
    y1=y1+h.*(35./384.*k1y+500./1113.*k3y+125./192.*k4y...
        -2187./6784.*k5y+11./84.*k6y);
    z1=z1+h.*(35./384.*k1z+500./1113.*k3z+125./192.*k4z...
        -2187./6784.*k5z+11./84.*k6z);
end
y1(isnan(y1))=1; % to fix z = 0, where y1 = NaN.