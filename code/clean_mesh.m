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
function [v,f] = clean_mesh(v,f,arg3)
%% remove unreferenced vertices
unref_v = setdiff(1:length(v),unique(reshape(f,1,3*length(f))));
if ~isempty(unref_v)
    for i = length(unref_v):-1:1
        f(f>=unref_v(i)) = f(f>=unref_v(i))-1;
    end
    v = v(setdiff(1:length(v),unref_v),1:3);
end
%% remove extra faces
[row,~] = find(f >length(v));
f(row,:) = [];
%% remove repeated faces
[~,unirow] = unique(sort(f,2), 'rows');
f = f(unirow,1:3);
%% remove zero area surfaces
e1 = v(f(:,3),:) - v(f(:,2),:);
e2 = v(f(:,1),:) - v(f(:,3),:);
a = cross(e1,e2);
area = ((a(:,1).^2 + a(:,2).^2 + a(:,3).^2).^(1/2))/2;
f = f(area>1e-15,1:3);

if nargin == 3
    for i = length(v):-1:1
        count = sum(sum((f == i)));
        if count == 1
            %% face with 1 connection 
            [row, ~] = find(f == i);
            f(row,:) = [];
            v(i,:) = [];
            f(f>=i) = f(f>=i)-1;
        end
    end
end