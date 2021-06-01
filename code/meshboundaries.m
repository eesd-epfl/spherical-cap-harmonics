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
function [bds, outer] = meshboundaries(f)

%  MESHBOUNDARIES Mesh Boundaries
%     [bds, n] = MESHBOUNDARIES(f) extracts boundaries of mesh into cells

%% this version is a little bit faster
% g = sparse(f, f(:, [2, 3, 1]), 1);
% [ii, jj, ~] = find(g+g' == 1);
% h = digraph(ii, jj);
% bds = {};
% visited = true(max(f(:)), 1);
% visited(ii) = false;
% n = 0;
% while 1
%   kk = find(~visited);
%   if isempty(kk)
%     break
%   end
%   n = n+1;
%   bds{n} = dfsearch(h, kk(1));
%   visited(bds{n}) = true;
% end
%% this version is older but works
nv = max(f(:));
v = zeros(nv, 2);
tr = triangulation(f, v);
fe = tr.freeBoundary;
if isempty(fe)
  bds = {};
  outer = [];
  return;
end
[~, order] = sort(fe(:, 1));
fe = fe(order, :);
ex = zeros(nv, 1);
vs = false(nv, 1);
nfe = size(fe, 1);
for i = 1 : nfe
  ex(fe(i, 1)) = fe(i, 2);
end
n = 0;
for i = unique(fe(:))'
  if ~vs(i)
    n = n + 1;
    [bds{n}, vs] = dfs(ex, vs, i);
  end
end
%% 
m = zeros(n, 1);
for i = 1 : n
  m(i) = size(bds{i}, 1);
end
[~, id] = sort(m, 'descend');
bds = bds(id);
outer = bds{1};

function [bd, vs] = dfs(ex, vs, i)
vs(i) = 1;
bd = i;
if ~vs(ex(i))
  [bd, vs] = dfs(ex, vs, ex(i));
  bd = [i; bd];
end
