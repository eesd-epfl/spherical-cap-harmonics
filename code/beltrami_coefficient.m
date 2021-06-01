function mu = beltrami_coefficient(v, f, map)
% Compute the Beltrami coefficient of a mapping.
% 
% If you use this code in your own work, please cite the following paper:
% [1] P. T. Choi and L. M. Lui, 
%     "Fast Disk Conformal Parameterization of Simply-Connected Open Surfaces."
%     Journal of Scientific Computing, 65(3), pp. 1065-1090, 2015.
%
% Copyright (c) 2014-2018, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi


nf = length(f);
Mi = reshape([1:nf;1:nf;1:nf], [1,3*nf]);
Mj = reshape(f', [1,3*nf]);

e1 = v(f(:,3),1:2) - v(f(:,2),1:2);
e2 = v(f(:,1),1:2) - v(f(:,3),1:2);
e3 = v(f(:,2),1:2) - v(f(:,1),1:2);

area = (-e2(:,1).*e1(:,2) + e1(:,1).*e2(:,2))'/2;
area = [area;area;area];

Mx = reshape([e1(:,2),e2(:,2),e3(:,2)]'./area /2 , [1, 3*nf]);
My = -reshape([e1(:,1),e2(:,1),e3(:,1)]'./area /2 , [1, 3*nf]);

Dx = sparse(Mi,Mj,Mx);
Dy = sparse(Mi,Mj,My);

dXdu = Dx*map(:,1);
dXdv = Dy*map(:,1);
dYdu = Dx*map(:,2);
dYdv = Dy*map(:,2);
dZdu = Dx*map(:,3);
dZdv = Dy*map(:,3);

E = dXdu.^2 + dYdu.^2 + dZdu.^2;
G = dXdv.^2 + dYdv.^2 + dZdv.^2;
F = dXdu.*dXdv + dYdu.*dYdv + dZdu.*dZdv;
mu = (E - G + 2 * 1i * F) ./ (E + G + 2*sqrt(E.*G - F.^2));
end
