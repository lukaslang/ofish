% Copyright 2014 Lukas Lang
%
% This file is part of OFISH.
%
%    OFISH is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    OFISH is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with OFISH.  If not, see <http://www.gnu.org/licenses/>.
function U = ofish(N, Ns, c, F, V, f1, f2, h, deg, alpha)
%OFISH Computes the optical flow on a sphere-like surface.
%
%   U = OFISH(N, Ns, c, F, V, f1, f2, h, alpha) takes a triangulation F, V 
%   of the unit sphere and a surface specified by Ns and c, and images f1, 
%   f2 evaluated at the nodal points of the triangulation and returns the 
%   optical flow U. Scalar h > 0 is a spacing parameter, deg the degree of
%   numerical quadrature, and alpha > 0 is the regularisation parameter.
%
%   U is defined on the barycentric centroids of the faces F and is of 
%   size [size(F, 1), 3].

m = size(F, 1);
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(size(f1, 1) == m);
assert(size(f2, 1) == m);
assert(size(f1, 2) == 6);
assert(size(f2, 2) == 6);
assert(alpha >= 0);
assert(h > 0);
assert(N > 0);
assert(isscalar(N));

% Compute linear system for optical flow.
[~, A, D, b] = surflinearsystem(F, V, Ns, c, 1:N, f1, f2, h, deg, 1e-6);

% Solve linear system.
u = gmres(A + alpha * D, b, [], 1e-6, 30);

% TODO: Allow specification of evaluation points.
xi = repmat([1/3, 1/3], size(F, 1), 1);

% TODO: Write a clean implementation of the synthesis.
% Recover vector field on the surface.
idx = 1;
for k=1:N
    [Y1c, Y2c] = trivspharmcoeff(k, F, V, xi);
    Y1(:, idx:(idx+2*k), :) = Y1c;
    Y2(:, idx:(idx+2*k), :) = Y2c;
    idx = idx + 2*k + 1;
end
Y = cat(2, Y1, Y2);
C = zeros(m, 2);
for k=1:m
    C(k, 1) = dot(squeeze(Y(k, :, 1)), u);
    C(k, 2) = dot(squeeze(Y(k, :, 2)), u);
end

% Compute rho at nodal points.
Vn = normalise(trinodalpts2(F, V));
[~, rho] = surfsynth(Ns, Vn, c);

% Compute tangent basis on surface.
[Dx, Dy] = surftanBasis(F, V, rho, xi);

% Recover vector field.
U = bsxfun(@times, C(:, 1), Dx) + bsxfun(@times, C(:, 2), Dy);

end