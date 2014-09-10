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
function [Dx, Dy] = surftanBasis(F, V, rho, xi)
%SURFTANBASIS Computes tangential basis of a sphere-like surface.
%
%   [Dx, Dy] = SURFTANBASIS(F, V, rho, xi) takes a triangulation F, V of
%   the unit sphere, a function rho evaluated at nodal points, and 
%   barycentric coordinates xi and returns the tangential basis Dx, Dy.
%
%   Note that Dx and Dy are of size m-by-3, where m is the number of
%   triangular faces F.

assert(size(F, 1) == size(rho, 1));
assert(size(F, 1) == size(xi, 1));
assert(size(rho, 2) == 6);
assert(size(xi, 2) == 2);

% Compute triangulated surface properties.
[Dx, Dy] = tritanBasis(F, V);

% Compute derivatives of polynomials at xi.
[DxA, DyA, Q] = tripoly2deriv(xi);

% Compute interpolation of derivatives of rho.
Dxrho = dot(rho, (DxA * Q)', 2);
Dyrho = dot(rho, (DyA * Q)', 2);

% Compute x.
x = trimap(F, V, xi);

% Compute interpolation of rho.
g = triinterp2(rho, xi);

Dx = bsxfun(@times, Dxrho, x) + bsxfun(@times, Dx, g);
Dy = bsxfun(@times, Dyrho, x) + bsxfun(@times, Dy, g);

end