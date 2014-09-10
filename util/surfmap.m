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
function x = surfmap(F, V, rho, xi)
%SURFMAP Maps barycentric coordinates to a sphere-like surface which is 
%interpolated with piecewise quadratic polynomials.
%
%   x = SURFMAP(F, V, rho, xi) takes a triangulation F, V, a function rho 
%   at nodal points, and barycentric coordinates xi and returns points x on
%   the sphere-like surface.
%
%   Vertices are assumed in clockwise order.
%
%   Note that xi must be of size m-by-3, where m is the number of
%   triangular faces. rho must be of size m-by-6.
%
%   x is a matrix of size m-by-3.

% Compute intermediate points on triangulated surface.
x = trimap(F, V, xi);

% Compute interpolation of rho.
g = triinterp2(rho, xi);

% Compute points on surface.
x = bsxfun(@times, g, x);

end