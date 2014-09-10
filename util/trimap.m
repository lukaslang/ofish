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
function x = trimap(F, V, xi)
%TRIMAP Maps barycentric coordinates to a triangulated surface.
%
%   x = TRIMAP(F, V, xi) takes a triangulation F, V and barycentric 
%   coordinates xi and returns points x on the triangulated surface.
%
%   Vertices are assumed in clockwise order.
%
%   Note that xi must be of size m-by-3, where m is the number of
%   triangular faces.
%
%   x is a matrix of size m-by-3.

% Compute vectors.
u = bsxfun(@times, xi(:, 1), V(F(:, 3), :) - V(F(:, 1), :));
v = bsxfun(@times, xi(:, 2), V(F(:, 2), :) - V(F(:, 1), :));

% Compute points.
x = V(F(:, 1), :) + u + v;

end