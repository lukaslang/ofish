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
function [Y1, Y2] = trivspharm(N, F, V, xi)
%TRIVSPHARM Generates fully normalised vector spherical harmonics on a
%triangulated surface with scalar spherical harmonics interpolated 
%piecewise quadratically.
%
%   [Y1, Y2] = TRIVSPHARM(N, F, V, xi) takes a triangulation F, V of the 
%   unit sphere and returns fully normalised vector spherical harmonics 
%   Y_Nj of degree N > 0 and j=-N,...,N at barycentric coordinates xi.
%
%   Note that size(Yi) = [n, 2*N + 1, 3] for i={1, 2}, where n is the
%   number of faces F. xi must be of size [n, 2].

assert(isscalar(N));
assert(N > 0);
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(size(xi, 1) == size(F, 1));
assert(size(xi, 2) == 2);

n = size(F, 1);

% Compute face normals.
Fn = facenormals(F, V);

% Project scalar spherical harmonics at nodal points.
Vn = trinodalpts2(F, V);
Ynj = zeros(n, 2*N+1, size(Vn, 3));
for k=1:size(Vn, 3)
    Ynj(:, :, k) = spharm(N, normalise(Vn(:, :, k)));
end

% Compute vector spherical harmonics.
Y1 = zeros(n, 2*N + 1, 3);
Y2 = zeros(n, 2*N + 1, 3);
for k=1:2*N+1
    Y1(:, k, :) = trigradp2(F, V, squeeze(Ynj(:, k, :)), xi);
    Y2(:, k, :) = cross(squeeze(Y1(:, k, :)), Fn);
end

% Normalise.
Y1 = Y1 ./ sqrt(N*(N+1));
Y2 = Y2 ./ sqrt(N*(N+1));

end