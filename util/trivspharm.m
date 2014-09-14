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
%   Note that size(Yi) = [m, 2*N + 1, 3] for i={1, 2}, where m is the
%   number of faces F. xi must be of size [m, 2, nq].

assert(isscalar(N));
assert(N > 0);
assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(size(xi, 1) == size(F, 1));
assert(size(xi, 2) == 2);

m = size(F, 1);
nq = size(xi, 3);

% Compute face normals.
Fn = facenormals(F, V);

% Project scalar spherical harmonics at nodal points.
Vn = trinodalpts2(F, V);
Ynj = zeros(m, 2*N+1, size(Vn, 3));
for k=1:size(Vn, 3)
    Ynj(:, :, k) = spharm(N, normalise(Vn(:, :, k)));
end

% Compute vector spherical harmonics.
Y1 = zeros(m, 2*N + 1, 3, nq);
Y2 = zeros(m, 2*N + 1, 3, nq);
for q=1:nq
    for k=1:2*N+1
        Y1(:, k, :, q) = trigradp2(F, V, squeeze(Ynj(:, k, :)), xi(:, :, q));
        Y2(:, k, :, q) = cross(squeeze(Y1(:, k, :, q)), Fn);
    end
end

% Normalise.
Y1 = Y1 ./ sqrt(N*(N+1));
Y2 = Y2 ./ sqrt(N*(N+1));

end