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
function g = triinterp2(f, xi)
%TRIINTERP2 Computes the piecewise quadratic interpolation of a function on
%the reference triangle.
%
%   g = TRIINTERP2(f, xi) takes a matrix f, which are the values at the 
%   nodal points, and returns the interpolation at points xi which are 
%   given in barycentric coordinates.
%
%   Note that f must be of size m-by-6 where each row corresponds to the
%   six quadrature points and xi is of size m-by-2-by-nq, where m is the 
%   number of triangles and nq the number of quadrature points.
%
%   Note that g is a matrix of size [m, nq].

m = size(xi, 1);
nq = size(xi, 3);
assert(size(f, 1) == m);
assert(size(f, 2) == 6);
assert(size(xi, 2) == 2);

g = zeros(m, nq);
parfor k=1:nq
    % Compute polynomials at xi.
    [A, Q] = tripoly2(xi(:, :, k));
    % Compute interpolation.
    g(:, k) = dot(f, (A * Q)', 2);
end
end