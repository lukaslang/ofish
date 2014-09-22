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
function [Y, DY] = trivspharmncoeff(N, F, V, xi)
%TRIVSPHARMNCOEFF Generates the coefficients and derivatives thereof for 
%fully normalised vector spherical harmonics for several degrees in the 
%tangent basis created by TRITANBASIS.
%
%   [Y, DY] = TRIVSPHARMNCOEFF(N, F, V, xi) calls TRIVSPHARMCOEFF for each
%   degree in N.
%
%   Note that size of Y is [m, dim, 2, nq], where m is the number of faces 
%   F and dim the degree of vector spherical harmonics (of both types) 
%   specified by N. xi must be of size [m, 2, nq], where nq is the number 
%   of quadrature points on each triangle.
%
%   The third index of Y corresponds to Dx and Dy of tritanBasis. The
%   second index of DY refers to the direction of the derivation.
%
%   N must be a vector of consecutive positive integers.

% Check if N is an interval of positive consecutive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(size(xi, 1) == size(F, 1));
assert(size(xi, 2) == 2);

m = size(F, 1);
nq = size(xi, 3);

% Compute dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

% Compute offset.
offset = (N(1)-1)^2 + 2*(N(1)-1);

Y = zeros(m, dim, 2, nq);
DY = zeros(m, 2, dim, 2, nq);
for k=N
    % Compute coefficients.
    [Y1, Y2, DxY1, DyY1, DxY2, DyY2] = trivspharmcoeff(k, F, V, xi);
    % Compute index.
    ord = 2*k;
    idx = k^2 - offset;
    
    Y(:, idx:idx+ord, :, :) = Y1;
    DY(:, 1, idx:idx+ord, :, :) = DxY1;
    DY(:, 2, idx:idx+ord, :, :) = DyY1;
    
    % Create indices.
    idx = idx + dim/2;
    
    Y(:, idx:idx+ord, :, :) = Y2;
    DY(:, 1, idx:idx+ord, :, :) = DxY2;
    DY(:, 2, idx:idx+ord, :, :) = DyY2;
end
end