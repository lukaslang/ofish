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
function Z = trivspharmdot(v, F, V, N, xi)
%TRIVSPHARMDOT Computes the dot product between a vector field and vector 
%spherical harmonics of multiple degrees.
%
%   Z = TRIVSPHARMDOT(v, F, V, N, xi) takes a vector field v at xi, a 
%   triangulation F, V of the unit sphere, degrees N, barycentric 
%   coordinates xi and returns the dot product between v and all vector 
%   spherical harmonics of degrees N.
%
%   Note that xi must be a matrix of size [m, 2, nq], where m is the number
%   of triangular faces size(F, 1) and nq is the number of quadrature
%   points.
%   
%   Note that v must be of size [m, 3, nq].
%
%   Z is of size [m, dim, nq], where dim is the dimension induced by the 
%   given degrees N, which must be a vector of consecutive positive 
%   integers.
%
%   Note that parallelisation was chosen this way since for each k the
%   number of dot products to calculate is 2*(2*k+1). Each of the 2*k+1
%   computations can be done with no communcation overhead.

assert(size(F, 2) == 3);
assert(size(V, 2) == 3);
assert(size(v, 2) == 3);
assert(size(v, 1) == size(F, 1));
assert(size(v, 3) == size(xi, 3));

% Check if N is an interval of positive consecutive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

nq = size(xi, 3);

% Compute dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

% Compute offset.
offset = (N(1)-1)^2 + 2*(N(1)-1);

% Compute dot products.
Z = zeros(size(F, 1), dim, nq);
for k=N
    % Create vector spherical harmonics of degree k.
    [Y1, Y2] = trivspharm(k, F, V, xi);
    % Run through quadrature dimension.
    for q=1:nq
        vq = v(:, :, q);        
        % Compute index.
        idx = k^2 - offset - 1;
        % Run through all orders.
        parfor l=1:2*k+1
            Z(:, idx + l, q) = dot(vq, squeeze(Y1(:, l, :, q)), 2);
        end
        % Create indices.
        idx = idx + dim/2;
        parfor l=1:2*k+1
            Z(:, idx + l, q) = dot(vq, squeeze(Y2(:, l, :, q)), 2);
        end
    end
end
end