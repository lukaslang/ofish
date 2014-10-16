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
function test_suite = matrixDTest
    initTestSuite;
end

function resultTest

% Generate icosahedron.
[F, V] = sphTriang(3);
m = size(F, 1);

% Generate constant function rho.
rho = ones(m, 6);

% Pick coordinates.
nq = 2;
xi = repmat([1/3, 1/3], [m, 1, nq]);

% Create a transformation matrix for orthonormal basis.
A(1, 1, :, :) = ones(m, 6);
A(1, 2, :, :) = zeros(m, 6);
A(2, 1, :, :) = zeros(m, 6);
A(2, 2, :, :) = ones(m, 6);

% Compute Christoffel symbols.
G = surfchristoffel(F, V, rho, xi);

% Compute metric.
parfor q=1:nq
    [Dx, Dy] = surftanBasis(F, V, rho, xi(:, :, q));
    g(:, :, :, q) = metricprops(Dx, Dy);
end

% Create vector fields U.
N = 5;
[Y, DY] = trivspharmncoeff(N, F, V, xi);
dim = size(Y, 2);

% Compute covariant derivative.
[Z11, Z12, Z21, Z22] = surfcovderiv(G, g, A, Y, DY, xi);
assertEqual(size(Z11), [m, dim, nq]);
assertEqual(size(Z12), [m, dim, nq]);
assertEqual(size(Z21), [m, dim, nq]);
assertEqual(size(Z22), [m, dim, nq]);

detphi = ones(m, nq);
w = ones(nq, 1);
a = triangArea(F, V);

D = matrixD(dim, Z11, Z12, Z21, Z22, detphi, w, a);
assertEqual(size(D), [dim, dim]);
assertAlmostEqual(D, D');

end