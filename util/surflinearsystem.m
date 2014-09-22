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
function [dim, A, D, b] = surflinearsystem(F, V, Ns, c, N, f1, f2, h, deg, tol)
%SURFLINEARSYSTEM Computes the linear system used in optical flow on a
%sphere-like surface.
%
%   [dim, A, D, b] = SURFLINEARSYSTEM(F, V, Ns, c, N, f1, f2, h, deg, tol)
%   takes a triangulation F, V, a sphere-like surface specified by Ns and 
%   c, solution degrees N, images f1, f2 evaluted at the nodal points of 
%   the triangulation, a finite-difference parameter h for the time 
%   derivative of f1, f2, a degree deg of quadrature, a scalar tol which is
%   is a tolerance parameter for numerical integration, and returns the 
%   linear system needed for optical flow computation.
%
%   Note that degrees N must be a vector of consecutive positive integers!
%   
%   The linear system returned is
%
%   (A + alpha * D) * x = b.
%
%   dim is the dimension of the linear system.
%
%   A and D are matrices of size dim-by-dim.
%
%   b is the right hand side and is of length dim.

% Check if N is an interval of consecutive positive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

m = size(F, 1);
assert(h > 0);
assert(size(f1, 2) == 6);
assert(size(f2, 2) == 6);
assert(size(f1, 1) == m);
assert(size(f2, 1) == m);

% TODO: Implement tol restriction for integration.

% Get quadrature rule.
[xiq, w] = triquadrature(deg);
xi = repmat(permute(xiq, [3, 2, 1]), [m, 1, 1]);
nq = length(w);

% Compute triangle areas to be used in integration.
a = triangArea(F, V);

% Compute dimension of vector spherical harmonics.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

% Compute approximate time derivative at quadrature points.
dfdt = triinterp2((f2 - f1) ./ h, xi);

% Compute surface gradient of first image.
gradf = trigradp2(F, V, f1, xi);

% Evaluate dot product between grad f1 and vspharm.
Z = trivspharmdot(gradf, F, V, N, xi);

% Compute rho at nodal points.
Vn = normalise(trinodalpts2(F, V));
[~, rho] = surfsynth(Ns, Vn, c);

% Compute determinant of pushforward at quadrature points.
detphi = detpushforward(F, V, rho, xi);

% Create matrix A.
A = matrixA(dim, Z, detphi, w, a);

% Compute Christoffel symbols.
G = surfchristoffel(F, V, rho, xi);

% Compute metric and orthonormal basis.
g = zeros(m, 2, 2, nq);
for q=1:nq
    [Dx, Dy] = surftanBasis(F, V, rho, xi(:, :, q));
    g(:, :, :, q) = metricprops(Dx, Dy);
end

% Compute orthonormal basis at nodal points.
xin(:, :, 1) = repmat([0, 0], [m, 1]);
xin(:, :, 2) = repmat([0, 1], [m, 1]);
xin(:, :, 3) = repmat([1, 0], [m, 1]);
xin(:, :, 4) = repmat([0, 1/2], [m, 1]);
xin(:, :, 5) = repmat([1/2, 1/2], [m, 1]);
xin(:, :, 6) = repmat([1/2, 0], [m, 1]);
ONB = zeros(2, 2, m, 6);
for k=1:6
    [Dx, Dy] = surftanBasis(F, V, rho, xin(:, :, k));
    ONB(:, :, :, k) = orthonormalise(Dx, Dy);
end

% Compute coefficients and derivatives of spherical harmonics.
[Y, DY] = trivspharmncoeff(N, F, V, xi);

% Compute covariant derivatives.
[Z11, Z12, Z21, Z22] = surfcovderiv(G, g, ONB, Y, DY, xi);

% Compute regulariser.
D = matrixD(dim, Z11, Z12, Z21, Z22, detphi, w, a);

% Create vector b.
b = zeros(dim, 1);
parfor k=1:dim
    b(k) = - surfintegral(dfdt .* squeeze(Z(:, k, :)) .* sqrt(abs(detphi)), w, a);
end

end

function A = matrixA(dim, Z, detphi, w, a)
    A = zeros(dim, dim);
    for p=1:dim
        for q=1:p
            A(p, q) = surfintegral(squeeze(Z(:, p, :) .* Z(:, q, :)) .* sqrt(abs(detphi)), w, a);
            A(q, p) = A(p, q);
        end
    end
end

function D = matrixD(dim, Z11, Z12, Z21, Z22, detphi, w, a)
    D = zeros(dim, dim);
    for p=1:dim
        for q=1:p
            f = Z11(:, p, :) .* Z11(:, q, :) + Z12(:, p, :) .* Z12(:, q, :) + Z21(:, p, :) .* Z21(:, q, :) + Z22(:, p, :) .* Z22(:, q, :);
            D(p, q) = surfintegral(squeeze(f) .* sqrt(abs(detphi)), w, a);
            D(q, p) = D(p, q);
        end
    end

    % Add diagonal.
    d = zeros(dim, 1);
    for p=1:dim
        f = Z11(:, p, :) .^2 + Z12(:, p, :) .^2 + Z21(:, p, :) .^2 + Z22(:, p, :) .^2;
        d(p) = surfintegral(squeeze(f) .* sqrt(abs(detphi)), w, a);
    end
    D = diag(d) + D / 2;
end