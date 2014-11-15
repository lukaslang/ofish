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
function [dim, A, D, b] = surflinearsystem(F, V, Ns, c, N, f1, f2, h, deg, tol, mem)
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
%   [dim, A, D, b] = SURFLINEARSYSTEM(F, V, Ns, c, N, f1, f2, h, deg, tol, mem)
%   the parameter mem additionally allows allows to specify the available 
%   memory in bytes.
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

% Compute rho at nodal points.
Vn = normalise(trinodalpts2(F, V));
[~, rho] = surfsynth(Ns, Vn, c);

% Compute determinant of pushforward at quadrature points.
detphi = detpushforward(F, V, rho, xi);

% Find indices where grad f is greater than epsilon. In areas where the
% length of the image gradient is almost zero inner products and thus
% surface integrals will be alomost zeros and thus can be excluded from
% numerical integration.
idx = any(sqrt(sum(gradf.^2, 2)) > tol, 3);

% Return all zero matrices if no data given.
if(numel(find(idx)) == 0)
    A = sparse([], [], [], dim, dim);
    D = A;
    b = sparse([], [], [], dim, 1);
    return;
end

% Constrain data.
Fc = F(idx, :);
gradfc = gradf(idx, :, :);
ac = a(idx);
dfdtc = dfdt(idx, :);
xic = xi(idx, :, :);
detphic = detphi(idx, :);

% Evaluate dot product between grad f1 and vspharm.
Z = trivspharmdot(gradfc, Fc, V, N, xic);

% Create matrix A.
A = matrixA(dim, Z, detphic, w, ac);

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

% Compute covariant derivatives.
if(nargin == 10)
    % Compute coefficients and derivatives of spherical harmonics.
    [Y, DY] = trivspharmncoeff(N, F, V, xi);
    [Z11, Z12, Z21, Z22] = surfcovderiv(G, g, ONB, Y, DY, xi);
else
    [Z11, Z12, Z21, Z22] = surfcovderivn(N, G, g, ONB, F, V, xi, mem);
end

% Compute regulariser.
D = matrixD(dim, Z11, Z12, Z21, Z22, detphi, w, a);

% Create vector b.
b = zeros(dim, 1);
parfor k=1:dim
    b(k) = - surfintegral(dfdtc .* squeeze(Z(:, k, :)) .* sqrt(abs(detphic)), w, ac);
end

end