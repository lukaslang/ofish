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
function G = surfchristoffel(F, V, rho, xi)
%SURFCHRISTOFFEL Computes the Christoffel symbols of a sphere-like surface.
%
%   G = SURFCHRISTOFFEL(F, V, rho, xi) takes a triangulation F, V of the 
%   unit sphere, a function rho evaluated at nodal points, and 
%   barycentric coordinates xi and returns the Christoffel symbols G.
%
%   rho must be of size [m, 6].
%
%   xi must be of size [m, 2, nq], where m is the number of triangular 
%   faces F and nq the quadrature dimension.
%
%   Note that G is of size [m, 2, 2, 2, nq], where G(:, i, j, k, :)
%   corresponds to \Gamma_{jk}^{i}.

m = size(F, 1);
nq = size(xi, 3);
assert(m == size(rho, 1));
assert(m == size(xi, 1));
assert(size(rho, 2) == 6);
assert(size(xi, 2) == 2);

% Compute triangulated surface properties.
[Dxt, Dyt] = tritanBasis(F, V);
Dxt = repmat(Dxt, [1, 1, nq]);
Dyt = repmat(Dyt, [1, 1, nq]);

% Compute first partial derivatives of polynomials at xi.
[DxA, DyA, Q] = tripoly2deriv(xi);

% Compute interpolation of derivatives of rho.
Dxrho = dot(rho, (DxA * Q)', 2);
Dyrho = dot(rho, (DyA * Q)', 2);

% Compute second partial derivatives of polynomials at xi.
[DxDxA, DyDxA, DxDyA, DyDyA, ~] = tripoly2deriv2(xi);

% Compute interpolation of second derivatives of rho.
DxDxrho = dot(rho, (DxDxA * Q)', 2);
DyDxrho = dot(rho, (DyDxA * Q)', 2);
DxDyrho = dot(rho, (DxDyA * Q)', 2);
DyDyrho = dot(rho, (DyDyA * Q)', 2);

% Compute x.
x = trimap(F, V, xi);

% Compute second partial derivatives of surface parameterisation.
DxDx = bsxfun(@times, DxDxrho, x) + bsxfun(@times, Dxrho, Dxt) + bsxfun(@times, Dxrho, Dxt);
DyDx = bsxfun(@times, DyDxrho, x) + bsxfun(@times, Dxrho, Dyt) + bsxfun(@times, Dyrho, Dxt);
DxDy = bsxfun(@times, DxDyrho, x) + bsxfun(@times, Dyrho, Dxt) + bsxfun(@times, Dxrho, Dyt);
DyDy = bsxfun(@times, DyDyrho, x) + bsxfun(@times, Dyrho, Dyt) + bsxfun(@times, Dyrho, Dyt);

% TODO: Improve dealing with quadrature points.
% Compute surface tangent basis.
Dx = zeros(m, 3, nq);
Dy = zeros(m, 3, nq);
ginv11 = zeros(m, nq);
ginv12 = zeros(m, nq);
ginv21 = zeros(m, nq);
ginv22 = zeros(m, nq);
for q=1:nq
    [u, v] = surftanBasis(F, V, rho, xi(:, :, q));
    Dx(:, :, q) = u;
    Dy(:, :, q) = v;
    
    % Compute inverse of metric.
    [~, ~, ginv] = metricprops(u, v);
    ginv11(:, q) = ginv(:, 1, 1);
    ginv12(:, q) = ginv(:, 1, 2);
    ginv21(:, q) = ginv(:, 2, 1);
    ginv22(:, q) = ginv(:, 2, 2);
end

% Compute derivatives of the metric.
Dxg11 = 2 * squeeze(dot(DxDx, Dx, 2));
Dyg11 = 2 * squeeze(dot(DyDx, Dx, 2));
Dxg12 = squeeze(dot(DxDx, Dy, 2)) + squeeze(dot(Dx, DxDy, 2));
Dyg12 = squeeze(dot(DyDx, Dy, 2)) + squeeze(dot(Dx, DyDy, 2));
Dxg21 = Dxg12;
Dyg21 = Dyg12;
Dxg22 = 2 * squeeze(dot(DxDy, Dy, 2));
Dyg22 = 2 * squeeze(dot(DyDy, Dy, 2));

% Compute Christoffel symbols. General formula is:
% G(I, K, L) = 0.5 * (ginvI1 .* (DLg1K + DKg1L - DxgKL) + ginvI2 .* (DLg2K + DKg2L - DygKL));
G = zeros(m, 2, 2, 2, nq);
G(:, 1, 1, 1, :) = 0.5 * (ginv11 .* (Dxg11 + Dxg11 - Dxg11) + ginv12 .* (Dxg21 + Dxg21 - Dyg11));
G(:, 1, 1, 2, :) = 0.5 * (ginv11 .* (Dyg11 + Dxg12 - Dxg12) + ginv12 .* (Dyg21 + Dxg22 - Dyg12));
G(:, 1, 2, 1, :) = 0.5 * (ginv11 .* (Dxg12 + Dyg11 - Dxg21) + ginv12 .* (Dxg22 + Dyg21 - Dyg21));
G(:, 1, 2, 2, :) = 0.5 * (ginv11 .* (Dyg12 + Dyg12 - Dxg22) + ginv12 .* (Dyg22 + Dyg22 - Dyg22));
G(:, 2, 1, 1, :) = 0.5 * (ginv21 .* (Dxg11 + Dxg11 - Dxg11) + ginv22 .* (Dxg21 + Dxg21 - Dyg11));
G(:, 2, 1, 2, :) = 0.5 * (ginv21 .* (Dyg11 + Dxg12 - Dxg12) + ginv22 .* (Dyg21 + Dxg22 - Dyg12));
G(:, 2, 2, 1, :) = 0.5 * (ginv21 .* (Dxg12 + Dyg11 - Dxg21) + ginv22 .* (Dxg22 + Dyg21 - Dyg21));
G(:, 2, 2, 2, :) = 0.5 * (ginv21 .* (Dyg12 + Dyg12 - Dxg22) + ginv22 .* (Dyg22 + Dyg22 - Dyg22));

end