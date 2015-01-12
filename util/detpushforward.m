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
function detphi = detpushforward(F, V, rho, xi)
%DETPUSHFORWARD Computes the determinant of the pushforward from a
%triangulated unit sphere to a sphere-like surface.
%
%   detphi = DETPUSHFORWARD(F, V, rho, xi) takes a triangulation F, V, a 
%   function rho defining a sphere-like surface and returns the determinant
%   if the pushforward phi at given barycentric coordinates xi.
%
%   Note that xi must be of size [m, 2, nq], where m is the number of
%   triangular faces and nq the number of quadrature points.
%
%   v is a vector of length m.
%
%   Note that F, V must be a triangulation of the unit sphere as dphi
%   relies on an extension to R3!.

m = size(F, 1);
nq = size(xi, 3);

% Compute points on triangles.
x = trimap(F, V, xi);

% Compute interpolation of rho at xi.
rhoi = triinterp2(rho, xi);

% Compute gradient of rho at xi.
[gradr, ~, ~] = trigradp2(F, V, rho, xi);

% Compute matrix dphi at each xi.
dphi = zeros(m, 3, 3, nq);
dphi(:, 1, 1, :) = rhoi + squeeze(x(:, 1, :) .* gradr(:, 1, :));
dphi(:, 1, 2, :) = x(:, 1, :) .* gradr(:, 2, :);
dphi(:, 1, 3, :) = x(:, 1, :) .* gradr(:, 3, :);
dphi(:, 2, 1, :) = x(:, 2, :) .* gradr(:, 1, :);
dphi(:, 2, 2, :) = rhoi + squeeze(x(:, 2, :) .* gradr(:, 2, :));
dphi(:, 2, 3, :) = x(:, 2, :) .* gradr(:, 3, :);
dphi(:, 3, 1, :) = x(:, 3, :) .* gradr(:, 1, :);
dphi(:, 3, 2, :) = x(:, 3, :) .* gradr(:, 2, :);
dphi(:, 3, 3, :) = rhoi + squeeze(x(:, 3, :) .* gradr(:, 3, :));

% Compute spherical surface normals.
normls = bsxfun(@rdivide, x, sqrt(sum(x .^ 2, 2)));

% Compute normals tensor product.
nnt(:, 1, 1, :) = normls(:, 1, :) .* normls(:, 1, :);
nnt(:, 1, 2, :) = normls(:, 1, :) .* normls(:, 2, :);
nnt(:, 1, 3, :) = normls(:, 1, :) .* normls(:, 3, :);
nnt(:, 2, 1, :) = normls(:, 2, :) .* normls(:, 1, :);
nnt(:, 2, 2, :) = normls(:, 2, :) .* normls(:, 2, :);
nnt(:, 2, 3, :) = normls(:, 2, :) .* normls(:, 3, :);
nnt(:, 3, 1, :) = normls(:, 3, :) .* normls(:, 1, :);
nnt(:, 3, 2, :) = normls(:, 3, :) .* normls(:, 2, :);
nnt(:, 3, 3, :) = normls(:, 3, :) .* normls(:, 3, :);

% Compute surface normals.
surfnormls = zeros(m, 3, nq);
for q=1:nq
    [Dx, Dy] = surftanBasis(F, V, rho, xi(:, :, q));
    surfnormls(:, :, q) = cross(Dx, Dy, 2);
end
% Normalise.
surfnormls = bsxfun(@rdivide, surfnormls, sqrt(sum(surfnormls .^ 2, 2)));

% Compute normals tensor product.
snnt(:, 1, 1, :) = surfnormls(:, 1, :) .* normls(:, 1, :);
snnt(:, 1, 2, :) = surfnormls(:, 1, :) .* normls(:, 2, :);
snnt(:, 1, 3, :) = surfnormls(:, 1, :) .* normls(:, 3, :);
snnt(:, 2, 1, :) = surfnormls(:, 2, :) .* normls(:, 1, :);
snnt(:, 2, 2, :) = surfnormls(:, 2, :) .* normls(:, 2, :);
snnt(:, 2, 3, :) = surfnormls(:, 2, :) .* normls(:, 3, :);
snnt(:, 3, 1, :) = surfnormls(:, 3, :) .* normls(:, 1, :);
snnt(:, 3, 2, :) = surfnormls(:, 3, :) .* normls(:, 2, :);
snnt(:, 3, 3, :) = surfnormls(:, 3, :) .* normls(:, 3, :);

% Compute matrix product between dphi and nnt.
dphinnt(:, 1, 1, :) = dphi(:, 1, 1, :) .* nnt(:, 1, 1, :) + dphi(:, 1, 2, :) .* nnt(:, 2, 1, :) + dphi(:, 1, 3, :) .* nnt(:, 3, 1, :);
dphinnt(:, 1, 2, :) = dphi(:, 1, 1, :) .* nnt(:, 1, 2, :) + dphi(:, 1, 2, :) .* nnt(:, 2, 2, :) + dphi(:, 1, 3, :) .* nnt(:, 3, 2, :);
dphinnt(:, 1, 3, :) = dphi(:, 1, 1, :) .* nnt(:, 1, 3, :) + dphi(:, 1, 2, :) .* nnt(:, 2, 3, :) + dphi(:, 1, 3, :) .* nnt(:, 3, 3, :);
dphinnt(:, 2, 1, :) = dphi(:, 2, 1, :) .* nnt(:, 1, 1, :) + dphi(:, 2, 2, :) .* nnt(:, 2, 1, :) + dphi(:, 2, 3, :) .* nnt(:, 3, 1, :);
dphinnt(:, 2, 2, :) = dphi(:, 2, 1, :) .* nnt(:, 1, 2, :) + dphi(:, 2, 2, :) .* nnt(:, 2, 2, :) + dphi(:, 2, 3, :) .* nnt(:, 3, 2, :);
dphinnt(:, 2, 3, :) = dphi(:, 2, 1, :) .* nnt(:, 1, 3, :) + dphi(:, 2, 2, :) .* nnt(:, 2, 3, :) + dphi(:, 2, 3, :) .* nnt(:, 3, 3, :);
dphinnt(:, 3, 1, :) = dphi(:, 3, 1, :) .* nnt(:, 1, 1, :) + dphi(:, 3, 2, :) .* nnt(:, 2, 1, :) + dphi(:, 3, 3, :) .* nnt(:, 3, 1, :);
dphinnt(:, 3, 2, :) = dphi(:, 3, 1, :) .* nnt(:, 1, 2, :) + dphi(:, 3, 2, :) .* nnt(:, 2, 2, :) + dphi(:, 3, 3, :) .* nnt(:, 3, 2, :);
dphinnt(:, 3, 3, :) = dphi(:, 3, 1, :) .* nnt(:, 1, 3, :) + dphi(:, 3, 2, :) .* nnt(:, 2, 3, :) + dphi(:, 3, 3, :) .* nnt(:, 3, 3, :);

% Compute pushforward operator.
pfwd = dphi - dphinnt + snnt;

% Compute determinant.
detphi = squeeze(pfwd(:, 1, 1, :) .* (pfwd(:, 2, 2, :) .* pfwd(:, 3, 3, :) - pfwd(:, 2, 3, :) .* pfwd(:, 3, 2, :)));
detphi = detphi - squeeze(pfwd(:, 1, 2, :) .* (pfwd(:, 2, 1, :) .* pfwd(:, 3, 3, :) - pfwd(:, 2, 3, :) .* pfwd(:, 3, 1, :)));
detphi = detphi + squeeze(pfwd(:, 1, 3, :) .* (pfwd(:, 2, 1, :) .* pfwd(:, 3, 2, :) - pfwd(:, 2, 2, :) .* pfwd(:, 3, 1, :)));

end