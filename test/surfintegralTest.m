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
function test_suite = surfintegralTest
    initTestSuite;
end

function refTriConstFunTest

% Create reference triangle.
F = [1, 2, 3];
V = [0, 0, 0; 0, 1, 0; 1, 0, 0];
m = size(F, 1);

% Generate constant function rho.
rho = ones(m, 6);

% Get quadrature rule.
[xi, w] = triquadrature(1);
xi = repmat(permute(xi, [3, 2, 1]), [m, 1, 1]);

% Evaluate constant function at quadrature points.
detphi = detpushforward(F, V, rho, xi);
f = ones(m, length(w)) .* sqrt(abs(detphi));

% Compute triangle area.
a = triangArea(F, V);

% Compute integral.
v = surfintegral(f, w, a);

assertTrue(isscalar(v));
assertAlmostEqual(v, 0.5);

end

function higherQuadratureTest

% Create reference triangle.
F = [1, 2, 3; 2, 4, 3];
V = [0, 0, 0; 0, 1, 0; 1, 0, 0; 1, 1, 0];
m = size(F, 1);

% Generate constant function rho.
rho = ones(m, 6);

% Get quadrature rule.
[xi, w] = triquadrature(7);
xi = repmat(permute(xi, [3, 2, 1]), [m, 1, 1]);

% Evaluate constant function at quadrature points.
detphi = detpushforward(F, V, rho, xi);
f = ones(m, length(w)) .* sqrt(abs(detphi));

% Compute triangle area.
a = triangArea(F, V);

% Compute integral.
v = surfintegral(f, w, a);

assertTrue(isscalar(v));
assertAlmostEqual(v, 1, 1e-6);

end

function triConstFunTest

% Create bigger triangle.
F = [1, 2, 3];
V = [0, 0, 0; 0, 2, 0; 2, 0, 0];
m = size(F, 1);

% Generate constant function rho.
rho = ones(m, 6);

% Get quadrature rule.
[xi, w] = triquadrature(1);
xi = repmat(permute(xi, [3, 2, 1]), [m, 1, 1]);

% Evaluate constant function at quadrature points.
detphi = detpushforward(F, V, rho, xi);
f = ones(m, length(w)) .* sqrt(abs(detphi));

% Compute triangle area.
a = triangArea(F, V);

% Compute integral.
v = surfintegral(f, w, a);

assertTrue(isscalar(v));
assertAlmostEqual(v, 2);

end

function constantFunctionTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(4);
m = size(F, 1);

% Generate constant function rho.
rho = ones(m, 6);

% Get quadrature rule.
[xi, w] = triquadrature(1);
xi = repmat(permute(xi, [3, 2, 1]), [m, 1, 1]);

% Evaluate constant function at quadrature points.
detphi = detpushforward(F, V, rho, xi);
f = ones(m, length(w)) .* sqrt(abs(detphi));

% Compute triangle area.
a = triangArea(F, V);

% Compute surface integral over function which is constant one.
v = surfintegral(f, w, a);

assertTrue(isscalar(v));
assertAlmostEqual(v, 4*pi, 1e-2);

end

function constantFunctionOnInterpolatedSphereTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Scale nodes.
radius = 2;
V = radius .* V;

Vn = trinodalpts2(F, V);
for k=1:6
    % Generate constant function rho.
    rho(:, k) = radius ./ sqrt(sum(Vn(:, :, k).^2, 2));
end

% Get quadrature rule.
[xi, w] = triquadrature(1);
xi = repmat(permute(xi, [3, 2, 1]), [m, 1, 1]);

% Evaluate constant function at quadrature points.
detphi = detpushforward(F, V, rho, xi);
f = ones(m, length(w)) .* sqrt(abs(detphi));

% Compute triangle area.
a = triangArea(F, V);

% Compute surface integral over function which is constant one.
v = surfintegral(f, w, a);

assertTrue(isscalar(v));
assertAlmostEqual(v, 4*pi*radius^2, 1e-2);

end

function compareTest

% Create triangulation of unit sphere.
[F, V] = sphTriang;
m = size(F, 1);

% Scale nodes.
radius = 1;
V = radius .* V;

Vn = trinodalpts2(F, V);
for k=1:6
    % Generate constant function rho.
    rho(:, k) = radius ./ sqrt(sum(Vn(:, :, k).^2, 2));
end

% Get quadrature rule.
[xi, w] = triquadrature(1);
xi = repmat(permute(xi, [3, 2, 1]), [m, 1, 1]);

% Evaluate constant function at quadrature points.
detphi = detpushforward(F, V, rho, xi);
f = ones(m, length(w)) .* sqrt(abs(detphi));

% Compute triangle area.
a = triangArea(F, V);

% Compute surface integral over function which is constant one.
vq = surfintegral(f, w, a);
assertTrue(isscalar(vq));

vl = triangIntegral(F, V, ones(m, 1), a);
assertTrue(isscalar(vl));

fprintf('\nDifference between linear and quadratic interpolation:\n');
fprintf('abs(vq-4pi*r^2): %f\n', abs(vq - 4*pi*radius^2));
fprintf('abs(vl-4pi*r^2): %f\n', abs(vl - 4*pi*radius^2));

end