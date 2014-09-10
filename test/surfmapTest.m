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
function test_suite = surfmapTest
    initTestSuite;
end

function resultTest

% Generate icosahedron.
[F, V] = sphTriang;
m = size(F, 1);

% Check with constant function rho.
rho = ones(m, 6);

x = surfmap(F, V, rho, zeros(m, 2));
assertFalse(isempty(x));
assertEqual(size(x), [m, 3]);
assertAlmostEqual(x, V(F(:, 1), :));

x = surfmap(F, V, rho, [ones(m, 1), zeros(m, 1)]);
assertFalse(isempty(x));
assertEqual(size(x), [m, 3]);
assertAlmostEqual(x, V(F(:, 3), :));

x = surfmap(F, V, rho, [zeros(m, 1), ones(m, 1)]);
assertFalse(isempty(x));
assertEqual(size(x), [m, 3]);
assertAlmostEqual(x, V(F(:, 2), :));

end

function sphereWithRadiusTwoTest

% Generate icosahedron.
[F, V] = sphTriang;
m = size(F, 1);

% Check with constant function rho.
rho = 2 * ones(m, 6);

x = surfmap(F, V, rho, zeros(m, 2));
assertFalse(isempty(x));
assertEqual(size(x), [m, 3]);
assertAlmostEqual(x, 2 * V(F(:, 1), :));

x = surfmap(F, V, rho, [ones(m, 1), zeros(m, 1)]);
assertFalse(isempty(x));
assertEqual(size(x), [m, 3]);
assertAlmostEqual(x, 2 * V(F(:, 3), :));

x = surfmap(F, V, rho, [zeros(m, 1), ones(m, 1)]);
assertFalse(isempty(x));
assertEqual(size(x), [m, 3]);
assertAlmostEqual(x, 2 * V(F(:, 2), :));

end

function sphereWithRadiusTwoVisualiseTest

% Generate icosahedron.
[F, V] = sphTriang(3);
m = size(F, 1);

% Check with constant function rho.
rho = 2 * ones(m, 6);

x = surfmap(F, V, rho, ones(m, 2) ./ 3);
assertFalse(isempty(x));
assertEqual(size(x), [m, 3]);
assertAlmostEqual(sqrt(sum(x.^2, 2)), 2*ones(m, 1), 0.01);

% Plot surface.
figure;
hold on;
trimesh(TriRep(F, V));
shading interp;
daspect([1, 1, 1]);
view(3);
quiver3(zeros(m, 1), zeros(m, 1), zeros(m, 1), x(:, 1), x(:, 2), x(:, 3));

end