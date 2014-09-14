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
function test_suite = triinterp2Test
    initTestSuite;
end

function resultTest

% Generate icosahedron.
[F, V] = sphTriang;
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [20, 3]);
assertEqual(size(V), [12, 3]);

% Create sample function.
f = ones(20, 6);

% Define centroids in barycentric coordinates.
xi = repmat([1/3, 1/3], 20, 1);

% Compute interpolation at xi.
g = triinterp2(f, xi);
assertFalse(isempty(g));
assertEqual(size(g), [20, 1]);
assertAlmostEqual(g, ones(20, 1));

% Check all vertices.
xi = repmat([0, 0], 20, 1);
% Compute interpolation at xi.
g = triinterp2(f, xi);
assertAlmostEqual(g, ones(20, 1));

xi = repmat([1, 0], 20, 1);
% Compute interpolation at xi.
g = triinterp2(f, xi);
assertAlmostEqual(g, ones(20, 1));

xi = repmat([0, 1], 20, 1);
% Compute interpolation at xi.
g = triinterp2(f, xi);
assertAlmostEqual(g, ones(20, 1));

end

function thirdDimensionTest

% Generate icosahedron.
[F, V] = sphTriang;
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [20, 3]);
assertEqual(size(V), [12, 3]);

% Create sample function.
f = ones(20, 6);

% Define centroids in barycentric coordinates.
xi = repmat([1/3, 1/3], [20, 1, 1]);

% Compute interpolation at xi.
g = triinterp2(f, xi);
assertFalse(isempty(g));
assertEqual(size(g), [20, 1]);
assertAlmostEqual(g, ones(20, 1));

xi(:, :, 1) = repmat([0, 0], 20, 1);
xi(:, :, 2) = repmat([1, 0], 20, 1);
xi(:, :, 3) = repmat([0, 1], 20, 1);

% Compute interpolation at xi.
g = triinterp2(f, xi);
assertAlmostEqual(g(:, 1), ones(20, 1));
assertAlmostEqual(g(:, 2), ones(20, 1));
assertAlmostEqual(g(:, 3), ones(20, 1));

end