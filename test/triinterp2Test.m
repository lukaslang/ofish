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
function tests = triinterp2Test
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Generate icosahedron.
[F, V] = sphTriang;
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
verifyEqual(testCase, size(F), [20, 3]);
verifyEqual(testCase, size(V), [12, 3]);

% Create sample function.
f = ones(20, 6);

% Define centroids in barycentric coordinates.
xi = repmat([1/3, 1/3], 20, 1);

% Compute interpolation at xi.
g = triinterp2(f, xi);
verifyFalse(testCase, isempty(g));
verifyEqual(testCase, size(g), [20, 1]);
verifyEqual(testCase, g, ones(20, 1), 'AbsTol', 1e-15);

% Check all vertices.
xi = repmat([0, 0], 20, 1);
% Compute interpolation at xi.
g = triinterp2(f, xi);
verifyEqual(testCase, g, ones(20, 1));

xi = repmat([1, 0], 20, 1);
% Compute interpolation at xi.
g = triinterp2(f, xi);
verifyEqual(testCase, g, ones(20, 1));

xi = repmat([0, 1], 20, 1);
% Compute interpolation at xi.
g = triinterp2(f, xi);
verifyEqual(testCase, g, ones(20, 1));

end

function thirdDimensionTest(testCase)

% Generate icosahedron.
[F, V] = sphTriang;
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));
verifyEqual(testCase, size(F), [20, 3]);
verifyEqual(testCase, size(V), [12, 3]);

% Create sample function.
f = ones(20, 6);

% Define centroids in barycentric coordinates.
xi = repmat([1/3, 1/3], [20, 1, 1]);

% Compute interpolation at xi.
g = triinterp2(f, xi);
verifyFalse(testCase, isempty(g));
verifyEqual(testCase, size(g), [20, 1]);
verifyEqual(testCase, g, ones(20, 1), 'AbsTol', 1e-15);

xi(:, :, 1) = repmat([0, 0], 20, 1);
xi(:, :, 2) = repmat([1, 0], 20, 1);
xi(:, :, 3) = repmat([0, 1], 20, 1);

% Compute interpolation at xi.
g = triinterp2(f, xi);
verifyEqual(testCase, g(:, 1), ones(20, 1));
verifyEqual(testCase, g(:, 2), ones(20, 1));
verifyEqual(testCase, g(:, 3), ones(20, 1));

end

function singleFaceTest(testCase)

% Create sample function.
f = ones(1, 6);

% Define centroids in barycentric coordinates.
xi = repmat([1/3, 1/3], [1, 1, 6]);

% Compute interpolation at xi.
g = triinterp2(f, xi);
verifyFalse(testCase, isempty(g));
verifyEqual(testCase, size(g), [1, 6]);
verifyEqual(testCase, g, ones(1, 6), 'AbsTol', 1e-15);

xi(:, :, 1) = repmat([0, 0], 1, 1);
xi(:, :, 2) = repmat([1, 0], 1, 1);
xi(:, :, 3) = repmat([0, 1], 1, 1);

% Compute interpolation at xi.
g = triinterp2(f, xi);
verifyEqual(testCase, g(:, 1), ones(1, 1));
verifyEqual(testCase, g(:, 2), ones(1, 1));
verifyEqual(testCase, g(:, 3), ones(1, 1));

end