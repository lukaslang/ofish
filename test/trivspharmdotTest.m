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
function tests = trivspharmdotTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Compute triangle areas.
a = triangArea(F, V);

% Pick coordinates.
xi = repmat([1/3, 1/3], m, 1);

% Arbitrary vector spherical harmonics.
deg = 1;
ord = 1;
[Y1, Y2] = trivspharm(deg, F, V, xi);
Y1nj = squeeze(Y1(:, ord, :));
Y2nj = squeeze(Y2(:, ord, :));

% Compute dot product with vector spherical harmonics of degrees N.
N = 1;
Z = trivspharmdot(Y1nj, F, V, N, xi);
verifyFalse(testCase, isempty(Z));
verifyEqual(testCase, size(Z), [m, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1)]);
verifyEqual(testCase, triangIntegral(F, V, Z(:, 1), a), 1, 'AbsTol', 1e-2);
for k=2:size(Z, 2)
    verifyEqual(testCase, triangIntegral(F, V, Z(:, k), a), 0, 'AbsTol', 1e-2);
end

% Compute dot product with vector spherical harmonics of degrees N.
Z = trivspharmdot(Y2nj, F, V, N, xi);
verifyFalse(testCase, isempty(Z));
verifyEqual(testCase, size(Z), [m, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1)]);
verifyEqual(testCase, triangIntegral(F, V, Z(:, 4), a), 1, 'AbsTol', 1e-2);
for k=[1:3, 5:6]
    verifyEqual(testCase, triangIntegral(F, V, Z(:, k), a), 0, 'AbsTol', 1e-2);
end

end

function quadratureTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Compute triangle areas.
a = triangArea(F, V);

% Pick coordinates.
nq = 6;
xi = repmat([1/3, 1/3], [m, 1, nq]);

% Arbitrary vector spherical harmonics.
deg = 1;
ord = 1;
[Y1, Y2] = trivspharm(deg, F, V, xi);
Y1nj = squeeze(Y1(:, ord, :, :));
Y2nj = squeeze(Y2(:, ord, :, :));

% Compute dot product with vector spherical harmonics of degrees N.
N = 1;
Z = trivspharmdot(Y1nj, F, V, N, xi);
verifyFalse(testCase, isempty(Z));
verifyEqual(testCase, size(Z), [m, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1), nq]);
verifyEqual(testCase, triangIntegral(F, V, Z(:, 1), a), 1, 'AbsTol', 1e-2);
for k=2:size(Z, 2)
    for q=1:nq
        verifyEqual(testCase, triangIntegral(F, V, Z(:, k, q), a), 0, 'AbsTol', 1e-2);
    end
end

% Compute dot product with vector spherical harmonics of degrees N.
Z = trivspharmdot(Y2nj, F, V, N, xi);
verifyFalse(testCase, isempty(Z));
verifyEqual(testCase, size(Z), [m, 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1), nq]);
verifyEqual(testCase, triangIntegral(F, V, Z(:, 4), a), 1, 'AbsTol', 1e-2);
for k=[1:3, 5:6]
    for q=1:nq
        verifyEqual(testCase, triangIntegral(F, V, Z(:, k, q), a), 0, 'AbsTol', 1e-2);
    end
end

end