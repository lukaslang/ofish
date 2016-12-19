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
function tests = trivspharmTest
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
n = size(F, 1);

% Pick coordinates.
xi = repmat([1/3, 1/3], n, 1);

% Create spherical harmonics.
N = 5;
[Y1, Y2] = trivspharm(N, F, V, xi);
verifyFalse(testCase, isempty(Y1));
verifyEqual(testCase, size(Y1), [n, 2*N + 1, 3]);
verifyFalse(testCase, isempty(Y2));
verifyEqual(testCase, size(Y2), [n, 2*N + 1, 3]);

% Pick coordinates.
nq = 3;
xi = repmat([1/3, 1/3], [n, 1, nq]);

% Create spherical harmonics.
N = 5;
[Y1, Y2] = trivspharm(N, F, V, xi);
verifyFalse(testCase, isempty(Y1));
verifyEqual(testCase, size(Y1), [n, 2*N + 1, 3, nq]);
verifyFalse(testCase, isempty(Y2));
verifyEqual(testCase, size(Y2), [n, 2*N + 1, 3, nq]);

end

function orthogonalityTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
n = size(F, 1);

% Pick coordinates.
xi = repmat([1/3, 1/3], n, 1);

% Create spherical harmonics.
N = 5;
[Y1, Y2] = trivspharm(N, F, V, xi);

% Compute R3 inner product.
ip = dot(Y1, Y2, 3);
verifyEqual(testCase, ip, zeros(n, 2*N + 1), 'AbsTol', 1e-15);

end

function visualiseTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
n = size(F, 1);

% Pick coordinates.
xi = repmat([1/3, 1/3], n, 1);
x = trimap(F, V, xi);

% Create spherical harmonics.
N = 3;
[Y1, Y2] = trivspharm(N, F, V, xi);
verifyFalse(testCase, isempty(Y1));
verifyEqual(testCase, size(Y1), [n, 2*N + 1, 3]);
verifyFalse(testCase, isempty(Y2));
verifyEqual(testCase, size(Y2), [n, 2*N + 1, 3]);

% Create spherical harmonics for visualisation.
Ynj = spharm(N, V);

figure;
for k=1:2*N+1
    f = Ynj(:, k);
    subplot(1, 2*N+1, k);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f);
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    % Plot vector field.
    quiver3(x(:, 1), x(:, 2), x(:, 3), Y1(:, k, 1), Y1(:, k, 2), Y1(:, k, 3), 1, 'r');
    quiver3(x(:, 1), x(:, 2), x(:, 3), Y2(:, k, 1), Y2(:, k, 2), Y2(:, k, 3), 1, 'b');
end
end