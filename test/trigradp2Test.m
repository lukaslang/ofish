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
function test_suite = trigradp2Test
    initTestSuite;
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
assertFalse(isempty(F));
assertFalse(isempty(V));
assertEqual(size(F), [20, 3]);
assertEqual(size(V), [12, 3]);

% Create sample function.
f = ones(20, 6);

% Define centroids.
xi = repmat([1/3, 1/3], 20, 1);

% Compute gradient.
g = trigradp2(F, V, f, xi);
assertFalse(isempty(g));
assertEqual(size(g), [20, 3]);

% Define centroids.
nq = 3;
xi = repmat([1/3, 1/3], [20, 1, nq]);

% Compute gradient.
g = trigradp2(F, V, f, xi);
assertFalse(isempty(g));
assertEqual(size(g), [20, 3, nq]);

end

function visualiseTest(testCase)

% Generate spherical triangulation.
[F, V] = sphTriang(3);
assertFalse(isempty(F));
assertFalse(isempty(V));

% Create quadrature points for each face.
Vn = trinodalpts2(F, V);

% Create spherical harmonics as test functions and project to surface.
N = 3;
for k=1:6
    fq(:, :, k) = spharm(N, normalise(Vn(:, :, k)));
end

% Compute points on surface.
xi = repmat([1/3, 1/3], size(F, 1), 1);
x = trimap(F, V, xi);

f = spharm(N, V);
figure;
for k=1:2*N+1
    % Compute gradient of k-th harmonic.
    d = squeeze(fq(:, k, :));
    g = trigradp2(F, V, d, xi);
    assertFalse(isempty(g));
    
    % Plot spherical harmonics.
    subplot(1, 2 * N + 1, k);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f(:, k));
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    
    % Plot gradient field.
    quiver3(x(:, 1), x(:, 2), x(:, 3), g(:, 1), g(:, 2), g(:, 3), 1);
end
end