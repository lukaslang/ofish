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
function tests = trivspharmcoeffTest
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

% Pick coordinates.
xi = repmat([1/3, 1/3], m, 1);

% Create spherical harmonics.
N = 5;
[Y1, Y2] = trivspharmcoeff(N, F, V, xi);
verifyFalse(testCase, isempty(Y1));
verifyEqual(testCase, size(Y1), [m, 2*N + 1, 2]);
verifyFalse(testCase, isempty(Y2));
verifyEqual(testCase, size(Y2), [m, 2*N + 1, 2]);

end

function quadratureDimensionTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Pick coordinates.
nq = 5;
xi = repmat([1/3, 1/3], [m, 1, nq]);

% Create spherical harmonics.
N = 5;
[Y1, Y2] = trivspharmcoeff(N, F, V, xi);
verifyFalse(testCase, isempty(Y1));
verifyEqual(testCase, size(Y1), [m, 2*N + 1, 2, nq]);
verifyFalse(testCase, isempty(Y2));
verifyEqual(testCase, size(Y2), [m, 2*N + 1, 2, nq]);

end

function orthogonalityTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Pick coordinates.
xi = repmat([1/3, 1/3], m, 1);

% Create spherical harmonics.
N = 5;
[Y1c, Y2c] = trivspharmcoeff(N, F, V, xi);

% Compute tangent basis.
[Dx, Dy] = tritanBasis(F, V);

% Compute vectors in R3.
for k=1:size(Y1c, 2)
    Y1(:, k, :) = bsxfun(@times, Y1c(:, k, 1), Dx) + bsxfun(@times, Y1c(:, k, 2), Dy);
    Y2(:, k, :) = bsxfun(@times, Y2c(:, k, 1), Dx) + bsxfun(@times, Y2c(:, k, 2), Dy);
end

% Compute R3 inner product.
ip = dot(Y1, Y2, 3);
verifyEqual(testCase, ip, zeros(m, 2*N + 1), 'AbsTol', 1e-15);

end

function visualiseTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Pick coordinates.
xi = repmat([1/3, 1/3], m, 1);
x = trimap(F, V, xi);

% Create spherical harmonics.
N = 3;
[Y1c, Y2c] = trivspharmcoeff(N, F, V, xi);
verifyFalse(testCase, isempty(Y1c));
verifyEqual(testCase, size(Y1c), [m, 2*N + 1, 2]);
verifyFalse(testCase, isempty(Y2c));
verifyEqual(testCase, size(Y2c), [m, 2*N + 1, 2]);

% Compute tangent basis.
[Dx, Dy] = tritanBasis(F, V);

% Compute vectors in R3.
for k=1:size(Y1c, 2)
    Y1(:, k, :) = bsxfun(@times, Y1c(:, k, 1), Dx) + bsxfun(@times, Y1c(:, k, 2), Dy);
    Y2(:, k, :) = bsxfun(@times, Y2c(:, k, 1), Dx) + bsxfun(@times, Y2c(:, k, 2), Dy);
end

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