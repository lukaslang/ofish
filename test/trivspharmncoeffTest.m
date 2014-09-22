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
function test_suite = trivspharmncoeffTest
    initTestSuite;
end

function resultTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Pick coordinates.
nq = 5;
xi = repmat([1/3, 1/3], [m, 1, nq]);

N = 1:5;
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

% Compute coefficients.
[Y, DY] = trivspharmncoeff(N, F, V, xi);
assertEqual(size(Y), [m, dim, 2, nq]);
assertEqual(size(DY), [m, 2, dim, 2, nq]);

end

function orthogonalityTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Pick coordinates.
xi = repmat([1/3, 1/3], m, 1);

% Create spherical harmonics.
N = 5;
Yc = trivspharmncoeff(N, F, V, xi);

% Compute tangent basis.
[Dx, Dy] = tritanBasis(F, V);

% Compute vectors in R3.
for k=1:size(Yc, 2)
    Y(:, k, :) = bsxfun(@times, Yc(:, k, 1), Dx) + bsxfun(@times, Yc(:, k, 2), Dy);
end

% Compute dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

Y1 = Y(:, 1:dim/2, :);
Y2 = Y(:, dim/2+1:end, :);

% Compute R3 inner product.
ip = dot(Y1, Y2, 3);
assertAlmostEqual(ip, zeros(m, 2*N + 1));

end

function visualiseCoefficientsTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Pick coordinates.
xi = repmat([1/3, 1/3], m, 1);
x = trimap(F, V, xi);

% Create spherical harmonics.
N = 3;
Yc = trivspharmncoeff(N, F, V, xi);

% Compute tangent basis.
[Dx, Dy] = tritanBasis(F, V);

% Compute vectors in R3.
for k=1:size(Yc, 2)
    Y(:, k, :) = bsxfun(@times, Yc(:, k, 1), Dx) + bsxfun(@times, Yc(:, k, 2), Dy);
end

% Compute dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

Y1 = Y(:, 1:dim/2, :);
Y2 = Y(:, dim/2+1:end, :);

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

function visualiseDerivativesTest

% Create triangulation of unit sphere.
[F, V] = sphTriang(5);
m = size(F, 1);

% Pick coordinates.
xi = repmat([1/3, 1/3], m, 1);
x = trimap(F, V, xi);

% Create spherical harmonics.
N = 2;
[~, DYc] = trivspharmncoeff(N, F, V, xi);

% Compute dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

DY1 = DYc(:, :, 1:dim/2, :);
DY2 = DYc(:, :, dim/2+1:end, :);

% Create spherical harmonics for visualisation.
Ynj = spharm(N, V);

% Visualise derivatives.
for k=1:2*N+1
    figure;
    f = Ynj(:, k);
    subplot(2, 5, 1);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f);
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    subplot(2, 5, 2);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), DY1(:, 1, k, 1), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    subplot(2, 5, 3);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), DY1(:, 1, k, 2), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    subplot(2, 5, 4);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), DY1(:, 2, k, 1), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    subplot(2, 5, 5);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), DY1(:, 2, k, 2), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    
    subplot(2, 5, 6);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f);
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    subplot(2, 5, 7);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), DY2(:, 1, k, 1), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    subplot(2, 5, 8);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), DY2(:, 1, k, 2), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    subplot(2, 5, 9);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), DY2(:, 2, k, 1), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    subplot(2, 5, 10);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), DY2(:, 2, k, 2), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
end
end