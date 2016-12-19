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
function tests = surfcovderivTest
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
[F, V] = sphTriang(3);
m = size(F, 1);

% Generate constant function rho.
rho = ones(m, 6);

% Pick coordinates.
nq = 16;
xi = repmat([1/3, 1/3], [m, 1, nq]);

% Create a transformation matrix for orthonormal basis.
A(1, 1, :, :) = ones(m, 6);
A(1, 2, :, :) = zeros(m, 6);
A(2, 1, :, :) = zeros(m, 6);
A(2, 2, :, :) = ones(m, 6);

% Compute Christoffel symbols.
G = surfchristoffel(F, V, rho, xi);

% Compute metric.
parfor q=1:nq
    [Dx, Dy] = surftanBasis(F, V, rho, xi(:, :, q));
    g(:, :, :, q) = metricprops(Dx, Dy);
end

% Create vector fields U.
N = 5;
[Y, DY] = trivspharmncoeff(N, F, V, xi);
dim = size(Y, 2);

% Compute covariant derivative.
[Z11, Z12, Z21, Z22] = surfcovderiv(G, g, A, Y, DY, xi);
verifyEqual(testCase, size(Z11), [m, dim, nq]);
verifyEqual(testCase, size(Z12), [m, dim, nq]);
verifyEqual(testCase, size(Z21), [m, dim, nq]);
verifyEqual(testCase, size(Z22), [m, dim, nq]);

end

function visualiseTest(testCase)

% Generate icosahedron.
[F, V] = sphTriang(4);
m = size(F, 1);

Vn = trinodalpts2(F, V);
for k=1:6
    % Generate constant function rho.
    rho(:, k) = 1 ./ sqrt(sum(Vn(:, :, k).^2, 2));
end

% Get quadrature rule.
deg = 1;
[xiq, w] = triquadrature(deg);
xi = repmat(permute(xiq, [3, 2, 1]), [m, 1, 1]);
nq = length(w);

% Compute orthonormal basis at nodal points.
xin(:, :, 1) = repmat([0, 0], [m, 1]);
xin(:, :, 2) = repmat([0, 1], [m, 1]);
xin(:, :, 3) = repmat([1, 0], [m, 1]);
xin(:, :, 4) = repmat([0, 1/2], [m, 1]);
xin(:, :, 5) = repmat([1/2, 1/2], [m, 1]);
xin(:, :, 6) = repmat([1/2, 0], [m, 1]);
for k=1:6
    [Dx, Dy] = surftanBasis(F, V, rho, xin(:, :, k));
    A(:, :, :, k) = orthonormalise(Dx, Dy);
end

% Compute Christoffel symbols.
G = surfchristoffel(F, V, rho, xi);

% Compute metric.
parfor q=1:nq
    [Dx, Dy] = surftanBasis(F, V, rho, xi(:, :, q));
    g(:, :, :, q) = metricprops(Dx, Dy);
end

% Create vector fields U.
N = 2;
[Y, DY] = trivspharmncoeff(N, F, V, xi);

% Compute covariant derivative.
[Z11, Z12, Z21, Z22] = surfcovderiv(G, g, A, Y, DY, xi);

Ynj = spharm(N, V);
for k=1:2*N+1
    figure;
    % Plot spherical harmonics.
    subplot(1, 2, 1);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), Ynj(:, k), 'EdgeColor', 'none');
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    title('Scalar and vector spherical harmonics.');
    
    % Plot vspharm.
    x = surfmap(F, V, rho, xi);
    [Dx, Dy] = surftanBasis(F, V, rho, xi);
    v = bsxfun(@times, Y(:, k, 1), Dx) + bsxfun(@times, Y(:, k, 2), Dy);
    quiver3(x(:, 1), x(:, 2), x(:, 3), v(:, 1), v(:, 2), v(:, 3), 1, 'c');
    
    % Plot norm of covariant derivative.
    subplot(1, 2, 2);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), sqrt(Z11(:, k, :).^2 + Z12(:, k, :).^2 + Z21(:, k, :).^2 + Z22(:, k, :).^2), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    colorbar;
    title('Norm of covariant derivative');
        
    % Plot covariant derivatives.
    figure;
    subplot(1, 4, 1);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), Z11(:, k, :), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    colorbar;

    subplot(1, 4, 2);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), Z12(:, k, :), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    colorbar;

    subplot(1, 4, 3);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), Z21(:, k, :), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    colorbar;
    
    subplot(1, 4, 4);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), Z22(:, k, :), 'EdgeColor', 'none');
    daspect([1, 1, 1]);
    view(3);
    colorbar;
end
end

function visualiseCovDerivTest(testCase)

% Generate icosahedron.
[F, V] = sphTriang(4);
m = size(F, 1);

Vn = trinodalpts2(F, V);
for k=1:6
    % Generate constant function rho.
    rho(:, k) = 1 ./ sqrt(sum(Vn(:, :, k).^2, 2));
end


% Compute orthonormal basis at nodal points.
xin(:, :, 1) = repmat([0, 0], [m, 1]);
xin(:, :, 2) = repmat([0, 1], [m, 1]);
xin(:, :, 3) = repmat([1, 0], [m, 1]);
xin(:, :, 4) = repmat([0, 1/2], [m, 1]);
xin(:, :, 5) = repmat([1/2, 1/2], [m, 1]);
xin(:, :, 6) = repmat([1/2, 0], [m, 1]);
for k=1:6
    [Dx, Dy] = surftanBasis(F, V, rho, xin(:, :, k));
    [a, e, f] = orthonormalise(Dx, Dy);
    A(:, :, :, k) = a;
end

% Pick coordinates.
nq = 1;
xi = repmat([1/3, 1/3], [m, 1, nq]);

% Compute Christoffel symbols.
G = surfchristoffel(F, V, rho, xi);

% Compute metric.
parfor q=1:nq
    [Dx, Dy] = surftanBasis(F, V, rho, xi(:, :, q));
    g(:, :, :, q) = metricprops(Dx, Dy);
end

% Create vector fields U.
N = 2;
[Y, DY] = trivspharmncoeff(N, F, V, xi);

% Compute covariant derivative.
[Z11, Z12, Z21, Z22] = surfcovderiv(G, g, A, Y, DY, xi);

for k=1:2*N+1
    x = surfmap(F, V, rho, xi);
    figure;
    % Plot spherical harmonics.
    subplot(1, 2, 1);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), sqrt(Z11(:, k).^2 + Z12(:, k).^2), 'EdgeColor', 'black');
    %shading interp;
    daspect([1, 1, 1]);
    view(3);
    
    % Compute cov. derivative in ONB.
    [Dx, Dy] = surftanBasis(F, V, rho, xi);
    [~, e, f] = orthonormalise(Dx, Dy);
    u = bsxfun(@times, Z11(:, k), e) + bsxfun(@times, Z12(:, k), f);
    quiver3(x(:, 1), x(:, 2), x(:, 3), u(:, 1), u(:, 2), u(:, 3), 1, 'k');
    
    % Plot tangent basis.
    quiver3(x(:, 1), x(:, 2), x(:, 3), e(:, 1), e(:, 2), e(:, 3), 1, 'r');
    quiver3(x(:, 1), x(:, 2), x(:, 3), f(:, 1), f(:, 2), f(:, 3), 1, 'b');
    
    % Plot vspharm.
    v = bsxfun(@times, Y(:, k, 1), Dx) + bsxfun(@times, Y(:, k, 2), Dy);
    quiver3(x(:, 1), x(:, 2), x(:, 3), v(:, 1), v(:, 2), v(:, 3), 1, 'c');
    
    
    subplot(1, 2, 2);
    axis([-1, 1, -1, 1, -1, 1]);
    hold on;
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), sqrt(Z21(:, k).^2 + Z22(:, k).^2), 'EdgeColor', 'black');
    %shading interp;
    daspect([1, 1, 1]);
    view(3);
    
    % Compute cov. derivative in ONB.
    u = bsxfun(@times, Z21(:, k), e) + bsxfun(@times, Z22(:, k), f);
    quiver3(x(:, 1), x(:, 2), x(:, 3), u(:, 1), u(:, 2), u(:, 3), 1, 'k');
    
    % Plot tangent basis.
    quiver3(x(:, 1), x(:, 2), x(:, 3), e(:, 1), e(:, 2), e(:, 3), 1, 'r');
    quiver3(x(:, 1), x(:, 2), x(:, 3), f(:, 1), f(:, 2), f(:, 3), 1, 'b');
    
    % Plot vspharm.
    quiver3(x(:, 1), x(:, 2), x(:, 3), v(:, 1), v(:, 2), v(:, 3), 1, 'c');
    
end
end