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
function test_suite = surftanBasisTest
    initTestSuite;
end

function resultTest

% Generate icosahedron.
[F, V] = sphTriang;

% Pick one triangle.
t = 1;
F = F(t, :);
m = size(F, 1);

% Generate constant function rho.
rho = ones(1, 6);

% Pick coordinates.
xi = repmat([0, 0], size(F, 1), 1);

% Compute triangulated surface properties.
[Dx, Dy] = surftanBasis(F, V, rho, xi);
assertFalse(isempty(Dx));
assertFalse(isempty(Dy));
assertEqual(size(Dx), [m, 3]);
assertEqual(size(Dy), [m, 3]);

end

function visualiseTest

% Generate triangle.
F = [1, 2, 3];
V = [0, 0, 1; 0, 1, 1; 1, 0, 1];

% Generate function rho.
rho = [1.9, 2, 2, 2.1, 1.9, 2.1];

% Plot surface.
figure;
hold on;
trimesh(TriRep(F, V));
daspect([1, 1, 1]);
view(3);

I = linspace(0, 1, 5);
for k=I
    I2 = linspace(0, 1-k, 5);
    for l=I2
    % Compute points on surface.
    xi = repmat([k, l], size(F, 1), 1);
    x = surfmap(F, V, rho, xi);

    % Compute triangulated surface properties.
    [Dx, Dy] = surftanBasis(F, V, rho, xi);

    % Plot tangent basis.
    quiver3(x(:, 1), x(:, 2), x(:, 3), Dx(:, 1), Dx(:, 2), Dx(:, 3), 0.1, 'r');
    quiver3(x(:, 1), x(:, 2), x(:, 3), Dy(:, 1), Dy(:, 2), Dy(:, 3), 0.1, 'b');
    end
end
end