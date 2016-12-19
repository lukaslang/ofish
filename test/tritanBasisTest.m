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
function tests = tritanBasisTest
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
m = size(F, 1);

% Compute triangulated surface properties.
[Dx, Dy] = tritanBasis(F, V);
verifyFalse(testCase, isempty(Dx));
verifyFalse(testCase, isempty(Dy));
verifyEqual(testCase, size(Dx), [m, 3]);
verifyEqual(testCase, size(Dy), [m, 3]);

end

function visualiseTest(testCase)

% Generate spherical triangulation.
[F, V] = sphTriang(3);
verifyFalse(testCase, isempty(F));
verifyFalse(testCase, isempty(V));

% Compute triangulated surface properties.
[Dx, Dy] = tritanBasis(F, V);

% Compute points on surface.
xi = repmat([1/6, 1/6], size(F, 1), 1);
x = trimap(F, V, xi);

% Plot surface.
figure;
hold on;
trimesh(TriRep(F, V));
daspect([1, 1, 1]);
view(3);

% Plot tangent basis.
quiver3(x(:, 1), x(:, 2), x(:, 3), Dx(:, 1), Dx(:, 2), Dx(:, 3), 1);
quiver3(x(:, 1), x(:, 2), x(:, 3), Dy(:, 1), Dy(:, 2), Dy(:, 3), 1);

end