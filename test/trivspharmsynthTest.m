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
function tests = trivspharmsynthTest
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
xi = repmat([1/3, 1/3], [m, 1]);

N = 10;
% Create coefficients.
u = zeros(2*(N^2 + 2*N), 1);
% Compute vector spherical harmonics synthesis.
[Yx, Yy] = trivspharmsynth(1:N, F, V, u, xi);

% Check results.
verifyFalse(testCase, isempty(Yx));
verifyFalse(testCase, isempty(Yy));
verifyEqual(testCase, size(Yx), [m, 1]);
verifyEqual(testCase, size(Yy), [m, 1]);
verifyEqual(testCase, Yx, zeros(m, 1));
verifyEqual(testCase, Yy, zeros(m, 1));

end

function intervalTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Pick coordinates.
xi = repmat([1/3, 1/3], [m, 1]);

N = 3:10;
% Create coefficients.
u = zeros(2*(N(end)^2 + 2*N(end) - N(1)^2 + 1), 1);
% Compute vector spherical harmonics synthesis.
[Yx, Yy] = trivspharmsynth(N, F, V, u, xi);

% Check results.
verifyFalse(testCase, isempty(Yx));
verifyFalse(testCase, isempty(Yy));
verifyEqual(testCase, size(Yx), [m, 1]);
verifyEqual(testCase, size(Yy), [m, 1]);
verifyEqual(testCase, Yx, zeros(m, 1));
verifyEqual(testCase, Yy, zeros(m, 1));

N = 1:10;
% Create coefficients.
u = zeros(2*(N(end)^2 + 2*N(end) - N(1)^2 + 1), 1);
% Compute vector spherical harmonics synthesis.
[Yx, Yy] = trivspharmsynth(N, F, V, u, xi);

% Check results.
verifyFalse(testCase, isempty(Yx));
verifyFalse(testCase, isempty(Yy));
verifyEqual(testCase, size(Yx), [m, 1]);
verifyEqual(testCase, size(Yy), [m, 1]);
verifyEqual(testCase, Yx, zeros(m, 1));
verifyEqual(testCase, Yy, zeros(m, 1));

end

function memConstraintTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Pick coordinates.
xi = repmat([1/3, 1/3], [m, 1]);

N = 30;
% Create coefficients.
u = zeros(2*(N^2 + 2*N), 1);
% Compute vector spherical harmonics synthesis.
mem = 2e6;
[Yx, Yy] = trivspharmsynth(1:N, F, V, u, xi, mem);

% Check results.
verifyFalse(testCase, isempty(Yx));
verifyFalse(testCase, isempty(Yy));
verifyEqual(testCase, size(Yx), [m, 1]);
verifyEqual(testCase, size(Yy), [m, 1]);
verifyEqual(testCase, Yx, zeros(m, 1));
verifyEqual(testCase, Yy, zeros(m, 1));

end

function compareCoefficientsTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Pick coordinates.
xi = repmat([1/3, 1/3], [m, 1]);

N = 5;
dim = 2*(N^2 + 2*N);
% Create coefficients.
u = ones(dim, 1);

% Compute vector spherical harmonics synthesis.
[Yx, Yy] = trivspharmsynth(1:N, F, V, u, xi);

% Compute coefficients.
Y = trivspharmncoeff(1:N, F, V, xi);
verifyEqual(testCase, size(Y), [m, dim, 2]);

% Check for equality of synthesis.
verifyEqual(testCase, Yx, sum(Y(:, :, 1), 2), 'AbsTol', 1e-13);
verifyEqual(testCase, Yy, sum(Y(:, :, 2), 2), 'AbsTol', 1e-13);

end