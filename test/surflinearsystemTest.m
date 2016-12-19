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
function tests = surflinearsystemTest
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

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;

% Create random data on nodal points.
f1 = randi(255, [m, 6]);
f2 = randi(255, [m, 6]);

% Set parameters.
h = 1;
tol = 1e-6;

% Create linear system.
N = 5;
deg = 7;
[dim, A, D, b] = surflinearsystem(F, V, Ns, c, 1:N, f1, f2, h, deg, tol);

% Check results.
verifyEqual(testCase, dim, 2*(N^2 + 2*N));
verifyFalse(testCase, isempty(A));
verifyEqual(testCase, size(A), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
verifyFalse(testCase, isempty(D));
verifyEqual(testCase, size(D), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
verifyFalse(testCase, isempty(b));
verifyTrue(testCase, isvector(b));
verifyEqual(testCase, size(b), [2*(N^2 + 2*N), 1]);

% Check if matrices are symmetric.
verifyEqual(testCase, A, A');
verifyEqual(testCase, D, D');

end

function intervalTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;

% Create random data on nodal points.
f1 = randi(255, [m, 6]);
f2 = randi(255, [m, 6]);

% Set parameters.
h = 1;
tol = 1e-6;

% Create linear system.
N = 3:7;
deg = 1;
[dim, A, D, b] = surflinearsystem(F, V, Ns, c, N, f1, f2, h, deg, tol);

% Check results.
expDim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);
verifyEqual(testCase, dim, expDim);
verifyFalse(testCase, isempty(A));
verifyEqual(testCase, size(A), [expDim, expDim]);
verifyFalse(testCase, isempty(D));
verifyEqual(testCase, size(D), [expDim, expDim]);
verifyFalse(testCase, isempty(b));
verifyTrue(testCase, isvector(b));
verifyEqual(testCase, size(b), [expDim, 1]);

% Check if matrices are symmetric.
verifyEqual(testCase, A, A');
verifyEqual(testCase, D, D');

end

function memTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;

% Create random data on nodal points.
f1 = randi(255, [m, 6]);
f2 = randi(255, [m, 6]);

% Set parameters.
h = 1;
tol = 1e-6;

% Create linear system.
N = 3:7;
mem = 1e9;
deg = 1;
[dim, A, D, b] = surflinearsystem(F, V, Ns, c, N, f1, f2, h, deg, tol, mem);

% Check results.
expDim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);
verifyEqual(testCase, dim, expDim);
verifyFalse(testCase, isempty(A));
verifyEqual(testCase, size(A), [expDim, expDim]);
verifyFalse(testCase, isempty(D));
verifyEqual(testCase, size(D), [expDim, expDim]);
verifyFalse(testCase, isempty(b));
verifyTrue(testCase, isvector(b));
verifyEqual(testCase, size(b), [expDim, 1]);

% Check if matrices are symmetric.
verifyEqual(testCase, A, A');
verifyEqual(testCase, D, D');

end

function noDataTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;

% Create random data on vertices.
f1 = zeros(m, 6);
f2 = zeros(m, 6);

% Set parameters.
h = 1;
tol = 1e-6;

% Create linear system.
N = 5;
deg = 1;
[dim, A, D, b] = surflinearsystem(F, V, Ns, c, 1:N, f1, f2, h, deg, tol);

% Check results.
verifyEqual(testCase, dim, 2*(N^2 + 2*N));
verifyFalse(testCase, isempty(A));
verifyEqual(testCase, size(A), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
verifyFalse(testCase, isempty(D));
verifyEqual(testCase, size(D), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
verifyFalse(testCase, isempty(b));
verifyTrue(testCase, isvector(b));
verifyEqual(testCase, size(b), [2*(N^2 + 2*N), 1]);

% Check if matrices are symmetric.
verifyEqual(testCase, A, A');
verifyEqual(testCase, D, D');

end

function partiallyNoDataTest(testCase)

% Create triangulation of unit sphere.
[F, V] = sphTriang(3);
m = size(F, 1);

% Set parameters for unit sphere.
Ns = 0;
Y = spharm(Ns, [0, 0, 1]);
c = 1 / Y;

% Create random data on vertices.
f1 = [randi(255, [m/2, 6]); zeros(m/2, 6)];
f2 = zeros(m, 6);

% Set parameters.
h = 1;
tol = 1e-6;

% Create linear system.
N = 5;
deg = 1;
[dim, A, D, b] = surflinearsystem(F, V, Ns, c, 1:N, f1, f2, h, deg, tol);

% Check results.
verifyEqual(testCase, dim, 2*(N^2 + 2*N));
verifyFalse(testCase, isempty(A));
verifyEqual(testCase, size(A), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
verifyFalse(testCase, isempty(D));
verifyEqual(testCase, size(D), [2*(N^2 + 2*N), 2*(N^2 + 2*N)]);
verifyFalse(testCase, isempty(b));
verifyTrue(testCase, isvector(b));
verifyEqual(testCase, size(b), [2*(N^2 + 2*N), 1]);

% Check if matrices are symmetric.
verifyEqual(testCase, A, A');
verifyEqual(testCase, D, D');

end