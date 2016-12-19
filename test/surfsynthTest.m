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
function tests = surfsynthTest
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
[~, V] = sphTriang(5);
n = size(V, 1);

% Set parameters.
Ns = 0;
% Compute coefficient such that rho is identically one, i.e. unit sphere.
c = 1/spharm(0, [0, 1, 0]);

% Compute synthesis.
[S, rho] = surfsynth(Ns, V, c);

% Check if points are still on unit sphere.
len = sqrt(sum(S .^ 2, 2));
verifyEqual(testCase, len, ones(n, 1), 'AbsTol', 1e-15);

% Check if points equal.
verifyEqual(testCase, S, V);

% Check if rho is identically one.
verifyEqual(testCase, rho, ones(n, 1), 'AbsTol', 1e-15);

end

function quadratureDimensionTest(testCase)

% Create triangulation of unit sphere.
[~, V] = sphTriang(5);
n = size(V, 1);

% Set parameters.
Ns = 0;
% Compute coefficient such that rho is identically one, i.e. unit sphere.
c = 1/spharm(0, [0, 1, 0]);

% Compute synthesis.
nq = 3;
[S, rho] = surfsynth(Ns, repmat(V, [1, 1, nq]), c);

% Check if points are still on unit sphere.
len = sqrt(sum(S .^ 2, 2));
verifyEqual(testCase, len, ones(n, 1, nq), 'AbsTol', 1e-15);

% Check if points equal.
verifyEqual(testCase, S, repmat(V, [1, 1, nq]), 'AbsTol', 1e-15);

% Check if rho is identically one.
verifyEqual(testCase, rho, ones(n, nq), 'AbsTol', 1e-15);

end

function resultIntervalTest(testCase)

% Create triangulation of unit sphere.
[~, V] = sphTriang(5);
n = size(V, 1);

% Set parameters.
Ns = 0:5;
dim = Ns(end)^2 + 2*Ns(end) - Ns(1)^2 + 1;
% Compute coefficient such that rho is identically one, i.e. unit sphere.
c = [1/spharm(0, [0, 1, 0]); zeros(dim - 1, 1)];

% Compute synthesis.
[S, rho] = surfsynth(Ns, V, c);

% Check if points are still on unit sphere.
len = sqrt(sum(S .^ 2, 2));
verifyEqual(testCase, len, ones(n, 1), 'AbsTol', 1e-15);

% Check if points equal.
verifyEqual(testCase, S, V);

% Check if rho is identically one.
verifyEqual(testCase, rho, ones(n, 1), 'AbsTol', 1e-15);

end

function sphereTest(testCase)

% Create triangulation of unit sphere.
[~, V] = sphTriang(5);
n = size(V, 1);

% Set parameters.
Ns = 0:5;
s = 1;
alpha = 1;

% Fit surface to sphere.
[c, Y] = surffit(Ns, V, alpha, s);
verifyFalse(testCase, isempty(c));
verifyFalse(testCase, isempty(Y));

% Compute synthesis from scratch.
[S, rho] = surfsynth(Ns, V, c);

% Check if points are still on unit sphere.
len = sqrt(sum(S .^ 2, 2));
verifyEqual(testCase, len, ones(n, 1), 'AbsTol', 1e-12);

% Check if points equal.
verifyEqual(testCase, S, V, 'AbsTol', 1e-12);

% Check if functions rho are identically one.
verifyEqual(testCase, rho, ones(n, 1), 'AbsTol', 1e-12);

end