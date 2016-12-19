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
function tests = trimapTest
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

x = trimap(F, V, zeros(m, 2));
verifyFalse(testCase, isempty(x));
verifyEqual(testCase, size(x), [m, 3]);
verifyEqual(testCase, x, V(F(:, 1), :));

x = trimap(F, V, [ones(m, 1), zeros(m, 1)]);
verifyFalse(testCase, isempty(x));
verifyEqual(testCase, size(x), [m, 3]);
verifyEqual(testCase, x, V(F(:, 3), :));

x = trimap(F, V, [zeros(m, 1), ones(m, 1)]);
verifyFalse(testCase, isempty(x));
verifyEqual(testCase, size(x), [m, 3]);
verifyEqual(testCase, x, V(F(:, 2), :));

end

function quadratureDimTest(testCase)

% Generate icosahedron.
[F, V] = sphTriang;
m = size(F, 1);

xi(:, :, 1) = zeros(m, 2);
xi(:, :, 2) = [ones(m, 1), zeros(m, 1)];
xi(:, :, 3) = [zeros(m, 1), ones(m, 1)];
x = trimap(F, V, xi);
verifyFalse(testCase, isempty(x));
verifyEqual(testCase, size(x), [m, 3, 3]);
verifyEqual(testCase, x(:, :, 1), V(F(:, 1), :));
verifyEqual(testCase, x(:, :, 2), V(F(:, 3), :));
verifyEqual(testCase, x(:, :, 3), V(F(:, 2), :));

end