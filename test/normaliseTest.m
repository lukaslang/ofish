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
function tests = normaliseTest
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
[~, V] = sphTriang;

Vn = normalise(V);
verifyEqual(testCase, Vn, V);

Vn = normalise(bsxfun(@times, V, rand(size(V, 1), 1)));
verifyEqual(testCase, Vn, V, 'AbsTol', 1e-15);

end

function quadratureDimensionTest(testCase)

% Generate icosahedron.
[~, V] = sphTriang;

nq = 5;
Vn = normalise(repmat(V, [1, 1, nq]));
verifyEqual(testCase, Vn, repmat(V, [1, 1, nq]));

end