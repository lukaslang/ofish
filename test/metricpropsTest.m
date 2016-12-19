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
function tests = metricpropsTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function triangulationTest(testCase)

% Generate icosahedron.
[F, V] = sphTriang;
m = size(F, 1);

% Compute tangent basis.
[Dx, Dy] = tritanBasis(F, V);
verifyFalse(testCase, isempty(Dx));
verifyFalse(testCase, isempty(Dy));

% Compute triangulated surface properties.
[g, detg, ginv] = metricprops(Dx, Dy);
verifyFalse(testCase, isempty(g));
verifyFalse(testCase, isempty(detg));
verifyFalse(testCase, isempty(ginv));
verifyEqual(testCase, size(g), [m, 2, 2]);
verifyEqual(testCase, size(detg), [m, 1]);
verifyEqual(testCase, size(ginv), [m, 2, 2]);

% Check if determinant is correct.
verifyEqual(testCase, detg, 4 * triangArea(F, V) .^ 2, 'AbsTol', 1e-15);

end