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
function test_suite = metricpropsTest
    initTestSuite;
end

function triangulationTest

% Generate icosahedron.
[F, V] = sphTriang;
m = size(F, 1);

% Compute tangent basis.
[Dx, Dy] = tritanBasis(F, V);
assertFalse(isempty(Dx));
assertFalse(isempty(Dy));

% Compute triangulated surface properties.
[g, detg, ginv] = metricprops(Dx, Dy);
assertFalse(isempty(g));
assertFalse(isempty(detg));
assertFalse(isempty(ginv));
assertEqual(size(g), [m, 2, 2]);
assertEqual(size(detg), [m, 1]);
assertEqual(size(ginv), [m, 2, 2]);

% Check if determinant is correct.
assertAlmostEqual(detg, 4 * triangArea(F, V) .^ 2);

end