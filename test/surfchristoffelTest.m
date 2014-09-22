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
function test_suite = surfchristoffelTest
    initTestSuite;
end

function resultTest

% Generate icosahedron.
[F, V] = sphTriang;
m = size(F, 1);

% Generate constant function rho.
rho = ones(m, 6);

% Pick coordinates.
nq = 5;
xi = repmat([1/3, 1/3], [m, 1, nq]);

% Compute Christoffel symbols.
G = surfchristoffel(F, V, rho, xi);
assertFalse(isempty(G));
assertEqual(size(G), [m, 2, 2, 2, nq]);
assertAlmostEqual(G, zeros(m, 2, 2, 2, nq));

end