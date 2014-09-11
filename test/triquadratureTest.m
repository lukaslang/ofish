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
function test_suite = triquadratureTest
    initTestSuite;
end

function resultTest

% Get quadrature points.
[xi, w] = triquadrature(1);
assertFalse(isempty(xi));
assertFalse(isempty(w));
assertEqual(size(xi), [1, 2]);
assertEqual(size(w), [1, 1]);

% Get quadrature points.
[xi, w] = triquadrature(7);
assertFalse(isempty(xi));
assertFalse(isempty(w));
assertEqual(size(xi), [16, 2]);
assertEqual(size(w), [16, 1]);

% Get quadrature points.
[xi, w] = triquadrature(15);
assertFalse(isempty(xi));
assertFalse(isempty(w));
assertEqual(size(xi), [64, 2]);
assertEqual(size(w), [64, 1]);

end