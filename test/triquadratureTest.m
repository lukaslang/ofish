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
function tests = triquadratureTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    cd('../');
end

function teardownOnce(testCase)
    cd('test');
end

function resultTest(testCase)

% Get quadrature points.
[xi, w] = triquadrature(1);
verifyFalse(testCase, isempty(xi));
verifyFalse(testCase, isempty(w));
verifyEqual(testCase, size(xi), [1, 2]);
verifyEqual(testCase, size(w), [1, 1]);

% Get quadrature points.
[xi, w] = triquadrature(7);
verifyFalse(testCase, isempty(xi));
verifyFalse(testCase, isempty(w));
verifyEqual(testCase, size(xi), [16, 2]);
verifyEqual(testCase, size(w), [16, 1]);

% Get quadrature points.
[xi, w] = triquadrature(15);
verifyFalse(testCase, isempty(xi));
verifyFalse(testCase, isempty(w));
verifyEqual(testCase, size(xi), [64, 2]);
verifyEqual(testCase, size(w), [64, 1]);

end