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
function tests = detpushforwardTest
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

xi = repmat([1/3, 1/3], m ,1);

% Check with constant function rho.
rho = ones(m, 6);
detphi = detpushforward(F, V, rho, xi);
verifyFalse(testCase, isempty(detphi));
verifyEqual(testCase, size(detphi), [m, 1]);
verifyEqual(testCase, detphi, ones(m, 1), 'AbsTol', 1e-15);

end

function quadratureDimTest(testCase)

% Generate icosahedron.
[F, V] = sphTriang;
m = size(F, 1);

nq = 4;
xi = repmat([1/3, 1/3], [m, 1, nq]);

% Check with constant function rho.
rho = ones(m, 6);
detphi = detpushforward(F, V, rho, xi);
verifyFalse(testCase, all(isempty(detphi)));
verifyEqual(testCase, size(detphi), [m, nq]);
verifyEqual(testCase, detphi, ones(m, nq), 'AbsTol', 1e-15);

end