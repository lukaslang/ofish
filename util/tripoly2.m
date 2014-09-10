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
function [A, Q] = tripoly2(xi)
%TRIPOLY2 Computes piecewise quadratic functions on triangles.
%
%   [A, Q] = TRIPOLY2(xi) takes barycentric coordinates xi and returns 
%   coefficients A of quadratic polynomials and all monomials Q at xi.
%
%   Note that xi must be a matrix of size m-by-2, where m is the number of
%   triangular faces.
%
%   A is a matrix of size 6-by-6 and Q is a matrix of size m-by-6.
m = size(xi, 1);

% Create coefficients of basis functions.
A = [1,-3,-3, 2, 4, 2;
     0, 0,-1, 0, 0, 2;
     0,-1, 0, 2, 0, 0;
     0, 0, 4, 0,-4,-4;
     0, 0, 0, 0, 4, 0;
     0, 4, 0,-4,-4, 0];
 
% Define all monomials up to degree two at xi.
Q = [ones(m, 1), xi(:, 1), xi(:, 2), xi(:, 1).^2, xi(:, 1).*xi(:, 2), xi(:, 2).^2]';

end