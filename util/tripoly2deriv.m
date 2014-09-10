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
function [DxA, DyA, Q] = tripoly2deriv(xi)
%TRIPOLY2DERIV Computes partial derivatives of piecewise quadratic 
%functions on triangles.
%
%   [DxA, DyA, Q] = TRIPOLY2DERIV(xi) takes barycentric coordinates xi and 
%   returns coefficients DxA, DyA of partial derivatives and all monomials 
%   Q at xi.
%
%   Note that xi must be a matrix of size m-by-2, where m is the number of
%   triangular faces.
%
%   DxA, DyA are matrices of size 6-by-6 and Q a matrix of size m-by-6.
m = size(xi, 1);

% Create coefficients of derivatives wrt. x.
DxA = [-3, 4, 4, 0, 0, 0;
        0, 0, 0, 0, 0, 0;
       -1, 4, 0, 0, 0, 0;
        0, 0,-4, 0, 0, 0;
        0, 0, 4, 0, 0, 0;
        4,-8,-4, 0, 0, 0];

% Create coefficients of derivatives wrt. y.
DyA = [-3, 4, 4, 0, 0, 0;
       -1, 0, 4, 0, 0, 0;
        0, 0, 0, 0, 0, 0;
        4,-4,-8, 0, 0, 0;
        0, 4, 0, 0, 0, 0;
        0,-4, 0, 0, 0, 0];
 
% Define all monomials up to degree two at xi.
Q = [ones(m, 1), xi(:, 1), xi(:, 2), xi(:, 1).^2, xi(:, 1).*xi(:, 2), xi(:, 2).^2]';

end