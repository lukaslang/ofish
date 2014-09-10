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
function [g, detg, ginv] = metricprops(Dx, Dy)
%METRICPROPS Computes metric, its inverse, and the determinant.
%
%   [g, detg, ginv] = METRICPROPS(Dx, Dy) takes a tangential basis Dx, Dy 
%   and returns the metric tensor g, the determinant detg, and the inverse 
%   ginv of the metric.
%
%   Note that Dx and Dy must be of size m-by-3, where m is the number of
%   triangular faces. 
%
%   g is of size m-by-2-by-2, detg a vector of length m, and ginv a matrix 
%   of size m-by-2-by-2.

% Compute metric.
g(:, 1, 1) = dot(Dx, Dx, 2);
g(:, 1, 2) = dot(Dx, Dy, 2);
g(:, 2, 1) = g(:, 1, 2);
g(:, 2, 2) = dot(Dy, Dy, 2);

% Compute its determinant.
detg = g(:, 1, 1) .* g(:, 2, 2) - g(:, 1, 2) .^2;

% Compute inverse.
ginv(:, 1, 1) = g(:, 2, 2);
ginv(:, 1, 2) = -g(:, 1, 2);
ginv(:, 2, 1) = ginv(:, 1, 2);
ginv(:, 2, 2) = g(:, 1, 1);
ginv = bsxfun(@rdivide, ginv, detg);

end