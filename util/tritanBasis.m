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
function [Dx, Dy] = tritanBasis(F, V)
%TRITANBASIS Computes tangential basis of a triangulated surface.
%
%   [Dx, Dy] = TRITANBASIS(F, V) takes a triangulation F, V and returns the
%   tangential basis Dx, Dy.
%
%   Note that Dx and Dy are of size m-by-3, where m is the number of
%   triangular faces F.

Dx = V(F(:, 3), :) - V(F(:, 1), :);
Dy = V(F(:, 2), :) - V(F(:, 1), :);

end