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
function Vn = trinodalpts2(F, V)
%TRINODALPTS2 Returns nodal points for quadratic interpolation on
%triangles.
%
%   Vn = TRINODALPTS2(F, V) takes a triangulation F, V and returns nodal
%   points Vn for quadratic interpolation.
%
%   Note that Vn is a matrix of size m-by-3-by-6.

Vn = zeros(size(F, 1), 3, 6);
Vn(:, :, 1) = V(F(:, 1), :);
Vn(:, :, 2) = V(F(:, 2), :);
Vn(:, :, 3) = V(F(:, 3), :);
Vn(:, :, 4) = 0.5 * (V(F(:, 2), :) + V(F(:, 1), :));
Vn(:, :, 5) = 0.5 * (V(F(:, 3), :) + V(F(:, 2), :));
Vn(:, :, 6) = 0.5 * (V(F(:, 1), :) + V(F(:, 3), :));

end