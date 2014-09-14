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
function V = normalise(V)
%NORMALISE Normalises row vectors.
%
%   V = normalise(V) takes an [n, 3, nq] matrix V and returns normalises 
%   row vectors for each V(:, :, q).

len = sqrt(sum(V.^2, 2));
V = bsxfun(@rdivide, V, len);

end