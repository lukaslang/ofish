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
function A = matrixA(dim, Z, detphi, w, a)
%MATRIXA Computes the data matrix.

nq = length(w);
m = length(a);
A = zeros(dim, dim);
for p=1:dim
    Zpq = bsxfun(@times, Z(:, p, :), Z(:, 1:p, :));
    f = bsxfun(@times, Zpq, permute(sqrt(abs(detphi)), [1, 3, 2]));    
    f = sum(bsxfun(@times, sum(bsxfun(@times, f, repmat(reshape(w, [1, 1, nq]), [m, p, 1])), 3), a), 1);
    A(p, 1:p) = f;
    A(1:p, p) = f;
end
end