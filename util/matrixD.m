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
function D = matrixD(dim, Z11, Z12, Z21, Z22, detphi, w, a)
%MATRIXD Computes the regularisation matrix.

nq = length(w);
D = zeros(dim, dim);
for p=1:dim
    Z11pq = bsxfun(@times, Z11(:, p, :), Z11(:, 1:p, :));
    Z12pq = bsxfun(@times, Z12(:, p, :), Z12(:, 1:p, :));
    Z21pq = bsxfun(@times, Z21(:, p, :), Z21(:, 1:p, :));
    Z22pq = bsxfun(@times, Z22(:, p, :), Z22(:, 1:p, :));
    f = bsxfun(@times, Z11pq + Z12pq + Z21pq + Z22pq, permute(sqrt(abs(detphi)), [1, 3, 2]));
    f = 0.5 * sum(bsxfun(@times, sum(bsxfun(@times, f, reshape(w, [1, 1, nq])), 3), a), 1);
    D(p, 1:p) = f;
    D(1:p, p) = f;
end

% Add diagonal.
D(1:(dim+1):end) = 2 * diag(D);

end