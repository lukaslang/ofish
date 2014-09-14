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
function [gradf, u, v] = trigradp2(F, V, f, xi)
%TRIGRADP2 Computes the surface gradient of a function which is interpolated 
%piecewise quadraticly.
%
%   [gradf, u, v] = TRIGRADP2(F, V, f, xi) takes a triangulation F, V of a 
%   surface, a matrix f, which are the values at the nodal points, and 
%   returns the surface gradient of the interpolated function f at given at
%   points xi, which are given in barycentric coordinates. In addition, the
%   coefficients u, v are returned.
%
%   Note that f must be of size [m, 6] where each row corresponds to the
%   six quadrature points and xi is of size [m, 2, nq]. m is the number of
%   triangles and nq the quadrature dimension.
%
%   Note that size(gradf) = [m, 3, nq].

assert(size(f, 1) == size(F, 1));
assert(size(f, 1) == size(xi, 1));
assert(size(f, 2) == 6);
assert(size(xi, 2) == 2);
nq = size(xi, 3);
m = size(F, 1);

gradf = zeros(m, 3, nq);
u = zeros(m, nq);
v = zeros(m, nq);
parfor q=1:nq
    % Compute derivatives of polynomials at xi.
    [DxA, DyA, Q] = tripoly2deriv(xi);

    % Compute triangulated surface properties.
    [Dx, Dy] = tritanBasis(F, V);

    % Compute metric properties.
    [~, ~, ginv] = metricprops(Dx, Dy);

    % Compute partial derivative of interpolation.
    Dxf = dot(f, (DxA * Q)', 2);
    Dyf = dot(f, (DyA * Q)', 2);

    % Compute coefficients of gradient.
    u(:, q) = (ginv(:, 1, 1) .* Dxf + ginv(:, 1, 2) .* Dyf);
    v(:, q) = (ginv(:, 2, 1) .* Dxf + ginv(:, 2, 2) .* Dyf);
    
    % Compute surface gradient.
    gradf(:, :, q) = bsxfun(@times, u(:, q), Dx) + bsxfun(@times, v(:, q), Dy);
end
end