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
function v = surfintegral(F, V, rho, f, xi, w, a)
%SURFINTEGRAL Computes the surface integral on a sphere-like surface.
%
%   v = SURFINTEGRAL(F, V, rho, f, xi, w, a) takes a triangulation F, V, a
%   function rho evaluated at nodal points, a function f evaluated at the
%   quadrature points xi at every triangle, and weights w and returns the
%   value of the integral. a is a vector containing the areas of the 
%   triangles F.
%
%   Note that xi must be of size k-by-2 and w of length k, where k is the 
%   number of quadrature points.
%
%   rho must be of size m-by-6, where m is the number of triangular faces.
%   f must be of size m-by-k. a must be a vector of length m.
m = size(F, 1);
assert(size(rho, 1) == m);
assert(size(rho, 2) == 6);
assert(size(f, 1) == m);
assert(size(f, 2) == length(w));
assert(size(xi, 1) == length(w));
assert(size(xi, 2) == 2);
assert(isvector(a));
assert(length(a) == m);

% Values of integral at quadrature points.
fq = zeros(m, length(w));

parfor k=1:length(w)
    % Compute determinant of pushforward at quadrature points.
    detphi = detPushforward(F, V, rho, repmat(xi(k, :), m, 1));
    % Evaluate function.
    fq(:, k) = f(:, k) .* sqrt(abs(detphi));
end

% Compute integral.
v = a' * (fq * w);

end

function detphi = detPushforward(F, V, rho, xiq)
    % Compute quadrature points on all triangles.
    xq = trimap(F, V, xiq);
    
    % Compute interpolation of rho at quadrature points.
    rhoq = triinterp2(rho, xiq);
    
    % Compute gradient of rho at quadrature points.
    [gradr, ~, ~] = trigradp2(F, V, rho, xiq);
    
    % Compute determinant of pushforward.
    detphi1 = (rhoq + xq(:, 1).*gradr(:, 1)).*((rhoq + xq(:, 2).*(gradr(:, 2))).*(rhoq + xq(:, 3).*gradr(:, 3)) - xq(:, 2).*xq(:, 3).*gradr(:, 2).*gradr(:, 3));
    detphi2 = xq(:, 1).*gradr(:, 2).*(xq(:, 2).*gradr(:, 1).*(rhoq + xq(:, 3).*gradr(:, 3)) - xq(:, 2).*xq(:, 3).*gradr(:, 1).*gradr(:, 2));
    detphi3 = xq(:, 1).*gradr(:, 3).*(xq(:, 2).*xq(:, 3).*gradr(:, 1).*gradr(:, 3) - xq(:, 3).*gradr(:, 1).*(rhoq + xq(:, 2).*gradr(:, 2)));
    detphi = detphi1 - detphi2 + detphi3;
end