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
function detphi = detpushforward(F, V, rho, xi)
%DETPUSHFORWARD Computes the determinant of the pushforward from a
%triangulated unit sphere to a sphere-like surface.
%
%   detphi = DETPUSHFORWARD(F, V, rho, xi) takes a triangulation F, V, a 
%   function rho defining a sphere-like surface and returns the determinant
%   if the pushforward phi at given barycentric coordinates xi.
%
%   Note that xi must be of size m-by-3, where m is the number of
%   triangular faces.
%
%   v is a vector of length m.

nq = size(xi, 3);

% Compute points on triangles.
x = trimap(F, V, xi);

% Compute interpolation of rho at xi.
rhoq = triinterp2(rho, xi);

% Compute gradient of rho at xi.
[gradr, ~, ~] = trigradp2(F, V, rho, xi);

% Compute determinant of pushforward.
detphi = rhoq.^3 + rhoq.^2 .* squeeze(dot(x, gradr, 2));

end