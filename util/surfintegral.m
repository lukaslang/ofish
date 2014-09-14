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
function v = surfintegral(f, w, a)
%SURFINTEGRAL Computes the surface integral on a sphere-like surface.
%
%   v = SURFINTEGRAL(f, w, a) takes a function evaluated at quadrature 
%   points, weights w, and triangular area a and returns the value of the 
%   integral.
%
%   Note that f must be of size [m, nq], where nq is the length of vector 
%   w. a must be a vector of length m.
%
%   Note that assertions have been turned off for efficiency!

% assert(size(f, 1) == length(a));
% assert(size(f, 2) == length(w));
% assert(isvector(a));
% assert(isvector(w));

% Compute integral.
v = a' * (f * w);

end