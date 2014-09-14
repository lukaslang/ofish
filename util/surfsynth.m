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
function [S, rho] = surfsynth(Ns, X, c, Y)
%SURFSYNTH Computes synthesis of a sphere-like surface for given points on 
%the unit sphere.
%
%   [S, rho] = surfsynth(Ns, X, c) takes coefficients c for scalar sphercial 
%   harmonics of degrees Ns and points X on the unit sphere and returns the
%   projection V of these to the sphere-like surface. rho is the radial
%   function evaluated at X.
%
%   [S, rho] = surfsynth(N, X, c, Y) additionally takes scalar spherical 
%   harmonics Y which have been evaluated at X. This is useful e.g. after 
%   calling surffit, which returns Y.
%
%   Note that Ns must be a vector of non-negative consecutive integers. 
%   X is an [n, 3, nq] matrix of points on the unit sphere. c is a vector of 
%   size dim, where dim is the number of scalar spherical harmonics of 
%   degrees specified by Ns.
%
%   S is a matrix of size [n, 3, nq]. rho is a matrix of length [n, nq].

% Check if Ns is an interval of consecutive positive integers.
assert(isvector(Ns));
assert(all(Ns >= 0));
assert(length(Ns) == Ns(end) - Ns(1) + 1);
assert(all((Ns == (Ns(1):Ns(end)))));

n = size(X, 1);
nq = size(X, 3);

% Compute and check dimension.
dim = Ns(end)^2 + 2*Ns(end) - Ns(1)^2 + 1;
assert(isvector(c));
assert(length(c) == dim);

if(nargin == 4)
    % Check if dimensions comply.
    assert(size(Y, 1) == n);
    assert(size(Y, 2) == dim);
    assert(size(Y, 3) == nq);
else
    % Evaluate scalar spherical harmonics at points X.
    Y = zeros(n, dim, nq);
    parfor q=1:nq
        Y(:, :, q) = spharmn(Ns, X(:, :, q));
    end
end

rho = zeros(n, nq);
S = zeros(n, 3, nq);
parfor q=1:nq
    % Recover surface function at X.
    rho(:, q) = Y(:, :, q) * c;

    % Compute coordinates of surface at points X.
    S(:, :, q) = bsxfun(@times, X(:, :, q), rho(:, q));
end