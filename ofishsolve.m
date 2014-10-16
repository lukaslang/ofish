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
function [u, L] = ofishsolve(dim, A, D, b, alpha, tol, maxit)
%OFISHSOLVE Solves the linear system.
%
%   [u, L] = OFISHSOLVE(dim, A, D, b, alpha, maxit) takes precomputed 
%   functions and solves the linear system (A + alpha*D)*u = b. tol is the 
%   desired relative residual. maxit is the maximum number of iterations 
%   during the linear system solve.
%
%   u is a vector of coefficients of the vector spherical harmonics basis 
%   and is of length dim.
%
%   L is a struct containing information about the linear system solve.

assert(alpha >= 0);
assert(dim > 0);
assert(all(size(D, 1) == [dim, dim]));
assert(all(size(A, 1) == [dim, dim]));
assert(isscalar(maxit));
assert(maxit > 0);

% Store norm of rhs.
L.rhs = norm(b, 2);

% Solve linear system.
ticId = tic;
[u, flag, relres, iter, resvec] = gmres(A + alpha * D, b, [], tol, maxit);

% Store solver information.
L.time = toc(ticId);
L.flag = flag;
L.relres = relres;
L.iter = iter;
L.resvec = resvec;
L.tol = tol;
L.maxit = maxit;
L.solver = 'gmres';
L.restart = 0;

end