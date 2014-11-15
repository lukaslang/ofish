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
function [Z11, Z12, Z21, Z22] = surfcovderivn(N, G, g, A, F, V, xi, mem)
%SURFCOVDERIVN Computes the covariant derivatives of the pushforward of
%vector spherical harmonics on a sphere-like surface.
%
%   Z = SURFCOVDERIV(N, G, g, A, xi) takes a vector field specified by 
%   a transformation matrix A from the tangent basis at points xi, an 
%   interval of degrees N of vector spherical harmonics and returns the 
%   coefficients Zij of the covariant derivative of vector spherical 
%   harmonics of degrees N. In addition, Christoffel symbols G and metric g
%   are needed.
%
%   Note that xi must be a matrix of size [m, 2, nq], where m is the number
%   of triangular faces and nq is the number of quadrature points.
%
%   A must be of size [2, 2, m, 6] and is interpolated quadratically.
%
%   N must be an interval of consecutive positive integers.
%
%   dim is the dimension of both tangential type vspharm!
%
%   Zij is of size [m, dim, nq] and denotes D_i u^j, which is the
%   derivation with respect to e_i (the ONB given by matrix A).

% Check if N is an interval of consecutive positive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

m = size(A, 3);
assert(size(A, 1) == 2);
assert(size(A, 2) == 2);
assert(size(A, 4) == 6);

nq = size(xi, 3);
assert(size(xi, 1) == m);
assert(size(xi, 2) == 2);

assert(size(G, 1) == m);
assert(size(G, 2) == 2);
assert(size(G, 3) == 2);
assert(size(G, 4) == 2);
assert(size(G, 5) == nq);

assert(size(g, 1) == m);
assert(size(g, 2) == 2);
assert(size(g, 3) == 2);
assert(size(g, 4) == nq);

% Compute and check dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);

% Compute offset for interval.
offset = (N(1)-1)^2 + 2*(N(1)-1);

% Compute intervals I that fit into mem.
I = {};
if(nargin == 8)
    % Determine memory required to store Zij.
    reqmem = 4*m*dim*nq*8;
    % Check if specified memory is sufficient.
    assert(mem > reqmem);
    % Compute intervals with remaining memory.
    s = N(1);
    while(s <= N(end))
        e = max(s, interval(N, s, m, nq, mem - reqmem));
        I{end+1} = s:e;
        s = e + 1;
    end
else
    I = mat2cell(N', ones(length(N), 1));
end

Z11 = zeros(m, dim, nq);
Z12 = zeros(m, dim, nq);
Z21 = zeros(m, dim, nq);
Z22 = zeros(m, dim, nq);
for k=1:length(I)
    % Generate interval.
    M = I{k};
    
    % Create indices.
    idx = [M(1)^2, M(end)^2+2*M(end)] - offset;
    len = idx(2) - idx(1) + 1;

    % Compute coefficients.
    [Y, DY] = trivspharmncoeff(M, F, V, xi);

    % Compute covariant derivatives.
    [Z11i, Z12i, Z21i, Z22i] = surfcovderiv(G, g, A, Y, DY, xi);

    % Explicitly free memory.
    Y = [];
    DY = [];
    
    % Store.
    Z11(:, idx(1):idx(2), :) = Z11i(:, 1:len, :);    
    Z11(:, idx(1)+dim/2:idx(2)+dim/2, :) = Z11i(:, len+1:2*len, :);
    
    Z12(:, idx(1):idx(2), :) = Z12i(:, 1:len, :);    
    Z12(:, idx(1)+dim/2:idx(2)+dim/2, :) = Z12i(:, len+1:2*len, :);
    
    Z21(:, idx(1):idx(2), :) = Z21i(:, 1:len, :);    
    Z21(:, idx(1)+dim/2:idx(2)+dim/2, :) = Z21i(:, len+1:2*len, :);
    
    Z22(:, idx(1):idx(2), :) = Z22i(:, 1:len, :);    
    Z22(:, idx(1)+dim/2:idx(2)+dim/2, :) = Z22i(:, len+1:2*len, :);
    
    % Explicitly free memory.
    Z11i = [];
    Z12i = [];
    Z21i = [];
    Z22i = [];
end

end


% TODO: Compute actual memory needed (e.g. add nq).
function e = interval(N, s, m, nq, mem)
    % Compute interval ending such that
    % (e^2 + 2*e^2 - s^2 + 1)*6*m*nq*8 = mem.
    % The first term is the number of spherical harmonics of degrees s:t.
    % The second term is the number of faces m times word length 8 bytes.
    % The number 6 comes from the fact that the inner for-loop needs to
    % store Y and DY, which are 6 times each dim*m*nq*8 bytes.
    e = min(N(end), floor(sqrt(2*(s^2 + mem/(6*m*nq*8) - 1)/3)));
end