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
function [Yx, Yy] = trivspharmcoeffsynth(N, F, V, u, xi, mem)
%TRIVSPHARMCOEFFSYNTH Computes vector spherical harmonics synthesis on a
%triangulated sphere.
%
%   [Yx, Yy] = TRIVSPHARMCOEFFSYNTH(N, F, V, u, xi) takes coefficients u 
%   for vector sphercial harmonics of degrees N and a triangulation F, V of
%   the unit sphere and returns the coefficients of a field Yi*Di at
%   barycentric coordinates xi. The coefficients are with respect to the
%   basis Dx, Dy computed by TRITANBASIS.
%
%   Note that Yx and Yy form a Helmholtz decomposition on the sphere.
%
%   Note that u must be a vector of length dim, where dim is the dimension
%   of vector spherical harmonics specified by the degrees N.
%
%   [Yx, Yy] = TRIVSPHARMCOEFFSYNTH(N, F, V, u, xi, mem) additionally takes
%   a memory constraint in bytes and allows TRIVSPHARMSYNTH to use up to
%   mem bytes.
%
%   Note that if mem is specified, TRIVSPHARMCOEFFSYNTH then creates 
%   several degrees of vector spherical harmonics at once using mem bytes.
%   However, at least one degree is generated in each iteration, possibly 
%   exceeding the specified memory! However, even if mem is specified, some
%   extra memory is consumed by Matlab for multiplication and addition.
%
%   Yx, Yy are vectors of length size(F, 1).
%
%   Note that N must be a vector of positive consecutive integers!

% Compute and check dimension.
dim = 2*(N(end)^2 + 2*N(end) - N(1)^2 + 1);
assert(isvector(u));
assert(length(u) == dim);

% Check if N is an interval of consecutive positive integers.
assert(isvector(N));
assert(all(N > 0));
assert(length(N) == N(end) - N(1) + 1);
assert(all((N == (N(1):N(end)))));

% Get number of faces.
m = size(F, 1);

% Compute offset for interval.
offset = (N(1)-1)^2 + 2*(N(1)-1);

% Compute intervals I that fit into mem.
I = {};
if(nargin == 6)
    assert(mem > 0);
    s = N(1);
    while(s <= N(end))
        e = max(s, interval(N, s, m, mem));
        I{end+1} = s:e;
        s = e + 1;
    end
else
    I = mat2cell(N', ones(length(N), 1));
end

% Compute triangle areas.
a = triangArea(F, V);

% Project scalar spherical harmonics at nodal points.
Vn = normalise(trinodalpts2(F, V));

% Compute triangulated surface properties.
[Dx, Dy] = tritanBasis(F, V);

% Compute metric properties.
[~, ~, ginv] = metricprops(Dx, Dy);
ginv11 = ginv(:, 1, 1);
ginv12 = ginv(:, 1, 2);
ginv21 = ginv(:, 2, 1);
ginv22 = ginv(:, 2, 2);

% Compute partial derivatives of polynomials at xi.
[DxA, DyA, Q] = tripoly2deriv(xi);
%DxA = permute(DxA * Q, [2, 3, 1]);
%DyA = permute(DyA * Q, [2, 3, 1]);
DxA = (DxA * Q)';
DyA = (DyA * Q)';

% Compute vector spherical harmonics coefficients.
Yx = zeros(m, 1);
Yy = zeros(m, 1);
% Note that for the memory constraints this loop must not run in parallel
% since it would otherwise use more resources! However, the inner loop
% might be subject to parallelisation.
for k=1:length(I)
    % Generate interval.
    M = I{k};
    
    % Create indices.
    idx = [M(1)^2, M(end)^2+2*M(end)] - offset;
    
    % Create temporary variables with coefficients.
    u1 = u(idx(1):idx(2));
    u2 = u(idx(1)+dim/2:idx(2)+dim/2);

    % Compute gradient of interpolation of spherical harmonics.
    DxY = zeros(m, idx(2) - idx(1) + 1);
    DyY = zeros(m, idx(2) - idx(1) + 1);
    for q=1:size(Vn, 3)
        Y = spharmn(M, Vn(:, :, q));
        DxY = DxY + bsxfun(@times, Y, DxA(:, q));
        DyY = DyY + bsxfun(@times, Y, DyA(:, q));
        Y = [];
    end
    
    % Normalise.
    d = sqrt(spharmeigs(M));
    u1 = u1 ./ d;
    u2 = u2 ./ d;
    
    % Compute synthesis of type 2 and type 3, first coefficient.
    Yx = Yx + bsxfun(@times, ginv11, DxY) * u1 + bsxfun(@times, ginv21, DyY) * u1 + bsxfun(@rdivide, DyY, 2 * a) * u2;
    % Compute synthesis of type 2 and type 3, second coefficient.
    Yy = Yy + bsxfun(@times, ginv12, DxY) * u1 + bsxfun(@times, ginv22, DyY) * u1 + bsxfun(@rdivide, -DxY, 2 * a) * u2;
    
    % Explicitly free memory.
    DxY = [];
    DyY = [];
end
end

function e = interval(N, s, m, mem)
    % Compute interval ending such that
    % (e^2 + 2*e^2 - s^2 + 1)*3*m*8 = mem.
    % The first term is the number of spherical harmonics of degrees s:t.
    % The second term is the number of faces m times word length 8 bytes.
    % The number 3 comes from the fact that the inner for-loop needs to
    % store Y, DxY, and DyY, which are each of size dim*m*8 bytes.
    e = min(N(end), floor(sqrt(s^2 + mem/(3*m*8)) - 1));
end