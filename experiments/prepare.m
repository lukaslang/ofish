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
clear;
close all;
clc;

% Define dataset.
name = 'cxcr4aMO2_290112';
%filename = 'frames-114-116-unfiltered.mat';
%filename = 'frames-120-122-unfiltered.mat';
filename = 'frames-140-142-unfiltered.mat';

% Set working directory.
path = fullfile('./', 'data', name);

% Import data.
disp('Loading image data.');
load(fullfile(path, filename));

% Import cell centres.
disp('Loading cell centres.');
load(fullfile(path, 'thresholdedcenters.mat'));

% Define cell centres.
%frame = 114;
%frame = 120;
frame = 140;

% Scaling in z-direction.
xscale = 1.6774;
yscale = 1.6774;
zscale = 7.1847;

% Set degrees of vector spherical harmonics basis.
N = 1:50;

% Finite difference time parameter.
h = 1;

% Numerical integration tolerance.
tol = 1e-6;

% Degree of numerical quadrature.
deg = 1;

% Parameters for radial projection of the data.
bandwidth = [0.8, 1.2];
layers = 80;

% Set surface fitting parameters.
Ns = 0:50;
beta = 0.25;
s = 1;

% Number of mesh refinements of the unit sphere triangulation.
ref = 7;

% Prepare cell centres.
X = xscale * F{frame}.X;
Y = yscale * F{frame}.Y;
Z = -zscale * F{frame}.Z;
shift = -min(Z);

% Fit sphere.
sc = mean([X, Y, Z + shift]);
sr = 300;
[sc, sr] = spherefit([X, Y, Z + shift], sc, sr);

% Center data.
X = X - sc(1);
Y = Y - sc(2);
Z = Z - sc(3) + shift;

% Fit sphere-like surface.
c = surffit(Ns, [X, Y, Z], beta, s);

% Create triangulation of the unit sphere.
[F, V] = sphTriang(ref);

% Get nodal points on surface.
Vn = normalise(trinodalpts2(F, V));

% Compute synthesis at nodal points.
[Vn, rhon] = surfsynth(Ns, Vn, c);

% Compute synthesis at vertices.
[Vs, rho] = surfsynth(Ns, V, c);

f = cell(2);
for k=1:2
    % Prepare data.
    u = flipdim(U{k}.u, 3);

    % Project data.
    [um, un, uo] = size(u);
    [X, Y, Z] = ndgrid(1:um, 1:un, 1:uo);

    % Compute radial maximum intensity projection at nodal points.
    rs = linspace(bandwidth(1), bandwidth(2), layers);
    fq = zeros(size(F, 1), 6);
    for q=1:6
        VB = kron(rs', Vn(:, :, q));
        fb = dataFromCube(sc(1) + VB(:, 1), sc(2) + VB(:, 2), sc(3) + VB(:, 3), xscale * X, yscale * Y, zscale * Z, u);
        fq(:, q) = max(reshape(fb, size(Vn, 1), length(rs)), [], 2);
    end
    f{k} = fq;
end

% Free memory.
clear VB;
clear U;
clear u;
clear X;
clear Y;
clear Z;
clear fb;
clear fq;

% Create output directory.
path = fullfile(path, 'generated');
mkdir(path);

[~, file, ~] = fileparts(filename);
disp('Computing linear system.');
tic;
[dim, A, D, b] = surflinearsystem(F, V, Ns, c, N, f{1}, f{2}, h, deg, tol);
elapsed = toc;
fprintf('Elapsed time is %.6f seconds.\n', elapsed);

% Define output file.
genFile = fullfile(path, sprintf('gen-%s-%i-%i-%i.mat', file, N(1), N(end), ref));
datFile = fullfile(path, sprintf('dat-%s-%i-%i-%i.mat', file, N(1), N(end), ref));

% Write output.
disp('Saving generated data.');
tic;
save(datFile, 'F', 'V', 'Vs', 'rho', 'rhon', 'f', 'Ns', 'c', 'N', 'ref', 'sc', 'sr', 's', 'beta', 'h', 'name', 'tol', 'deg', 'elapsed', 'bandwidth', 'layers', 'frame', '-v7.3');
save(genFile, 'dim', 'A', 'D', 'b', '-v7.3');
toc;