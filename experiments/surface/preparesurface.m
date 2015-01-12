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
%filename = 'frames-136-138-unfiltered.mat';
%filename = 'frames-138-140-unfiltered.mat';
%filename = 'frames-140-142-unfiltered.mat';
%filename = 'frames-142-144-unfiltered.mat';
%filename = 'frames-144-146-unfiltered.mat';
%filename = 'frames-146-148-unfiltered.mat';
filename = 'frames-148-150-unfiltered.mat';

% Set working directory.
path = fullfile('./', 'data', name);

% Import data.
disp('Loading image data.');
load(fullfile(path, filename));

% Import cell centres.
disp('Loading cell centres.');
C = load(fullfile(path, 'thresholdedcenters.mat'));

% Define cell centres.
frames = [148, 150];

% Scaling of data.
xscale = 1.6774;
yscale = 1.6774;
zscale = 7.1847;

% Parameters for radial projection of the data.
bandwidth = [0.8, 1.2];
layers = 80;

% Set surface fitting parameters.
Ns = 0:50;
beta = 1e-4;
s = 3+eps;

% Number of mesh refinements of the unit sphere triangulation.
ref = 7;

% Create triangulation of the unit sphere.
[F, V] = sphTriang(ref);

S = cell(2);
for k=1:2
    frame = frames(k);

    % Prepare cell centres.
    X = xscale * C.F{frame}.Y;
    Y = yscale * C.F{frame}.X;
    Z = -zscale * C.F{frame}.Z;
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

    % Get nodal points on surface.
    Vn = normalise(trinodalpts2(F, V));

    % Compute synthesis at nodal points.
    [Vn, rhon] = surfsynth(Ns, Vn, c);

    % Compute synthesis at vertices.
    [Vs, rho] = surfsynth(Ns, V, c);

    % Prepare data.
    u = flipdim(U{k}.u, 3);

    % Project data.
    [um, un, uo] = size(u);
    [X, Y, Z] = ndgrid(1:um, 1:un, 1:uo);

    % Compute radial maximum intensity projection at nodal points.
    rs = linspace(bandwidth(1), bandwidth(2), layers);
    fn = zeros(size(F, 1), 6);
    for q=1:6
        VB = kron(rs', Vn(:, :, q));
        fb = dataFromCube(sc(1) + VB(:, 1), sc(2) + VB(:, 2), sc(3) + VB(:, 3), xscale * X, yscale * Y, zscale * Z, u);
        fn(:, q) = max(reshape(fb, size(Vn, 1), length(rs)), [], 2);
    end
    
    % Store surface data.
    S{k}.F = F;
    S{k}.V = V;
    S{k}.Vs = Vs;
    S{k}.rho = rho;
    S{k}.rhon = rhon;
    S{k}.fn = fn;
    S{k}.c = c;    
    S{k}.f = fn;
end

% Free memory.
clear C;
clear F;
clear V;
clear VB;
clear U;
clear u;
clear X;
clear Y;
clear Z;
clear fb;
clear fn;

% Create output directory.
path = fullfile(path, 'generated', 'surface');
mkdir(path);

% Define output file.
[~, file, ~] = fileparts(filename);
datFile = fullfile(path, sprintf('surf-%s-%i-%i-%i.mat', file, Ns(1), Ns(end), ref));

% Write output.
disp('Saving generated data.');
tic;
save(datFile, 'S', 'Ns', 'c', 'ref', 'sc', 'sr', 's', 'beta', 'name', 'bandwidth', 'layers', '-v7.3');
toc;