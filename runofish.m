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
% Set working directory.
path = fullfile('./', 'data', name);

% Import data.
disp('Loading image data.');
%load(fullfile(path, 'frames-114-116-unfiltered.mat'));
%load(fullfile(path, 'frames-120-122-unfiltered.mat'));
load(fullfile(path, 'frames-140-142-unfiltered.mat'));

% Import cell centres.
disp('Loading cell centres.');
C = load(fullfile(path, 'thresholdedcenters.mat'));

% Load colormap for proper visualisation.
load(fullfile(path, 'cmapblue.mat'));

% Define cell centres.
frames = [140, 142];

% Scaling of data.
xscale = 1.6774;
yscale = 1.6774;
zscale = 7.1847;

% Set maximum degree of vector spherical harmonics basis.
N = 10;

% Set regularisation parameter for optical flow.
alpha = 1;

% Finite difference time parameter.
h = 1;

% Degree of numerical quadrature.
deg = 1;

% Parameters for radial projection of the data.
bandwidth = [0.8, 1.2];
layers = 80;

% Set surface fitting parameters.
Ns = 0:10;
beta = 0.5;
s = 1;

% Number of mesh refinements of the unit sphere triangulation.
ref = 7;

% Create triangulation of the unit sphere.
[F, V] = sphTriang(ref);

f = cell(2);
fq = cell(2);
for k=[2, 1]
    frame = frames(k);

    % Prepare cell centres (and flip X and Y axis).
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

    % Compute radial maximum intensity projection at vertices.
    rs = linspace(bandwidth(1), bandwidth(2), layers);
    VB = kron(rs', Vs);
    fb = dataFromCube(sc(1) + VB(:, 1), sc(2) + VB(:, 2), sc(3) + VB(:, 3), xscale * X, yscale * Y, zscale * Z, u);
    f{k} = max(reshape(fb, size(Vs, 1), length(rs)), [], 2);
    
    % Compute radial maximum intensity projection at nodal points.
    fq{k} = zeros(size(F, 1), 6);
    for q=1:6
        VB = kron(rs', Vn(:, :, q));
        fb = dataFromCube(sc(1) + VB(:, 1), sc(2) + VB(:, 2), sc(3) + VB(:, 3), xscale * X, yscale * Y, zscale * Z, u);
        fq{k}(:, q) = max(reshape(fb, size(Vn, 1), length(rs)), [], 2);
    end
    
    % Plot function rho on the unit sphere.
    figure;
    hold on;
    axis([-1, 1, -1, 1, -1, 1]);
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), rho, 'EdgeColor', 'none');
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    colorbar;

    % Plot surface.
    figure;
    hold on;
    trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3));
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    
    % Plot data.
    figure;
    colormap(cmap);
    subplot(1, 2, 1);
    trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), f{k}, 'EdgeColor', 'none');
    shading interp;
    daspect([1, 1, 1]);
    view(3);
    
    subplot(1, 2, 2);
    axis([-1, 1, -1, 1, -1, 1]);
    trisurf(F, V(:, 1), V(:, 2), V(:, 3), f{k}, 'EdgeColor', 'none');
    shading interp;
    daspect([1, 1, 1]);
    view(3);
end

% Free memory.
clear C;
clear VB;
clear U;
clear u;
clear X;
clear Y;
clear Z;
clear fb;

disp('Computing optical flow...');
ticId = tic;
u = ofish(N, Ns, c, F, V, fq{1}, fq{2}, h, deg, alpha);
elapsedTime = toc(ticId);
fprintf('Elapsed time %.6f seconds.\n', elapsedTime);

% Get incenters of triangles.
TR = TriRep(F, Vs);
IC = TR.incenters;

% Plot result.
figure;
hold on;
trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), f{1});
shading interp;
colormap(cmap);
daspect([1, 1, 1]);
view(3);
quiver3(IC(:, 1), IC(:, 2), IC(:, 3), u(:, 1), u(:, 2), u(:, 3), 0, 'r');

% Project and scale flow.
up = projecttoplane(u);

% Compute colour space scaling.
nmax = max(sqrt(sum(up.^2, 2)));

% Compute colour of projection.
col = double(squeeze(computeColour(up(:, 1)/nmax, up(:, 2)/nmax))) ./ 255;

figure;
trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', col, 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);

% Plot colourwheel.
figure;
cw = colourWheel;
surf(1:200, 1:200, zeros(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
daspect([1, 1, 1]);
view(3);