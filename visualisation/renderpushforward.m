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

% This is used to render and save vector spherical harmonics for
% demonstration purpose.
clear;
close all;
clc;

mkdir(fullfile('./', 'renderings', 'vspharm'));

clear;
close all;
clc;

% Define dataset.
name = 'cxcr4aMO2_290112';
% Set working directory.
path = fullfile('./', 'data', name);

% Import data.
disp('Loading image data.');
%load(fullfile(path, 'frames-114-116-filtered.mat'));
%load(fullfile(path, 'frames-114-116-unfiltered.mat'));
load(fullfile(path, 'frames-140-142-unfiltered.mat'));

% Import cell centres.
disp('Loading cell centres.');
load(fullfile(path, 'thresholdedcenters.mat'));

% Load colormap for proper visualisation.
load(fullfile(path, 'cmapblue.mat'));

%frame = 114;
frame = 140;

% Set parameters.
Ns = 0:10;
alpha = 0.25;
s = 1;

% Prepare cell centres.
DX = F{frame}.X;
DY = F{frame}.Y;
DZ = -4.2832 * F{frame}.Z;
shift = -min(DZ);

% Fit sphere.
sc = mean([DX, DY, DZ + shift]);
sr = 300;
[sc, sr] = spherefit([DX, DY, DZ + shift], sc, sr);
DX = DX - sc(1);
DY = DY - sc(2);
DZ = DZ - sc(3) + shift;

% Fit spherical surface.
c = surffit(Ns, [DX, DY, DZ], alpha, s);

% Create triangulation.
[F, V] = sphTriang(3);

% Get quadrature points.
Vn = trinodalpts2(F, V);
for k=1:6
    Vn(:, :, k) = normalise(Vn(:, :, k));
end

% Copmpute rho at quadrature points.
for k=1:6
    [~, rhoq] = surfsynth(Ns, Vn(:, :, k), c);
    rho(:, k) = rhoq;
end

% Compute Jacobian at centroids.
xi = repmat([1/3, 1/3], size(F, 1), 1);
[Dx, Dy] = surftanBasis(F, V, rho, xi);

% Compute synthesis at vertices.
[Vs, ~] = surfsynth(Ns, V, c);

% Compute points on surface.
x = surfmap(F, V, rho, xi);

% Run through all degrees and all orders.
N = 2;
for k=N
    % Create vector spherical harmonics.
    [Y1, Y2] = trivspharmcoeff(k, F, V, xi);
    % Create scalar spherical harmonics for visualisation.
    Ynj = spharm(k, V);
    fh = createFigure3;
    for l=1:2*k+1
        cla;
        % Compute pushforward.
        u = bsxfun(@times, Y1(:, l, 1), Dx) + bsxfun(@times, Y1(:, l, 2), Dy);
        v = bsxfun(@times, Y2(:, l, 1), Dx) + bsxfun(@times, Y2(:, l, 2), Dy);
        
        trisurf(F, Vs(:, 1), Vs(:, 2), Vs(:, 3), Ynj(:, l));
        shading interp;
        view(3);
        % Plot vector fields.
        quiver3(x(:, 1), x(:, 2), x(:, 3), u(:, 1), u(:, 2), u(:, 3), 1, 'r');
        quiver3(x(:, 1), x(:, 2), x(:, 3), v(:, 1), v(:, 2), v(:, 3), 1, 'b');
        adjustFigure3;
        %set(gca, 'ZTick', -1:0.5:1);
        %set(gca, 'CLim', [-1, 1]);
        % Save image.
        savefigure(fh, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-600dpi.png', k, l)), '-png', '-r600');
        %savefigure(fh, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-1200dpi.png', k, l)), '-png', '-r1200');
        %savefigure(fh, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-600dpi.jpg', k, l)), '-jpg', '-r600', '-q100');
        %savefigure(fh, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-1200dpi.jpg', k, l)), '-jpg', '-r1200', '-q100');
    end
end
% Save last figure with colorbar.
colorbar;
adjustFigure3;
%set(gca, 'ZTick', -1:0.5:1);
cbar = findobj(fh, 'tag', 'Colorbar');
set(cbar, 'YTick', -1:0.25:1);
set(cbar, 'TickLength', [.02 .02], 'YColor', [0 0 0]);
savefigure(fh, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-colourbar-600dpi.png', k, l)), '-png', '-r600');
%savefigure(fh, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-colourbar-1200dpi.png', k, l)), '-png', '-r1200');
%savefigure(fh, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-colourbar-600dpi.jpg', k, l)), '-jpg', '-r600', '-q100');
%savefigure(fh, fullfile('./', 'renderings', 'vspharm', sprintf('vspharm-deg-%i-ord-%i-colourbar-1200dpi.jpg', k, l)), '-jpg', '-r1200', '-q100');