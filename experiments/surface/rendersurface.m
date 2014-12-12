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

% Define dataset and get result files.
name = 'cxcr4aMO2_290112';

% Import data.
disp('Loading precomputed data.');
path = fullfile('./', 'data', name, 'generated', 'surface');
filename = 'frames-140-142-unfiltered-0-50-7';
D = load(fullfile(path, sprintf('surf-%s.mat', filename)));

% Load colormap for proper visualisation.
load(fullfile('./', 'data', name, 'cmapblue.mat'));

% Define renderings path.
renderPath = fullfile('./', 'renderings', name, 'surface');
mkdir(renderPath);
mkdir(fullfile(renderPath, 'data2'));
mkdir(fullfile(renderPath, 'data3'));
mkdir(fullfile(renderPath, 'rho2'));
mkdir(fullfile(renderPath, 'rho3'));

% Define scaling factor of spherical mesh.
gridscale = 1.01;

% Find min and max values of rho.
rhomin = min(min(D.S{1}.rho, D.S{2}.rho));
rhomax = max(max(D.S{1}.rho, D.S{2}.rho));

% Plot data.
fn = zeros(size(D.S{1}.V, 1), 1);
for l=1:2
    % Choose surface.
    S = D.S{l};

    % Create spherical mesh.
    nsphere = 31;
    [Xm, Ym, Zm] = sphere(nsphere-1);

    % Remove poles.
    Xm = Xm(2:end-1, :);
    Ym = Ym(2:end-1, :);
    Zm = Zm(2:end-1, :);
    
    % Compute synthesis at mesh vertices.
    [Vm, ~] = surfsynth(D.Ns, [Xm(:), Ym(:), Zm(:)], D.S{l}.c);
    Vm = Vm .* gridscale;
    Xm = reshape(Vm(:, 1), nsphere-2, nsphere);
    Ym = reshape(Vm(:, 2), nsphere-2, nsphere);
    Zm = reshape(Vm(:, 3), nsphere-2, nsphere);
    
    % Plot function rho on the unit sphere.
    F = createFigure3;
    caxis([rhomin rhomax]);
    trisurf(S.F, S.V(:, 1), S.V(:, 2), S.V(:, 3), S.rho, 'EdgeColor', 'none');
    shading interp;
    view(3);
    set(gca, 'ZLim', [-1, 1]);
    set(gca, 'XLim', [-1, 1]);
    set(gca, 'YLim', [-1, 1]);
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], 'XMinorTick', 'on', 'YMinorTick', 'on', 'ZMinorTick', 'on', 'YGrid', 'off');
    adjustFigure3;
    set(gca, 'XTick', -1:0.5:1);
    set(gca, 'YTick', -1:0.5:1);
    set(gca, 'ZTick', -1:0.5:1);
    savefigure(F, fullfile(renderPath, 'rho3', sprintf('rho3-%s-%i-600dpi.png', filename, l)), '-png', '-r600');
    savefigure(F, fullfile(renderPath, 'rho3', sprintf('rho3-%s-%i-1200dpi.png', filename, l)), '-png', '-r1200');
    savefigure(F, fullfile(renderPath, 'rho3', sprintf('rho3-%s-%i-600dpi.jpg', filename, l)), '-jpg', '-r600', '-q100');
    savefigure(F, fullfile(renderPath, 'rho3', sprintf('rho3-%s-%i-1200dpi.jpg', filename, l)), '-jpg', '-r1200', '-q100');

    % Rotate by pi.
    [az, el] = view;
    view(az + 180, el);
    savefigure(F, fullfile(renderPath, 'rho3', sprintf('rho3-%s-%i-rotated-600dpi.png', filename, l)), '-png', '-r600');
    savefigure(F, fullfile(renderPath, 'rho3', sprintf('rho3-%s-%i-rotated-1200dpi.png', filename, l)), '-png', '-r1200');
    savefigure(F, fullfile(renderPath, 'rho3', sprintf('rho3-%s-%i-rotated-600dpi.jpg', filename, l)), '-jpg', '-r600', '-q100');
    savefigure(F, fullfile(renderPath, 'rho3', sprintf('rho3-%s-%i-rotated-1200dpi.jpg', filename, l)), '-jpg', '-r1200', '-q100');

    colorbar;
    cbar = findobj(F, 'tag', 'Colorbar');
    set(cbar, 'YTick', 280:10:400);
    set(cbar, 'TickLength', [.02 .02], 'YColor', [0 0 0]);
    view(2);
    adjustFigure3;
    set(gca, 'XTick', -1:0.5:1);
    set(gca, 'YTick', -1:0.5:1);
    set(gca, 'ZTick', -1:0.5:1);
    savefigure(F, fullfile(renderPath, 'rho2', sprintf('rho2-%s-%i-600dpi.png', filename, l)), '-png', '-r600');
    savefigure(F, fullfile(renderPath, 'rho2', sprintf('rho2-%s-%i-1200dpi.png', filename, l)), '-png', '-r1200');
    savefigure(F, fullfile(renderPath, 'rho2', sprintf('rho2-%s-%i-600dpi.jpg', filename, l)), '-jpg', '-r600', '-q100');
    savefigure(F, fullfile(renderPath, 'rho2', sprintf('rho2-%s-%i-1200dpi.jpg', filename, l)), '-jpg', '-r1200', '-q100');

    % Convert data from nodal points to a data at vertices.
    fv = S.f(:, 1:3);
    fn(S.F) = fv;
    
    F = createFigure3(cmap);
    trisurf(S.F, S.Vs(:, 1), S.Vs(:, 2), S.Vs(:, 3), fn, 'EdgeColor', 'none');
    shading interp;
    view(3);
    set(gca, 'XTick', -300:150:300);
    set(gca, 'YTick', -300:150:300);
    set(gca, 'ZTick', -300:150:300);
    adjustFigure3;
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-600dpi.png', filename, l)), '-png', '-r600');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-1200dpi.png', filename, l)), '-png', '-r1200');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-600dpi.jpg', filename, l)), '-jpg', '-r600', '-q100');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-1200dpi.jpg', filename, l)), '-jpg', '-r1200', '-q100');
    
    % Rotate by pi.
    [az, el] = view;
    view(az + 180, el);
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-rotated-600dpi.png', filename, l)), '-png', '-r600');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-rotated-1200dpi.png', filename, l)), '-png', '-r1200');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-rotated-600dpi.jpg', filename, l)), '-jpg', '-r600', '-q100');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-rotated-1200dpi.jpg', filename, l)), '-jpg', '-r1200', '-q100');
    
    view(2);
    savefigure(F, fullfile(renderPath, 'data2', sprintf('data2-%s-%i-600dpi.png', filename, l)), '-png', '-r600');
    savefigure(F, fullfile(renderPath, 'data2', sprintf('data2-%s-%i-1200dpi.png', filename, l)), '-png', '-r1200');
    savefigure(F, fullfile(renderPath, 'data2', sprintf('data2-%s-%i-600dpi.jpg', filename, l)), '-jpg', '-r600', '-q100');
    savefigure(F, fullfile(renderPath, 'data2', sprintf('data2-%s-%i-1200dpi.jpg', filename, l)), '-jpg', '-r1200', '-q100');
    
    % Create a plot with grid.
    view(3);
    surf(Xm, Ym, Zm, ones(nsphere-2, nsphere), 'EdgeColor', [0.7, 0.7, 0.7], 'FaceColor', 'none', 'FaceAlpha', 1, 'EdgeAlpha', 1);
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-grid-600dpi.png', filename, l)), '-png', '-r600');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-grid-1200dpi.png', filename, l)), '-png', '-r1200');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-grid-600dpi.jpg', filename, l)), '-jpg', '-r600', '-q100');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-grid-1200dpi.jpg', filename, l)), '-jpg', '-r1200', '-q100');
    
    % Rotate by pi.
    [az, el] = view;
    view(az + 180, el);
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-grid-rotated-600dpi.png', filename, l)), '-png', '-r600');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-grid-rotated-1200dpi.png', filename, l)), '-png', '-r1200');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-grid-rotated-600dpi.jpg', filename, l)), '-jpg', '-r600', '-q100');
    savefigure(F, fullfile(renderPath, 'data3', sprintf('data3-%s-%i-grid-rotated-1200dpi.jpg', filename, l)), '-jpg', '-r1200', '-q100');

    view(2);
    savefigure(F, fullfile(renderPath, 'data2', sprintf('data2-%s-%i-grid-600dpi.png', filename, l)), '-png', '-r600');
    savefigure(F, fullfile(renderPath, 'data2', sprintf('data2-%s-%i-grid-1200dpi.png', filename, l)), '-png', '-r1200');
    savefigure(F, fullfile(renderPath, 'data2', sprintf('data2-%s-%i-grid-600dpi.jpg', filename, l)), '-jpg', '-r600', '-q100');
    savefigure(F, fullfile(renderPath, 'data2', sprintf('data2-%s-%i-grid-1200dpi.jpg', filename, l)), '-jpg', '-r1200', '-q100');
end