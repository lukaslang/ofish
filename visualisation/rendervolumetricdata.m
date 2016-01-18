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

% The purpose of this script is to render and save volumetric microscopy
% data used during experiments.
clear;
close all;
clc;

% Define dataset.
name = 'cxcr4aMO2_290112';
% Set working directory.
path = fullfile('./', 'data', name);

% Import data.
disp('Loading image data.');
filename = 'frames-140-141-unfiltered';
load(fullfile(path, sprintf('%s.mat', filename)));

% Load colormap for proper visualisation.
load(fullfile(path, 'cmapblue.mat'));

% Scaling of data.
xscale = 1.6774;
yscale = 1.6774;
zscale = 7.1847;

% Create directory.
renderPath = fullfile('./', 'renderings', 'volumetric');
mkdir(renderPath);

% Run through both frames.
for k=1:2
    % Prepare data (flip z-axis).
    u = flipdim(U{k}.u, 3);
    % Switch x and y axis.
    u = permute(u, [2, 1, 3]);
    
    % Create grid data.
    [um, un, uo] = size(u);
    [X, Y, Z] = ndgrid(1:um, 1:un, 1:uo);
    
    % Create volumetric plot.
    F = createFigure3(cmap, xscale, um * xscale, yscale, un * yscale, zscale, uo * zscale);
    set(F, 'renderer', 'opengl');
    vol3d('CData', u, 'XData', X * xscale, 'YData', Y * yscale, 'ZData', Z * zscale);
    box on;
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 14);
    set(gca, 'TickDir', 'out')
    set(gca, 'TickLength', [.02 .02]);
    set(gca, 'XTick', 0:200:900);
    set(gca, 'YTick', 0:200:900);
    set(gca, 'ZTick', 0:100:400);    
    set(gca, 'XMinorTick', 'off');
    set(gca, 'YMinorTick', 'off');
    set(gca, 'ZMinorTick', 'off');
    view(3);
    export_fig(fullfile(renderPath, sprintf('raw3-%s-%i-600dpi.png', filename, k)), '-png', '-r600', '-opengl', '-transparent', F); 
    export_fig(fullfile(renderPath, sprintf('raw3-%s-%i-1200dpi.png', filename, k)), '-png', '-r1200', '-opengl', '-transparent', F);
    export_fig(fullfile(renderPath, sprintf('raw3-%s-%i-600dpi.jpg', filename, k)), '-jpg', '-r600', '-opengl', '-transparent', F); 
    export_fig(fullfile(renderPath, sprintf('raw3-%s-%i-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-opengl', '-transparent', F); 
    view(2);
    % Save one 2D image for control reasons.
    export_fig(fullfile(renderPath, sprintf('raw2-%s-%i-600dpi.png', filename, k)), '-png', '-r600', '-opengl', '-transparent', F); 
end