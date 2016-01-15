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
resultsPath = fullfile('./', 'results', name, 'ofish');
resultsname = '2016-01-10-15-53-14-frames-140-141-unfiltered-1-50-7';
%resultsname = '2016-01-14-10-08-43-frames-140-141-unfiltered-1-50-7-sphere';
load(fullfile(resultsPath, sprintf('%s.mat', resultsname)));

% Import data.
disp('Loading precomputed data.');
path = fullfile('./', 'data', name, 'generated');
filename = 'frames-140-141-unfiltered-1-50-7';
%filename = 'frames-140-141-unfiltered-1-50-7-sphere';
D = load(fullfile(path, sprintf('dat-%s.mat', filename)));

% Load colormap for proper visualisation.
load(fullfile('./', 'data', name, 'cmapblue.mat'));

% Define renderings path.
renderPath = fullfile('./', 'renderings', name, 'ofish', resultsname);
mkdir(renderPath);
mkdir(fullfile(renderPath, 'residual'));
mkdir(fullfile(renderPath, 'coefficients'));
mkdir(fullfile(renderPath, 'data2'));
mkdir(fullfile(renderPath, 'data3'));
mkdir(fullfile(renderPath, 'flow2'));
mkdir(fullfile(renderPath, 'flow3'));
mkdir(fullfile(renderPath, 'rho2'));
mkdir(fullfile(renderPath, 'rho3'));
mkdir(fullfile(renderPath, 'surfvel2'));
mkdir(fullfile(renderPath, 'surfvel3'));
mkdir(fullfile(renderPath, 'signnormsurfvel2'));
mkdir(fullfile(renderPath, 'signnormsurfvel3'));
mkdir(fullfile(renderPath, 'motion2'));
mkdir(fullfile(renderPath, 'motion3'));

% Restriction allows to search among the results.
%e = cell2mat(E);
%idx = find([e.s] == 1);

% Select results to render.
%idx = [1:8];
idx = 1:length(E);

% Since the data for every experiment is the same we just pick one.
k = 1;

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
    if(rhomax - rhomin > eps)
        caxis([rhomin rhomax]);
    end
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
    set(cbar, 'YTick', 280:20:460);
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
    set(gca, 'XTick', -450:150:450);
    set(gca, 'YTick', -450:150:450);
    set(gca, 'ZTick', -450:150:450);
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

% Open a file to save flow scaling.
fid = fopen(fullfile(renderPath, 'radii.txt'), 'w');

% Compute time derivative of rho at nodal points.
veln = (D.S{2}.rhon - D.S{1}.rhon) ./ D.h;

% Specify evaluation points.
xi = repmat([1/3, 1/3], size(S.F, 1), 1);

% Evaluate time derivative of rho at evaluation points.
velev = triinterp2(veln, xi);

% Compute surface velocity.
vel = bsxfun(@times, velev, trimap(D.S{1}.F, D.S{1}.V, xi));

% Recover vector field.
velp = projecttoplane(vel);

% Compute colour space scaling.
nmax = max(sqrt(sum(velp.^2, 2)));

% Save flow scaling.
fprintf(fid, 'Surface velocity: %.4f\n', nmax);

% Scale velocity.
c = double(squeeze(computeColour(velp(:, 1)/nmax, velp(:, 2)/nmax))) ./ 255;

% Create colourwheel.
cw = colourwheelbg;

% Plot surface velocity.
F = createFigure3(cmap);
trisurf(S.F, S.Vs(:, 1), S.Vs(:, 2), S.Vs(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
C = surf(-350:-151, -350:-151, -380*ones(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
view(3);
set(gca, 'XTick', -450:150:450);
set(gca, 'YTick', -450:150:450);
set(gca, 'ZTick', -450:150:450);
adjustFigure3;
savefigure(F, fullfile(renderPath, 'surfvel3', sprintf('surfvel3-%s-%i-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'surfvel3', sprintf('surfvel3-%s-%i-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'surfvel3', sprintf('surfvel3-%s-%i-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'surfvel3', sprintf('surfvel3-%s-%i-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');
delete(C);
% Rotate by pi.
[az, el] = view;
view(az + 180, el);
C = surf(151:350, 151:350, -380*ones(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
savefigure(F, fullfile(renderPath, 'surfvel3', sprintf('surfvel3-%s-%i-rotated-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'surfvel3', sprintf('surfvel3-%s-%i-rotated-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'surfvel3', sprintf('surfvel3-%s-%i-rotated-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'surfvel3', sprintf('surfvel3-%s-%i-rotated-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');
delete(C);
view(2);
C = surf(281:380, -330:-231, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
savefigure(F, fullfile(renderPath, 'surfvel2', sprintf('surfvel2-%s-%i-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'surfvel2', sprintf('surfvel2-%s-%i-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'surfvel2', sprintf('surfvel2-%s-%i-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'surfvel2', sprintf('surfvel2-%s-%i-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');

% Compute signed norm.
velnorm = sqrt(sum(vel.^2, 2));
c = sign(velev) .* velnorm;

% Plot signed norm of surface velocity.
F = createFigure3(jet);
trisurf(S.F, S.Vs(:, 1), S.Vs(:, 2), S.Vs(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
view(3);
set(gca, 'XTick', -450:150:450);
set(gca, 'YTick', -450:150:450);
set(gca, 'ZTick', -450:150:450);
adjustFigure3;
savefigure(F, fullfile(renderPath, 'signnormsurfvel3', sprintf('signnormsurfvel3-%s-%i-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'signnormsurfvel3', sprintf('signnormsurfvel3-%s-%i-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'signnormsurfvel3', sprintf('signnormsurfvel3-%s-%i-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'signnormsurfvel3', sprintf('signnormsurfvel3-%s-%i-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');
% Rotate by pi.
[az, el] = view;
view(az + 180, el);
savefigure(F, fullfile(renderPath, 'signnormsurfvel3', sprintf('signnormsurfvel3-%s-%i-rotated-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'signnormsurfvel3', sprintf('signnormsurfvel3-%s-%i-rotated-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'signnormsurfvel3', sprintf('signnormsurfvel3-%s-%i-rotated-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'signnormsurfvel3', sprintf('signnormsurfvel3-%s-%i-rotated-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');

colorbar;
cbar = findobj(F, 'tag', 'Colorbar');
%set(cbar, 'YTick', -10:2:6);
set(cbar, 'TickLength', [.02 .02], 'YColor', [0 0 0]);
view(2);
adjustFigure3;
savefigure(F, fullfile(renderPath, 'signnormsurfvel2', sprintf('signnormsurfvel2-%s-%i-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'signnormsurfvel2', sprintf('signnormsurfvel2-%s-%i-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'signnormsurfvel2', sprintf('signnormsurfvel2-%s-%i-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'signnormsurfvel2', sprintf('signnormsurfvel2-%s-%i-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');

for k=idx

% Plot residual vector.
F = createFigure;
plot(0:length(E{k}.L.resvec)-1, E{k}.L.resvec/E{k}.L.rhs, 'b-');
if(isempty(E{k}.L.restart))
    pos = E{k}.L.iter(2);
else
    pos = (E{k}.L.iter(1)-1)*E{k}.L.restart+E{k}.L.iter(2);
end
plot(pos, E{k}.L.relres, 'rx');
text(pos, E{k}.L.relres,  sprintf('%0.5f', E{k}.L.relres), 'horizontal', 'right', 'vertical', 'bottom');
axis on;
adjustFigure;
savefigure(F, fullfile(renderPath, 'residual', sprintf('%s-%i.png', filename, k)));

% Plot coefficients.
F = createFigure;
bar(E{k}.u);
axis square;
adjustFigure;
savefigure(F, fullfile(renderPath, 'coefficients', sprintf('%s-%i.png', filename, k)));

% Recover vector field.
U = projecttoplane(E{k}.U);

% Compute colour space scaling.
nmax = max(sqrt(sum(U.^2, 2)));

% Save flow scaling.
fprintf(fid, 'Experiment %i: %.4f\n', k, nmax);

% Choose surface.
S = D.S{1};

% Create colourwheel.
cw = colourwheelbg;

% Plot optical flow.
F = createFigure3(cmap);
c = double(squeeze(computeColour(U(:, 1)/nmax, U(:, 2)/nmax))) ./ 255;
trisurf(S.F, S.Vs(:, 1), S.Vs(:, 2), S.Vs(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
C = surf(-350:-151, -350:-151, -380*ones(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
view(3);
set(gca, 'XTick', -450:150:450);
set(gca, 'YTick', -450:150:450);
set(gca, 'ZTick', -450:150:450);
adjustFigure3;
savefigure(F, fullfile(renderPath, 'flow3', sprintf('flow3-%s-%i-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'flow3', sprintf('flow3-%s-%i-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'flow3', sprintf('flow3-%s-%i-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'flow3', sprintf('flow3-%s-%i-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');
delete(C);
% Rotate by pi.
[az, el] = view;
view(az + 180, el);
C = surf(151:350, 151:350, -380*ones(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
savefigure(F, fullfile(renderPath, 'flow3', sprintf('flow3-%s-%i-rotated-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'flow3', sprintf('flow3-%s-%i-rotated-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'flow3', sprintf('flow3-%s-%i-rotated-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'flow3', sprintf('flow3-%s-%i-rotated-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');
delete(C);
view(2);
C = surf(281:380, -330:-231, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
savefigure(F, fullfile(renderPath, 'flow2', sprintf('flow2-%s-%i-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'flow2', sprintf('flow2-%s-%i-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'flow2', sprintf('flow2-%s-%i-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'flow2', sprintf('flow2-%s-%i-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');

% Recover total velocity.
M = projecttoplane(E{k}.U + vel);

% Compute colour space scaling.
nmax = max(sqrt(sum(M.^2, 2)));

% Save flow scaling.
fprintf(fid, 'Experiment %i total motion: %.4f\n', k, nmax);

% Plot total velocity.
F = createFigure3(cmap);
c = double(squeeze(computeColour(M(:, 1)/nmax, M(:, 2)/nmax))) ./ 255;
trisurf(S.F, S.Vs(:, 1), S.Vs(:, 2), S.Vs(:, 3), 'FaceColor', 'flat', 'FaceVertexCData', c, 'EdgeColor', 'none');
C = surf(-350:-151, -350:-151, -380*ones(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
view(3);
set(gca, 'XTick', -450:150:450);
set(gca, 'YTick', -450:150:450);
set(gca, 'ZTick', -450:150:450);
adjustFigure3;
savefigure(F, fullfile(renderPath, 'motion3', sprintf('motion3-%s-%i-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'motion3', sprintf('motion3-%s-%i-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'motion3', sprintf('motion3-%s-%i-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'motion3', sprintf('motion3-%s-%i-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');
delete(C);
% Rotate by pi.
[az, el] = view;
view(az + 180, el);
C = surf(151:350, 151:350, -380*ones(200, 200), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
savefigure(F, fullfile(renderPath, 'motion3', sprintf('motion3-%s-%i-rotated-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'motion3', sprintf('motion3-%s-%i-rotated-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'motion3', sprintf('motion3-%s-%i-rotated-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'motion3', sprintf('motion3-%s-%i-rotated-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');
delete(C);
view(2);
C = surf(281:380, -330:-231, zeros(100, 100), cw, 'FaceColor','texturemap', 'EdgeColor', 'none');
savefigure(F, fullfile(renderPath, 'motion2', sprintf('motion2-%s-%i-600dpi.png', filename, k)), '-png', '-r600');
savefigure(F, fullfile(renderPath, 'motion2', sprintf('motion2-%s-%i-1200dpi.png', filename, k)), '-png', '-r1200');
savefigure(F, fullfile(renderPath, 'motion2', sprintf('motion2-%s-%i-600dpi.jpg', filename, k)), '-jpg', '-r600', '-q100');
savefigure(F, fullfile(renderPath, 'motion2', sprintf('motion2-%s-%i-1200dpi.jpg', filename, k)), '-jpg', '-r1200', '-q100');

%close all;
end

% Close file.
fclose(fid);