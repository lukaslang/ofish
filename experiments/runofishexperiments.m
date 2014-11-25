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

% Import data.
disp('Loading precomputed data.');
name = 'cxcr4aMO2_290112';
path = fullfile('./', 'data', name, 'generated');

% Load file.
filename = 'frames-140-142-unfiltered-1-50-7';
D = load(fullfile(path, sprintf('dat-%s.mat', filename)));
G = load(fullfile(path, sprintf('gen-%s.mat', filename)));

% Create folder for results.
resultsPath = fullfile('./', 'results', name, 'ofish');
mkdir(resultsPath);

% Specify memory to use.
mem = 100e9;

% Set tolerance for relative residual.
tol = 1e-6;

% Set maximum number of iterations.
maxit = 100;

% Set range for alpha.
rng1 = [1e-3, 1e-2, 1e-1, 0, 1e0, 1e1, 1e2, 1e3];

% Run experiments.
run = 1;
runs = length(rng1);
E = cell(runs, 1);
for alpha=rng1
    fprintf('Computing flow %d/%d: alpha=%g.\n', run, runs, alpha);
    ticId = tic;
    [u, L] = ofishsolve(G.dim, G.A, G.D, G.b, alpha, tol, maxit);
    elapsedTime = toc(ticId);
    fprintf('Elapsed time %.6f seconds.\n', elapsedTime);

    % Store experiment.
    E{run}.u = u;
    E{run}.L = L;
    E{run}.alpha = alpha;

    run = run + 1;
end

disp('Recovering vector fields.');
ticId = tic;

% Choose surface.
S = D.S{1};

% Specify evaluation points.
xi = repmat([1/3, 1/3], size(S.F, 1), 1);

% Compute tangent basis on surface.
[Dx, Dy] = surftanBasis(S.F, S.V, S.rhon, xi);

% Save to experiments.
for k=1:runs
    % Compute synthesis.
    [Yx, Yy] = trivspharmsynth(D.N, S.F, S.V, E{k}.u, xi, mem);
    
    % Recover vector field on the surface.
    E{k}.U = bsxfun(@times, Yx, Dx) + bsxfun(@times, Yy, Dy);
    
    % Store evaluation points.
    E{k}.xi = xi;
end

elapsedTime = toc(ticId);
fprintf('Elapsed time %.6f seconds.\n', elapsedTime);

% Create filename.
wsFilename = sprintf('%s-%s.mat', datestr(now, 'yyyy-mm-dd-HH-MM-SS'), filename);
% Save workspace.
save(fullfile(resultsPath, wsFilename), 'E', '-v7.3');