% SLAM_essai

clear all; clc;

pathfun = 'C:\Users\Maxime\Documents\These\DATA images\2015-11-30_PPLN_SFG';

% L'utilisateur sélectionne le fichier d'images qu'il veut analyser
folder_name=uigetdir(pathfun,'Select your path containing your file(s) !');
if isa(folder_name, 'char') % a folder has been chosen
    cd(folder_name);
end

[fname, folder_name, FILTERINDEX] = uigetfile('*.tif','Select your gaussian image !','stack.tif');
if ~FILTERINDEX
    error('File not chosen ! Program ends.');
end

[fname2, folder_name2, FILTERINDEX] = uigetfile('*.tif','Select your donut image !','stack.tif');
if ~FILTERINDEX
    error('File not chosen ! Program ends.');
end

gauss_img00 = double(imread(fullfile(folder_name, fname)));

donut_img00 = double(imread(fullfile(folder_name2, fname2)));

max_abs = max(max(max(gauss_img00)), max(max(donut_img00)));
min_abs = min(min(min(gauss_img00)), min(min(donut_img00)));

gauss_img = (gauss_img00 - min(min(gauss_img00)))*(max_abs - min_abs)/(max(max(gauss_img00))-min(min(gauss_img00))) + min_abs;

donut_img = (donut_img00 - min(min(donut_img00)))*(max_abs - min_abs)/(max(max(donut_img00))-min(min(donut_img00))) + min_abs;

close all;

figure; subplot(1, 2, 1); imagesc(gauss_img); colormap gray; title('Gaussian'); colorbar
subplot(1, 2, 2); imagesc(donut_img); colormap gray; title('Donut'); colorbar

for r = 1.5:-0.3:0

results_img = gauss_img - r*donut_img;

results_img(results_img < 0) = 0;

figure; imagesc(results_img); colormap gray; colorbar; 
title(sprintf('Gauss - %g*donut', r));

end
