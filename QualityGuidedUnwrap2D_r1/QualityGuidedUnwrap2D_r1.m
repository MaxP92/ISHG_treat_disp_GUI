%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QualityGuidedUnwrap2D implements 2D quality guided path following phase
% unwrapping algorithm.
%
% Inputs: 1. Complex image in .mat double format
%         2. Binary mask (optional)          
% Outputs: 1. Unwrapped phase image
%          2. Phase quality map
%
% This code can easily be extended for 3D phase unwrapping.
% Technique adapted from:
% D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
% Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
%
% Posted by Bruce Spottiswoode on 22 December 2008
% 2010/07/23  Modified by Carey Smith
% 
% 2016/04/22 : adatapted by Maxime PINSARD 
%
% took 383 sec = 6 min for a 500x200 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars ; close all; clc

%% Replace with your images

path00 = 'C:\Users\pc\Documents\';
addpath(fullfile(path00, 'These\codes Matlab\PhaseUnwrapping2D\QualityGuidedUnwrap2D_r1'));
addpath(fullfile(path00, 'These\codes Matlab\reduction_nb_img'));

% choice = menu('Data', 'Complex image', 'amp + phase image');

[fname, folder_name, FILTERINDEX] = uigetfile('*.mat','Select your mat !','im.mat');
cd(folder_name);
if ~FILTERINDEX
    error('File not chosen ! Program ends.');
else
    fprintf('%s\n', fname);
end%Load complex image
ff = load(fname); % phase image
fff=struct2cell(ff);
IM = fff{1};
choice = 2; % not complex image

if choice == 1
    
    im_mag=abs(IM);                             %Magnitude image
    im_phase=angle(IM); %Phase image
else
    im_phase = double(IM); %Phase image
    
%     im_phase = pi*double(IM); %Phase image
%     fprintf(2, 'pi multiplied !\n');
    
    [fname, folder_name, FILTERINDEX] = uigetfile('*.tif','Select your amp image (can be in stack) !','amp.tif');
    cd(folder_name);
    if ~FILTERINDEX
        error('File not chosen ! Program ends.');
    else
        fprintf('%s\n', fname);
    end%Load complex image
    
    im_mag = imread(fname); % load 1st frame
    im_mag  = im2double( im_mag );
    
end

% dialog box to crop or multiply by pi etc.

[minR , maxR , minC, maxC, keep_cond, mult_cond] = choosedialog_unwrap(size(im_phase, 1), size(im_phase, 2));

if ~keep_cond
    im_phase = im_phase(minR:maxR, minC: maxC)*((pi-1)*mult_cond + 1); %round(size(im_phase,2)*5/8) % 625
    im_mag = im_mag(minR:maxR, minC: maxC)*((pi-1)*mult_cond + 1);
    IM = IM(minR:maxR, minC: maxC)*((pi-1)*mult_cond + 1);
    % menu('WARNING  : Size of image reduced !', 'Ok');
end

if (mult_cond && keep_cond) % multiply by pi
    im_phase = im_phase*pi;
end

%% Replace with your mask (if required)
mag_max = max(im_mag(:));
indx1 = find(im_mag < 0.1*mag_max);  %Intensity = mag^2, so this = .04 threshold on the intensity
im_mask = ones(size(IM));
im_mask(indx1) = 0;                  %Mask
if(~exist('im_mask','var'))
  im_mask = ones(size(IM));          %Mask (if applicable)
end
figure; imagesc(im_mag.*im_mask),   colormap(gray), axis image, axis off, title('Initial masked magnitude'); colorbar;
figure; imagesc(im_phase.*im_mask), colormap(gray), axis image, axis off, title('Initial masked phase'); colorbar;

im_unwrapped = nan(size(IM));        %Initialze the output unwrapped version of the phase
adjoin = zeros(size(IM));            %Zero starting matrix for adjoin matrix
unwrapped_binary = zeros(size(IM));  %Binary image to mark unwrapped pixels

%% Calculate phase quality map
im_phase_quality = PhaseDerivativeVariance_r1(im_phase);

%% Automatically (default) or manually identify starting seed point on a phase quality map 
minp = im_phase_quality(2:end-1, 2:end-1); minp = min(minp(:));
maxp = im_phase_quality(2:end-1, 2:end-1); maxp = max(maxp(:));
msgbox('You should verify that your phase is in [-pi, pi] !')
choice1 = menu('Start Point (seed)', 'Manual', 'auto')-1;

if ~choice1   % Chose starting point interactively
  figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), colorbar, axis image, axis off; title('Phase quality map'); 
  uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
  [xpoint,ypoint] = ginput(1);                %Select starting point for the guided floodfill algorithm
  colref = round(xpoint);
  rowref = round(ypoint);
  close;  % close the figure;
else   % Chose starting point = max. intensity, but avoid an edge pixel
  [rowrefn,colrefn] = find(im_mag(2:end-1, 2:end-1) >= 0.99*mag_max);
  rowref = rowrefn(1)+1; % choose the 1st point for a reference (known good value)
  colref = colrefn(1)+1; % choose the 1st point for a reference (known good value)
end

%% Unwrap
t1=tic;
disp('Unwrap has begun !')
im_unwrapped(rowref,colref) = im_phase(rowref,colref);                          %Save the unwrapped values
unwrapped_binary(rowref,colref,1) = 1;
if im_mask(rowref-1, colref, 1)==1;  adjoin(rowref-1, colref, 1) = 1; end       %Mark the pixels adjoining the selected point
if im_mask(rowref+1, colref, 1)==1;  adjoin(rowref+1, colref, 1) = 1; end
if im_mask(rowref, colref-1, 1)==1;  adjoin(rowref, colref-1, 1) = 1; end
if im_mask(rowref, colref+1, 1)==1;  adjoin(rowref, colref+1, 1) = 1; end
im_unwrapped = GuidedFloodFill_r1(im_phase, im_mag, im_unwrapped, unwrapped_binary, im_phase_quality, adjoin, im_mask);    %Unwrap
toc(t1);

im_unwrapped(1, :) =0; im_unwrapped(end, :)=0; % remove the edges
im_unwrapped(:, 1) =0; im_unwrapped(:, end)=0; % !!!

save(['im_unwrapped_' num2str(t1) '.mat'], 'im_unwrapped');

%% Plot images
%figure; imagesc(im_mag),       colormap(gray), colorbar, axis square, axis off; title('QG Magnitude image'); 
%figure; imagesc(im_phase),     colormap(gray), colorbar, axis square, axis off; title('QG Wrapped phase'); 

ch = menu('Cmap', 'jet', 'hot', 'hsv (to show no changes)')-1;

switch ch 
    
    case 1
    c_d = hot;

    case 0
    c_d = jet;
    case 2
      c_d = hsv;  
end


figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), axis square, axis off, title('QG Phase quality map'); colorbar;
screensize = get( 0, 'Screensize' ); % to get screen size
fact = 4/5;
left_offset_fig = 40;
top_offset_fig = 100;
h_f=figure(5); set(h_f,'Color', [1 1 1], 'outerposition',...
    [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
    screensize(3)*fact screensize(4)*fact]);
if size(im_phase,1)> size(im_phase,2)
    h1=subplot(1,2,1);
else
    h1=subplot(2,1,1);
end
imagesc(im_phase); colormap(h1, hsv), axis image, axis off; title('Wrapped phase');
hb1=colorbar; title(hb1, 'Phase (rad)')

if size(im_phase,1)> size(im_phase,2)
    h2 = subplot(1,2,2);
else
    h2 = subplot(2,1,2);
end
imagesc(im_unwrapped); colormap(h2, c_d), axis image, axis off; title('QG Unwrapped phase');
hb2=colorbar; title(hb2, 'Phase (rad)')

ch00 = menu('Part for advanced plot', 'Go on', 'Exit')-1;

if ch00
    return
end

%% corr

ch = menu('Cmap', 'jet', 'hot')-1;

if ch
    c_d = hot;
    cc = colormap(c_d); cc(end,:) = [0, 1, 0];
    cc(1,:) = [0, 0, 1];
else
    c_d = jet;
    cc = colormap(c_d); cc(end,:) = [0, 0, 0];
cc(1,:) = [1, 1, 1];
end

avg_nearest_neighbors = 2 - menu('avg_nearest_neighbors', 'Yes', 'No');

if avg_nearest_neighbors
    im_unwrapped = conv2(im_unwrapped, ones(3)/9, 'same');
end

figure;
h1=subplot(2,1,1);
imagesc(im_unwrapped); colormap(h1,c_d), axis image, axis off; title('QG Unwrapped phase');
hb2=colorbar; title(hb2, 'Phase (rad)')

lim_sup = 7;
lim_inf = -lim_sup;
im2 = im_unwrapped;
im2(im2 > lim_sup) = lim_sup + .05;
im2(im2 < lim_inf) = lim_inf - .05;

h2=subplot(2,1,2);
imagesc(im2); colormap(h2,cc), axis image, axis off; title('QG Unwrapped phase scaled');
hb2=colorbar; title(hb2, 'Phase (rad)')

%% try

h1=subplot(2,2,1); imagesc(im_phase(:, 625:end)); colormap(h1, jet), axis image, axis off; 
title(sprintf('1. QG Unwrapped phase \nzoomed'));
hb2=colorbar; title(hb2, 'Phase (rad)'); set(h1, 'FontSize', 16);
h2=subplot(2,2,3); imagesc(im_unwrapped(:, :)); colormap(h2, jet), axis image, axis off; 
title(sprintf('2. QG Unwrapped phase \n only on ROI'));
hb2=colorbar; title(hb2, 'Phase (rad)') % to avoid edges
set(h2, 'FontSize', 16);
h3=subplot(2,2,2); imagesc(im_unwrapped(10:end-10, 10:end-10)); colormap(h3, jet), axis image, axis off; 
title(sprintf('3. QG Unwrapped phase \n only on ROI, edges removed'));
hb3=colorbar; title(hb3, 'Phase (rad)') % to avoid edges
set(h3, 'FontSize', 16);

im2=im_unwrapped(10:end-10, 10:end-10);
lim_sup = 2.3;
lim_inf = -lim_sup;

im2(im2 > lim_sup) = lim_sup + .5;
im2(im2 < lim_inf) = lim_inf - .5;
c_d = jet;
fact0=1;
cc = colormap(c_d); cc(end-fact0+1:end,:) = repmat([0, 0, 0],fact0,1);
cc(1:fact0,:) = repmat([1, 1, 1], fact0,1);

h4=subplot(2,2,4); imagesc(im2); colormap(h4, cc), axis image, axis off; 
title(sprintf('4. QG Unwrapped phase \n only on ROI, edges \n removed, re-scaled'));
hb4=colorbar; title(hb4, 'Phase (rad)') % to avoid edges
set(h4, 'FontSize', 16);
