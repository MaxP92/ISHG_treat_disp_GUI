function im_unwrapped = unwrap_2D_MP( im_phase,  im_mag,  plot_now, auto_sead, ch_map, ch_adv_plt, save_mat)
% im_unwrapped = unwrap_2D_MP( im_phase,  im_mag)
%
%
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
% 2016/04/22 : adapted by Maxime PINSARD
%
% took 383 sec = 6 min for a 500x200 matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clearvars ; close all; clc


% dialog box to crop or multiply by pi etc.
addpath('C:\Users\pc\Documents\These\codes Matlab\Codes_I-SHG\MP\QualityGuidedUnwrap2D_r1');

im_phase = im_phase*pi; % has to be in [-pi, pi]

%% Replace with your mask (if required)
mag_max = max(im_mag(:));
indx1 = im_mag < 0.1*mag_max;  %Intensity = mag^2, so this = .04 threshold on the intensity
im_mask = ones(size(im_phase));
im_mask(indx1) = 0;                  %Mask
if(~exist('im_mask','var'))
    im_mask = ones(size(im_phase));          %Mask (if applicable)
end
% figure; imagesc(im_mag.*im_mask),   colormap(gray), axis image, axis off, title('Initial masked magnitude'); colorbar;
% figure; imagesc(im_phase.*im_mask), colormap(gray), axis image, axis off, title('Initial masked phase'); colorbar;

im_unwrapped = nan(size(im_phase));        %Initialze the output unwrapped version of the phase
adjoin = zeros(size(im_phase));            %Zero starting matrix for adjoin matrix
unwrapped_binary = zeros(size(im_phase));  %Binary image to mark unwrapped pixels

%% Calculate phase quality map
im_phase_quality = PhaseDerivativeVariance_r1(im_phase);

%% Automatically (default) or manually identify starting seed point on a phase quality map
minp = im_phase_quality(2:end-1, 2:end-1); minp = min(minp(:));
maxp = im_phase_quality(2:end-1, 2:end-1); maxp = max(maxp(:));
% msgbox('You should verify that your phase is in [-pi, pi] !')
% choice1 = menu('Start Point (seed)', 'Manual', 'auto')-1;

if ~auto_sead   % Chose starting point interactively
    figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), colorbar, axis image, axis off; title('Phase quality map');
    uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
    [xpoint,ypoint] = ginput(1);                %Select starting point for the guided floodfill algorithm
    colref = round(xpoint);
    rowref = round(ypoint);
    close;  % close the figure;
else   % Chose starting point = max. intensity, but avoid an edge pixel
    
    [rowrefn,colrefn] = find(im_mag(2:end-1, 2:end-1) >= 0.99*mag_max);
    if (isempty(rowrefn) || isempty(colrefn))
        mag_max =  max(max(im_mag(2:end-1, 2:end-1)));
        [rowrefn,colrefn] = find(im_mag(2:end-1, 2:end-1) >= 0.99*mag_max);
    end
    rowref = rowrefn(1)+1; % choose the 1st point for a reference (known good value)
    colref = colrefn(1)+1; % choose the 1st point for a reference (known good value)
end

%% Unwrap
t1=tic;
disp('Unwrap has begun !')
hb = msgbox('Check the console to see the progression of the unwrapping !');
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

close(hb)

im_phase = im_phase/pi; % has to be in [-pi, pi]
im_unwrapped = im_unwrapped/pi; % has to be in [-pi, pi]
%% corr

if save_mat
    
    save(['im_unwrapped_' num2str(t1) '.mat'], 'im_unwrapped');
end

switch plot_now 
    
    case 1
    %% Plot images
    %figure; imagesc(im_mag),       colormap(gray), colorbar, axis square, axis off; title('QG Magnitude image');
    %figure; imagesc(im_phase),     colormap(gray), colorbar, axis square, axis off; title('QG Wrapped phase');
        
    if ch_map
       
        c_d = 'jet';
        
    else
        c_d = 'hsv';
        
    end

    figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), axis image, axis off, title('QG Phase quality map'); colorbar;
    screensize = get( 0, 'Screensize' ); % to get screen size
    fact = 4/5;
    left_offset_fig = 40;
    top_offset_fig = 100;
    h_f=figure(5); set(h_f,'Color', [1 1 1], 'outerposition',...
        [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
        screensize(3)*fact screensize(4)*fact]);
    if size(im_phase,1)> 9/16*size(im_phase,2)
        h1=subplot(1,2,1);
    else
        h1=subplot(2,1,1);
    end
    imagesc(h1, im_phase); colormap(h1, hsv), axis image, axis off; title('Wrapped phase');
    hb1=colorbar; title(hb1, 'Phase (/pi)')
    
    if size(im_phase,1)> 9/16*size(im_phase,2)
        h2 = subplot(1,2,2);
    else
        h2 = subplot(2,1,2);
    end
%     im_unwrapped(isnan(im_unwrapped)) = 1.05;
    imagesc(h2, im_unwrapped); colormap(h2, c_d), axis image, axis off; title('QG Unwrapped phase');
    hb2=colorbar; title(hb2, 'Phase (/pi)')
    
%     ch_adv_plt = menu('Advanced plot : rescaling without extremes', 'Go on', 'Exit')-1;
    
    if ~ch_adv_plt
        return
    end
    
    if ch_map
        c_d = 'jet';
        cc = colormap(c_d); cc(end,:) = [0, 1, 0];
        cc(1,:) = [0, 0, 1];
    else
        c_d = 'hsv';
        cc = colormap(c_d); cc(end,:) = [0, 0, 0];
        cc(1,:) = [1, 1, 1];
    end
    
    avg_nearest_neighbors = 2 - menu('avg_nearest_neighbors', 'Yes', 'No');
    
    if avg_nearest_neighbors
        im_unwrapped = conv2(im_unwrapped, ones(3)/9, 'same');
    end
    
    figure;
    if size(im_unwrapped, 1)> 9/16*size(im_unwrapped, 2)
        h1 = subplot(1,2,1);
    else
        h1 = subplot(2,1,1);
    end
    imagesc(im_unwrapped); colormap(h1,c_d), axis image, axis off; title('QG Unwrapped phase');
    hb2=colorbar; title(hb2, 'Phase (/pi)')
    
    pas_content = 1;
    while pas_content == 1
        prompt = { 'Inf. bounds scaling (in /pi)', 'Sup. bounds scaling (in /pi)'};
        dlg_title = 'Paramètres';
        num_lines = 1;
        def = {'-7', '7'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        if isempty(answer)
            answer = def;
        end
        lim_inf = str2double(answer{1}); 
        lim_sup = str2double(answer{2}); 


    %     lim_sup = 7;
    %     lim_inf = -lim_sup;
        im2 = im_unwrapped;
        im2(im2 > lim_sup) = lim_sup + .05;
        im2(im2 < lim_inf) = lim_inf - .05;

        if size(im2, 1)> 9/16*size(im2,2)
            h2 = subplot(1,2,2);
        else
            h2 = subplot(2,1,2);
        end
        imagesc(im2); colormap(h2,cc), axis image, axis off; title('QG Unwrapped phase scaled');
        hb2=colorbar; title(hb2, 'Phase (/pi)')

        pas_content = menu('Scaling ok?', 'Yes', 'No')-1;
% %         if pas_content
    end
end

end

