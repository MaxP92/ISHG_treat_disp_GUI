% I_SHG_2.m : script for MAIN code, alternative to the GUI
%
% CHARLES-ANDRÉ COUTURE
% 2014-02-24
%
% Edited by Maxime Rivard
% 2015-06-01
%
% Re-shaped, completed, edited by Maxime PINSARD
% 2015-10-23
%
% Analyse stack of iSHG images to extract the relative phase for each pixel
% (optionally averaged on a 3x3 pixel region)
%
% Needs fit_I_SHG.m to work
%
% V14 2016.6.9

clearvars; close all; clc

%% Initialisation et entrée des paramètres

% pathfun = 'C:\Users\Stephane\Documents\Postdoc INRS\Data\I-SHG galvos\'; % main path, you can change it to browse faster (optional)
pathfun ='C:\Users\Maxime\Documents\These\codes Matlab\Codes_I-SHG\MP'; addpath(pathfun);

% Valeur max (et min) de la colorbar pour l'affichage des images de
% soustraction et du stack qui sert à faire le fit pour trouver la phase
contrast = 0.7;

% Mémoriser le folder matlab de départ
current_folder = pwd;

%% Plot of the 1st image, image reading

% L'utilisateur sélectionne le fichier d'images qu'il veut analyser
folder_name=uigetdir(pathfun,'Select your path containing your file(s) !');
if isa(folder_name, 'char') % a folder has been chosen
    cd(folder_name);
end

[fname, folder_name, FILTERINDEX] = uigetfile('*.tif','Select your stack of images !','stack.tif');
cd(folder_name);
if ~FILTERINDEX
    error('File not chosen ! Program ends.');
else
    fprintf('%s\n', fname);
end

% Obtenir l'info sur le stack d'images OLD VERSION
% info = imfinfo(fname);
% num_images = numel(info);
% A1_first = imread(fname,1);
% to be removed

%% user choose parameters

prompt = {'Résolution en x (um):','Résolution en y (um):','Phase de départ (en degrés):',...
    '2nd phase (en degrés):', '3rd phase :', '4th phase :', '2nd step', 'Nombre de répétitions de chaque image'...
    'Nbins for hist in phase', 'Nbins for hist in interf. contrast', 'Mode contrast? (=1, otherwise raw)',...
'Use system inversion', 'method calc phase : 1 for classic, 0 for classic with 1D matrices, 2 for algo 3 phases, 3 for algo 3 phases with tilt corr., 4 for algo 3 phases with tilt + vib. corr.'};
dlg_title = 'Paramètres';
num_lines = 1;
% Valeurs par défaut
def = {'0.2','0.2','0','180','360','30','30','1', '250', '100', '1','1', '1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
% Les réponses en caractères sont converties en chiffres qui sont enregistrés dans des variables.
resx = str2double(cell2mat(answer(1)));
resy = str2double(cell2mat(answer(2)));
start_phase = str2double(cell2mat(answer(3)));
second_phase = str2double(cell2mat(answer(4)));
third_phase = str2double(cell2mat(answer(5)));
fourth_phase = str2double(cell2mat(answer(6)));
snd_step = str2double(cell2mat(answer(7)));
avg_img = str2double(cell2mat(answer(8)));
% Intervalle pour le nombre de bins en x et y pour les histogrammes -- IMPORTANT
int_x = 1/ str2double(cell2mat(answer(9)));
int_y = str2double(cell2mat(answer(10)));
contr_mode = str2double(cell2mat(answer(11)));
use_invA = 1-str2double(cell2mat(answer(12)));
method = str2double(cell2mat(answer(13)));
  % 2 for CA, 1 for Rivard, 4 for Rivard with 1D matrices, 3 for algo 3 phases

% Menus pour les différents choix
analyse = menu('Type d''analyse','Phase relative','Vérifier la variation en cos', 'Default choices (Rel. phase, Hist. interf. contrast, no filt., no obj., no save, eng.)');

if analyse~=3 % pas default
    shg = menu('Histo2D','Intensité SHG','Contraste interféro');
    crop = menu('Région à analyser','Toute l''image','Sélectionner une région')-1;
    
    filt=menu('fltrage de l''image','Oui','Non');
    
    obj_test = menu('Y-a-t''il des régions à enlever de l''image (e.g. des cellules)?','Oui','Non');
    save_test = menu('Faut-il sauver chaque image du stack de différences (Il faudra indiquer où les sauver) ?','Oui','Non');
    langue = menu('Langue du texte des figures','Eng','Fr');
    
else % default
    analyse =1;
    shg = 2;
    crop = 0;
    filt = 2;
    obj_test = 2;
    save_test = 2;
    langue = 1;
    
end

close(gcf);
undocked_fig = 1; h = 1; % plot fig outside

ch = menu ('As it is', 'Normalization', 'Nearest neighbor', 'Nearest neigh + norm.');

if ch > 1
    if ch == 2
        near_avg_ini = 0; norm_ini = 1;
    elseif ch == 3
        near_avg_ini = 1; norm_ini = 0;
    elseif ch == 4
        near_avg_ini = 1; norm_ini = 1;
    else
        
        near_avg_ini = 0; norm_ini = 0;
    end
end

%% Load stack and plot first image

[ img_3D, num_images, xTitle_dflt, yTitle_dflt, phi_mat_default, img_shg, ...
    screensize, fact, left_offset_fig, top_offset_fig, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, cmap_redgreen, cmap_blkredblkgrn, cmap_brbgsat, cmap_hsvsat, h,...
 order_frame, start_phase, diff_phase, x_phase, nonequal_spacing] = ...
    load_stack_plot_ISHG( fname, shg, undocked_fig, h, ...
 start_phase, second_phase, third_phase, fourth_phase, snd_step, ...
    contr_mode, avg_img, near_avg_ini, norm_ini );


%% Soustraction des images GSH interférométriques pour obtenir les images du contraste interférométrique

[ x, y, contr, rect, img_shg , x_phase ] = contrast_by_subtract_img_ISHG( img_3D, num_images, filt, crop, resx, resy, shg, img_shg, undocked_fig, h, ...
screensize, fact, left_offset_fig, top_offset_fig, order_frame, contr_mode, diff_phase, x_phase, use_range  );
% function to calculate contrast by subtraction of images

%% Exclusions de certaines portions de l'image pour le calcul

if obj_test == 1
    prompt = {'Nombre d''objets à enlever de l''analyse :'};
    dlg_title = 'Input';
    num_lines = 1;
    answer = inputdlg(prompt,dlg_title,num_lines);
    nbr_objects = str2double(answer{1});
else
    nbr_objects =0;
end

[ contr, img_shg ] = exclude_area_ISHG( contr, obj_test, crop, shg, img_3D, rect, img_shg, nbr_objects, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig  );
% function to exclude the zones

%% Réorganisation des images du contraste interférométrique dans un ordre croissant + plot

if contr_mode
    h = diff_phase_stack_ISHG( contr, x, y, start_phase, diff_phase, contrast, xTitle_dflt, yTitle_dflt, ...
        screensize, fact, left_offset_fig, top_offset_fig, undocked_fig, h, axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, cmap_redgreen, 1, 4, x_phase, nonequal_spacing );
end
% function to calculate and plot the difference for the different angles

%% Sauvegarde du stack d'images du contraste interférométrique

% Pour sauver le stack de différences

if save_test==1
    save_results_stack_ISHG( contr, x, y, contrast, cmap_redgreen, xTitle_dflt, yTitle_dflt, screensize, fact, left_offset_fig, top_offset_fig,...
        axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );
    % function to save the results image of inferferometric contrast in tiff images
end

%% Analyse du contraste sinusoidal des images du contraste interférométrique

% Si on veut simplement analyser le comportement en cos
if analyse==2
    
    short_version = 0;
    [ x, y, complete, test, model, r2 ] = analyse_cos_ISHG( contr, img_3D, resx, resy, rect, start_phase, diff_phase, crop, xTitle_dflt, yTitle_dflt, screensize, fact, left_offset_fig, top_offset_fig, ...
        axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, short_version, x_phase );
    % function to analyze the behavior in cosinus
    
else
    %% Analyse complète des images du contraste interférométrique, pour trouver la phase
    
    % On revient à l'ancien dossier - voir au tout début du code
    cd(current_folder);
    
    test_temps = menu('Phase analysis ?','Stop here','Calc. phase as it is', ...
'Calc. phase with avg neighbour averaging', 'Calc. phase with median neighbour averaging');
    
    test_avg = test_temps - 2;
    
    if test_temps==1
        disp('User chose to stop here.');
        return;
    end
    
    [ phase1, amp, err, phase_test_hist, stringb1, stringc1, stringb2, stringc2, ...
    img, title_hist1, yaxis_hist1, Titre4, Counts, Legendehisto, titl_box, Titre1] = ...
        phase_calc_ISHG( contr, start_phase, diff_phase, langue, obj_test, x, y, xTitle_dflt, yTitle_dflt, phi_mat_default, shg, img_shg,...
        screensize, fact, left_offset_fig, top_offset_fig, ...
        axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, cmap_blkredblkgrn, cmap_brbgsat, test_avg,...
        x_phase, use_invA, method);
    %  To calculate the phase map, with complete analyze of interferometric
    %  contrast images
    
    %% Calcul pour les histogrammes
    
    [ ph_hist, amp_hist, data_tot, ylabel_hist3, title_hist3, ~ ] = hist_tot_ISHG( phase1, img, ...
        int_x, int_y, phi_mat_default, Counts, Titre4, title_hist1, yaxis_hist1,...
        axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig );
    %  To calculate and plot the histograms (2D and 3D)
    
    %% Analyse région spécifique
    
    max_amp = max(max(data_tot(:,2)));
    fprintf('Maximum of the amplitude calculated = %g \n', max_amp);
    
    [pb, ph, keep_cond, offset_pi2] = choosedialog(titl_box, max_amp);
    
    if (~keep_cond || pb ~= 0 || ph ~= max_amp || offset_pi2)
        
        offset = 0; % offset of 0 default
        
        [phase_test_hist, second_hist, phase_test, h] = map_spec_region_ISHG( x, y, phase1, amp, img_shg, shg, ...
            xTitle_dflt, yTitle_dflt, Titre1, phi_mat_default, pb, ph,...
            axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, ...
            undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig, cmap_brbgsat );
        % To calculate and plot the histogram 2D and 3D for a specific range of
        % intensity
        
        ph_hist_current00 = phase_test_hist;
        
        if offset_pi2
            offset = 0.5; % offset of pi/2
            phase_test_hist = Offset_ph_hist( phase_test_hist, ph_hist_current00, offset);
        end
        
        %% 3D hist
        
        [ ~, ~, ~, ~,h] = hist_3D_ISHG(  phase_test_hist, second_hist, int_x, int_y, title_hist3, ylabel_hist3,...
            Counts, phi_mat_default, axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, ...
            screensize, fact, left_offset_fig, top_offset_fig, offset_pi2 );
        
        %% 2D hist
        
        hhist1 = 0;
        
        [hhist1, hist2_xdata, hist2_ydata] = hist2D_ISHG( fact, left_offset_fig, ...
            top_offset_fig, phase_test_hist, int_x, phi_mat_default, Counts, Titre4, ...
            axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, screensize, offset_pi2, hhist1 );
        
        
        fprintf(2, '\nYou should implement neutral offset (as in GUI) in me !\n');
        fprintf(2, '\nYou should implement neutral offset (as in GUI) in me !\n');
        
        %% fit hist 2D
        M = 0; hfun2=100;
        loop_plz =1 ;
        while loop_plz
            fit_chosen = menu('Fit of phase distr.', 'Gaussian MANU', 'Lorentzian MANU', 'Product of both MANU', 'pseudo-Voigt (sum) MANU',  'Pearson VII (MANU)', 'Von Mises (MANU)',...
                '2 Gauss', '2 Lorentz','2 Prod. Gaussian/Lorentzian', '2 Sum Gaussian/Lorentzian (pseudo-Voigt)','2 Pearson VII', '2 Von Mises', 'STOP');
            
            if fit_chosen == 13
                loop_plz = 0;
            else
                if (fit_chosen == 5 || fit_chosen == 11)
                    
                    prompt = { 'M = (1 : Lorentz, Inf : Gauss, < 1 : super-Lorentz '};
                    dlg_title = 'M for Pearson VII';
                    num_lines = 1;
                    % Valeurs par défaut
                    def = {'0.5'};
                    answer = inputdlg(prompt,dlg_title,num_lines,def);
                    M = str2double(cell2mat(answer(1)));
                end
                [ stringb1, stringc1, stringb2, stringc2, cmap_brbg_fit ] = fit_hist2d_ishg( hist2_xdata, hist2_ydata, hhist1, Legendehisto, offset, fit_chosen, M );
                
                 fit_chosen = menu('Would you like the fit width to be applied on colormap of phase map ?', ...
                'Yes', 'No');
            hfun2 = pl_phasemap_ISHG(1, screensize, fact, left_offset_fig, ...
                top_offset_fig, hfun2, phase_test, x, y, ...
                xTitle_dflt, yTitle_dflt, Titre1, cmap_brbg_fit, phi_mat_default, ...
                axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz);
                
            end
        end
        
    end
    %end
    
end

% On revient à l'ancien dossier - voir au tout début du code
cd(current_folder);