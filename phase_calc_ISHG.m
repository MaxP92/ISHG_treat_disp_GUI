function [ phase1, amp, err, phase_test_hist, stringb1, stringc1, stringb2, stringc2, img, title_hist1, yaxis_hist1, Titre4, Counts, Legendehisto, titl_box, Titre1 ] = phase_calc_ISHG( contr, start_phase, diff_phase, ...
    langue, obj_test, x, y, xTitle_dflt, yTitle_dflt, phi_mat_default, shg, img_shg, screensize, fact, left_offset_fig, top_offset_fig, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, cmap_blkredblkgrn, cmap_phase_sat, test_avg, x_phase, use_invA, method )
% [ phase1, amp, err, phase_test, stringb1, stringc1, img, title_hist1, yaxis_hist1, Titre4, Counts, Legendehisto, titl_box, Titre1  ] = phase_calc_ISHG( contr, start_phase, diff_phase, langue, obj_test, x, y, xTitle_dflt, yTitle_dflt, phi_mat_default, shg, img_shg, screensize, fact, left_offset_fig, top_offset_fig, ...
% axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, cmap_blkredblkgrn, cmap_phase_sat, test_avg, x_phase, use_invA, method )
%
% 2015-10-27 edited by Maxime PINSARD
%
%  To calculate the phase map, with complete analyze of interferometric
%  contrast images

[ Titre1, Titre2, Titre3, Titre4, Titre5, Titre6, Counts, Yaxis1, Yaxis2, Legendehisto ] = string_axis_ISHG( langue );
%  The strings for xlabel, title, lengends etc. of axis

phase_test_hist = 0; % init for nothing
stringb1 = 0; stringc1 = 0;
stringb2 = 0; stringc2 = 0;

%% Calculation and display of the phase map

weight = 0; coef_centre = 5; coef_cross = 2; coef_diag = 1; sigma_med = 0;
    mean_phase_meth = 0; 

if test_avg == 1 % mean
    
    weight = 1;
    
    prompt = { 'Centre coeff', 'Cross coeff', 'Diag coeff'};
    dlg_title = 'Average on adjacent pixels';
    num_lines = 1;
    % Valeurs par défaut
    def = {'1', '0', '0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    % Les réponses en caractères sont converties en chiffres qui sont enregistrés dans des variables.
    coef_centre = str2double(cell2mat(answer(1)));
    coef_cross = str2double(cell2mat(answer(2)));
    coef_diag = str2double(cell2mat(answer(3)));
    mean_phase_meth = 1; % mean
    
elseif test_avg == 2 % median
        weight = 1;
    
    prompt = {'Sigma value for the exponential median filter (default 25.5, 0=unweighted filter)'};
    dlg_title = 'Weight';
    num_lines = 1;
    % Valeurs par défaut
    def = {'0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    sigma_med = str2double(cell2mat(answer(1)));
    mean_phase_meth = 0; % median
    
end

[ phase1, amp, err, ~ ] = phase_map_ISHG( contr, x, y, start_phase, diff_phase,...
    xTitle_dflt, yTitle_dflt, phi_mat_default, Titre1, Titre2, Titre3, obj_test, weight, ...
    coef_centre, coef_cross, coef_diag, screensize, fact, left_offset_fig, top_offset_fig, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h,...
 cmap_blkredblkgrn, cmap_phase_sat, x_phase, use_invA, method ,mean_phase_meth, sigma_med);
% On calcule la phase pour chaque pixel (moyenne sur 9 pixels) - voir fit_I_SHG.m

if shg == 1
    img = img_shg;
    title_hist1 = Titre5;
    yaxis_hist1 = Yaxis1;
    titl_box = 'SHG intensity';
else
    img = amp;
    title_hist1 = Titre6;
    yaxis_hist1 = Yaxis2;
    titl_box = 'interf. contrast';
end

% !!!!!!!!!!!!!! cut here

end

