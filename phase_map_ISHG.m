function [ phase, amp, err, h, divers_vect ] = phase_map_ISHG( contr, x, y, xTitle_dflt, yTitle_dflt, phi_mat_default, Titre1, Titre2, Titre3, ...
    obj_test, weight, coef_centre, coef_cross, coef_diag, screensize, fact, left_offset_fig, top_offset_fig, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, cmap_default, cmap_obj_test, ...
x_phase, use_invA, method, mean_meth, sigma_med, parallel_comput, decomp_LU, double_contrast, contr_mode, ...
skip_asking_coeff_median, k_max, epsilon, epsilon1, reit_prev_phase, prev_res, ramp_phshft, mask_val_sat, avg_on_ph )
% [ phase, amp, err, h ] = phase_map_ISHG( contr, x, y, start_phase, diff_phase, xTitle_dflt, yTitle_dflt, phi_mat_default, Titre1, Titre2, Titre3, ...
% obj_test, weight, coef_centre, coef_cross, coef_diag, screensize, fact, left_offset_fig, top_offset_fig, ...
% axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, cmap_blkredblkgrn, cmap_obj_test, x_phase, inv_A, method)
%
% 2015 - 10 edited by Maxime PINSARD
%
%  To calculate and plot the phase map

% On calcule la phase pour chaque pixel (moyenne sur 9 pixels) - voir fit_I_SHG.m

tic

%% averaging on nearest neighbors

% mean_meth = 1; % 0 for median, 1 for mean (mean is FAST)
% real_sigma=1;

if (weight && ~avg_on_ph) % % avg on contr !
    func_hdl = load_stack_plot_ISHG;
     [contr1, ~] = func_hdl.avg_nearest_ut(contr, skip_asking_coeff_median, coef_centre, coef_cross, coef_diag, mean_meth, sigma_med, sigma_med, 1); %real_sigma
     if ~isempty(contr1); contr = contr1; end % % empty if user canceled
end

%% phase

% phase=zeros(size(contr,1),size(contr,2));
% amp=zeros(size(contr,1),size(contr,2));
% err = zeros(size(contr,1),size(contr,2));

% !!
% figure; plot(squeeze(mean(mean(contr,1),2))); hold on;
% plot(squeeze(mean(mean(contr,1),2)),'rx','Markersize',18); title('Is that a cos behavior ??');
short_version = 1;
[ ~, ~, ~, ~, ~, ~ ] = analyse_cos_ISHG( contr, 0, 0, 0, 0, 0, xTitle_dflt, yTitle_dflt, screensize, fact, left_offset_fig, top_offset_fig, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, short_version, x_phase, double_contrast, contr_mode  );
title(gca, 'Is that a cos behavior ??');


% switch use_invA % 1 for inv. of matrix then multiply
%     case 1
%         A = [length(x_phase'), sum(cos(x_phase'/180*pi)), sum(sin(x_phase'/180*pi)); ...
%             sum(cos(x_phase'/180*pi)), sum(cos(x_phase'/180*pi).^2), sum(cos(x_phase'/180*pi).*sin(x_phase'/180*pi)); ...
%             sum(sin(x_phase'/180*pi)), sum(cos(x_phase'/180*pi).*sin(x_phase'/180*pi)), sum(sin(x_phase'/180*pi).^2)];
%         inv_A = A^(-1);
%     case 0 %  0 for inv. of system
%         inv_A = 0;
% end

verbose = 0;

switch method
    
    case 2 % algo 3 phases not tilted no vib
        
        method = 3;
        algo_tilted = 0;
        algo_vib = 0;
        fact_time = 0.000105; % secondes
        
        if verbose 
            htime = msgbox(sprintf('It will take around %.3g seconds/iteration (can be ~ 10 iterations)', ...
            fact_time*size(contr,1)*size(contr,2)));
        end
    case 3 % algo 3 phases not tilted, vib
        method = 3;
        algo_tilted = 0;
        algo_vib = 1;
        
        fact_time = 0.000250; % secondes
        if verbose
            htime = msgbox(sprintf('It will take around %.3g seconds/iteration (can be ~ 10 iterations)', fact_time*size(contr,1)*size(contr,2)));
        end
    case 4 % algo 3 phases % tilted, no vib
        method = 3;
        algo_tilted = 1;
        algo_vib = 0;
        
        fact_time = 0.000250; % secondes
        if verbose
            htime = msgbox(sprintf('It will take around %.3g seconds/iteration (can be ~ 10 iterations)', fact_time*size(contr,1)*size(contr,2)));
        end
    case 5 % tilted + vib.
        method = 3;
        algo_tilted = 1;
        algo_vib = 1;
        
        fact_time = 0.000270; % secondes
        if verbose
            htime = msgbox(sprintf('It will take around %.3g seconds/iteration (can be ~ 10 iterations)', fact_time*size(contr,1)*size(contr,2)));
        end
    case 0 % Pinsard : with 1D matrices
        % example 2 cornea : took 3.4 sec
       
        algo_tilted = 0;
        algo_vib = 0;
        if verbose
            htime =  msgbox(sprintf('It will take around %.3g seconds', 0.000120*size(contr,1)*size(contr,2)));
        end
    case 1 % Rivard method 2D matrices
        algo_tilted = 0;
        algo_vib = 0;
        if verbose
            htime = msgbox(sprintf('It will take around %.3g seconds', 0.000130*size(contr,1)*size(contr,2)));
        end
end

% for a 150*100 image, 6sec/it. with parallel and 3.8sec without (tilted no
% vib) !!
% 8.4 - 40 sec with parallel, 12.3sec without (tilted no
% vib) !!
% 15 - 38 sec with parallel, 44sec for 1000x200 without
% parallel_comput = 1;
% decomp_LU = 1; % 1 recommended, otherwise calculation of inverse matrix
corr = 4 ; % 0 for no corr, 1 for custom, 2 for unwrap with matlab function
%  3 for corr on pixels, at each phase shifts
% 4 for corr on phase shifts, at each pixels
% 5 for both
% corr_d = 0; % correction on x and y direction
% corr_90 = 0;

[ phase, amp, err, ~, ~, divers_vect ] = phase_algo_3frames3( contr, x_phase, use_invA, k_max, epsilon, epsilon1, ...
    algo_tilted, algo_vib, corr, h, method, parallel_comput, decomp_LU, reit_prev_phase, prev_res, ramp_phshft );

if (weight && avg_on_ph) % % avg on phase res !
    func_hdl = load_stack_plot_ISHG;
     [phase1, ~] = func_hdl.avg_nearest_ut(phase, skip_asking_coeff_median, coef_centre, coef_cross, coef_diag, mean_meth, sigma_med, sigma_med, 1); %real_sigma
     if ~isempty(phase1)
         phase = phase1;
         [amp, ~] = func_hdl.avg_nearest_ut(amp, skip_asking_coeff_median, coef_centre, coef_cross, coef_diag, mean_meth, sigma_med, sigma_med, 1); %real_sigma
         [err, ~] = func_hdl.avg_nearest_ut(err, skip_asking_coeff_median, coef_centre, coef_cross, coef_diag, mean_meth, sigma_med, sigma_med, 1); %real_sigma
     end % % empty if user canceled
end

phase = phase.*mask_val_sat; amp = amp.*mask_val_sat;  err = err.*mask_val_sat; 

mm = min(min(amp));
if mm < 0
    %     amp = amp - mm;
    im = amp(:);
    nb = round(length(im)/10);
    [N,edges_hist] =histcounts(log(abs(im)), nb);
    coef = 0.95;
    x1=find(cumsum(N)>=length(im)*(1-coef)); if ~isempty(x1); x1 = x1(1); else; x1 = 1; end
    x2=find(cumsum(N)>=length(im)*coef); if ~isempty(x2); x2 = x2(1); else; x2 = length(N); end
    if x1~= x2
        amp(abs(amp) < exp(edges_hist(x1))) = edges_hist(x1);  amp(abs(amp) > exp(edges_hist(x2))) = edges_hist(x2);
    end
end

try
    delete(htime );
    %     delete(hh);
catch
    disp('Couldn`t close bar')
end
% !!!! Does NOT take into account the edges pixels, cf C.A. code for
% this !!!!

toc

% Objects removed

% ******** Pour avoir une figure 3D de la phase **********
% figure;
% surf(x,y,phase)
% colormap([0 0 0;0.0625 0 0;0.125 0 0;0.1875 0 0;0.25 0 0;0.3125 0 0;0.375 0 0;0.4375 0 0;0.5 0 0;0.5625 0 0;0.625 0 0;0.6875 0 0;0.75 0 0;0.8125 0 0;0.875 0 0;0.9375 0 0;1 0 0;1 0.0625 0.0625;1 0.125 0.125;1 0.1875 0.1875;1 0.25 0.25;1 0.3125 0.3125;1 0.375 0.375;1 0.4375 0.4375;1 0.5 0.5;1 0.5625 0.5625;1 0.625 0.625;1 0.6875 0.6875;1 0.75 0.75;1 0.8125 0.8125;1 0.875 0.875;1 0.9375 0.9375;1 1 1;0.9375 1 0.9375;0.875 1 0.875;0.8125 1 0.8125;0.75 1 0.75;0.6875 1 0.6875;0.625 1 0.625;0.5625 1 0.5625;0.5 1 0.5;0.4375 1 0.4375;0.375 1 0.375;0.3125 1 0.3125;0.25 1 0.25;0.1875 1 0.1875;0.125 1 0.125;0.0625 1 0.0625;0 1 0;0 0.933333337306976 0;0 0.866666674613953 0;0 0.800000011920929 0;0 0.733333349227905 0;0 0.666666686534882 0;0 0.600000023841858 0;0 0.533333361148834 0;0 0.466666668653488 0;0 0.400000005960464 0;0 0.333333343267441 0;0 0.266666680574417 0;0 0.200000002980232 0;0 0.133333340287209 0;0 0.0666666701436043 0;0 0 0]);
% colorbar
% title('Phase relative (\pi\timesrad)','FontSize',30)
% xlabel('x (\mum)')
% ylabel('y (\mum)')
% axis image
% axis xy
% set(gca,'TickDir','out')
% set(gca,'CLim',[-1 1])
% ******** END Pour avoir une figure 3D de la phase **********

% Affichage graphique de la phase

h = pl_phasemap_ISHG(undocked_fig, screensize, fact, left_offset_fig, ...
    top_offset_fig, h, phase, x, y, ...
    xTitle_dflt, yTitle_dflt, Titre1, cmap_default, phi_mat_default, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz);

% in C.A.'s code
%colormap(hsv); % in Rivard's code

if obj_test==1 % Si on a enlevé des cellules, elles devraient avoir une phase fixe de 1.05
    set(gca,'CLim',[-1 1.05])
    %colormap([0 0 0;0.0625 0 0;0.125 0 0;0.1875 0 0;0.25 0 0;0.3125 0 0;0.375 0 0;0.4375 0 0;0.5 0 0;0.5625 0 0;0.625 0 0;0.6875 0 0;0.75 0 0;0.8125 0 0;0.875 0 0;0.9375 0 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0;0.4375 0 0;0.375 0 0;0.3125 0 0;0.25 0 0;0.1875 0 0;0.125 0 0;0.0625 0 0;0 0 0;0 0.0625 0;0 0.125 0;0 0.1875 0;0 0.25 0;0 0.3125 0;0 0.375 0;0 0.4375 0;0 0.5 0;0 0.5625 0;0 0.625 0;0 0.6875 0;0 0.75 0;0 0.8125 0;0 0.875 0;0 0.9375 0;0 1 0;0 0.933333337306976 0;0 0.866666674613953 0;0 0.800000011920929 0;0 0.733333349227905 0;0 0.666666686534882 0;0 0.600000023841858 0;0 0.533333361148834 0;0 0.466666668653488 0;0 0.400000005960464 0;0 0.333333343267441 0;0 0.266666680574417 0;0 0.200000002980232 0;0 0.133333340287209 0;0 0.0666666701436043 0;1 1 0]);
    % %     colormap ([1 0 0;1 0.09375 0;1 0.1875 0;1 0.2813 0;1 0.375 0;1 0.4688 0;1 0.5625 0;1 0.6563 0;1 0.75 0;1 0.8438 0;1 0.9375 0;0.9688 1 0;0.875 1 0;0.7813 1 0;0.6875 1 0;0.5938 1 0;0.5 1 0;0.4063 1 0;0.3125 1 0;0.2188 1 0;0.125 1 0;0.03125 1 0;0 1 0.0625;0 1 0.1563;0 1 0.25;0 1 0.3438;0 1 0.4375;0 1 0.5313;0 1 0.625;0 1 0.7188;0 1 0.8125;0 1 0.9063;0 1 1;0 0.9063 1;0 0.8125 1;0 0.7188 1;0 0.625 1;0 0.5313 1;0 0.4375 1;0 0.3438 1;0 0.25 1;0 0.1563 1;0 0.0625 1;0.03125 0 1;0.125 0 1;0.2188 0 1;0.3125 0 1;0.4063 0 1;0.5 0 1;0.5938 0 1;0.6863 0 1;0.7882 0 1;0.8902 0 1;0.9922 0 1;1 0 0.9059;1 0 0.8039;1 0 0.7059;1 0 0.6039;1 0 0.502;1 0 0.4;1 0 0.298;1 0 0.1961;1 0 0.09412;0 0 0]);
    % HSV with sat
    colormap (cmap_obj_test);
else
    set(gca,'CLim',[-1 1])
end

% export_fig phase.jpg

% *********** Affichage graphique Interf contrast + erreur ***********

undocked_fig = 1; % want to see these fig outside

if undocked_fig
    h_f=figure(18); set(h_f,'Color', [1 1 1], 'outerposition',...
        [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
        screensize(3)*fact screensize(4)*fact]);
    % set(gcf, 'Position', get(0,'Screensize'));
    
    if length(x) >= (4/3)^2*length(y)
        h3 = subplot(2,1,1);
        h4 = subplot(2,1,2);
    else
        
        h3 = subplot(1,2,1);
        h4 = subplot(1,2,2);
    end
    
end

% % disp interf. contrast
draw_plots_ISHG( 0, 0, amp, x, y, h3, ...
    xTitle_dflt, yTitle_dflt, Titre2, parula, 0, '', 0, 1, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );

% Affichage de l'erreur

draw_plots_ISHG( 0, 0, err, x, y, h4, ...
    xTitle_dflt, yTitle_dflt, Titre3, parula, 0, '', 0, 1, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );


end

