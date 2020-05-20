function [ phase_test_hist, phase_test, hist2_xdata, hist2_ydata, hhist1 ] = hist_spec_region_ISHG( x, y, phase, amp, img_shg, shg, int_x, int_y, xTitle_dflt, yTitle_dflt, Titre1, phi_mat_default, Counts, Titre4, Legendehisto, ylabel_hist3, title_hist3, nbinsx, pb, ph,...
axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig1, undocked_fig2, hfun2, hact2, screensize, fact, left_offset_fig, top_offset_fig, cmap_brbgsat, offset_pi2, offset  )
% [  phase_test_hist, phase_test, hist2_xdata, hist2_ydata, hhist1 ] = hist_spec_region_ISHG( x, y, phase, amp, img_shg, shg, int_x, int_y, xTitle_dflt, yTitle_dflt, Titre1, phi_mat_default, Counts, Titre4, Legendehisto, ylabel_hist3, title_hist3, nbinsx, pb, ph,...
%  axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig1, undocked_fig2, hfun2, hact2, screensize, fact, left_offset_fig, top_offset_fig, cmap_brbgsat, offset_pi2)
% 
% 2015 - 10 edited by Maxime PINSARD
% 
% To calculate and plot the histogram 2D and 3D for a specific range of
% intensity

    % nbinsy2 = ceil((ph-pb)/int_y);
    
    phase_test = phase;
    if shg==1 % SHG mode chosen
        shg_test=img_shg;
        AA = find(img_shg(1:size(phase_test,1), 1:size(phase_test,2))>ph | img_shg(1:size(phase_test,1), 1:size(phase_test,2))<pb);
        shg_test(AA) = max(max(img_shg))+1;
    else % interferometric contrast chosen
        amp_test=amp;
        AA = find(amp(1:size(phase_test,1), 1:size(phase_test,2))>ph | amp(1:size(phase_test,1), 1:size(phase_test,2))<pb);
        amp_test(AA) = 0;
    end
    
    % phase_test = phase(1:150,400:500);
    % amp_test = amp(1:150,400:500);
    
    % old code with for loops
    % %         for i = 1:size(phase_test,1)
    % %             for j = 1:size(phase_test,2)
    % %                 if shg==1
    % %                     if (img_shg(i,j)>ph || img_shg(i,j)<pb)
    % %
    % %                         phase_test(i,j) = 1.05;
    % %                         shg_test(i,j) = max(max(img_shg))+1;
    % %                     end
    % %                 else
    % %                     if (amp(i,j)>ph || amp(i,j)<pb)
    % %                         phase_test(i,j) = 1.05;
    % %                         amp_test(i,j) = 0;
    % %                     end
    % %                 end
    % %             end
    % %         end
    
    phase_test(AA) = 1.05;
    
    if shg==1
        shg_test_hist = reshape(shg_test,1,size(shg_test,1)*size(shg_test,2));
        shg_test_hist(shg_test==max(max(img_shg))+1) = [];
    else
        amp_test_hist = reshape(amp_test,1,size(amp_test,1)*size(amp_test,2));
        amp_test_hist(amp_test==0) = [];
    end
    phase_test_hist = reshape(phase_test,1,size(phase_test,1)*size(phase_test,2));
    phase_test_hist(phase_test_hist==1.05) = [];
    
    %% Offset of pi/2 
    
    if offset_pi2
       
         phase_test_hist = phase_test_hist + offset;
        
        phase_test_hist(phase_test_hist > 1) = phase_test_hist(phase_test_hist > 1) - 2;
    end
    
    %% 2D hist
    
    undocked_fig_2d = 1; % fig outside
    
    if undocked_fig_2d
        figure('Color', [1 1 1],'outerposition',...
            [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
            screensize(3)*fact screensize(4)*fact]);
        hhist1 = axes; % create axes in current figure
    end
    
    hist2_xdata = -1+int_x/2:int_x:1-int_x/2;
    
    draw_plots_ISHG( 2, 0, phase_test_hist, hist2_xdata, 0, hhist1, ...
        phi_mat_default, Counts, Titre4, 0, 0, '', 0, 0, ...
        16, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );
    
    if offset_pi2
        mbversion = version; % query the MatLab version
        axes(hhist1);
        if (str2double(mbversion(end-3)) > 0 && str2double(mbversion(end-2)) >= 4) % Matlab releases posterior to R2014a
            set(hhist1,'XTickLabel',{'\pi /2' '-\pi' '-\pi /2' '0' '\pi /2'}, 'TickLabelInterpreter', 'tex');
        else % because of Matalb f***ing update !
            set(hhist1,'XTickLabel','p /2 | -p | -p /2 | 0 | p /2', 'fontname', 'symbol');
        end
    end
     
    hist2_ydata = hist(phase_test_hist, hist2_xdata); % does not do the plot
    
%     
    %% Phase map uploaded
    
    if undocked_fig1
        figure('Color', [1 1 1],'outerposition',...
            [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
            screensize(3)*fact screensize(4)*fact]);
        hfun2 = axes; % create axes in current figure
    end
    
    if hfun2 == 0
        hfun2 = axes; % create axes in current figure
    end
    
    % set(gcf, 'Position', get(0,'Screensize'));
    
    
    % in C.A.'s code
    
    draw_plots_ISHG( 0, 1, phase_test, x, y, hfun2, ...
        xTitle_dflt, yTitle_dflt, Titre1, cmap_brbgsat, 0, phi_mat_default, 0, 1, ...
        axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );
    
    set(gca,'CLim',[-1 1.05]);
    
    %colormap([0 0 0;0.0625 0 0;0.125 0 0;0.1875 0 0;0.25 0 0;0.3125 0 0;0.375 0 0;0.4375 0 0;0.5 0 0;0.5625 0 0;0.625 0 0;0.6875 0 0;0.75 0 0;0.8125 0 0;0.875 0 0;0.9375 0 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0;0.4375 0 0;0.375 0 0;0.3125 0 0;0.25 0 0;0.1875 0 0;0.125 0 0;0.0625 0 0;0 0 0;0 0.0625 0;0 0.125 0;0 0.1875 0;0 0.25 0;0 0.3125 0;0 0.375 0;0 0.4375 0;0 0.5 0;0 0.5625 0;0 0.625 0;0 0.6875 0;0 0.75 0;0 0.8125 0;0 0.875 0;0 0.9375 0;0 1 0;0 0.933333337306976 0;0 0.866666674613953 0;0 0.800000011920929 0;0 0.733333349227905 0;0 0.666666686534882 0;0 0.600000023841858 0;0 0.533333361148834 0;0 0.466666668653488 0;0 0.400000005960464 0;0 0.333333343267441 0;0 0.266666680574417 0;0 0.200000002980232 0;0 0.133333340287209 0;0 0.0666666701436043 0;0 0 0]);
    %colormap([0 0 0;0.0625 0 0;0.125 0 0;0.1875 0 0;0.25 0 0;0.3125 0 0;0.375 0 0;0.4375 0 0;0.5 0 0;0.5625 0 0;0.625 0 0;0.6875 0 0;0.75 0 0;0.8125 0 0;0.875 0 0;0.9375 0 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0;0.4375 0 0;0.375 0 0;0.3125 0 0;0.25 0 0;0.1875 0 0;0.125 0 0;0.0625 0 0;0 0 0;0 0.0625 0;0 0.125 0;0 0.1875 0;0 0.25 0;0 0.3125 0;0 0.375 0;0 0.4375 0;0 0.5 0;0 0.5625 0;0 0.625 0;0 0.6875 0;0 0.75 0;0 0.8125 0;0 0.875 0;0 0.9375 0;0 1 0;0 0.933333337306976 0;0 0.866666674613953 0;0 0.800000011920929 0;0 0.733333349227905 0;0 0.666666686534882 0;0 0.600000023841858 0;0 0.533333361148834 0;0 0.466666668653488 0;0 0.400000005960464 0;0 0.333333343267441 0;0 0.266666680574417 0;0 0.200000002980232 0;0 0.133333340287209 0;0 0.0666666701436043 0;1 1 0]);
    %colormap([1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0;0.4375 0 0;0.375 0 0;0.3125 0 0;0.25 0 0;0.1875 0 0;0.125 0 0;0.0625 0 0;0 0 0;0 0.0625 0;0 0.125 0;0 0.1875 0;0 0.25 0;0 0.3125 0;0 0.375 0;0 0.4375 0;0 0.5 0;0 0.5625 0;0 0.625 0;0 0.6875 0;0 0.75 0;0 0.8125 0;0 0.875 0;0 0.9375 0;0 1 0;0 0.9375 0;0 0.875 0;0 0.8125 0;0 0.75 0;0 0.6875 0;0 0.625 0;0 0.5625 0;0 0.5 0;0 0.4375 0;0 0.375 0;0 0.3125 0;0 0.25 0;0 0.1875 0;0 0.125 0;0 0.0625 0;0 0 0;0.0714285746216774 0 0;0.142857149243355 0 0;0.214285716414452 0 0;0.28571429848671 0 0;0.357142865657806 0 0;0.428571432828903 0 0;0.5 0 0;0.571428596973419 0 0;0.642857134342194 0 0;0.714285731315613 0 0;0.785714268684387 0 0;0.857142865657806 0 0;0.928571403026581 0 0;1 0 0;1 1 0]);
    %colormap([0 0 0;0.06667 0 0;0.1333 0 0;0.2 0 0;0.2667 0 0;0.3333 0 0;0.4 0 0;0.4667 0 0;0.5333 0 0;0.6 0 0;0.6667 0 0;0.7333 0 0;0.8 0 0;0.8667 0 0;0.9333 0 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0;0.4375 0 0;0.375 0 0;0.3125 0 0;0.25 0 0;0.1875 0 0;0.125 0 0;0.0625 0 0;0 0 0;0 0.06667 0;0 0.1333 0;0 0.2 0;0 0.2667 0;0 0.3333 0;0 0.4 0;0 0.4667 0;0 0.5333 0;0 0.6 0;0 0.6667 0;0 0.7333 0;0 0.8 0;0 0.8667 0;0 0.9333 0;0 1 0;0 0.9375 0;0 0.875 0;0 0.8125 0;0 0.75 0;0 0.6875 0;0 0.625 0;0 0.5625 0;0 0.5 0;0 0.4375 0;0 0.375 0;0 0.3125 0;0 0.25 0;0 0.1875 0;0 0.125 0;0 0.0625 0;0 0 0;0 0 0]);
    %colormap ([1 0 0;1 0.09375 0;1 0.1875 0;1 0.2813 0;1 0.375 0;1 0.4688 0;1 0.5625 0;1 0.6563 0;1 0.75 0;1 0.8438 0;1 0.9375 0;0.9688 1 0;0.875 1 0;0.7813 1 0;0.6875 1 0;0.5938 1 0;0.5 1 0;0.4063 1 0;0.3125 1 0;0.2188 1 0;0.125 1 0;0.03125 1 0;0 1 0.0625;0 1 0.1563;0 1 0.25;0 1 0.3438;0 1 0.4375;0 1 0.5313;0 1 0.625;0 1 0.7188;0 1 0.8125;0 1 0.9063;0 1 1;0 0.9063 1;0 0.8125 1;0 0.7188 1;0 0.625 1;0 0.5313 1;0 0.4375 1;0 0.3438 1;0 0.25 1;0 0.1563 1;0 0.0625 1;0.03125 0 1;0.125 0 1;0.2188 0 1;0.3125 0 1;0.4063 0 1;0.5 0 1;0.5938 0 1;0.6863 0 1;0.7882 0 1;0.8902 0 1;0.9922 0 1;1 0 0.9059;1 0 0.8039;1 0 0.7059;1 0 0.6039;1 0 0.502;1 0 0.4;1 0 0.298;1 0 0.1961;1 0 0.09412;0 0 0]);
    % in Rivard's code
    
    data_tot2(:,1) = phase_test_hist;
    if shg==1
        data_tot2(:,2) = shg_test_hist;
        ylim_hist3 = [min(shg_test_hist) max(shg_test_hist)];
    else
        data_tot2(:,2) = amp_test_hist;
        ylim_hist3 = [min(amp_test_hist) max(amp_test_hist)];
    end
    
    %     nbinsy2 = ceil((max(data_tot2(:,2))-min(data_tot2(:,2)))/int_y);
    nbinsy2 = int_y;

    %% Hist 3D
    
    if undocked_fig2
         figure('Color', [1 1 1],'outerposition',...
        [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
            screensize(3)*fact screensize(4)*fact]);
        hact2 = axes; % create axes in current figure
    end
    
    if hact2 == 0
        hact2 = axes; % create axes in current figure
    end
    
    draw_plots_ISHG( 3, 0, data_tot2, [nbinsx nbinsy2], 0, hact2,...
        phi_mat_default, ylabel_hist3, title_hist3, 0, 0, Counts, ylim_hist3, 0, ...
        axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );
    
    if offset_pi2
        mbversion = version; % query the MatLab version
        axes(hact2);
        if (str2double(mbversion(end-3)) > 0 && str2double(mbversion(end-2)) >= 4) % Matlab releases posterior to R2014a
            set(hact2,'XTickLabel',{'\pi /2' '-\pi' '-\pi /2' '0' '\pi /2'}, 'TickLabelInterpreter', 'tex');
        else % because of Matalb f***ing update !
            set(hact2,'XTickLabel','p /2 | -p | -p /2 | 0 | p /2', 'fontname', 'symbol');
        end
    end
    
end

