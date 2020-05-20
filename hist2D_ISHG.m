function [h5, hist2_xdata, hist2_ydata] = hist2D_ISHG( fact, left_offset_fig, ...
    top_offset_fig, phasee, int_x, phi_mat_default, Counts, Titre4, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, screensize, offset_pi2, h5 )
%hist2D_ISHG : plot the repartition of the phase among the pixels of the
%image
%  
% edited 2016.3 Maxime Pinsard
% 
% [h5, hist2_xdata, hist2_ydata] = hist2D_ISHG( fact, left_offset_fig, ...
%     top_offset_fig, ph_hist, int_x, phi_mat_default, Counts, Titre4, ...
%     axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, screensize, offset_pi2, h5 )

%% init 

ph_hist = reshape(phasee,1,size(phasee,1)*size(phasee,2));

ph_hist(ph_hist==1.05) = [];

%% hist 2D

hist2_xdata = -1+int_x/2:int_x:1-int_x/2;
% hist2_xdata =length(-1:int_x:1);

undocked_fig_2d = 1; % we want this plot outside

try
    get(h5, 'BeingDeleted');  % if off, h5 has not been deleted
catch % h5 has been deleted
    
    if undocked_fig_2d
        figure(50); % % 50 is for 2D hist
        set(50, 'outerposition',...
            [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
            screensize(3)*fact screensize(4)*fact]);
        h5 = axes; % create axes in current figure
    end
    
end

draw_plots_ISHG( 2, 0, ph_hist, hist2_xdata, 0, h5,...
    phi_mat_default, Counts, Titre4, 0, 0, '', 0, 0, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );

colorbar('peer', h5, 'off');

if offset_pi2
    mbversion = version; % query the MatLab version
    axes(h5);
    if (str2double(mbversion(end-3)) > 0 && str2double(mbversion(end-2)) >= 4) % Matlab releases posterior to R2014a
        set(h5,'XTickLabel',{'\pi /2' '-\pi' '-\pi /2' '0' '\pi /2'}, 'TickLabelInterpreter', 'tex');
    else % because of Matalb f***ing update !
        set(h5,'XTickLabel','p /2 | -p | -p /2 | 0 | p /2', 'fontname', 'symbol');
    end
end

hist2_ydata = hist(ph_hist, hist2_xdata); % does not do the plot

end

