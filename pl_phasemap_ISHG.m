function [h, hfig] = pl_phasemap_ISHG(undocked_fig, screensize, fact, left_offset_fig, ...
top_offset_fig, h, phase, x, y, ...
    xTitle_dflt, yTitle_dflt, Titre1, cmap_blkredblkgrn, phi_mat_default, ...
axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz)
% pl_phasemap_ISHG 
% 
% edited 2013.3.16 Maxime PINSARD
% 
%   h = pl_phasemap_ISHG(undocked_fig, screensize, fact, left_offset_fig, ...
% top_offset_fig, h, phase, x, y, ...
%     xTitle_dflt, yTitle_dflt, Titre1, cmap_blkredblkgrn, phi_mat_default, ...
% axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz)

if undocked_fig
    try
        get(h,'Children'); % error if deleted
        hh=get(h, 'parent');
        if strcmp(get(hh, 'name'),'I_SHG_GUI'); hfig = 1; else; hfig = 0; end
    catch % if it has been deleted
        hfig = 1;
    end
    if hfig == 1
        hfig = figure('Color', [1 1 1], 'outerposition',...
            [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
            screensize(3)*fact screensize(4)*fact]);
        % set(gcf, 'Position', get(0,'Screensize'));
        h = axes; % create axes in current figure
    end
else
    hfig = 0;
end

if h == 0
    h = axes; % create axes in current figure
end

% % cmap = [0 0 0;0.0625 0 0;0.125 0 0;0.1875 0 0;0.25 0 0;0.3125 0 0;0.375 0 0;0.4375 0 0;0.5 0 0;0.5625 0 0;0.625 0 0;0.6875 0 0;0.75 0 0;0.8125 0 0;0.875 0 0;0.9375 0 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0;0.5 0 0;0.4375 0 0;0.375 0 0;0.3125 0 0;0.25 0 0;0.1875 0 0;0.125 0 0;0.0625 0 0;0 0 0;0 0.0625 0;0 0.125 0;0 0.1875 0;0 0.25 0;0 0.3125 0;0 0.375 0;0 0.4375 0;0 0.5 0;0 0.5625 0;0 0.625 0;0 0.6875 0;0 0.75 0;0 0.8125 0;0 0.875 0;0 0.9375 0;0 1 0;0 0.933333337306976 0;0 0.866666674613953 0;0 0.800000011920929 0;0 0.733333349227905 0;0 0.666666686534882 0;0 0.600000023841858 0;0 0.533333361148834 0;0 0.466666668653488 0;0 0.400000005960464 0;0 0.333333343267441 0;0 0.266666680574417 0;0 0.200000002980232 0;0 0.133333340287209 0;0 0.0666666701436043 0;0 0 0];

draw_plots_ISHG( 0, 1, phase, x, y, h, ...
    xTitle_dflt, yTitle_dflt, Titre1, cmap_blkredblkgrn, 0, phi_mat_default, 0, 1, ...
axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );

caxis(h, [-1, 1]);

end

