function [ data_tot, ylabel_hist3, title_hist3, nbinsx, h ] = hist_3D_ISHG( phasee, img, int_x, int_y, title_hist3, ylabel_hist3, Counts, phi_mat_default,...
axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig, offset_pi2, range_sat )
% [ data_tot, ylabel_hist3, title_hist3, nbinsx, h ] = hist_3D_ISHG( ph_hist, second_hist, int_x, int_y, title_hist3, ylabel_hist3, Counts, phi_mat_default,...
% axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig, offset_pi2 )
% 
% 2015-10 edited by Maxime PINSARD
% 
% : plot the repartition of the phase and the interf. contrast among the pixels of the
% image.
%  CAUTION : the count is lower than in 2D hist because it is also spread among the interf. contrast (in addition to the phase) 

% ph_hist(ph_hist==1.05) = [];

%% init

ph_hist = reshape(phasee,1,size(phasee,1)*size(phasee,2));

ph_hist(ph_hist==1.05) = [];

if range_sat == 0
    second_hist = reshape(img,1,size(img,1)*size(img,2));
    second_hist(second_hist==0) = [];
else
    amp_test=img;
    amp_test_hist = reshape(amp_test,1,size(amp_test,1)*size(amp_test,2));
    amp_test_hist(amp_test==0) = [];
    second_hist = amp_test_hist;
    %5 hist 3D
end
second_hist = second_hist(~isnan(second_hist));
ph_hist = ph_hist(~isnan(second_hist)); % concordance of sizes
if length(second_hist)< length(ph_hist)
    second_hist(end:end+length(ph_hist)-length(second_hist)) = 0;
end

data_tot(:,1) = ph_hist;
data_tot(:,2) = second_hist;
%data_tot(:,2) = shg_hist./2^15;

nbinsx = length(-1:int_x:1);
% nbinsy = ceil((max(data_tot(:,2))-min(data_tot(:,2)))/int_y);
nbinsy = int_y;

ylim_hist3 = [min(second_hist), (2-(max(second_hist)>min(second_hist)))*max(second_hist)]; % for the case where max(second_hist) == min(second_hist)

if undocked_fig
    try
        get(h,'Children'); % error if deleted
    catch % if it has been deleted
        figure('outerposition',...
            [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
            screensize(3)*fact screensize(4)*fact]);
        h = axes; % create axes in current figure
    end
end

if h == 0
    h = axes; % create axes in current figure
end

draw_plots_ISHG( 3, 0, [data_tot; [1.0,0;-1.0,0]], [nbinsx nbinsy], 0, h,...
    phi_mat_default, ylabel_hist3, title_hist3, 0, 0, Counts, ylim_hist3, 0, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );% just for a better display

zlabel(h, 'Number of pixels', 'FontSize', yaxis_sz);

hh = colorbar(h) ;
% title(hh, '');     
ylabel(hh, Counts,'Interpreter','latex','fontsize',clrbr_tl_sz);
% set(hL,'Rotation',90);
% poss = get (hh, 'Position');
% % set(hL,'Position', [(55 + (axes_font_size/12-1)*24) 90 0])
% set(hL,'Position', [((poss(1) + poss(3))*60 + 10) poss(2)*200 + 40 0])

% set(gcf, 'Position', get(0,'Screensize'));
% set(gcf,'renderer','opengl');

% put 'counts' as colrbar title

if offset_pi2
    mbversion = version; % query the MatLab version
    axes(h);
    if (str2double(mbversion(end-3)) > 0 && str2double(mbversion(end-2)) >= 4) % Matlab releases posterior to R2014a
        set(h,'XTickLabel',{'\pi /2' '-\pi' '-\pi /2' '0' '\pi /2'}, 'TickLabelInterpreter', 'tex');
    else % because of Matalb f***ing update !
        set(h,'XTickLabel','p /2 | -p | -p /2 | 0 | p /2', 'fontname', 'symbol');
    end
end

end

