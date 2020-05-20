function draw_plots_ISHG( hist_plot, phase_contrast, contr, x, y, hh, xTitle, yTitle, figtitle, cmap, caxis_img, clrbr_title, ylim_vect, set_axis_img, axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz )
% draw_plots_ISHG( hist_plot, phase_contrast, contr, x, y, hh, xTitle, yTitle, figtitle, cmap, caxis_img, clrbr_title, ylim_vect, set_axis_img, axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz )
%
% hist_plot : states if it's an histogram or image
% phase_contrast : states if it's an image in phase contrast or not
% contr : 2D matrix to plot
% x, y : limit of the plot (y optional in some case)
% hh : handle of the AXES, not figure (may be optional)
% xTitle, yTitle, figtitle : string for resp. xlabel, ylabel and title
% cmap : colormap
% caxis_img : limit of the colormap (optional in some case)
% clrbr_title : title of the colorbar
% ylim_vect : limit fir the y axis (optional in some case)
% set_axis_img : set axis to 'image'
% 
% 2015-10-27 edited by Maxime PINSARD
% 
%  To plot the figures of the ISHG program : hist or image figure


if (hist_plot == 2 || hist_plot == 3)  % histogram or not
    
    if hist_plot == 2   
        %% if it is an histogram 2D
        
        histo =  histogram(hh, contr, x); % does only the plot
        
        histo.FaceColor = [0 0 0.5];
        histo.EdgeColor = [0 0 0.5];
        % set(gcf, 'Position', get(0,'Screensize'));
        %         hist_xdata = -1+int_x/2:int_x:1-int_x/2;
        %         hist_ydata = hist(ph_hist, hist_xdata);
%         h = findobj(hh ,'Type','patch');
%         set(h,'FaceColor','b')
        
    elseif hist_plot == 3
        %% if it is an histogram 3D
        
        hist3(hh, contr, x); % contr is Nx2 vect
        % set(gcf, 'Position', get(0,'Screensize'));
        % set(gcf,'renderer','opengl');
        set(gcf, 'renderer', 'zbuffer')
        set(get(hh, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto');
    
        view(hh, [0,90]); % equivalent to 2D
        xlim(hh, [-1, 1]);
        ylim(hh, ylim_vect);
        
    end
    
    set(hh,'xtick',[-1 -0.5 0 0.5 1])
    mbversion = version; % query the MatLab version
    if (str2double(mbversion(1:3)) >= 8.3) % Matlab releases posterior to R2014a
        set(hh,'XTickLabel',{'-\pi' '-\pi /2' '0' '\pi /2' '\pi'}, 'TickLabelInterpreter', 'tex');
    else % because of Matalb f***ing update !
        set(hh,'XTickLabel','-p | -p /2 | 0 | p /2 | p', 'fontname', 'symbol');
    end
    xlabel(hh, '\boldmath$\phi_{mat}$','Interpreter','Latex','FontSize',xaxis_sz)
    colormap(hh, parula); % default colormap
    
    if hh ~= 0
        axes(hh);
    else
        axes( 'Color',[1 1 1]);
    end

else
    %% plot of the image with phase contrast
    
    axes(hh) ;
    hold off;
    imagesc(x, y, contr);

    set(hh,'TickDir','out');
    if (~isnan(sum(caxis_img)) && sum(caxis_img) ~= 0 ) % % contr
        caxis(hh, caxis_img);
    else % %  phase
        if phase_contrast
            caxis(hh, [-1, 1]);
        end
    end
    colormap(hh, cmap);
    
end

% % set(hh,'xlim',[min( y ) max( y )]);
% % set(hh,'ylim',[min( x ) max( x )]);

h = colorbar(hh) ;
title(h, clrbr_title,'Interpreter','latex','fontsize',clrbr_tl_sz);

if phase_contrast % in the case of the plot with a phase contrast
   if (~isnan(sum(caxis_img)) && sum(caxis_img) ~= 0 ) % contr
%         cc = contr(1:round(y(end)), 1:round(x(end)));
%         set(h,'ylim',[min(min( cc )), max(max( cc ))]);
%         set(h,'ytick',linspace(min(min( cc )), max(max( cc )), 5));
        set(h,'ylim',[caxis_img(1), caxis_img(2)]);
        set(h,'ytick',linspace(caxis_img(1), caxis_img(2), 5));

    else % phase
        set(h,'ylim',[-1, 1]);
        set(h,'ytick',[-1, -0.5, 0, 0.5, 1]);
            mbversion = version; % query the MatLab version
        if (str2double(mbversion(1:3)) >= 8.3) % Matlab releases posterior to R2014a
            set(h,'YTickLabel',{'-\pi' '-\pi /2' '0' '\pi /2' '\pi'}, 'TickLabelInterpreter', 'tex');
        else % because of Matalb f***ing update !
            set(h,'YTickLabel','-p | -p /2 | 0 | p /2 | p', 'fontname', 'symbol');
        end 
    end
    
end

if set_axis_img
    axis(hh, 'image');
end

set(hh,'FontSize', axes_font_size)
xlabel(hh, xTitle, 'Interpreter', 'Latex', 'FontSize', xaxis_sz);
ylabel(hh, yTitle, 'Interpreter', 'Latex', 'FontSize', yaxis_sz);
title(hh, figtitle, 'FontSize', title_sz);

end

