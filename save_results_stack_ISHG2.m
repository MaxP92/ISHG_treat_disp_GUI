function save_results_stack_ISHG2( plot_diff_menu, plot_phase, contr, in_gui, undk_p1_chck, screensize, fact, contrast, cmap, xTitle_dflt, yTitle_dflt, left_offset_fig, top_offset_fig, axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, x, y )
% save_results_stack_ISHG2
%
% 2017-10 by Maxime PINSARD
%
% WARNING : GUI-only option
%
% To save the results image of inferferometric contrast in tiff images
choice_save_meth = menu('Save method','matlab .fig', '.tif color');

folder_name2 = uigetdir('Select the folder where to put the images to save !');
cd(folder_name2);


title_fig_all = get(plot_diff_menu, 'String');

set(undk_p1_chck, 'Value', 0);


for K = 1:size(contr,3)
    % Sauvegarde des images rouges/vertes du contraste interférométrique
    
    set(plot_diff_menu, 'Value', 2+K);
    
    title_fig_current = title_fig_all{2+K};
    
    
    switch choice_save_meth
        
        case 1 % .fig
            
            switch in_gui
                
                case 1
                    
                    Fig2 = figure;
                    
                    
                    copyobj(plot_phase, Fig2);
                    
                case 0
                    
                    Fig2 = figure('outerposition',...
                        [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
                        screensize(3)*fact screensize(4)*fact]);
                    h2 = axes;
                    
                    % [left bottom width height]
                    
                    Maximum = max(max(abs(contr(:,:,K))));
                    
                    draw_plots_ISHG( 0, 0, contr(:,:,K), x, y, h2, ...
                        xTitle_dflt, yTitle_dflt, title_fig_current, cmap, [-1*contrast*Maximum contrast*Maximum], '', 0, 1, ...
                        axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );
            end
            savefig(Fig2, sprintf('%s.fig', title_fig_current));
            
            close(Fig2);
        case 2 % .mat
            
            contr_corr = contr(:,:,K) + abs(min(min(contr(:,:,K))));
            
            imwrite(uint16(contr_corr), sprintf('%s.tif', title_fig_current));
    end
    
    %         if hfig1 == 0
    %             hfig1 = figure('outerposition',...
    %                 [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
    %                 screensize(3)*fact screensize(4)*fact]);
    %             h2 = axes;
    %         else
    %             figure(hfig1);
    %         end
    %         % [left bottom width height]
    %
    %         Maximum = max(max(abs(contr(:,:,K))));
    %
    %         draw_plots_ISHG( 0, 0, contr(:,:,K), x, y, h2, ...
    %         xTitle_dflt, yTitle_dflt, '', cmap, [-1*contrast*Maximum contrast*Maximum], '', 0, 1, ...
    %         axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );
    %
    %         saveas(hfig1,['ColorStack' num2str(K)],'tif');
    
end

end

