function [h, list_contr_name] = diff_phase_stack_ISHG( contr, x, y, contrast, xTitle_dflt, yTitle_dflt, screensize, fact, ...
    left_offset_fig, top_offset_fig, undocked_fig, h, axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, cmap_redgreen, ind_begin, ind_end, x_phase )
% h = diff_phase_stack_ISHG( contr, x, y, start_phase, diff_phase,
% contrast, screensize, fact, left_offset_fig, top_offset_fig, undocked_fig, h, axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, cmap_redgreen, ind_begin, ind_end, x_phase )
%
% 2015-10-26 edited by Maxime PINSARD
%
% draw the plot for the different differences of angles
%
% old version in V13

% % %% Réorganisation des images du contraste interférométrique dans un ordre croissant
% %
% % % On replace les images dans un ordre croissant de phase
% % tmp = contr;
% % clear contr
% % contr = zeros(size(tmp));
% %
% % % On a par exemple 0-180 180-360 15-195 195-375 ...
% % % On veut 0-180 15-195 ... ... 180-360 195-375 ...
% % contr(:,:,1:size(tmp,3)/2) = tmp(:,:,1:2:size(tmp,3)-1);
% % contr(:,:,size(tmp,3)/2+1:size(contr,3)) = tmp(:,:,2:2:size(tmp,3));

%% Affichage du stack de différences

% cmap = parula; % C.A.'s code

Maximum = max(max(max((contr(:,:,:)))));
Minimum = min(min(min((contr(:,:,:)))));

list_contr_name = cell(1, ind_end - ind_begin + 1);

for ii = ind_begin:ind_end % 
    
%     try
        
        %         if nonequal_spacing
        %             pos = ii;
        deg1 = num2str(x_phase(ii));
        deg2 = num2str(x_phase(ii)+180);
        %         else % equal
        %             if diff_phase <= 90
        %
        %                 pos(ii) = 1 + 90*(ii-1)/diff_phase;
        %                 deg1 = num2str(start_phase + 90*(ii-1));
        %                 deg2 = num2str(start_phase + 90*(ii-1) + 180);
        %             end
        %         end
        %h_stack = subplot(2,2,1);

        if undocked_fig
            %             try
            %                 get(h,'Children'); % error if deleted
            %             catch % if it has been deleted
            figure('outerposition',...
                [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
                screensize(3)*fact screensize(4)*fact]);
            % [left bottom width height]
            h = axes; % create axes in current figure
            %             end
        end
        
        list_contr_name{ii} = [deg2 '° - ' deg1 '°'];
        
        draw_plots_ISHG( 0, 1, contr(:,:,ii), x, y, h, xTitle_dflt, yTitle_dflt, ...
            list_contr_name{ii}, cmap_redgreen, [contrast*Minimum contrast*Maximum], 'Interf. contrast', 0, 1, ...
            axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz ); % '\boldmath$\phi_{mat}$'
        
        
%     catch MEE
%         warning('Could not plot one diff phase');
%         fprintf(2, 'Plot # %d ; %s \n %s \n', ii, MEE.identifier, MEE.message);
%         
%     end
end

% step = 2;
% contr = contr(1:step:end);
% fprintf(2, 'warning : only few phases considered to calculate (diff_phase line 54)\n'); % !!!!

end

