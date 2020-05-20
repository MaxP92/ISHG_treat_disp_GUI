function array_frames = erase_img_avg_stck(img_3D, axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, hh_plot, nb_avg, nb_ps)
% erase_img_avg_stck : to erase some frame before averaging (if certain are
% shifted or ugly)
%
% nb_avg is handles.avg_img_edt 
%
% edited 2018.7 by Maxime PINSARD
% 

d = dialog('Position',[300 300 150 100],'Name', 'AH!'); % left bottom width height

% % ph = max_amp;
% % pb = 0;
% % offset_pi2 = 0;

% draw_plots_ISHG( 0, 0, img_3D(:,:,1), 1:size(img_3D(:,:,1), 2), 1:size(img_3D(:,:,1), 1), hh_plot, ...
%     'X (pixels)', 'Y (pixels)', 'First image of the stack', gray, 0, 'SHG int. (a. u.)', 0, 1, ...
%     axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );

% % txt_min = uicontrol('Parent',d,...
% %     'Style','text',...
% %     'Position',[50 170 25 25],...
% %     'String','Min');
% % 
% % txt_max = uicontrol('Parent',d,...
% %     'Style','text',...
% %     'Position',[50 145 25 25],...
% %     'String','Max');
ind_disp =1;
nb_avg_current = nb_avg;
array_frames = reshape(1:size(img_3D, 3), [ nb_avg, nb_ps]);

txt_1 = uicontrol('Parent',d,...
    'Style','text',...
    'Position',[10 50 100 50],...
    'String',sprintf('nb frame in this phase (# %d), fr # %d :', 1, 1)); % [10 150 100 50]

txt_nb_avg = uicontrol('Parent',d,...
    'Style','text',...
    'Position',[110 50 50 50],...
    'String',num2str(nb_avg));

btn_next = uicontrol('Parent',d,...
    'Position',[60 40 30 20],...
    'String','->',...
    'Callback', @nxt_cllbck);

btn_prec = uicontrol('Parent',d,...
    'Position',[30 40 30 20],...
    'String','<-',...
    'Callback', @prec_cllbck);

erase_btn = uicontrol('Parent',d,...
    'Position',[10 10 40 20],...
    'String','Del.',...
    'Callback', @erase_frm_cllbck);

end_btn = uicontrol('Parent',d,...
    'Position',[110 20 40 30],...
    'String','Close',...
    'Callback','delete(gcf)'); %#ok<*NASGU>
% hor, vert, w, h
% Wait for d to close before running to completion
uiwait(d);


    function nb_avg_update()
        num_ph = ceil(ind_disp/nb_avg);
        set(txt_1, 'String', sprintf('nb frame in this phase (# %d), fr. # %d :', num_ph, ind_disp));
        nb_avg_current = sum(~isnan(array_frames(:,num_ph)));
        set(txt_nb_avg, 'String', num2str(nb_avg_current));
    end

    function callbackdata = nxt_cllbck(~, callbackdata)
% %         disp(ind_disp)

        ind00 = min(max(array_frames(~isnan(array_frames))), ind_disp+1);
        arr = array_frames(ind00:end);
        arr = arr(~isnan(arr));
        if ~isempty(arr) % not empty
            ind_disp = arr(1);
        else
            ind00 = max(min(array_frames(~isnan(array_frames))), ind00);
            arr = array_frames(1:ind00);
            arr = arr(~isnan(arr));
            if ~isempty(arr) % not empty
                ind_disp = arr(1);
            else
                return
            end
        end
% %         disp(ind_disp)
        
        draw_plots_ISHG( 0, 0, img_3D(:,:,ind_disp), 1:size(img_3D(:,:,1), 2), 1:size(img_3D(:,:,1), 1), hh_plot, ...
    'X (pixels)', 'Y (pixels)', 'X image of the stack', 'gray', 0, 'SHG int. (a. u.)', 0, 1, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );
%         imagesc(hh_plot, img_3D(:,:,ind_disp))
        nb_avg_update();
    end

    function callbackdata = prec_cllbck(~, callbackdata)
        
% %         disp(ind_disp)

        ind00 = max(min(array_frames(~isnan(array_frames))), ind_disp-1);
        arr = array_frames(1:ind00);
        arr = arr(~isnan(arr));
        if ~isempty(arr) % not empty
            ind_disp = arr(end);
        else
            disp(min(array_frames(~isnan(array_frames))) )

            ind00 = min(max(array_frames(~isnan(array_frames))), ind00);
            arr = array_frames(ind00:end);
            arr = arr(~isnan(arr));
            if ~isempty(arr) % not empty
                ind_disp = arr(end);
            else
                return
            end
        end

        draw_plots_ISHG( 0, 0, img_3D(:,:,ind_disp), 1:size(img_3D(:,:,1), 2), 1:size(img_3D(:,:,1), 1), hh_plot, ...
    'X (pixels)', 'Y (pixels)', 'X image of the stack', 'gray', 0, 'SHG int. (a. u.)', 0, 1, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );
%         imagesc(hh_plot, img_3D(:,:,ind_disp))
        nb_avg_update();
    end

    function callbackdata = erase_frm_cllbck(~,callbackdata)
        if nb_avg_current <= 1
            errordlg('You cannot do that: there will be no frame left for this phase-shift !', 'Not enough frame');
        else
            array_frames(ind_disp) = nan;
            cla(hh_plot);
%             ind_disp = ind_disp+1;
           nb_avg_update();
        end
    end

end

