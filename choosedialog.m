function func_hdl = choosedialog

func_hdl.choosedialog_range = @choosedialog_range;
func_hdl.dialog_tiltcorr = @dialog_tiltcorr;
func_hdl.unwrap_prepare = @unwrap_prepare;

end

function [pb, ph, keep_cond, offset_pi2] = choosedialog_range(titl_box, max_amp)
% choosedialog : open dialog box for choices
%
% edited 2016.3.31 by Maxime PINSARD
%
% [pb, ph, keep_cond, offset_pi2] = choosedialog(titl_box, max_amp)

d = dialog('Position',[300 300 300 200],'Name',['Specific range in ' titl_box]); % left bottom width height

%         txt = uicontrol('Parent',d,...
%             'Style','text',...
%             'Position',[20 200 250 40],...
%             'String',['Choose your specific range in '
%             titl_box]);= 0;
ph = max_amp;
pb = 0;
offset_pi2 = 0;

txt_min = uicontrol('Parent',d,...
    'Style','text',...
    'Position',[50 170 25 25],...
    'String','Min');

txt_max = uicontrol('Parent',d,...
    'Style','text',...
    'Position',[50 145 25 25],...
    'String','Max');


minrange = uicontrol('Parent',d,...
    'Style','edit',...
    'Position',[90 170 100 25],...
    'String','0',...
    'Callback',@min_cllbck);

maxrange = uicontrol('Parent',d,...
    'Style','edit',...
    'Position',[90 145 100 25],...
    'String',num2str(max_amp),...
    'Callback',@max_cllbck);

keep_cond = 0;


btn2 = uicontrol('Parent',d,...
    'Style','check',...
    'Position',[50 100 150 40],...
    'String','Keep existing range',...
    'Callback',@chck_cllbck);

cgrange = uicontrol('Parent',d,...
    'Style','check',...
    'Position',[50 70 150 40],...
    'String','Apply a pi/2 offset',...
    'Callback',@cg_cllbck);

btn = uicontrol('Parent',d,...
    'Position',[50 30 130 40],...
    'String','OK',...
    'Callback','delete(gcf)'); %#ok<*NASGU>


% Wait for d to close before running to completion
uiwait(d);

%              keep_cond = get(btn2, 'Value');

    function chck_cllbck(btn2, callbackdata)
        keep_cond = get(btn2, 'Value');
        if keep_cond
            set(maxrange, 'Enable', 'off');
            set(minrange, 'Enable', 'off');
            set(maxrange, 'String', num2str(max_amp));
            set(minrange, 'String', 0);
        else
            set(maxrange, 'Enable', 'on');
            set(minrange, 'Enable', 'on');
        end
    end

    function cg_cllbck(cgrange, callbackdata)
        offset_pi2 = get(cgrange, 'Value');
        
    end

    function min_cllbck(minrange, callbackdata)
        pb =  str2double(get(minrange,'String'));
        
    end

    function max_cllbck(maxrange, callbackdata)
        ph =  str2double(get(maxrange,'String'));
        
    end

end

function [unwrap_flag, ch_detilt] = dialog_tiltcorr
% choosedialog : open dialog box for choices
%
% edited 2019.3.31 by Maxime PINSARD
%

unwrap_flag = 0;
ch_detilt = 1;

d = dialog('Position',[300 300 300 200],'Name','Params untilt'); % left bottom width height

btn2 = uicontrol('Parent',d,...
    'Style','check',...
    'Position',[50, 170, 200, 25],...
    'String','unwrap phase',...
    'Callback',@chck_unwrap_cllbck);

useunwr = uicontrol('Parent',d,...
    'Style','check',...
    'Position',[150, 170, 200, 25],...
    'String','use prev unwrap',...
    'Callback',@useunwr_cllbck);

surface_fit = uicontrol('Parent',d,...
    'Style','radio',...
    'Position',[50 115 150 25],...
    'String','Surface fit (best)',...
    'Callback',@surface_fitradio_cllbck);

coplanar_points = uicontrol('Parent',d,...
    'Style','radio',...
    'Position',[50 70 150 40],...
    'String','Coplanar points',...
    'Callback',@coplanar_pointsradio_cllbck);

load_mat = uicontrol('Parent',d,...
    'Style','radio',...
    'Position',[50 30 130 40],...
    'String','Load .mat',...
    'Callback',@load_matradio_cllbck);

% cgrange = uicontrol('Parent',d,...
%     'Style','check',...
%     'Position',[50 70 150 40],...
%     'String','Apply a pi/2 offset',...
%     'Callback',@cg_cllbck);
%
uicontrol('Parent',d,...
    'Position',[190 50 100 110],...
    'String','OK',...
    'Callback','delete(gcf)');


% Wait for d to close before running to completion
uiwait(d);

    function useunwr_cllbck(useunwr, callbackdata)
        if get(useunwr, 'Value') % check
            unwrap_flag = 2;
            set(btn2, 'Value', 0);
        elseif get(btn2, 'Value') % check
           unwrap_flag = 1;
        else
           unwrap_flag = 0;
        end
    end

    function chck_unwrap_cllbck(btn2, callbackdata)
       if get(btn2, 'Value') % check
           unwrap_flag = 1;
           set(useunwr, 'Value', 0);
       elseif get(useunwr, 'Value') % check
           unwrap_flag = 2;
       else
           unwrap_flag = 0;
       end
    end

    function surface_fitradio_cllbck(surface_fit, callbackdata)
        if get(surface_fit, 'Value')
            ch_detilt = 1;
            set(coplanar_points, 'Value', 0);
            set(load_mat, 'Value', 0);
        end
        
    end

    function coplanar_pointsradio_cllbck(coplanar_points, callbackdata) %#ok<*INUSD>
        if get(coplanar_points, 'Value')
            ch_detilt = 2;
            set(surface_fit, 'Value', 0);
            set(load_mat, 'Value', 0);
        end
        
    end

    function load_matradio_cllbck(load_mat, callbackdata)
        if get(load_mat, 'Value')
            ch_detilt = 3;
            set(coplanar_points, 'Value', 0);
            set(surface_fit, 'Value', 0);
        end
        
    end

end

function handles = unwrap_prepare(handles)

ans1 = get(handles.h_phase, 'Children'); phase1 = ans1.CData;

prompt = {'plot maps', 'Auto-sead', 'Use 1st frame for mag. ? (1, 0 use iterf. contr., >1 to load SHG image)', ...
'ADvanced plot, thresholding', 'Colormap (1 = jet, 0 = hsv)', 'Save in .mat auto', 'Unwrap where (0 for 2D, 1 for X, 2 for Y)?'};
dlg_title = 'Paramètres';
num_lines = 1;
def = {'1', '1', '0', '1', '0', '0', '0'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
plot_now = str2double(answer{1});
auto_sead = str2double(answer{2});
use_frst_frame = str2double(answer{3});
ch_adv_plt = str2double(answer{4});
ch_map = str2double(answer{5});
save_mat = str2double(answer{6});
unwrapwh = str2double(answer{7});

if use_frst_frame == 1
    magn = handles.img_3D(:,:,1);
elseif use_frst_frame > 1
    [fname1, folder_name, FILTERINDEX] = uigetfile('*.tif','Select your SHG .tif !','shg.tif');
    
    if (~FILTERINDEX)
        warning('File not chosen ! I`ll use interf. contr.');
        magn = handles.amp_cell{end};
    else
        magn = single(imread(fullfile(folder_name, fname1)));
    end
else
    magn = handles.amp_cell{end};
end
if unwrapwh == 0
    handles.im_unwrapped = unwrap_2D_MP( phase1, magn, plot_now, auto_sead, ch_map, ch_adv_plt, save_mat);
else
%     if unwrapwh == 1 % X
    handles.im_unwrapped = unwrap(phase1*pi,[],unwrapwh)/pi;
%     end
    figure(1); imagesc(handles.im_unwrapped); colormap(jet); colorbar; axis image;
end
end
