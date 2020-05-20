function varargout = I_SHG_GUI(varargin)
% I_SHG_GUI MATLAB code for I_SHG_GUI.fig
%      I_SHG_GUI, by itself, creates a new I_SHG_GUI or raises the existing
%      singleton*.
%
%      H = I_SHG_GUI returns the handlce to a new I_SHG_GUI or the handle to
%      the existing singleton*.
%
%      I_SHG_GUI('CALLBACK',hObject,~,handles,...) calls the local
%      function named CALLBACK in I_SHG_GUI.M with the given input arguments.
%
%      I_SHG_GUI('Property','Value',...) creates a new I_SHG_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before I_SHG_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to I_SHG_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help I_SHG_GUI

% Last Modified by GUIDE v2.5 04-Oct-2019 15:49:28
% Maxime Pinsard
clc
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @I_SHG_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @I_SHG_GUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before I_SHG_GUI is made visible.
function I_SHG_GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to I_SHG_GUI (see VARARGIN)


%% Initialisation et entrée des paramètres

handles.genpath1 = 'C:\Users\pc\Documents\These';%'C:\Users\arash1372\Desktop\phd';%
addpath(fullfile(handles.genpath1, 'codes Matlab\Codes_I-SHG'));

foder_files =  '\Phase_EOM\anlg_galvos_mirrorOK';%'C:\Users\arash1372\Desktop\phd'; %
handles.pathfun = fullfile(handles.genpath1, foder_files); %'DATA images\'; 
% main path, you can change it to browse faster (optional)
try %#ok<TRYNC>
    cd(handles.pathfun);
end
% Mémoriser le folder matlab de départ
handles.current_folder = pwd;

%% Buttons def.

handles.point_center_batch_plotting = 0;
set(handles.load_stck_button, 'Enable', 'off');
handles.crop = get(handles.crop_selec_button,'Value') + 1;
handles.filt = get(handles.filter_button,'Value');
handles.resx = str2double(handles.res_x_edt.String);
handles.resy = str2double(handles.res_y_edt.String);
handles.start_phase = str2double(handles.first_phase_edt.String);
handles.diff_phase = str2double(handles.step_ps_edt.String);
handles.langue = handles.language_menu.Value;
set(handles.int_contrast_button, 'Enable', 'off');
set(handles.saving_popup, 'Enable', 'off');
set(handles.cos_button, 'Enable', 'off');
set(handles.relative_phase_button, 'Enable', 'off');
set(handles.hist_range_button, 'Enable', 'off');
handles.init_hist_ok = 0;
% % set(handles.center_coeff_edt, 'Visible', 'off');
% % set(handles.center_coeff_txt, 'Visible', 'off');
% % set(handles.cross_coeff_edt, 'Visible', 'off');
% % set(handles.cross_coeff_txt, 'Visible', 'off');
% % set(handles.diag_coeff_edt, 'Visible', 'off');
% % set(handles.diag_coeff_txt, 'Visible', 'off');
handles.undocked_fig1 = 0; % no outside fig a priori
handles.undocked_fig2 = 0; % no outside fig a priori
handles.FILTERINDEX = 0;
handles.contrast_dflt = 0.75;
set(handles.contrast_edt, 'String', num2str(handles.contrast_dflt));
handles.nbinsx_dflt = 100;
handles.nbinsy_dflt = 100;
handles.saved_ref = 0;
handles.contr00=0;
set(handles.interf_pannel, 'Visible', 'off');
set(handles.analyze_panel, 'Visible', 'off');
set(handles.hist_button, 'Enable', 'off');
set(handles.two_d_hist_popup, 'Enable', 'off');
set(handles.nbinsX_edt, 'Enable', 'off');
set(handles.nbinsY_edt, 'Enable', 'off');
set(handles.hist_range_button, 'Visible', 'off');
set(handles.min_hist_edt, 'Visible', 'off');
set(handles.max_hist_edt, 'Visible', 'off');
set(handles.ampmin_txt, 'Visible', 'off');
set(handles.ampmax_txt, 'Visible', 'off');
set(handles.plot_diff_menu, 'Visible', 'off');
set(handles.sat_value_slider, 'Visible', 'off');
set(handles.sat_txt, 'Visible', 'off');
set(handles.fit_hist_button, 'Visible', 'off');
set(handles.baseline_hist_edt, 'Visible', 'off'); set(handles.baseline_hist_lbl, 'Visible', 'off');
set(handles.apply_offset_button, 'Enable', 'off')
handles.offset_pi2 = 0;
handles.offset = 0;
set(handles.plot_phase, 'Visible', 'off');
set(handles.cmap1_menu, 'String', {'Colormap left', 'Red/Black/Green', 'Black/Red/Black/Green', 'HSV', 'Parula', 'Cubehelix', 'Grey'});
set(handles.cmap2_menu, 'String', {'Colormap hist', 'Cubehelix', 'HSV', 'Parula', 'Grey', 'Red/Black/Green', 'Black/Red/Black/Green'});
handles.name_c_p_menu_str={'Colormap phase', 'Red/Black/Green', 'Black/Red/Black/Green', 'B/R/B/G with saturation','HSV','HSV with saturation','Parula','Cubehelix','Grey'};
set(handles.cmap_p_menu, 'String', handles.name_c_p_menu_str);
handles.file_not_good = 1;
handles.undocked_phase = 0;
set(handles.reset_range_chck, 'Enable', 'off');
handles.hhist1 = 0;
handles.order_frame_vect=0;
set(handles.fit_hist_button, 'String', ...
    {'Apply on hist fit function ...', ' Gaussian MANU', ' Lorentzian MANU', ' Prod. Gaussian/Lorentzian MANU',...
    ' Sum Gaussian/Lorentzian (pseudo-Voigt) MANU ', ' Pearson VII (MANU) ', ' Von Mises (MANU) ', ' 2 Gauss', ' 2 Lorentz', ' 2 Prod. Gaussian/Lorentzian ', ' 2 Sum Gaussian/Lorentzian (pseudo-Voigt)', ' 2 Pearson VII ', ' 2 Von Mises'});
set(handles.edt_pearson, 'Visible', 'off');
set(handles.text_pearson, 'Visible', 'off');
set(handles.action_map_menu, 'Visible', 'off');
set(handles.expboth2miji_chck, 'Visible', 'off');
set(handles.edt_pearson, 'String', 0.5);
handles.h_princ = handles.plot_princ;
h_rndm1 = figure(1);
handles.h_princ_out = axes;
close(h_rndm1); % fig deleted = axes non valid
handles.h_phase = handles.plot_phase;
h_rndm1 = figure(1);
handles.h_phase_out = axes;
close(h_rndm1); % fig deleted = axes non valid
handles.phase_calib = [];
handles.amp_calib = []; handles.do_realign_calib = 0;
handles.contrast = handles.contrast_dflt;
handles.h_hist = handles.plot_second;
h_rndm1 = figure(1);
handles.h_hist_out = axes;
close(h_rndm1); % fig deleted = axes non valid
handles.exp_phase_before_ictr = 1;
handles.fit_clrwheelBRBG = 0;
set(handles.contr_radio, 'Visible', 'on');
set(handles.raw_radio, 'Visible', 'on');
handles.dipimage = 0;
handles.sigma_med = 0;
handles.weight = 0;
handles.coef_centre = 1;
handles.coef_cross = 1;
handles.coef_diag = 1;
handles.mean_phase_meth = 1; % mean
handles.max_c = [];
handles.phase_cell = {};
handles.amp_cell = {};
handles.err_cell = {};
handles.phase_cell00 = {};
handles.amp_cell00 = {};
handles.range_sat = 0;
handles.x_phase = 0;
handles.flag_save_fig_out = 0;
handles.expboth = 0;
handles.surf_fit_tilt_auto_bacth = 0;
handles.int_y = str2double(handles.nbinsY_edt.String);
% set(handles.double_contr_chck, 'Enable', 'off');
set(handles.saving_popup, 'Enable', 'off');

handles.screensize = get( 0, 'Screensize' ); % to get screen size
handles.fact = 4/5;
handles.left_offset_fig = 40;
handles.top_offset_fig = 100;
handles.predef_foldr_put_mijph = 0;
handles.predef_filenm_put_mijph = 0;

handles.range_interf_keep = 0;
handles.choose_file_tiff_legacy = '';
set(handles.it_max_algo_edt, 'Visible', 'off');
set(handles.max_it_algo_txt, 'Visible', 'off');
set(handles.epsilon_ph_th_edt, 'Visible', 'off');
set(handles.eps_threshold_ph_txt, 'Visible', 'off');
set(handles.eps_threshold_tilt_txt, 'Visible', 'off');
set(handles.epsilon_tilt_th_edt, 'Visible', 'off');
setappdata(0 , 'param_algo_changed', 0); 
handles.cmap_default1 = hsv;
handles.str_ratiof_sigma = 'ratio f sigma, shg';
handles.str_hist_ctr_shg = 'Corr. hist ctr/shg';
handles.double_contrast = str2double(handles.double_contr_chck.String(end)); % if order 0, 180, 360, ... use (0-180) - (180-360) to double the contrast
% % handles.double_contr_chck.Enable = 'on';
double_contr_chck_Callback(handles.double_contr_chck, 0, handles);

% handles.mij_obj = false;
global mij_obj
mij_obj = false;

[handles.xTitle_dflt, handles.yTitle_dflt, handles.phi_mat_default, handles.cmap_brbgsat, ...
    handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz, ...
    handles.cmap_cubehelix, handles.cmap_redgreen, handles.cmap_blkredblkgrn, handles.cmap_hsvsat] = load_axis_param;

handles.cmap_brbgsat_dflt = handles.cmap_brbgsat;
handles.cmap_hsvsat_dflt = handles.cmap_hsvsat;
handles.cmap_sat_dflt = handles.cmap_brbgsat;
handles.cmap_sat_dflt = handles.cmap_hsvsat;
handles.Titre1 = ''; handles.title_hist1 = '';
 first_phase_edt_Callback(hObject, 0, handles);
 [ handles.Titre1, handles.Titre2, handles.Titre3, handles.Titre4, handles.Titre5, handles.Titre6, handles.Counts, handles.Yaxis1, handles.Yaxis2, handles.Legendehisto ] = string_axis_ISHG( handles.langue );
handles.Titre4_modified = handles.Titre4;

handles.title_hist1 = handles.Titre6;
handles.yaxis_hist1 = handles.Yaxis2;

% Choose default command line output for I_SHG_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes I_SHG_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = I_SHG_GUI_OutputFcn(hObject, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in open_dir_main.
function open_dir_main_Callback(hObject, ~, handles)
% hObject    handle to open_dir_main (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% L'utilisateur sélectionne le fichier d'images qu'il veut analyser
folder_name=uigetdir(handles.pathfun,'Select your path containing your file(s) !');
if isa(folder_name, 'char') % a folder has been chosen
    cd(folder_name);
end
handles.folder_name = folder_name;
fprintf('Path chosen : %s \n', handles.folder_name);
set(hObject, 'String', ['...' handles.folder_name(end-min(38, length(handles.folder_name)-1):end)]);
set(hObject, 'BackgroundColor', [1 1 204/255]); set(hObject, 'ForegroundColor', [0 0 0]); % turn the color yellow

% set(handles.gen_param,'Enable','on') ;
% set(handles.open_param_button,'Enable','on') ;

guidata(hObject, handles);

function [handles, file_not_good] = numimg_fname_func(hObject,handles, fname0, file_not_good)

if isa(fname0, 'cell')% many single files
    num_images=length(handles.fname);
    fname = handles.fname{1};
    set(handles.ind_files_chck, 'Value', 1);
else % one stack file
    fname = fullfile(handles.folder_name, fname0);
    warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
    warning('off','MATLAB:imagesci:tifftagsread:badTagValueDivisionByZero');
    InfoImage=imfinfo(fname); % disp a warning 'Division by zero when processing YResolution'
    warning( 'on','MATLAB:imagesci:tiffmexutils:libtiffWarning');
    warning('on','MATLAB:imagesci:tifftagsread:badTagValueDivisionByZero');
    warning('ON', 'all');
%     mImage=InfoImage(1).Width;
%     nImage=InfoImage(1).Height;
    num_images=length(InfoImage);
end
fprintf('File chosen : %s \n', fname);
set(hObject, 'String', ['...' fname(end-min(57, length(fname)-1):end)]);
set(hObject, 'BackgroundColor', [0.8, 1, 0.8]); set(hObject, 'ForegroundColor', [0 0 0]); % turn the color green
cd(handles.folder_name);
fprintf('Path chosen : %s \n', handles.folder_name);
set(handles.open_dir_main, 'String', ['...' handles.folder_name(end-min(40, length(handles.folder_name)-1):end)]);
set(handles.open_dir_main, 'BackgroundColor', [1 1 204/255]); set(handles.open_dir_main, 'ForegroundColor', [0 0 0]); % turn the color yellow

if ~strcmp(fname(end-2:end), 'tif')
    file_not_good = 1;
else
    set( handles.nb_img_in_stack_edt, 'String', num_images);
    moy = str2double(get(handles.avg_img_edt,'String'));
    slice_per_step = str2double(handles.slice_per_step_edt.String);
    if ~isnan(slice_per_step)
        set(handles.step_ps_edt, 'String', num2str(moy*360/num_images*slice_per_step/2));
        set(handles.snd_step_edt, 'String', num2str(moy*360/num_images*slice_per_step/2));
    end
end

% --- Executes on button press in choose_file_tiff.
function choose_file_tiff_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to choose_file_tiff (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.choose_file_tiff_legacy = get(handles.choose_file_tiff,'String');
if isempty(handles.choose_file_tiff_legacy); handles.choose_file_tiff_legacy= '0'; end
if strcmp(handles.choose_file_tiff_legacy(1), '&') % special for direct call
    pp='';ll=regexp(pwd,filesep,'split'); for i=1:length(ll); pp=fullfile(pp,ll{i});if strcmp(ll{i}, 'tmp'); break; end; end
    cd(pp);
    if strcmp(handles.choose_file_tiff_legacy(end-3:end), '.tif') % a priori, a file stack not folder
        handles.folder_name =  handles.open_dir_main.String;
        handles.fname = handles.choose_file_tiff_legacy(2:end); % remove '&'
    else % many single files
        handles.folder_name = handles.choose_file_tiff_legacy(2:end); % remove '&'
        [handles.fname, ~,~,~]=list_fname_func(handles.folder_name,[], 0);
    end
    if ~strcmp(handles.folder_name(2:3), ':\') % not a full path
        handles.folder_name = fullfile(pp, handles.folder_name);
    end
    FILTERINDEX =1;
else% normal
    set(handles.choose_file_tiff,'String', '');

    if (get(handles.batch_chck, 'Value') && get(handles.ind_files_chck, 'Value'))
        handles.dirnames = uigetdirMultiSelect('', 'Select many folders'); % custom func
        if isempty(handles.dirnames); FILTERINDEX = 0;
        else; FILTERINDEX = 1;
        end
        handles.fname = 0;
        handles.folder_name = handles.dirnames{1};
    else
        if ~get(handles.ind_files_chck, 'Value') % stack
            msg_open = 'Select your stack of images !';
            dflt = 'stack.tif';
        else
            msg_open = 'Select your single images !';
            dflt = 'img.tif';
        end
        [handles.fname, handles.folder_name, FILTERINDEX] = uigetfile('*.tif', msg_open, dflt,'MultiSelect','on');
    end
end
if length(handles.folder_name)>1; try cd(handles.folder_name);end; end %#ok<TRYNC>
file_not_good = 0; % a priori
batch_mode = 0;
selfiles_before = 0; % for single files in each folder bacth, select by hand each files if 1

if (~FILTERINDEX && (~get(handles.batch_chck,'Value') || ~handles.FILTERINDEX))
    set(handles.choose_file_tiff,'String', handles.choose_file_tiff_legacy);
    error('File not chosen ! Program ends.');
else
   
    if (~get(handles.ind_files_chck, 'Value') && get(handles.batch_chck, 'Value')) % a number of file > 1 (stacks) has been chosen, batch
        
        batch_mode = 1;
        
    elseif ~get(handles.batch_chck, 'Value') %isa(handles.fname, 'char') % a correct file has been chosen no batch

        [handles, file_not_good] = numimg_fname_func(hObject,handles, handles.fname, file_not_good);
            
    elseif (get(handles.batch_chck, 'Value') && get(handles.ind_files_chck, 'Value')) % batch of folders with single tif files 
        if selfiles_before
            for i=1:length(handles.dirnames) %#ok<UNRCH>
                [handles.fname{i}, handles.folder_name, FILTERINDEX] = uigetfile('*.tif', 'select your files in this dir.', fullfile(handles.dirnames{i}, 'img.tif'),'MultiSelect','on');
            end
            % % cell of cells
        end
        batch_mode = 1;
    else
        file_not_good = 1;
    end
    
    if batch_mode
       try
            handles = batch_util(hObject, handles, '1', 'treat', '1', '1', '1', '0', '0', '0', '1','raw', '0'); % % standard_notload, nm1 ,put_peaks ,do_mijiph,do_clrwheel,saved_ref ,fit_batch ,fit_clrwheelBRBG,range_interf_keep,prefix ,sel_files
       catch ME
            handles.predef_foldr_put_mijph = 0;
            handles.predef_filenm_put_mijph = 0;
            handles.range_interf_keep = 0;
            handles.fit_clrwheelBRBG = 0;
            handles.expboth = 0;
            disp(ME); 
            for i = length(ME.stack):-1:1
                try %#ok<TRYNC>
                    fprintf(2, ME.stack(i).name);
%                     fprintf('\n');
%                     fprintf(2, ME.stack(i).message);
                    fprintf('\n');
                    disp(ME.stack(i).line);  
                end
                fprintf('\n');
            end
       end
    end

    if file_not_good
        disp('File not good !');
        set(hObject, 'String', 'File not good !');
        set(hObject, 'BackgroundColor', [1 153/255 153/255]); set(hObject, 'ForegroundColor', [0 0 0]); % turn the color green
    end
end

handles.file_not_good = file_not_good;
handles.FILTERINDEX = FILTERINDEX;

guidata(hObject, handles);

function [handles,standard_notload,put_peaks ,do_mijiph,do_clrwheel, sel_files,num_fit, prefix] = prompt_batch_util(handles, standard_notload,nm1 ,put_peaks ,do_mijiph,do_clrwheel,saved_ref ,fit_batch ,fit_clrwheelBRBG,range_interf_keep,prefixdf ,sel_files, use_prev_str)
    
    def = {standard_notload, nm1 , put_peaks, do_mijiph, do_clrwheel, saved_ref , '0', fit_batch , fit_clrwheelBRBG, '1', range_interf_keep, '0', prefixdf , sel_files};
    prefix  = inputdlg( {'You chose batch mode (auto-actions): Standard(1) or load treated (0)', ...
    'Folder for treat (will be created)', 'put peaks at -0.5 and 0.5\pi ?',...
    'Do Miji phase save (1 for auto, 2 for asking save location)?', 'Clr wheel hist save ?', ...
    sprintf('ask for ld calib %s', use_prev_str), 'save imgs uncorrected (if calib)', ...
    'fit hist2D save_res_array?', 'fit hist2D with clrwheelBRBG', ...
    'func for fit (1: Gauss, 2 Lorentz, 3 Prod. Gaussian/Lorentzian, 4 Sum Gaussian/Lorentzian (pseudo-Voigt),5 Pearson VII ,6  Von Mises (+6 for auto 2)',...
     'keep interf. imposed range', ...
    'do surf fit tilt (1:lin2D, 2:parabola, 0 no and >2 ask)' , 'prefix of saved filename:', 'man. sel. files in each folder'},...
     'Batch', 1, def); 

    if isempty( prefix)
    %              prefix = def;
        handles.predef_foldr_put_mijph=0;
        handles.predef_filenm_put_mijph=0;
        error('user canceled.')
    end
    standard_notload =str2double(prefix{1});
    handles.nm1 = prefix{2}; put_peaks = str2double(prefix{3}); 
    do_mijiph= str2double(prefix{4}); do_clrwheel= str2double(prefix{5}); 
    handles.saved_ref = str2double(prefix{6}); handles.save_raw_befcalib = str2double(prefix{7}); handles.fit_batch = str2double(prefix{8}); 
    handles.fit_clrwheelBRBG = str2double(prefix{9}); num_fit = str2double(prefix{10}); num_fit=max(num_fit, 1 ); % min(
    handles.range_interf_keep = str2double(prefix{11});
    handles.surf_fit_tilt_auto_bacth = str2double(prefix{12});
    sel_files= str2double(prefix{end});
    prefix =prefix{end-1};

function [fname, dirnames,vect_loop,max_l]=list_fname_func(folder_name, dirnames,vect_loop)
    listing = dir(folder_name);
    fname = {};
    max_l = 0;
    for i3 =3:length(listing) % '.' '..' % loop of files in folder
        nm=listing(i3).name;
        if ~listing(i3).isdir % not a folder 
            if (length(nm) >= 5 && strcmp(nm(end-3:end), '.tif') && ~strcmp(nm(end-10:end-4), 'SUMISHG') ...
            && (~isnan(str2double(nm(end-4))) || strcmp(nm(end-4), '°')) && ~isnan(str2double(nm(1))))
               fname{end+1} = nm; %#ok<AGROW>
               if length(nm) > max_l; max_l = length(nm); end
            end
        else % folder
            if (~strcmp(nm, 'treat') && (length(nm) < 5 || ~strcmp(nm(end-3:end), '.npy')))
                dirnames{end+1} = fullfile(folder_name, nm); %#ok<AGROW>
                vect_loop = vect_loop+1;
            end
        end
    end

function handles = batch_util(hObject, handles, standard_notload,nm1 ,put_peaks ,do_mijiph,do_clrwheel,saved_ref ,fit_batch ,fit_clrwheelBRBG,range_interf_keep,prefixdf ,sel_files)

use_prev_str = ''; % dflt
if (~isempty(handles.phase_calib) && sum(handles.phase_calib(isnan(handles.phase_calib))) ~= 0); use_prev_str = '(-1 for use prev)'; end

[handles,standard_notload,put_peaks ,do_mijiph,do_clrwheel, sel_files,num_fit, prefix] = ...
    prompt_batch_util(handles, standard_notload,nm1 ,put_peaks ,do_mijiph,do_clrwheel,saved_ref ,fit_batch ,fit_clrwheelBRBG,range_interf_keep,prefixdf ,sel_files, use_prev_str);

if ((isempty(handles.phase_calib) || sum(handles.phase_calib(isnan(handles.phase_calib))) == 0) && handles.saved_ref); handles.saved_ref = 1; end
if handles.fit_batch
    disp( 'You chose several files : it will treat them as batch and save max of the phase in a cell file. Choose prefix of saved filename');
end
if (sel_files && ~selfiles_before && standard_notload)
    for i=1:length(handles.dirnames)
        [handles.fname{i}, handles.folder_name, ~] = uigetfile('*.tif', 'select your files in this dir.', fullfile(handles.dirnames{i}, 'img.tif'),'MultiSelect','on');
    end
end

%         callback_contrast = get(handles.int_contrast_button,'Callback'); % get the callback of the button that we want to activate
%         callback_load = get(handles.load_stck_button,'Callback');
%         callback_phase = get(handles.relative_phase_button,'Callback');
%         callback_corr_phase = get(handles.correct_phase_subtract,'Callback');
% %         callback_hist2d = get(handles.two_d_hist_popup,'Callback');
%         callback_map_range = get(handles.hist_range_button, 'Callback');
set(handles.histmode_menu, 'Value', 3); handles.shg = 0;
% % callback_action_map_menu = get(handles.action_map_menu,'Callback');

if handles.saved_ref == 1
    calib.do_realign = 0; % init struct
%         handles.saved_ref = 2-menu('Subtract current ref ?', 'Yes', 'No');
% %     [fname, folder_name, FILTERINDEX] = uigetfile('*.mat', 'Select your .mat file of quartz (calib) phase (cancel for no corr)!', 'calib.mat');
   func_hdl = load_stack_plot_ISHG; 
    [calib, ~, ~, ~, ~, ~, ~, ~, folder_name] = correct_phase_ref_subtract( 0, ...
0, 0, 0, 0, -1, calib, 0, 0, func_hdl.scaling_img, 'fig'); 
    FILTERINDEX = calib.FILTERINDEX;
    handles.phase_calib = calib.phase;
    handles.amp_calib=calib.amp ;
    handles.do_realign_calib=calib.do_realign ;
    handles.choose_file_tiff_legacy = ['...' folder_name(end-min(57, length(folder_name)-1):end)];
    set(handles.choose_file_tiff, 'String', handles.choose_file_tiff_legacy);

elseif handles.saved_ref == -1
    if ~isfield(handles, 'phase_calib')
        FILTERINDEX = 0;
    else % ok
        FILTERINDEX = 1;
    end
else; FILTERINDEX = 0;
end

if ~FILTERINDEX
    handles.saved_ref = 0;
    disp('User cancelled the corr file loading : no correction will be applied')
else
% %     if isa(fname, 'char') % a correct file has been chosen
% %         if strcmp(fname(end-2:end), 'mat')
% %             m=load(fullfile(folder_name, fname));
% %             m2=struct2cell(m);
% %             handles.phase_calib=m2{1};
% %             handles.saved_ref = 1;
% %         end
% %     end
    handles.saved_ref = 1;
end
% %         handles.crop_batch = 2-menu('Do you want to crop (see code of load stack for parameters?)', 'Yes', 'No');

handles.nm_clrwheel = sprintf('clrwheel%s.fig',prefix);
handles.nm_clrwheelBRBG = sprintf('clrwheel%s_corrref.fig',prefix);

if handles.saved_ref % % calib ok

    handles.predef_filenm_put_mijph = [handles.predef_filenm_put_mijph, '_corrref'];
    
% %         else '_phmp'
    handles.nm_clrwheel = 'clrwheel_corrref.fig';
    handles.nm_clrwheelBRBG = [handles.nm_clrwheelBRBG(1:end-4), 'BRBG', '.fig'];
end

if standard_notload
    if get(handles.ind_files_chck, 'Value') % % many singles files
        vect_loop = length(handles.dirnames);
    else % % stacks
        vect_loop = length(handles.fname);
        cd (handles.folder_name);
    end
else % load treated
    vect_loop = 100; % high values, user will cancel to stop
end
plot_hist_off = 0;
ii = 0;
handles.expboth = 1;
while ii < vect_loop
    ii = ii+1;
    set(handles.cmap_p_menu, 'Value',5); % HSV
     handles.cmap_default1 = hsv;
    try, close (1);close (50);close (51); catch ME; disp('nofig2close'); disp(ME); end  %#ok<NOCOM>
    handles.predef_foldr_put_mijph = 0;
    handles.predef_filenm_put_mijph = 0;
    safemode = 0;
    if standard_notload
        handles.phase_cell = {};
        handles.amp_cell = {};
        handles.err_cell = {};
        handles.x_cell = {};
        handles.y_cell = {};
        str0={'View plot phase ...'};
        set(handles.plot_phase_popup, 'Value', 1);
        set(handles.plot_phase_popup, 'String', str0);
        guidata(hObject, handles);
        if get(handles.ind_files_chck, 'Value') % % many singles files
            nmdsp = handles.dirnames{ii};
            handles.folder_name=handles.dirnames{ii} ; % % will be cd by load
            [handles.fname, handles.dirnames,vect_loop,max_l]=list_fname_func(handles.folder_name, handles.dirnames,vect_loop);
            if length(handles.fname) < 3 % need minimum 3 ps
            if ~isempty(handles.fname) % no files
               disp('WARN: not enough ps in this folder'); disp(handles.folder_name);
            end
            continue;
            end
            if length(handles.fname) < 360/str2double(get(handles.step_ps_edt, 'String'))*str2double(get(handles.slice_per_step_edt, 'String'))/2
                disp('WARN: not enough ps'); safemode = 1;
            end
            fname_cell = handles.fname;
            for k =1:length(handles.fname)
                if length(handles.fname{k}) < max_l
                    handles.fname{k} = [handles.fname{k}(1:end-1-4), '0', handles.fname{k}(end-4:end)];
                end
            end
            [~,ind]=sort(handles.fname);
            handles.fname = fname_cell(ind);
            set( handles.nb_img_in_stack_edt, 'String', num2str(length(handles.fname)));
            
        else  % stack
            if isa(handles.fname, 'cell'); nmdsp = handles.fname{ii}; else; nmdsp = handles.fname; end
        end
        handles.count_batch = ii;
        fprintf('file # %d ... %s\n', handles.count_batch, nmdsp);
        [handles, file_not_good] = numimg_fname_func(hObject,handles, handles.fname, 0);
        if file_not_good; fprintf(2,'no tiff ??\n'); return; end
        mkdir(handles.folder_name, handles.nm1); % treat
        guidata(hObject, handles);
        if safemode; try, handles = ld_stck3d_util(hObject, handles, ''); end %#ok<NOCOM,TRYNC>
        else; handles = ld_stck3d_util(hObject, handles, ''); end
        guidata(hObject, handles);
        handles = int_contrast_util(hObject, handles);
        guidata(hObject, handles);
        handles = rel_phase_util(hObject, handles);
        guidata(hObject, handles);
% %             hgfeval({callback_load, hObject,eventdata});
% %             hgfeval({callback_contrast, hObject,eventdata});
% %             relative_phase_button_Callback(   hObject,eventdata , handles);
% %         handles.phase_cell{end+1} = handles.phase1;
% %         handles.amp_cell(end+1) = {handles.amp1};
    if handles.min_hist_edt.Value == 0
        handles.max_hist_edt.Value = 1e9; % a very large value to be sure to take every value of phase into consideration
    end
%             hgfeval({callback_map_range,hObject,eventdata});
%             hist_range_button_Callback
% %             if handles.saved_ref
% %                 fprintf('saved ref \n');
% %                 hgfeval({callback_corr_phase,hObject,eventdata});
% %             end
    else % load
        disp('LOAD MODE BACTH : load your treated files !')
        [handles,FILTERINDEX] = load_phase_util(hObject, handles, '.tif'); % 3 for .tiff
        if ~FILTERINDEX; break; end
    end
    
    batch_furthertreat_core(hObject, handles, prefix, do_mijiph, put_peaks,do_clrwheel, sel_files,num_fit, plot_hist_off);

    %             max_c{ii} = max(handles.hist2_ydata);
    try
        close('2');close('3'); % close('18');close('1'); close('19');
    catch ME
        disp(ME)
    end
%             input('wait') % !!!
end

for ii = 1:length(handles.fname)
    fprintf('%s\n', handles.fname{ii});
end
%         max_c = handles.max_c;
%         save('batch_max_hist.mat', 'max_c');

function batch_furthertreat_core(hObject, handles, prefix, do_mijiph, put_peaks,do_clrwheel, sel_files,num_fit, plot_hist_off)

if do_mijiph ~= 2 % 2 =  ask for location of save
    handles.predef_filenm_put_mijph = sprintf('%s_ctr×%d' , prefix, 2^(handles.double_contrast+1));
end
    
set(handles.open_dir_main, 'String', ['...' handles.folder_name(end-min(40, length(handles.folder_name)-1):end)]);    
if do_mijiph ~= 2 % % 2 for asking save location
    handles.predef_foldr_put_mijph =handles.folder_name;
    if (~strcmp(handles.folder_name(end-4:end), handles.nm1) && ~strcmp(handles.folder_name(end-5:end-1), handles.nm1))
        handles.predef_foldr_put_mijph = fullfile(handles.folder_name, handles.nm1);
    end
end
guidata(hObject, handles);

if (do_mijiph && handles.save_raw_befcalib) 
% %         hgfeval({callback_action_map_menu, hObject,0});
    if handles.save_raw_befcalib
        handles.predef_filenm_put_mijph_raw =  handles.predef_filenm_put_mijph;
    end
    nm = handles.predef_filenm_put_mijph;
    handles.predef_filenm_put_mijph = handles.predef_filenm_put_mijph_raw;
    handles =  expmiji_phase_util(handles, get(handles.expboth2miji_chck, 'Value'), 'phmp', sprintf('%s_ctr×%d' , prefix, 2^(handles.double_contrast+1)), handles.predef_foldr_put_mijph, handles.predef_filenm_put_mijph);
    handles.predef_filenm_put_mijph = nm; % restore
    guidata(hObject, handles);
end
if handles.saved_ref
    if handles.save_raw_befcalib
        handles.predef_filenm_put_mijph_raw =  handles.predef_filenm_put_mijph;
    end
    handles.predef_filenm_put_mijph = [handles.predef_filenm_put_mijph, '_corrref'];
end

if handles.saved_ref % % corr by ref, licit
%         a = get(handles.plot_phase, 'Children'); phase1 =  a.CData;
%          [calib, handles.amp_cell, handles.phase_cell, ...
% handles.x, handles.y, handles.x_cell, handles.y_cell, ~, folder_name] = correct_phase_ref_subtract( handles.amp_cell, ...
% handles.phase_cell, handles.resx, handles.x_cell, handles.y_cell, 0, calib, phase1, length(handles.plot_diff_menu.String)- handles.plot_diff_menu.Value );

    [handles, ~] = corr_ph_ref_util(hObject, handles);
    guidata(hObject, handles);
end
set(handles.action_map_menu, 'Value', 4); % miji phase 
if (put_peaks) % % corr position of peaks to have -0.5pi 0.5 pi
    switch plot_hist_off
        case 1
            if handles.fit_clrwheelBRBG; fit_clrwheelBRBG=handles.fit_clrwheelBRBG; handles.fit_clrwheelBRBG = 0;end
            if handles.fit_batch; fit_batch=handles.fit_batch; handles.fit_batch=0; end % just need 2D plot for now
            handles = two_d_hist(hObject, handles); % simple func
            handles.fit_batch = fit_batch; handles.fit_clrwheelBRBG = fit_clrwheelBRBG; % restore
        case 0
            int_x =  single(1./str2double(handles.nbinsX_edt.String));
            ph_hist = handles.phase_cell{end}(:);
            ph_hist(ph_hist==1.05) = [];
            handles.hist2_xdata = -1+int_x/2:int_x:1-int_x/2;
            handles.hist2_ydata = hist(ph_hist, handles.hist2_xdata); % %#ok<HIST> % does not do the plot
    end
    [~, ind] = max(handles.hist2_ydata); 
    cond_offset=abs(abs(handles.hist2_xdata(ind))-0.5) > 0.1;
    if cond_offset % % offset too high
        if handles.hist2_xdata(ind) <0
            off = -0.5-handles.hist2_xdata(ind);
        else
            off = 0.5-handles.hist2_xdata(ind);
        end
        set(handles.offset_neutral, 'String', num2str(off) ); 
        apply_offset_button_Callback(hObject, 0, handles);
        im=get(handles.plot_phase, 'Children'); handles.phase_cell{end+1}= im.CData; 
    end
     handles = add_data_in_list(hObject, handles, 1, 1, 1, 0, [1, Inf], [1, Inf]);
     guidata(hObject, handles);
end

if (do_mijiph) % && ~handles.fit_clrwheelBRBG) % in handles.predef_foldr_put_mijph
% %         hgfeval({callback_action_map_menu, hObject,0});
    handles =  expmiji_phase_util(handles, get(handles.expboth2miji_chck, 'Value'), 'phmp', sprintf('%s_ctr×%d' , prefix, 2^(handles.double_contrast+1)), handles.predef_foldr_put_mijph, handles.predef_filenm_put_mijph); % sprintf('BRBG%s', handles.predef_filenm_put_mijph)
end
guidata(hObject, handles);
if do_clrwheel % clrwheel normal HSV
    handles = clrwheel_twod_hist_util(handles);
    try, savefig(handles.hclrwh, fullfile(handles.predef_foldr_put_mijph, handles.nm_clrwheel) ,'compact'); %#ok<NOCOM>
    catch ME; disp(ME);
    end
    guidata(hObject, handles);
end
if (handles.fit_clrwheelBRBG || handles.fit_batch)
    if (plot_hist_off && put_peaks && cond_offset) % % offset too high
        close(50); % 2D hist
    end
    if handles.fit_batch
        handles = two_d_hist(hObject, handles); % simple func, will do the fit if asked
    else
        num_cell = get(handles.plot_phase_popup, 'Value')-1;
        [handles.hhist1, handles.hist2_xdata, handles.hist2_ydata] = hist2D_ISHG(  handles.fact,  handles.left_offset_fig, ...
        handles.top_offset_fig, handles.phase_cell{num_cell}, single(1./str2double(handles.nbinsX_edt.String)),  handles.phi_mat_default,  handles.Counts,  handles.Titre4_modified, ...
        handles.axes_font_size,  handles.xaxis_sz,  handles.yaxis_sz,  handles.title_sz,  handles.clrbr_tl_sz, handles.screensize, handles.offset_pi2, handles.hhist1 );

        set(handles.fit_hist_button, 'Value', 1+num_fit); % 9 = 2 lorentz
        handles = fit_hist_util(handles);
    end
    if handles.fit_clrwheelBRBG
        if strcmp(get(handles.plot_phase, 'Visible'), 'off')
            set(handles.plot_phase, 'Visible', 'on');
            set(handles.plot_princ, 'Visible', 'off');
            % colorbar(handles.plot_princ,'off')
% %                 set(handles.cmap_p_menu, 'Value', 3);
            set(handles.plot_phase_popup, 'Value', length(handles.plot_phase_popup.String));
            set(handles.cmap_p_menu, 'Value', length(handles.cmap_p_menu.String));
            plot_phase_popup_Callback(hObject, 0, handles);

        end
        handles = clrwheel_twod_hist_util(handles); % BRBG
        try, savefig(handles.hclrwh, fullfile(handles.predef_foldr_put_mijph, handles.nm_clrwheelBRBG) ,'compact'); %#ok<NOCOM>
        catch ME; disp(ME);
        end
        if (do_mijiph && handles.fit_clrwheelBRBG) % in handles.predef_foldr_put_mijph
%                 hgfeval({callback_action_map_menu, hObject,0});
            handles =  expmiji_phase_util(handles, get(handles.expboth2miji_chck, 'Value'), 'phmp', sprintf('%s_ctr×%d' , prefix, 2^(handles.double_contrast+1)), handles.predef_foldr_put_mijph, handles.predef_filenm_put_mijph);
        end
    end
end
if handles.surf_fit_tilt_auto_bacth % 1 or 2
    [~, ~, ~]=subplane_mp(handles.phase_cell{end},0,0, 0, handles.h_phase, ...
    0, handles.surf_fit_tilt_auto_bacth, 0);
end
guidata(hObject, handles);
% %             hgfeval({callback_hist2d, hObject,0});


% --- Executes on selection change in cos_button.
function cos_button_Callback(hObject, ~, handles)
% hObject    handle to cos_button (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cos_button contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cos_button

handles.crop = handles.crop_selec_button.Value + 1;
handles.short_version = 0;

[ handles.x, handles.y, handles.complete, handles.test, handles.model, handles.r2 ] = ...
    analyse_cos_ISHG( handles.contr, handles.img_3D, handles.resx, handles.resy, handles.rect, ...
    handles.crop, handles.xTitle_dflt, handles.yTitle_dflt, handles.screensize, handles.fact, ...
    handles.left_offset_fig, handles.top_offset_fig, handles.axes_font_size, handles.xaxis_sz,...
    handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz, handles.short_version, handles.x_phase, handles.double_contrast, handles.contr_mode, handles.cmap_redgreen);

handles.y_cell = {handles.y}; handles.x_cell = {handles.x};

[ handles.Titre1, handles.Titre2, handles.Titre3, handles.Titre4, handles.Titre5, handles.Titre6, handles.Counts, handles.Yaxis1, handles.Yaxis2, handles.Legendehisto ] = string_axis_ISHG( handles.langue );
handles.Titre4_modified = handles.Titre4;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cos_button_CreateFcn(hObject, ~, handles)
% hObject    handle to cos_button (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in histmode_menu.
function histmode_menu_Callback(hObject, ~, handles)
% hObject    handle to histmode_menu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns histmode_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from histmode_menu

state = get(hObject,'Value');

switch state
    case 1 % no choice position
        set(handles.load_stck_button, 'Enable', 'off');
        handles.shg = 0;
    case 2 % SHG intensity
        if ~handles.file_not_good
            set(handles.load_stck_button, 'Enable', 'on');
        end
%         save_img_chck = 1;
        handles.shg = 1;
    case 3 % Interferometric contrast
        if ~handles.file_not_good
            set(handles.load_stck_button, 'Enable', 'on');
        end
        handles.shg = 0;
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function histmode_menu_CreateFcn(hObject, ~, handles)
% hObject    handle to histmode_menu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in crop_selec_button.
function crop_selec_button_Callback(hObject, ~, handles)
% hObject    handle to crop_selec_button (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of crop_selec_button

% --- Executes on button press in filter_button.
function filter_button_Callback(hObject, ~, handles)
% hObject    handle to filter_button (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filter_button

handles.filt = get(hObject,'Value');

guidata(hObject, handles);

function nb_exclu_edt_Callback(hObject, ~, handles)
% hObject    handle to nb_exclu_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nb_exclu_edt as text
%        str2double(get(hObject,'String')) returns contents of nb_exclu_edt as a double

handles.nbr_objects = str2double(get(hObject, 'String'));

if handles.nbr_objects >=1
    handles.obj_test = 1;
else
    handles.obj_test = 0;
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nb_exclu_edt_CreateFcn(hObject, ~, handles)
% hObject    handle to nb_exclu_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in language_menu.
function language_menu_Callback(hObject, ~, handles)
% hObject    handle to language_menu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns language_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from language_menu


% --- Executes during object creation, after setting all properties.
function language_menu_CreateFcn(hObject, ~, handles)
% hObject    handle to language_menu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function res_x_edt_Callback(hObject, ~, handles)
% hObject    handle to res_x_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of res_x_edt as text
%        str2double(get(hObject,'String')) returns contents of res_x_edt as a double

handles.resx = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function res_x_edt_CreateFcn(hObject, ~, handles)
% hObject    handle to res_x_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function res_y_edt_Callback(hObject, ~, handles)
% hObject    handle to res_y_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of res_y_edt as text
%        str2double(get(hObject,'String')) returns contents of res_y_edt as a double


handles.resy = str2double(get(hObject,'String'));

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function res_y_edt_CreateFcn(hObject, ~, handles)
% hObject    handle to res_y_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slice_per_step_edt_Callback(hObject, ~, handles)
% hObject    handle to slice_per_step_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slice_per_step_edt as text
%        str2double(get(hObject,'String')) returns contents of slice_per_step_edt as a double

handles.start_phase = str2double(get(hObject,'String'));

set(handles.max_ps_indic_edt, 'String', num2str(str2double(handles.first_phase_edt.String) + 360/2*max(2, str2double(handles.slice_per_step_edt.String)) - str2double(handles.step_ps_edt.String)));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slice_per_step_edt_CreateFcn(hObject, ~, handles)
% hObject    handle to slice_per_step_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function two_phase_edt_Callback(hObject, ~, handles)
% hObject    handle to two_phase_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of two_phase_edt as text
%        str2double(get(hObject,'String')) returns contents of two_phase_edt as a double

handles.step_phase = str2double(get(hObject,'String'));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function two_phase_edt_CreateFcn(hObject, ~, handles)
% hObject    handle to two_phase_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function avg_img_edt_Callback(hObject, ~, handles)
% hObject    handle to avg_img_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of avg_img_edt as text
%        str2double(get(hObject,'String')) returns contents of avg_img_edt as a double
avg = str2double(get(hObject,'String'));
stp = str2double(get(handles.step_ps_edt, 'String'));

if avg > 1
    set(handles.sel_only_few_fr_per_ps, 'Enable', 'on');
else
    set(handles.sel_only_few_fr_per_ps, 'Value', 0);
    set(handles.sel_only_few_fr_per_ps, 'Enable', 'off');
end

new_stp = stp*avg/handles.avg_img;
handles.avg_img = avg;

set(handles.step_ps_edt, 'String', num2str(new_stp));
set(handles.snd_step_edt, 'String', num2str(new_stp));

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function avg_img_edt_CreateFcn(hObject, ~, handles)
% hObject    handle to avg_img_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.avg_img = 1;

guidata(hObject, handles);

% --- Executes on button press in upd_graph_button.
function upd_graph_button_Callback(hObject, ~, handles)
% hObject    handle to upd_graph_button (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in deflt_param_chck.
function deflt_param_chck_Callback(hObject, ~, handles)
% hObject    handle to deflt_param_chck (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deflt_param_chck

dflt_chck = get(hObject,'Value');

switch dflt_chck
    case 1 % default config is checked
        
        set(handles.language_menu, 'Value', 1);
        set(handles.histmode_menu, 'Value', 3);
        set(handles.load_stck_button, 'Enable', 'on');
        handles.shg = 0;
        set(handles.avg_img_edt, 'String', '1');
        set(handles.crop_selec_button, 'Value', 0);
        set(handles.filter_button, 'Value', 0);
        set(handles.crop_selec_button, 'Value', 0);
        set(handles.nb_exclu_edt, 'String', '0');
        set(handles.contrast_edt, 'String', num2str(handles.contrast_dflt));
        set(handles.nbinsX_edt, 'String', num2str(handles.nbinsx_dflt));
        set(handles.nbinsY_edt, 'String', num2str(handles.nbinsy_dflt));
        set(handles.weight_filtr_chck, 'Value', 0);
        %         set(handles.center_coeff_edt, 'Visible', 'off');
        %         set(handles.center_coeff_txt, 'Visible', 'off');
        %         set(handles.cross_coeff_edt, 'Visible', 'off');
        %         set(handles.cross_coeff_txt, 'Visible', 'off');
        %         set(handles.diag_coeff_edt, 'Visible', 'off');
        %         set(handles.diag_coeff_txt, 'Visible', 'off');
        %         set(handles.center_coeff_edt, 'String', '1');
        %         set(handles.cross_coeff_edt, 'String', '0');
        %         set(handles.diag_coeff_edt, 'String', '0');
        
    case 0 % default config not checked
end

guidata(hObject, handles);

function handles = int_contrast_util(hObject, handles)

handles.crop = handles.crop_selec_button.Value;
handles.filt = handles.filter_button.Value;

handles.contr_mode = get(handles.contr_radio, 'Value');

handles.x_phase = handles.x_phase0;

if handles.use_range % str2double(get( handles.useonly_end_edt, 'String'))
    handles.use_range(1) =  str2double(get( handles.useonly_strt_edt, 'String'));
    handles.use_range(2) =  str2double(get( handles.useonly_stp_edt, 'String'));
    handles.use_range(3) =  str2double(get( handles.useonly_end_edt, 'String'));
    
    if handles.use_range(1) <= 0 % array begin at 0
        handles.use_range(1) =1;
    end
    if handles.use_range(2) < 1 % step > = 1
        handles.use_range(2) =1;
    end
    if handles.use_range(3) < handles.use_range(1) % 
        handles.use_range(3) = handles.use_range(1);
    end
% % else
% %     handles.use_range = 0;
end

%% Soustraction des images GSH interférométriques pour obtenir les images du contraste interférométrique

handles.avg_median = handles.med_contr_chck.Value;

if (isa(handles.fname, 'cell') && get(handles.batch_chck, 'Value')) % batch files
    skip_asking_coeff_median = 1;
    handles.crop = 0; % crop is performed in image reading in batch
    set(handles.crop_selec_button, 'Value', 0);
else
    skip_asking_coeff_median = 0;
end

% % if handles.slice_per_ps_step == 1
% %     msgbox('If you want increasing order with no contrast mode, you should choose 0 slices per step instead and re-load your stack')
% % end
min_val = str2double(get(handles.val_min_edt, 'String'));
max_val = str2double(get(handles.val_maxsat_edt, 'String'));
handles.slice_per_ps_step =str2double(handles.slice_per_step_edt.String); % keep this one for auto !!
handles.double_contrast = max(0, floor(log(str2double(handles.double_contr_chck.String(end)))/log(2)-1));

[ handles.x, handles.y, handles.contr00, handles.rect, handles.img_shg, handles.x_phase, handles.img_3D , handles.mask_val_sat, handles.order_frame_vect ] = ...
    contrast_by_subtract_img_ISHG( handles.img_3D, handles.num_images, handles.filt, ...
    handles.crop, handles.resx, handles.resy, handles.shg, handles.img_shg, handles.undocked_fig1, handles.plot_princ,...
    handles.screensize, handles.fact, handles.left_offset_fig, handles.top_offset_fig, handles.slice_per_ps_step, ...
    handles.contr_mode, handles.diff_phase, handles.x_phase, handles.use_range, handles.avg_median, ...
    handles.double_contrast, skip_asking_coeff_median, handles.force_incr_order_chck.Value, ...
    handles.disp_vect00, get(handles.norm_interframes_ctr_chk, 'Value'), [min_val, max_val]);
% function to calculate contrast by subtraction of images
handles.img_shg_max = max(handles.img_shg(:));
handles.y_cell = {handles.y}; handles.x_cell = {handles.x};

if length(handles.rect) ~= 1
    
    fprintf('%f, %f, %fx%f \n', handles.rect(1), handles.rect(2), handles.rect(3), handles.rect(4));% [xmin ymin width height]
end

%% Exclusions de certaines portions de l'image pour le calcul

handles.nbr_objects = str2double(handles.nb_exclu_edt.String);

if handles.nbr_objects == 0
    handles.obj_test = 0;
else
    handles.obj_test = 1;
end

[ handles.contr00, handles.img_shg ] = exclude_area_ISHG( handles.contr00, handles.obj_test, handles.crop, handles.shg, handles.img_3D, handles.rect, handles.img_shg, ...
    handles.nbr_objects, 1, handles.plot_princ, handles.screensize, handles.fact, handles.left_offset_fig, handles.top_offset_fig  );
% % function to exclude the zones

%% Réorganisation des images du contraste interférométrique dans un ordre croissant + plot

% Valeur max (et min) de la colorbar pour l'affichage des images de
% soustraction et du stack qui sert à faire le fit pour trouver la phase
handles.contrast = str2double(handles.contrast_edt.String);
handles.list_contr_name = [];

if get(handles.disp_interm_plots_chck, 'Value') % disp interm. plots
    if handles.contr_mode == 1 % contrast mode

        [handles.h, handles.list_contr_name] = diff_phase_stack_ISHG( handles.contr00, handles.x, handles.y, handles.contrast, handles.xTitle_dflt, handles.yTitle_dflt, ...
            handles.screensize, handles.fact, handles.left_offset_fig, handles.top_offset_fig, handles.undocked_fig1, handles.h_princ, handles.axes_font_size, ...
            handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz, handles.cmap_redgreen, 1, size(handles.contr00, 3), handles.x_phase );
        % function to calculate and plot the difference for the different angles

    else
        htime = msgbox('Images re-ordered in increasing phase-shift');
        handles.list_contr_name = textscan(num2str(handles.x_phase),'%s');
        handles.list_contr_name = handles.list_contr_name{1}';
    end
end
handles.contr =  handles.contr00;

set(handles.saving_popup, 'Enable', 'on');
set(handles.cos_button, 'Enable', 'on');
set(handles.relative_phase_button, 'Enable', 'on');
set(handles.analyze_panel, 'Visible', 'on');
set(handles.plot_diff_menu, 'Visible', 'on');
% aa = handles.plot_diff_menu.String;

% if length(aa) < 3
aa =  [{'View Plot ...'; '1st image of stack'}; handles.list_contr_name'];
set(handles.plot_diff_menu, 'String', aa);
set(handles.plot_diff_menu, 'Value', length(aa));
% end

if handles.undocked_fig1
    handles.h_princ_out = handles.h_princ;
end

set(handles.plot_phase_popup, 'String', {'View plot phase ...'});
set(handles.plot_phase_popup, 'Value', 1);

if ((isa(handles.fname, 'cell') && get(handles.batch_chck, 'Value')) && handles.contr_mode == 0) % a number of file > 1 has been chosen
    delete(htime );
end

% % handles.cmap_princ_curr = handles.cmap_redgreen;

% --- Executes on button press in int_contrast_button.
function int_contrast_button_Callback(hObject, ~, handles)
% hObject    handle to int_contrast_button (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 handles = int_contrast_util(hObject, handles);
 
guidata(hObject, handles);

function handles = ld_stck3d_util(hObject, handles, clickType)

cd(handles.folder_name);

num_images= 0;
if isa(handles.fname, 'cell')
    if (get(handles.batch_chck, 'Value') &&  ~get(handles.ind_files_chck, 'Value')) % many stacks ,  batch
        fname = handles.fname{handles.count_batch};
    else % many files not batch
        fname = handles.fname{1};
        num_images = length(handles.fname);
    end
    
else % one stck file
    fname = handles.fname;
end

warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('off','MATLAB:imagesci:tifftagsread:badTagValueDivisionByZero');
InfoImage=imfinfo(fname); % disp a warning 'Division by zero when processing YResolution'
warning( 'on','MATLAB:imagesci:tiffmexutils:libtiffWarning');
warning('on','MATLAB:imagesci:tifftagsread:badTagValueDivisionByZero');
% warning('ON', 'all');
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
if num_images == 0
    num_images = length(InfoImage);
end

if ((isa(handles.fname, 'cell') && ~get(handles.batch_chck, 'Value')) ...
|| (get(handles.batch_chck, 'Value') &&  get(handles.ind_files_chck, 'Value'))) % many single files not batch
    fname = handles.fname;
end

% set(handles.step_ps_edt, 'String', num2str(360/num_images*str2double(handles.slice_per_step_edt.String)/2));
% set(handles.snd_step_edt, 'String', num2str(360/num_images*str2double(handles.slice_per_step_edt.String)/2));
if sum(strcmp(clickType, {'alt', 'open'}))% right-click
    %rightclick
else
    cla(handles.plot_phase, 'reset');
    set(handles.plot_phase, 'Visible', 'off');
    set(handles.plot_princ, 'Visible', 'on');
end

handles.resx = str2double(handles.res_x_edt.String);
handles.resy = str2double(handles.res_y_edt.String);
handles.init_phase = str2double(handles.first_phase_edt.String);
% % handles.second_phase = str2double(handles.two_phase_edt.String);
% % handles.third_phase = str2double(handles.three_phase_edt.String);
handles.first_stp = str2double(handles.step_ps_edt.String);
handles.snd_step = str2double(handles.snd_step_edt.String);
handles.slice_per_ps_step = str2double(handles.slice_per_step_edt.String);

handles.near_avg_ini = get(handles.near_chck, 'Value');
handles.norm_ini = get(handles.norm_chck, 'Value');

handles.contr_mode = get(handles.contr_radio, 'Value');

get(handles.specify_order, 'Value');

if (isa(handles.fname, 'cell') && get(handles.batch_chck, 'Value')) % batch files
    skip_asking_coeff_median = 1;
else
    skip_asking_coeff_median = 0;
end

scl_X = str2double(get(handles.scaling_X_edt, 'String'));
scl_Y = str2double(get(handles.scaling_Y_edt, 'String'));
if (scl_X <= 1 && scl_Y <= 1); scaling = [];
else; scaling = [scl_X, scl_Y];
end

handles.avg_img = str2double(get(handles.avg_img_edt,'String'));

[handles.xTitle_dflt, handles.yTitle_dflt, handles.phi_mat_default, handles.cmap_brbgsat, ...
    handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz, ...
    handles.cmap_cubehelix, handles.cmap_redgreen, handles.cmap_blkredblkgrn, handles.cmap_hsvsat] = load_axis_param; 

if ~isvalid(handles.h_princ)
    handles.h_princ = handles.plot_phase;
end

func_hdl = load_stack_plot_ISHG;
[ handles.img_3D, handles.type_im, handles.num_images, handles.img_shg, ...
    handles.h_princ, handles.slice_per_ps_step, handles.start_phase, handles.diff_phase, ...
    handles.x_phase, handles.nonequal_spacing, x_center, y_center, handles.disp_vect00, contr_mode ] = ...
    func_hdl.load_stack_plot_ISHG_func( fname, handles.shg, handles.undocked_fig1, handles.h_princ, ...
    handles.init_phase, handles.first_stp, handles.snd_step, handles.slice_per_ps_step, ...
    handles.contr_mode, handles.avg_img, handles.near_avg_ini, handles.norm_ini, skip_asking_coeff_median, get(handles.point_ROI_batch_chck,'Value'),...
    nImage, mImage, num_images, handles.force_incr_order_chck.Value, handles.axes_font_size, ...
    handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz, handles.screensize,...
    handles.fact, handles.left_offset_fig, handles.top_offset_fig, get(handles.sel_only_few_fr_per_ps, 'Value'), ...
    get(handles.invphshft_order_chk, 'Value'), get(handles.norm_interframes_chk, 'Value') , scaling);

% % handles.cmap_princ_curr = handles.cmap_redgreen;
if length(contr_mode)>1
    new_steperps = contr_mode(2);
    set(handles.slice_per_step_edt, 'String',num2str(new_steperps));
    contr_mode=contr_mode(1);
    handles.slice_per_ps_step=new_steperps;
end
if ~contr_mode
    set(handles.raw_radio, 'Value', 1);
end
handles.crop_batch = get(handles.crop_chck_box, 'Value');
if handles.crop_batch
    [handles.img_3D, ~, ~, ~] = plot_crop_util(handles, handles.img_3D, x_center, y_center);
end

if ~sum(strcmp(clickType, {'alt', 'open'}))% left-click
    set(handles.int_contrast_button, 'Enable', 'on');

    handles.x_phase0 = handles.x_phase;

    handles.langue = get(handles.language_menu, 'Value');

    [ handles.Titre1, handles.Titre2, handles.Titre3, handles.Titre4, handles.Titre5, handles.Titre6, handles.Counts, handles.Yaxis1, handles.Yaxis2, handles.Legendehisto ] = string_axis_ISHG( handles.langue );

    set(handles.interf_pannel, 'Visible', 'on');

    aa =  {'View plot ...', '1st image of stack'};
    set(handles.plot_diff_menu, 'String', aa);
    set(handles.plot_diff_menu, 'Value', 2);

    if handles.undocked_fig1
        handles.h_princ_out = handles.h_princ;
    end

    if (handles.slice_per_ps_step > 1 && handles.contr_radio.Value == 1)
        set(handles.double_contr_chck, 'Enable', 'on');
    else
        set(handles.double_contr_chck, 'Enable', 'off');
    end

    handles.phase_cell = {}; % reset
    handles.amp_cell = {}; % reset
    handles.err_cell = {}; % reset
    handles.x_cell = {}; % reset
    handles.y_cell = {}; % reset
    % disp(size(handles.img_3D, 3))
    % disp(str2double(handles.slice_per_step_edt.String))

    set(handles.saving_popup, 'Enable', 'on');
    set(handles.action_map_menu, 'Visible', 'on');

    handles.use_range = 0;

    nb =str2double(get(handles.useonly_end_edt, 'String'));
    if ((~handles.range_interf_keep || nb <= 0 ) || nb > length(handles.x_phase))
        set(handles.useonly_end_edt, 'String', length(handles.x_phase)); % if contrast mode or not, the number of 'contr' frames should be the same if the user took the right number
    end
    set(handles.useonly_stp_edt, 'String', '1');
    set(handles.useonly_strt_edt, 'String', '1');
    handles.order_frame_vect = 0; % reinit
    % % handles.cmap_princ_curr =  gray;
end

function [im, ph, vy, vx] = plot_crop_util(handles, im, x_center, y_center)

state = get(handles.plot_princ, 'Visible');
switch state
    case 'off' % phase
        hh= handles.plot_phase;
        cmp = handles.cmap_default1;
        ph = 1;
    case 'on' % plot left
        hh= handles.plot_princ;
        cmp = handles.cmap_princ_curr ;
        ph = 0;
end
if isempty(im)
    a=get(handles.plot_phase, 'Children'); im=a.CData; 
end
xmin = str2double(handles.xmin_crop_edit.String); ymin = str2double(handles.ymin_crop_edit.String);
width = str2double(handles.width_crop_edit.String); height =str2double(handles.height_crop_edit.String);
% % 36.478836, 52.087302, 122.751323x122.751323
% % mask=logical(im*0);
if ~get(handles.point_ROI_batch_chck,'Value')
    vy = max(1, round(ymin)): min(size(im, 1), round(ymin + height - 1));
    vx = max(1,round(xmin)): min(size(im, 2), round(xmin + width - 1));
else % point center at every frame
    vy = max(1, y_center - round(height/2)): min(size(im, 1), y_center + round(height/2));
    vx =  max(1, x_center - round(width/2)): min(size(im, 2), x_center + round(width/2));
end
im = im(vy, vx, :);
imagesc(hh, im(:,:,1)); axis image; colormap(cmp); colorbar(hh);

% --- Executes on button press in load_stck_button.
function load_stck_button_Callback(hObject, ~, handles)
% hObject    handle to load_stck_button (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clickType = get(ancestor(hObject, 'figure'), 'SelectionType'); % normal
handles = ld_stck3d_util(hObject, handles,clickType);

guidata(hObject, handles);


function handles = rel_phase_util(hObject, handles)

cla(handles.plot_princ, 'reset');
set(handles.plot_princ, 'Visible', 'off');
set(handles.plot_phase, 'Visible', 'on');
handles.weight = get(handles.weight_filtr_chck, 'Value');

handles.use_invA = 1-get(handles.notuse_A_chck, 'Value');

handles.method_calc_phase = get(handles.three_ph_algo_chck, 'Value') - 2;

if handles.method_calc_phase < 0
    handles.method_calc_phase = 0;
    set(handles.three_ph_algo_chck, 'Value', 2);
end

handles.parallel_comput = get(handles.parallel_comp_chck, 'Value');
handles.decomp_LU = get(handles.LU_inv_chck, 'Value');
%     handles.method_calc_phase = 3; % 3 phases algo
% else
%     if get(handles.vect_chck, 'Value'); % on vector
%         % 2 for CA, 1 for Rivard, 4 for Rivard with 1D matrices, 3 for algo 3 phases
%         handles.method_calc_phase = 4; % 3 phases algo
%     else
%         handles.method_calc_phase = 1; % Rivard
%     end
% end

% if handles.weight
%     handles.coef_centre = str2double(get(handles.center_coeff_edt, 'String'));
%     handles.coef_cross = str2double(get(handles.cross_coeff_edt, 'String'));
%     handles.coef_diag = str2double(get(handles.diag_coeff_edt, 'String'));
% else
%     handles.coef_centre = 1;
%     handles.coef_cross = 0;
%     handles.coef_diag = 0;
% end

% handles.cmap_default1 = handles.cmap_blkredblkgrn;
if (isa(handles.fname, 'cell') && get(handles.batch_chck, 'Value')) % batch files
    skip_asking_coeff_median = 1;
    handles.crop = 0; % crop is performed in image reading in batch
    set(handles.crop_selec_button, 'Value', 0);
else
    skip_asking_coeff_median = 0;
end
handles.cmap_default1 = hsv;

k_max = str2double(get(handles.it_max_algo_edt, 'String'));
epsilon = str2double(get(handles.epsilon_ph_th_edt, 'String')); % 1e-4 in article
epsilon1 = str2double(get(handles.epsilon_tilt_th_edt, 'String')); % rad, 1e-2; % 1e-5 in article
epsilon_tilt_th_edt_Callback(hObject, 0, handles);
epsilon_ph_th_edt_Callback(hObject, 0, handles);
it_max_algo_edt_Callback(hObject, 0, handles);

ramp_phshft = get(handles.ramps_phshft_flag, 'Value');
reit_prev_phase = get(handles.reuse_curr_ph_chck, 'Value');
prev_res ={};
if reit_prev_phase
    num = get(handles.plot_phase_popup, 'Value')-1;
    prev_res = {handles.phase_cell{num}, handles.amp_cell{num}, handles.err_cell{num}, handles.divers_vect};
end

[ handles.phase1, handles.amp1, handles.err, handles.h_phase, handles.divers_vect] = phase_map_ISHG( handles.contr, handles.x, handles.y,...
    handles.xTitle_dflt, handles.yTitle_dflt, handles.phi_mat_default, handles.Titre1, handles.Titre2, ...
    handles.Titre3, handles.obj_test, handles.weight, handles.coef_centre, handles.coef_cross, handles.coef_diag,...
    handles.screensize, handles.fact, handles.left_offset_fig, handles.top_offset_fig, ...
    handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz, ...
    handles.undocked_phase, handles.h_phase, handles.cmap_default1, handles.cmap_brbgsat, handles.x_phase,...
    handles.use_invA, handles.method_calc_phase, handles.mean_phase_meth, handles.sigma_med,...
    handles.parallel_comput, handles.decomp_LU, handles.double_contrast, handles.contr_mode, ...
    skip_asking_coeff_median, k_max, epsilon, epsilon1, reit_prev_phase, prev_res, ramp_phshft, handles.mask_val_sat, get(handles.avg_on_ph_chck, 'Value') );
% On calcule la phase pour chaque pixel - voir fit_I_SHG.m

handles.max_amp00 = round(max(max(handles.amp1)));
handles.phase_cell{end+1} = handles.phase1;
handles.phase_cell00{end+1} = handles.phase1;

handles.x_cell{end+1} = handles.x_cell{end};
handles.y_cell{end+1} = handles.y_cell{end};

handles.err_cell{end+1} = handles.err;

set(handles.action_map_menu, 'Visible', 'on');
set(handles.expboth2miji_chck, 'Visible', 'on');

if handles.shg == 1 % shg intensity
% %     handles.amp_cell(end+1) = {handles.img_shg};
    handles.title_hist1 = handles.Titre5;
    handles.yaxis_hist1 = handles.Yaxis1;
else % interf contrast
    
    handles.title_hist1 = handles.Titre6;
    handles.yaxis_hist1 = handles.Yaxis2;
end
handles.amp_cell{end+1} = handles.amp1;
if handles.amp1 ~= 0 % shg img loaded
    str1=get(handles.two_d_hist_popup, 'String');
    if (~strcmp(str1{end} , handles.str_hist_ctr_shg) && ~strcmp(str1{end-1} , handles.str_hist_ctr_shg))
        str1{end+1}= handles.str_hist_ctr_shg;
        set(handles.two_d_hist_popup, 'String', str1);
    end
end

set(handles.cmap_p_menu, 'Value',5); % HSV
if ~isa(handles.amp_cell00, 'cell'); handles.amp_cell00 = {handles.amp_cell00}; end
handles.amp_cell00(end+1) = handles.amp_cell(end);

% colorbar(handles.plot_princ,'off')
% aa = handles.plot_phase_popup.String;
% aa = [aa, 'Phase map total'];
% set(handles.plot_phase_popup, 'String', aa);
% set(handles.plot_phase_popup, 'Value', 2);

set(handles.hist_button, 'Enable', 'on');
set(handles.two_d_hist_popup, 'Enable', 'on');
set(handles.nbinsX_edt, 'Enable', 'on');
set(handles.nbinsY_edt, 'Enable', 'on');

handles.init_hist_ok = 0; % flag for init success

str = get(handles.plot_phase_popup, 'String' );

if isa(str, 'char')
    str = {str, 'Phase map total'};
    set(handles.plot_phase_popup, 'String', str)
else
    %     if (~strcmp(str{end}, 'Phase map total'))
    str{end+1} = 'Phase map total';
    set(handles.plot_phase_popup, 'String', str)
    %     end
end

set(handles.plot_phase_popup,'Value', length(str) );

handles.Titre4_modified = handles.Titre4;

aa = handles.plot_diff_menu.String;

aa = [aa; {'Interf. contrast map'; 'Rel. Error'}];
set(handles.plot_diff_menu, 'String', aa) ;

% --- Executes on button press in relative_phase_button.
function relative_phase_button_Callback(hObject, ~, handles)
% hObject    handle to relative_phase_button (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = rel_phase_util(hObject, handles);

guidata(hObject, handles);

% --- Executes on button press in hist_button.
function hist_button_Callback(hObject, eventdata, handles)
% hObject    handle to hist_button (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.int_x = 1./str2double(handles.nbinsX_edt.String); % int_x is inverse of the nbins

handles.int_y = str2double(handles.nbinsY_edt.String);%% Histogramme 3D

num_cell = get(handles.plot_phase_popup, 'Value')-1;
if (handles.shg && length(handles.img_shg)>1);amp = handles.img_shg;
else;amp = handles.amp_cell{num_cell};
end
[ handles.data_tot, handles.ylabel_hist3, handles.title_hist3, handles.nbinsx, handles.h_hist ] = hist_3D_ISHG(handles.phase_cell{num_cell}, amp, handles.int_x, handles.int_y, ...
    handles.title_hist1, handles.yaxis_hist1, handles.Counts, handles.phi_mat_default,...
    handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz, handles.undocked_fig2, ...
    handles.h_hist, handles.screensize, handles.fact, handles.left_offset_fig, handles.top_offset_fig, handles.offset_pi2, handles.range_sat );
% Histogram 3D function

handles.max_amp = round(max(max(handles.data_tot(:,2))));
set(handles.max_hist_edt, 'String', num2str(handles.max_amp));
min_amp = round(min(min(handles.data_tot(:,2))));
set(handles.min_hist_edt, 'String', num2str(min_amp));

% % fprintf('Maximum of the amplitude calculated = %g \n', handles.max_amp);

if handles.undocked_fig2
    handles.h_hist_out = handles.h_hist;
else
    handles.h_hist_out = 1;
end

% % set(handles. ,'Enable', 'on');

guidata(hObject, handles);

function handles = two_d_hist(hObject, handles)
% % not a callback

handles.int_x = 1./str2double(handles.nbinsX_edt.String); % int_x is inverse of the nbins

%% hist 2D
% iii=handles.plot_diff_menu.Value-1;
num_cell = get(handles.plot_phase_popup, 'Value')-1;
[handles.hhist1, handles.hist2_xdata, handles.hist2_ydata] = hist2D_ISHG(  handles.fact,  handles.left_offset_fig, ...
    handles.top_offset_fig, handles.phase_cell{num_cell},  handles.int_x,  handles.phi_mat_default,  handles.Counts,  handles.Titre4_modified, ...
    handles.axes_font_size,  handles.xaxis_sz,  handles.yaxis_sz,  handles.title_sz,  handles.clrbr_tl_sz, handles.screensize, handles.offset_pi2, handles.hhist1 );

set(handles.fit_hist_button, 'Visible', 'on');
set(handles.baseline_hist_edt, 'Visible', 'on'); set(handles.baseline_hist_lbl, 'Visible', 'on');
set(handles.edt_pearson, 'Visible', 'on');
set(handles.text_pearson, 'Visible', 'on');

% WARNING : you can't use the same var to store some data because it is
% reinitialised each time
if (get(handles.batch_chck, 'Value') && (handles.fit_batch || handles.fit_clrwheelBRBG)) % a number of file > 1 has been chosen
        
    handles.fit_chosen=8; handles.M =0.5;
    func_hdl = fit_hist2d_funcs;
    [ handles.stringb1, handles.stringc1, handles.stringb2, handles.stringc2, handles.cmap_brbg_fit, handles.vect_fit_struct, fit_warning ] = ...
        func_hdl.fit_hist2d_ishg( handles.hist2_xdata, handles.hist2_ydata, handles.hhist1, handles.Legendehisto, handles.offset, handles.fit_chosen, handles.M, ...
    str2double(get(handles.norm_fact_hist, 'String')),  str2double(get(handles.baseline_hist_edt, 'String')));
    
    if handles.fit_batch
        name = 'cell_hist_rename.mat';

        [~, I]  = max(handles.hist2_ydata);
        handles.max_c = handles.hist2_xdata(I);
        %     handles.fID = fopen('max_hist_batch_rename-me.txt', 'a');% %s.txt', handles.fname{1}), 'a');
        %     fprintf(handles.fID, sprintf('%f; \n', handles.max_c));
        %     fclose(handles.fID);
        %     fprintf('%f; \n', handles.max_c);
        %     figure(1); hold on;

        if handles.count_batch == 1
            handles.cell_hist_save{1, 1} = 'name';
            handles.cell_hist_save{1, 2} = 'hist_vect';
            handles.cell_hist_save{1, 3} = 'vect_fit';
            handles.cell_hist_save{1, 4} = 'max_hist';
            handles.cell_hist_save{1, 5} = 'center fit';
            handles.cell_hist_save{1, 6} = 'width fit';
            handles.cell_hist_save{1, 7} = 'good fit ?';
            handles.cell_hist_save{1, 8} = 'Phase map';
            handles.cell_hist_save{1, 9} = 'Avg intensity';

        else
            handles.cell_hist_save=load(name);
            handles.cell_hist_save = handles.cell_hist_save.cell_hist_save; % because a struct is loaded
        end
        handles.cell_hist_save{handles.count_batch+1, 1} = handles.fname{handles.count_batch};
        handles.cell_hist_save{handles.count_batch+1, 2} = handles.hist2_ydata;
        handles.cell_hist_save{handles.count_batch+1, 3} = handles.vect_fit_struct;
        handles.cell_hist_save{handles.count_batch+1, 4} = handles.max_c;
        handles.cell_hist_save{handles.count_batch+1, 5} = handles.stringb1;
        handles.cell_hist_save{handles.count_batch+1, 6} = handles.stringc1;
        handles.cell_hist_save{handles.count_batch+1, 7} = fit_warning;
        handles.cell_hist_save{handles.count_batch+1, 8} = handles.phase1;
        handles.cell_hist_save{handles.count_batch+1, 9} = mean(mean(mean(handles.img_3D)));
        cell_hist_save = handles.cell_hist_save;  %#ok<NASGU> % saved
        save(name, 'cell_hist_save');

        %     if handles.count_batch == length(handles.fname) % the end
        %
        %         cell_hist_save = handles.cell_hist_save;
        %         save(sprintf('cell_hist_%d_rename.mat', handles.count_batch), 'cell_hist_save');
        %     end
    end
end

set_histrange_vis(handles);

guidata(hObject, handles);

% --- Executes on button press in hist_range_button.
function hist_range_button_Callback(hObject, eventdata, handles)
% hObject    handle to hist_range_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.pb = str2double(handles.min_hist_edt.String); % min lambda_temp
handles.ph = str2double(handles.max_hist_edt.String); % max lambda_temp

if (handles.ph == handles.pb || handles.ph == 0) % if range not set
% %     prompt = {'Valeur minimale de l''amplitude:','Valeur maximale de l''amplitude:'};
% %     dlg_title = 'Input';
% %     num_lines = 1;
% %     answer = inputdlg(prompt,dlg_title,num_lines);
    fhandle = choosedialog; % compulsory
    if handles.shg; str = 'shg amp'; else; str = 'ctr amp'; end
    [pb, ph, keep_cond, offset_pi2] = fhandle.choosedialog_range(str, handles.ph);
% %     handles.pb = str2double(answer{1});
% %     handles.ph = str2double(answer{2});
    if ~keep_cond
        handles.pb = pb;
        handles.ph = ph;
        set(handles.min_hist_edt, 'String', num2str(handles.pb));
        set(handles.max_hist_edt, 'String', num2str(handles.ph));
        if offset_pi2~=handles.offset_pi2
            handles.offset_pi2 = offset_pi2;
            hist_button_Callback(hObject, eventdata, handles)
        end
    end
end

bx = str2double(handles.nbinsX_edt.String);

if bx > 100
    bx=100;
    disp('Nbins on phase too high for a disp of 3D hist ! Set to 100');
end

handles.int_x = 1./bx; % int_x is inverse of the nbins
handles.int_y = str2double(handles.nbinsY_edt.String);
if ~isa(handles.amp_cell00, 'cell'); handles.amp_cell00 = {handles.amp_cell00}; end

[phase_test, img, handles.h_phase, handles.range_sat] = map_spec_region_ISHG( handles.x_cell{end}, handles.y_cell{end},...
 handles.phase_cell00{end}, handles.amp_cell00{end}, ...
    handles.xTitle_dflt, handles.yTitle_dflt, handles.Titre1, handles.phi_mat_default, handles.pb, handles.ph,...
    handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz, ...
    handles.undocked_phase, handles.h_phase, handles.screensize, handles.fact, handles.left_offset_fig, ...
    handles.top_offset_fig, handles.cmap_sat_dflt);
% To calculate and plot the histogram 2D and 3D for a specific range of
% intensity

% cc = 2;
% while 1 % all the time true
%     if isfield(handles, sprintf('phase%d', cc)) % phasecc already exists
%         cc = cc+1;
%     else
%         handles.phase1 = phase_test;
%     end
% end
% %     handles.init_hist_ok = 1; % flag for init success

set(handles.plot_phase, 'Visible', 'on');
set(handles.plot_princ, 'Visible', 'off');
% colorbar(handles.plot_princ,'off')
set(handles.cmap_p_menu, 'Value', 3);

handles.cmap_default1 = handles.cmap_sat_dflt;

set(handles.sat_value_slider, 'Visible', 'on');
set(handles.sat_txt, 'Visible', 'on');
set(handles.sat_value_slider, 'Value', 63/64);

set(handles.reset_range_chck, 'Enable', 'on');

str = get(handles.plot_phase_popup, 'String' );
if ~strcmp(str{end}, 'Phase map sat')
    str{end+1} = 'Phase map sat';
    set(handles.plot_phase_popup, 'String', str);
else % there was a phase map sat just before
    handles.phase_cell(end) = [];
    handles.amp_cell(end) = [];
    handles.x_cell(end) = [];
    handles.y_cell(end) = [];
end

handles.phase_cell(end+1) = {phase_test};
handles.amp_cell(end+1) = {img};
handles.x_cell{end+1} = handles.x_cell{end};
handles.y_cell{end+1} = handles.y_cell{end};

handles.err_cell(end+1) = handles.err_cell(end);

set(handles.plot_phase_popup,'Value', length(str) );

handles.Titre4_modified = sprintf('%s in [%.1f, %.1f] interf. contrast', handles.Titre4, handles.pb, handles.ph);

if handles.undocked_fig1
    handles.h_phase_out = handles.h_phase;
    % else
    %     handles.h_phase_out = 1;
end

guidata(hObject, handles);

function min_hist_edt_Callback(hObject, ~, handles)
% hObject    handle to min_hist_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_hist_edt as text
%        str2double(get(hObject,'String')) returns contents of min_hist_edt as a double

% if str2double(get(hObject,'String')) > 0 % negative values are not good
%     set(handles.hist_range_button, 'Enable', 'on');
% else
%     set(handles.hist_range_button, 'Enable', 'off');
% end

set(handles.reset_range_chck, 'Value', 0);

if str2double(get(hObject,'String')) < 0
    
    set(hObject,'String', 0)
end

if str2double(get(hObject,'String')) > 0 % value above 0 are accepted
    set(handles.hist_range_button, 'Enable', 'on');
else
    set(handles.hist_range_button, 'Enable', 'off');
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function min_hist_edt_CreateFcn(hObject, ~, handles)
% hObject    handle to min_hist_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function max_hist_edt_Callback(hObject, ~, handles)
% hObject    handle to max_hist_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_hist_edt as text
%        str2double(get(hObject,'String')) returns contents of max_hist_edt as a double

set(handles.reset_range_chck, 'Value', 0);

if str2double(get(hObject,'String')) > 0 % value above 0 are accepted
    set(handles.hist_range_button, 'Enable', 'on');
else
    set(handles.hist_range_button, 'Enable', 'off');
end

if isfield(handles, 'max_amp00')
    if str2double(get(hObject,'String')) > handles.amp_cell{ get(handles.plot_phase_popup, 'Value')-1}% %handles.max_amp00

        set(hObject,'String', num2str(handles.max_amp00))
    end
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function max_hist_edt_CreateFcn(hObject, ~, handles)
% hObject    handle to max_hist_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in weight_filtr_chck.
function weight_filtr_chck_Callback(hObject, ~, handles)
% hObject    handle to weight_filtr_chck (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of weight_filtr_chck

avg_weight_filtr = get(hObject,'Value');

switch avg_weight_filtr
    case 1 % filter checked
        %         set(handles.center_coeff_edt, 'Visible', 'on');
        %         set(handles.center_coeff_txt, 'Visible', 'on');
        %         set(handles.cross_coeff_edt, 'Visible', 'on');
        %         set(handles.cross_coeff_txt, 'Visible', 'on');
        %         set(handles.diag_coeff_edt, 'Visible', 'on');
        %         set(handles.diag_coeff_txt, 'Visible', 'on');
        %         set(handles.center_coeff_edt, 'String', '5');
        %         set(handles.cross_coeff_edt, 'String', '2');
        %         set(handles.diag_coeff_edt, 'String', '1');
        prompt = {'Method (1- = mean, 2+ = median)', 'Center  coeff', 'Cross  coeff', 'Diag.  coeff'};
        dlg_title = 'Weight';
        num_lines = 1;
        % Valeurs par défaut
        def = {'2', '1', '1', '1'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        meth_avg = str2double(answer{1});
        handles.coef_centre = str2double(answer{2});
        handles.coef_cross = str2double(answer{3});
        handles.coef_diag = str2double(answer{4});
        
        if (handles.coef_cross == 0 && handles.coef_diag == 0)
            handles.mean_phase_meth = 0;
            handles.weight = 0;
            msgbox('You chose 0 value weight for adjacent pixels. The avg filter won`t be performed.');
            set(handles.no_fltr_radio, 'Value', 1);
            
        else % do a filter
             handles.weight = 1;
            
            
            if meth_avg <= 1 % mean
                handles.mean_phase_meth = 1;
            else % median
                handles.mean_phase_meth = 0; % median
                
                if (handles.coef_cross == handles.coef_diag == handles.coef_centre)
                    handles.sigma_med = 0; % not jointWMF, unweighted fast method
                else % weighted, slow
                    handles.sigma_med = -1; 
                end
            end
           
        end
    case 0 % filter NOT checked
        handles.mean_phase_meth = 0;
        %         set(handles.center_coeff_edt, 'Visible', 'off');
        %         set(handles.center_coeff_txt, 'Visible', 'off');
        %         set(handles.cross_coeff_edt, 'Visible', 'off');
        %         set(handles.cross_coeff_txt, 'Visible', 'off');
        %         set(handles.diag_coeff_edt, 'Visible', 'off');
        %         set(handles.diag_coeff_txt, 'Visible', 'off');
        %         set(handles.center_coeff_edt, 'String', '1');
        %         set(handles.cross_coeff_edt, 'String', '0');
        %         set(handles.diag_coeff_edt, 'String', '0');
end

guidata(hObject, handles);


% --- Executes on selection change in cmap1_menu.
function cmap1_menu_Callback(hObject, ~, handles)
% hObject    handle to cmap1_menu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cmap1_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cmap1_menu

% try
%     get(handles.h_princ, 'Children');
%     try
%         get(handles.h_princ,'CData')
%     catch
%         hh_fun0 = get(0, 'Children'); % current fig  not GUI
%         hh_fun = hh_fun0(1);% current fig  not GUI
%         if strcmp(hh_fun.Name, 'I_SHG_GUI')
%            hh_fun = hh_fun0(2); % current fig  not GUI
%         end
%         if get(hh_fun,'Number') == 18
%             hh_fun = hh_fun0(3);
%         end
%         hh_fig_fun=get(hh_fun, 'Children'); % current axes not GUI
%
%         handles.h_princ = hh_fig_fun(2);
% %         disp('ok')
%
%     end
% catch % fig deleted
%     handles.h_princ = gca;
% end
% h1=gca;
% axes(handles.h_princ); % disp gca

if ~handles.undocked_fig1
    handles.h_princ =  handles.plot_princ;
else
    hh_fun0 = get(0, 'Children'); % list fig
    hh_fun = hh_fun0(1);% current fig  not GUI
    
    if length(hh_fun0) <= 1 % only GUI open
        handles.h_princ =  handles.plot_princ;
    else
        
        if strcmp(hh_fun.Name, 'I_SHG_GUI')
            hh_fun = hh_fun0(2);
            if (get(hh_fun,'Number') == 18 || get(hh_fun,'Number') == 19)
                if length(hh_fun0) <= 2
                    handles.h_princ =  handles.plot_princ;
                else
                    hh_fun = hh_fun0(3);
                    h1 = get(hh_fun, 'Children'); % current axes not GUI
                    handles.h_princ=h1(end);
                    axes(handles.h_princ);
                end
            else
                h1 = get(hh_fun, 'Children'); % current axes not GUI
                handles.h_princ=h1(end);
                axes(handles.h_princ);
            end
            
            %              try
            %                  get(handles.h_hist, 'Children'); % if deleted, fails
            %                  try
            %                      get(handles.h_hist,'CData');
            %                  catch
            %                  end
            %              end
        else
            if (get(hh_fun,'Number') == 18 || get(hh_fun,'Number') == 19)
                hh_fun = hh_fun0(2);
                
            end
            h1 = get(hh_fun, 'Children'); % current axes not GUI
            handles.h_princ=h1(end);
            axes(handles.h_princ);
        end
        
    end
    
end

cmap_content = cellstr(get(hObject,'String'));
cmap_choice = cmap_content{get(hObject,'Value')};

ch = 0;

if strcmp(cmap_choice, 'Red/Black/Green')
    handles.cmap_princ_curr = handles.cmap_redgreen;
    ch = 1;
elseif strcmp(cmap_choice, 'HSV')
    handles.cmap_princ_curr = hsv;
    ch = 1;
elseif strcmp(cmap_choice, 'Parula')
    handles.cmap_princ_curr = parula;
    ch = 1;
elseif strcmp(cmap_choice, 'Cubehelix')
    handles.cmap_princ_curr = handles.cmap_cubehelix;
    ch = 1;
elseif strcmp(cmap_choice, 'Black/Red/Black/Green')
    handles.cmap_princ_curr = handles.cmap_blkredblkgrn;
    ch = 1;
    % elseif strcmp(cmap_choice, 'B/R/B/G with saturation')
    %     colormap(handles.plot_princ, handles.cmap_brbgsat);
    % elseif strcmp(cmap_choice, 'HSV with saturation')
    %     colormap(handles.plot_princ, handles.cmap_hsvsat);
elseif strcmp(cmap_choice, 'Grey')
    handles.cmap_princ_curr =  gray;
    ch = 1;
else
    disp('Colormap unchanged on axes 1.')
end

if ch
    colormap(handles.h_princ, handles.cmap_princ_curr);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cmap1_menu_CreateFcn(hObject, ~, handles)
% hObject    handle to cmap1_menu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cmap2_menu.
function cmap2_menu_Callback(hObject, ~, handles)
% hObject    handle to cmap2_menu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cmap2_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cmap2_menu

if ~handles.undocked_fig2
    handles.h_hist =  handles.plot_second;
else
    hh_fun0 = get(0, 'Children'); % list fig
    hh_fun = hh_fun0(1);% current fig  not GUI
    
    if length(hh_fun0) <= 1 % only GUI open
        handles.h_hist =  handles.plot_second;
    else
        
        if strcmp(hh_fun.Name, 'I_SHG_GUI')
            hh_fun = hh_fun0(2);
            if (get(hh_fun,'Number') == 18 || get(hh_fun,'Number') == 19)
                if length(hh_fun0) <= 2
                    handles.h_hist =  handles.plot_second;
                else
                    hh_fun = hh_fun0(3);
                    h1 = get(hh_fun, 'Children'); % current axes not GUI
                    handles.h_hist=h1(end);
                    axes(handles.h_hist);
                end
            else
                h1 = get(hh_fun, 'Children'); % current axes not GUI
                handles.h_hist=h1(end);
                axes(handles.h_hist);
            end
            
            %              try
            %                  get(handles.h_hist, 'Children'); % if deleted, fails
            %                  try
            %                      get(handles.h_hist,'CData');
            %                  catch
            %                  end
            %              end
        else
            if (get(hh_fun,'Number') == 18 || get(hh_fun,'Number') == 19)
                hh_fun = hh_fun0(2);
                
            end
            h1 = get(hh_fun, 'Children'); % current axes not GUI
            handles.h_hist=h1(end);
            axes(handles.h_hist);
        end
        
    end
    
end
% try
%     get(handles.h_hist, 'Children');
%     try
%         get(handles.h_hist,'CData');
%     catch
%
% %         if strcmp(hh_fun.Name, 'I_SHG_GUI')
% %             hh_fun = hh_fun0(2); % current fig  not GUI
% %         end
%
%         hh_fig_fun=get(hh_fun, 'Children'); % current axes not GUI
%
%         handles.h_hist = hh_fig_fun(2);
%         disp('ok')
%
%     end
% catch % fig deleted
%     handles.h_hist = gca;
% end
%
% g=handles.h_hist
% end

cmap_content = cellstr(get(hObject,'String'));
cmap_choice = cmap_content{get(hObject,'Value')};

if strcmp(cmap_choice, 'Red/Black/Green')
    colormap(handles.h_hist, handles.cmap_redgreen);
elseif strcmp(cmap_choice, 'HSV')
    colormap(handles.h_hist, hsv);
elseif strcmp(cmap_choice, 'Parula')
    colormap(handles.h_hist, parula);
elseif strcmp(cmap_choice, 'Cubehelix')
    colormap(handles.h_hist, handles.cmap_cubehelix);
elseif strcmp(cmap_choice, 'Grey')
    colormap(handles.h_hist, gray);
elseif strcmp(cmap_choice, 'Black/Red/Black/Green')
    colormap(handles.h_hist, handles.cmap_blkredblkgrn);
else
    disp('Colormap unchanged on axes 2.')
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cmap2_menu_CreateFcn(hObject, ~, handles)
% hObject    handle to cmap2_menu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in plot_diff_menu.
function plot_diff_menu_Callback(hObject, eventdata, handles) %#ok<*INUSL>
% hObject    handle to plot_diff_menu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns histmode_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from histmode_menu

% if isa(eventdata, 'matlab.ui.eventdata.ActionData') 
hhere = handles.plot_diff_menu; % % hObject val will be of the calling function  !!
% else % direct call
%    hhere = hObject;
% end
ii = get(hhere,'Value')-1;
str_l=get(hhere,'String');

if ii <= 0
    ii = get(handles.plot_diff_menu,'Value') - 1; % don't use hObject here ! it does not works with an explicit call to the callback
    if ii == 0
        ii = 1;
    end
end

% if ii < 6

cla(handles.plot_phase, 'reset');
cla(handles.plot_princ, 'reset');
set(handles.plot_princ, 'Visible', 'on');
set(handles.plot_phase, 'Visible', 'off');

%     colorbar(handles.plot_phase,'off')
axes(handles.plot_princ);
%     colorbar(handles.plot_princ)
cond1 = 1; 
if (strcmp(str_l{ii+1}, 'mapavg') || strcmp(str_l{ii+1}(1:9), 'Interf. c'))
    cond1 = 0; 
end
if (strcmp(str_l{ii+1}, 'mapavg'))
    str_l{ii+1} = ['Interf. c', str_l{ii+1}];
end

if (ii >= 2 && ii <= size(handles.contr00, 3)+1 && cond1)
    
    [handles.h_princ, ~] = diff_phase_stack_ISHG( handles.contr00, handles.x, handles.y, handles.contrast, handles.xTitle_dflt, handles.yTitle_dflt, ...
        handles.screensize, handles.fact, handles.left_offset_fig, handles.top_offset_fig, handles.undocked_fig1, handles.plot_princ, handles.axes_font_size, ...
        handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz, handles.cmap_redgreen, ii-1, ii-1, handles.x_phase );
    
else
    if handles.undocked_fig1
        
        try
            get(handles.h_princ,'Children'); % error if deleted
            %                 disp('ok')
        catch % if it has been deleted
       
            figure('outerposition',...
                [min(handles.screensize(3)*(1-handles.fact), handles.left_offset_fig) min(handles.screensize(4)*(1-handles.fact), handles.top_offset_fig) ...
                handles.screensize(3)*handles.fact handles.screensize(4)*handles.fact]);
            % [left bottom width height]
            handles.h_princ = axes; % create axes in current figure
        end
    end
    
    if (ii == 1 &&  ~strcmp(str_l{ii+1}(1:9), 'Interf. c'))  % i = 1, 1st image
        
        draw_plots_ISHG( 0, 0, handles.img_3D(:,:,1), 1:size(handles.img_3D(:,:,1), 2), 1:size(handles.img_3D(:,:,1), 1), handles.h_princ, ...
            'X (pixels)', 'Y (pixels)', 'First image of the stack', gray, 0, 'SHG int. (a. u.)', 0, 1, ...
            handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz );
        
    elseif (~mod(ii - size(handles.contr00, 3), 2) ||  strcmp(str_l{ii+1}(1:9), 'Interf. c'))  % contraste interfero (even number after contr)
        % ii == size(handles.contr00, 3) + 2
        %         parula0=parula;
        %cat(1,[1 1 1],parula0(2:end,:))
% %         handles.cmap_princ_curr = parula;
%         num_cell = max(1,(ii - size(handles.contr00, 3))/2);
           if ii >1
            n=find(contains(str_l, 'Interf. c'), 1, 'last' )-1-ii;
            num_cell = max(1, length(handles.amp_cell) - n);
           else; num_cell = 1; 
           end
        draw_plots_ISHG( 0, 0, handles.amp_cell{num_cell}, handles.x_cell{num_cell}, handles.y_cell{num_cell}, handles.h_princ, ...
            handles.xTitle_dflt, handles.yTitle_dflt, handles.Titre2, parula , 0, '', 0, 1, ...
            handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz );
        
    elseif (mod(ii - size(handles.contr00, 3), 2) ||  strcmp(str_l{ii+1}(1:10), 'Rel. Error'))  % err. (odd number after contr)
        % ii == size(handles.contr00, 3) + 3
        % Affichage de l'erreur
% %         handles.cmap_princ_curr = parula;
        num_cell = (ii - size(handles.contr00, 3)-1)/2;
        draw_plots_ISHG( 0, 0, handles.err_cell{num_cell}, handles.x_cell{num_cell}, handles.y_cell{num_cell}, handles.h_princ, ...
            handles.xTitle_dflt, handles.yTitle_dflt, handles.Titre3, parula, 0, '', 0, 1, ...
            handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz );
    
    end
end
% pause(0.1);

if handles.undocked_fig1
    handles.h_princ_out = handles.h_princ;
end


% hf=figure;
% ha=get(hf, 'Children'); ha(1) = colorbar, ha(2) = axes
% haa=get(ha(2), 'Children'); haa = Image handle
% hf2 = figure(2); ha2 = axes; copyobj(haa, ha2); %copy image object into
% new figure (does not conserve axes)
% hf2 = figure(2); set(ha(2), 'Parent', hf2); % copy all axes + img into
% new figure, but let the original one empty (and colormap not copied)

guidata(hObject, handles);

% --- Executes on button press in undk_p1_chck.
function undk_p1_chck_Callback(hObject, eventdata, handles)
% hObject    handle to undk_p1_chck (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of undk_p1_chck

state = get(hObject,'Value');

switch state
    case 1 % undock fig
        handles.undocked_fig1 = 1;
        handles.h_princ = handles.h_princ_out;
    case 0
        handles.undocked_fig1 = 0;
        handles.h_princ = handles.plot_princ ;
end

guidata(hObject, handles);

% --- Executes on button press in undk_p2_chck.
function undk_p2_chck_Callback(hObject, ~, handles)
% hObject    handle to undk_p2_chck (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of undk_p2_chck

state = get(hObject,'Value');

switch state
    case 1 % undock
        handles.undocked_fig2 = 1;
        handles.h_hist = handles.h_hist_out;
        %         get(handles.h_hist, 'Children') % !!!!
    case 0
        handles.undocked_fig2 = 0;
        handles.h_hist = handles.plot_second;
end

guidata(hObject, handles);



function contrast_edt_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of contrast_edt as text
%        str2double(get(hObject,'String')) returns contents of contrast_edt as a double

handles.contrast = str2double(get(hObject,'String'));
guidata(hObject, handles);

if (isnan(handles.contrast) || handles.contrast < 0)
    disp('Wrong value of contrast entered !');
    set(handles.contrast_edt, 'String', num2str(handles.contrast_dflt));
end

plot_diff_menu_Callback(hObject, 0, handles)

% guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function contrast_edt_CreateFcn(hObject, ~, handles)
% hObject    handle to contrast_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function nbinsX_edt_Callback(hObject, ~, handles)
% hObject    handle to nbinsX_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nbinsX_edt as text
%        str2double(get(hObject,'String')) returns contents of nbinsX_edt as a double

nbinsX = str2double(get(hObject,'String'));

if nbinsX < 1
    disp('Wrong value for nbins in X !');
    set(hObject, 'String', num2str(handles.nbinsx_dflt));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nbinsX_edt_CreateFcn(hObject, ~, handles)
% hObject    handle to nbinsX_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nbinsY_edt_Callback(hObject, ~, handles)
% hObject    handle to nbinsY_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nbinsY_edt as text
%        str2double(get(hObject,'String')) returns contents of nbinsY_edt as a double

nbinsY = str2double(get(hObject,'String'));

if nbinsY < 1
    disp('Wrong value for nbins in Y !');
    set(hObject, 'String', num2str(handles.nbinsy_dflt));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nbinsY_edt_CreateFcn(hObject, ~, handles)
% hObject    handle to nbinsY_edt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sat_value_slider_Callback(hObject, ~, handles)
% hObject    handle to sat_value_slider (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% between 0 and 1

taille = 64;
val = 1 - get(hObject,'Value');

handles.cmap_brbgsat =  [handles.cmap_brbgsat_dflt(1:end-round(val*taille), :); repmat([1 1 0], round(val*taille) , 1)]; % for yellow

handles.cmap_hsvsat = [handles.cmap_hsvsat_dflt(1:end-round(val*taille), :); repmat([0 0 0], round(val*taille) , 1)]; % for black


cmap_content = cellstr(get(handles.cmap1_menu,'String'));
cmap_choice = cmap_content{get(handles.cmap1_menu,'Value')};
if (strcmp(cmap_choice, 'HSV with saturation') || strcmp(cmap_choice, 'HSV'))
    cmap = handles.cmap_hsvsat;
else
    cmap = handles.cmap_brbgsat;
end

draw_plots_ISHG( 0, 1,  handles.phase_test,  handles.x,  handles.y, handles.plot_princ, ...
    handles.xTitle_dflt,  handles.yTitle_dflt,  handles.Titre1,  cmap, 0,  handles.phi_mat_default, 0, 1, ...
    handles.axes_font_size,  handles.xaxis_sz,  handles.yaxis_sz,  handles.title_sz,  handles.clrbr_tl_sz );

set(handles.plot_princ,'CLim',[-1 1.05]);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sat_value_slider_CreateFcn(hObject, ~, handles)
% hObject    handle to sat_value_slider (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function handles = fit_hist_util(handles)

handles.fit_chosen = get(handles.fit_hist_button,'Value') - 1;
% % disp(num2str(handles.fit_chosen))

if (handles.fit_chosen == 5 || handles.fit_chosen == 10)
    handles.M = str2double(get(handles.edt_pearson, 'String'));
else
    handles.M = 0.5;
end

% switch state
%
%     case 1 % gaussian
%     case 2 % lorentzian
%
%     case 3 % prod lorentzian/ gaussian
%
%     case 4 % sum lorentzian/ gaussian (Voigt)
%
% end
try
    ax=get(50, 'Children');hist1=get(ax, 'Children');
    if isa(hist1, 'cell'); hist1 = hist1{end}; end
    if isa(hist1, 'matlab.graphics.primitive.Data'); hist1 = hist1(end); end
    hist2_ydata = hist1.BinCounts; hist2_xdata = hist1.BinEdges(1:end-1);
catch
    hist2_xdata=handles.hist2_xdata(1:end-1); hist2_ydata=handles.hist2_ydata(1:end-1);
end
func_hdl = fit_hist2d_funcs;
[ handles.stringb1, handles.stringc1, handles.stringb2, handles.stringc2, handles.cmap_brbg_fit, handles.fit2dstruct, ~  ] = ...
    func_hdl.fit_hist2d_ishg(hist2_xdata, hist2_ydata, handles.hhist1, handles.Legendehisto, handles.offset, handles.fit_chosen, handles.M, ...
str2double(get(handles.norm_fact_hist, 'String')),  str2double(get(handles.baseline_hist_edt, 'String')));

% code to change colormap on correct figure
if ~handles.undocked_phase
    handles.h_phase =  handles.plot_phase;
else
    hh_fun0 = get(0, 'Children'); % list fig
    hh_fun = hh_fun0(1);% current fig  not GUI
    
    if length(hh_fun0) <= 1 % only GUI open
        handles.h_phase =  handles.plot_phase;
    else
        
        if strcmp(hh_fun.Name, 'I_SHG_GUI')
            hh_fun = hh_fun0(2);
            if (get(hh_fun,'Number') == 18 || get(hh_fun,'Number') == 19)
                if length(hh_fun0) <= 2
                    handles.h_phase =  handles.plot_phase;
                else
                    hh_fun = hh_fun0(3);
                    h1 = get(hh_fun, 'Children'); % current axes not GUI
                    handles.h_phase=h1(end);
                    axes(handles.h_phase);
                end
            else
                h1 = get(hh_fun, 'Children'); % current axes not GUI
                handles.h_phase=h1(end);
                axes(handles.h_phase);
            end
            
        else
            if (get(hh_fun,'Number') == 18 || get(hh_fun,'Number') == 19)
                hh_fun = hh_fun0(2);
                
            end
            h1 = get(hh_fun, 'Children'); % current axes not GUI
            handles.h_phase=h1(end);
            axes(handles.h_phase);
        end
        
    end
end

handles.cmap_brbg_fit(handles.cmap_brbg_fit <0) = 0;
handles.cmap_brbg_fit(handles.cmap_brbg_fit >1) = 1;

colormap(handles.h_phase, handles.cmap_brbg_fit);

handles.cmap_default1 = handles.cmap_brbg_fit;

str = get(handles.cmap_p_menu, 'String' );
% % whos str
% name = 'B/R/B/G of fit';
if ~strcmp(str{end}, 'B/R/B/G of fit')
    handles.name_c_p_menu_str{end+1} = 'B/R/B/G of fit';
    str{end+1} = 'B/R/B/G of fit';
    set(handles.cmap_p_menu, 'String', str);
    %     disp('ok')
end
set(handles.cmap_p_menu, 'Value', size(str,1));

% --- Executes on button press in fit_hist_button.
function fit_hist_button_Callback(hObject, ~, handles)
% hObject    handle to fit_hist_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = fit_hist_util(handles);

str1=get(handles.two_d_hist_popup, 'String');
if (length(str1)< 4 && ~strcmp(str1{end} , handles.str_ratiof_sigma) && ~strcmp(str1{end-1} , handles.str_ratiof_sigma))
    str1{end+1}= handles.str_ratiof_sigma;
    set(handles.two_d_hist_popup, 'String', str1);
end

guidata(hObject, handles);


% --- Executes on selection change in cmap_p_menu.
function cmap_p_menu_Callback(hObject, ~, handles)
% hObject    handle to cmap_p_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cmap_p_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cmap_p_menu

% try
%     get(handles.h_phase, 'Children');
%     try
%         get(handles.h_phase,'CData')
%     catch
%         hh_fun0 = get(0, 'Children'); % current fig  not GUI
%         hh_fun = hh_fun0(1);% current fig  not GUI
%         if strcmp(hh_fun.Name, 'I_SHG_GUI')
%             hh_fun = hh_fun0(2); % current fig  not GUI
%         end
%         if get(hh_fun,'Number') == 18
%             hh_fun = hh_fun0(3);
%         end
%         hh_fig_fun=get(hh_fun, 'Children'); % current axes not GUI
%
%         handles.h_phase = hh_fig_fun(2);
%         %         disp('ok')
%
%     end
% catch % fig deleted
%     handles.h_phase = gca;
% end
%
% axes(handles.h_phase);

if ~handles.undocked_phase
    handles.h_phase =  handles.plot_phase;
else
    hh_fun0 = get(0, 'Children'); % list fig
    hh_fun = hh_fun0(1);% current fig  not GUI
    
    if length(hh_fun0) <= 1 % only GUI open
        handles.h_phase =  handles.plot_phase;
    else
        
        if strcmp(hh_fun.Name, 'I_SHG_GUI')
            hh_fun = hh_fun0(2);
            if (get(hh_fun,'Number') == 18 || get(hh_fun,'Number') == 19)
                if length(hh_fun0) <= 2
                    handles.h_phase =  handles.plot_phase;
                else
                    hh_fun = hh_fun0(3);
                    h1 = get(hh_fun, 'Children'); % current axes not GUI
                    handles.h_phase=h1(end);
                    axes(handles.h_phase);
                end
            else
                h1 = get(hh_fun, 'Children'); % current axes not GUI
                handles.h_phase=h1(end);
                axes(handles.h_phase);
            end
            
        else
            if (get(hh_fun,'Number') == 18 || get(hh_fun,'Number') == 19)
                hh_fun = hh_fun0(2);
                
            end
            h1 = get(hh_fun, 'Children'); % current axes not GUI
            handles.h_phase=h1(end);
            axes(handles.h_phase);
        end
        
    end
end

cmap_content = cellstr(get(hObject,'String'));
cmap_choice = cmap_content{get(hObject,'Value')};

ch = 0;

if strcmp(cmap_choice, 'Red/Black/Green')
    handles.cmap_default1 = handles.cmap_redgreen;
    ch = 1;
    
elseif strcmp(cmap_choice, 'HSV')
    handles.cmap_default1 = hsv;
    ch = 1;
    %     colormap(handles.h_phase, hsv);
elseif strcmp(cmap_choice, 'Parula')
    handles.cmap_default1 = parula;
    ch = 1;
    %     colormap(handles.h_phase, parula);
elseif strcmp(cmap_choice, 'Cubehelix')
    handles.cmap_default1 = handles.cmap_cubehelix;
    ch = 1;
    %     colormap(handles.h_phase, handles.cmap_cubehelix);
elseif strcmp(cmap_choice, 'Black/Red/Black/Green')
    handles.cmap_default1 = handles.cmap_blkredblkgrn;
    ch = 1;
    %     colormap(handles.h_phase, handles.cmap_blkredblkgrn);
elseif strcmp(cmap_choice, 'B/R/B/G with saturation')
    handles.cmap_default1 = handles.cmap_brbgsat;
    ch = 1;
    %     colormap(handles.h_phase, handles.cmap_brbgsat);
elseif strcmp(cmap_choice, 'HSV with saturation')
    handles.cmap_default1 = handles.cmap_hsvsat;
    ch = 1;
    %     colormap(handles.h_phase, handles.cmap_hsvsat);
    
    set(handles.sat_value_slider, 'Value', 63/64);
elseif strcmp(cmap_choice, 'Grey')
    %     disp('ok')
    handles.cmap_default1 = gray;
    ch = 1;
    %     colormap(handles.h_phase, gray);
elseif strcmp(cmap_choice, 'B/R/B/G of fit')
    %     colormap(handles.h_phase, handles.cmap_brbg_fit);
    handles.cmap_default1 = handles.cmap_brbg_fit;
    ch = 1;
else
    disp('Colormap unchanged on axes phase.')
end

if ch
    colormap(handles.h_phase, handles.cmap_default1);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cmap_p_menu_CreateFcn(hObject, ~, handles)
% hObject    handle to cmap_p_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in undck_plotphase.
function undck_plotphase_Callback(hObject, eventdata, handles)
% hObject    handle to undck_plotphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of undck_plotphase

state = get(hObject,'Value');

switch state
    case 1 % undock
        handles.undocked_phase = 1;
        handles.h_phase =  handles.h_phase_out;
    case 0
        handles.undocked_phase = 0;
        handles.h_phase =  handles.plot_phase;
end

guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function plot_phase_DeleteFcn(hObject, ~, handles)
% hObject    handle to plot_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edt_pearson_Callback(hObject, eventdata, handles)
% hObject    handle to edt_pearson (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_pearson as text
%        str2double(get(hObject,'String')) returns contents of edt_pearson as a double


% --- Executes during object creation, after setting all properties.
function edt_pearson_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_pearson (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function fit_hist_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fit_hist_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function three_phase_edt_Callback(hObject, eventdata, handles)
% hObject    handle to three_phase_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of three_phase_edt as text
%        str2double(get(hObject,'String')) returns contents of three_phase_edt as a double


% --- Executes during object creation, after setting all properties.
function three_phase_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to three_phase_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function step_ps_edt_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to step_ps_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of step_ps_edt as text
%        str2double(get(hObject,'String')) returns contents of step_ps_edt as a double

set(handles.max_ps_indic_edt, 'String', num2str(str2double(handles.first_phase_edt.String) + 360/2*max(2, str2double(handles.slice_per_step_edt.String)) - str2double(handles.step_ps_edt.String)));
set(handles.snd_step_edt, 'String',  get(hObject, 'String'));

% --- Executes during object creation, after setting all properties.
function step_ps_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to step_ps_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in contr_radio.
function contr_radio_Callback(hObject, eventdata, handles)
% hObject    handle to contr_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of contr_radio

if get(hObject,'Value')
    set(handles.double_contr_chck, 'Enable', 'on');
else
    set(handles.double_contr_chck, 'Enable', 'off');
end

double_contr_util(handles.double_contr_chck, handles, 'normal');

% --- Executes on button press in raw_radio.
function raw_radio_Callback(hObject, eventdata, handles)
% hObject    handle to raw_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of raw_radio

if get(hObject,'Value')
    set(handles.double_contr_chck, 'Enable', 'off');
end
double_contr_util(handles.double_contr_chck, handles, 'normal');

function snd_step_edt_Callback(hObject, eventdata, handles)
% hObject    handle to snd_step_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snd_step_edt as text
%        str2double(get(hObject,'String')) returns contents of snd_step_edt as a double


% --- Executes during object creation, after setting all properties.
function snd_step_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snd_step_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in notuse_A_chck.
function notuse_A_chck_Callback(hObject, eventdata, handles)
% hObject    handle to notuse_A_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of notuse_A_chck


% --- Executes on button press in three_ph_algo_chck.
function three_ph_algo_chck_Callback(hObject, eventdata, handles)
% hObject    handle to three_ph_algo_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_phase_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_phase_popup

val = get(hObject,'Value');

if val >= 4 % % algo 3 frames
    set(handles.it_max_algo_edt, 'Visible', 'on');
    set(handles.max_it_algo_txt, 'Visible', 'on');
    set(handles.epsilon_ph_th_edt, 'Visible', 'on');
    set(handles.eps_threshold_ph_txt, 'Visible', 'on');
    if val >= 6 % % with tilt
        set(handles.epsilon_tilt_th_edt, 'Visible', 'on');
        set(handles.eps_threshold_tilt_txt, 'Visible', 'on');
        if val == 6
            set(handles.epsilon_ph_th_edt, 'Visible', 'off');
            set(handles.eps_threshold_ph_txt, 'Visible', 'off');
        elseif val == 7
            set(handles.epsilon_ph_th_edt, 'Visible', 'on');
            set(handles.eps_threshold_ph_txt, 'Visible', 'on');
        end
    else
        set(handles.epsilon_tilt_th_edt, 'Visible', 'off');
        set(handles.eps_threshold_tilt_txt, 'Visible', 'off');
    end
else
    set(handles.it_max_algo_edt, 'Visible', 'off');
    set(handles.max_it_algo_txt, 'Visible', 'off');
    set(handles.epsilon_ph_th_edt, 'Visible', 'off');
    set(handles.eps_threshold_ph_txt, 'Visible', 'off');
    set(handles.eps_threshold_tilt_txt, 'Visible', 'off');
    set(handles.epsilon_tilt_th_edt, 'Visible', 'off');

end

% --- Executes on button press in med_fltr_radio.
function med_fltr_radio_Callback(hObject, ~, handles)
% hObject    handle to med_fltr_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of med_fltr_radio

cond = get(hObject,'Value'); % if checked;

if cond
    prompt = {'Sigma value for the exponential median filter (default 25.5, 0=unweighted filter)'};
    dlg_title = 'Weight';
    num_lines = 1;
    % Valeurs par défaut
    def = {'0'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    if ~isempty(answer)
        handles.sigma_med = str2double(answer{1});
        handles.weight = 1;
        handles.mean_phase_meth = 0; % median
    else
        cond=0;
    end
end

if ~cond
    handles.weight = 0;
    handles.mean_phase_meth = 1; % mean
    
end

guidata(hObject, handles);


% --- Executes on button press in no_fltr_radio.
function no_fltr_radio_Callback(hObject, ~, handles)
% hObject    handle to no_fltr_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of no_fltr_radio

if get(hObject,'Value') % if checked
    handles.weight = 0;
else
    handles.weight = 1;
end

guidata(hObject, handles);



function offset_neutral_Callback(hObject, eventdata, handles)
% hObject    handle to offset_neutral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset_neutral as text
%        str2double(get(hObject,'String')) returns contents of offset_neutral as a double

switch str2double(get(hObject,'String'))
    case 0
        set(handles.apply_offset_button, 'Enable', 'off')
    otherwise
        set(handles.apply_offset_button, 'Enable', 'on')
end

% --- Executes during object creation, after setting all properties.
function offset_neutral_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset_neutral (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in apply_offset_button.
function apply_offset_button_Callback(hObject, eventdata, handles)
% hObject    handle to apply_offset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~strcmp(get(handles.plot_phase, 'Visible'), 'on'); plot_phase_popup_Callback(handles.plot_phase_popup, 0, handles); end

handles.off_neutral = str2double(get(handles.offset_neutral, 'String'));

im=get(handles.plot_phase, 'Children'); phase1 = im.CData; %handles.phase_cell{end} ;

str1=get(handles.plot_phase_popup, 'String'); ii=get(handles.plot_phase_popup, 'Value');
off = max(0, length(str1) - ii); 

phase1(phase1 == 1.05) = 1.05e6; % because otherwise it gets corrected also

phase_corr = phase1 + handles.off_neutral; % in frac of pi

phase_corr(phase_corr>1 & phase_corr~= (1.05e6+ handles.off_neutral)) = phase_corr(phase_corr>1 & phase_corr~= (1.05e6+ handles.off_neutral)) - 2;
phase_corr(phase_corr<-1 & phase_corr~= (1.05e6+ handles.off_neutral)) = phase_corr(phase_corr<-1 & phase_corr~= (1.05e6+ handles.off_neutral)) + 2;

phase_corr(phase_corr == (1.05e6+ handles.off_neutral)) = 1.05;

handles.phase_cell{end+1} = phase_corr;
handles.phase_cell00{end+1} = phase_corr;

handles = add_data_in_list(hObject, handles, 1, 1, 1, off, [1, Inf], [1, Inf]); % doxy, do_amp, do_err, off
if ~isa(handles.amp_cell00, 'cell'); handles.amp_cell00 = {handles.amp_cell00}; end
handles.amp_cell00(end+1) = handles.amp_cell(end-off);

add_plot_ph(hObject, handles, sprintf('Phase map %d', length(handles.phase_cell)));

guidata(hObject, handles);

plot_phase_popup_Callback(handles.plot_phase_popup, eventdata, handles);

guidata(hObject, handles);


% --- Executes on button press in specify_order.
function specify_order_Callback(hObject, eventdata, handles)
% hObject    handle to specify_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of specify_order



function useonly_strt_edt_Callback(hObject, eventdata, handles)
% hObject    handle to useonly_strt_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of useonly_strt_edt as text
%        str2double(get(hObject,'String')) returns contents of useonly_strt_edt as a double

handles.use_range = 1;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function useonly_strt_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to useonly_strt_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function useonly_stp_edt_Callback(hObject, eventdata, handles)
% hObject    handle to useonly_stp_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of useonly_stp_edt as text
%        str2double(get(hObject,'String')) returns contents of useonly_stp_edt as a double

handles.use_range = 1;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function useonly_stp_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to useonly_stp_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function useonly_end_edt_Callback(hObject, eventdata, handles)
% hObject    handle to useonly_end_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of useonly_end_edt as text
%        str2double(get(hObject,'String')) returns contents of useonly_end_edt as a double

handles.use_range = 1;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function useonly_end_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to useonly_end_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in near_chck.
function near_chck_Callback(hObject, eventdata, handles)
% hObject    handle to near_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of near_chck


% --- Executes on button press in norm_chck.
function norm_chck_Callback(hObject, eventdata, handles)
% hObject    handle to norm_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of norm_chck


% --- Executes on button press in LU_inv_chck.
function LU_inv_chck_Callback(hObject, eventdata, handles)
% hObject    handle to LU_inv_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LU_inv_chck


% --- Executes on button press in parallel_comp_chck.
function parallel_comp_chck_Callback(hObject, ~, handles)
% hObject    handle to parallel_comp_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parallel_comp_chck

fprintf(2, ' Parallel computing did not really improved for number of pixels app. < 1e5\n');
fprintf(2, 'Use only for matrices more than like 300x300 if long algorithms\n');

if get(hObject,'Value')
    p = gcp('nocreate'); % get current pool
    if ~size(p,1) % if no pool
        parpool; % start parallel
    end
end


function xmin_crop_edit_Callback(hObject, eventdata, handles)
% hObject    handle to xmin_crop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xmin_crop_edit as text
%        str2double(get(hObject,'String')) returns contents of xmin_crop_edit as a double


% --- Executes during object creation, after setting all properties.
function xmin_crop_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xmin_crop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ymin_crop_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ymin_crop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ymin_crop_edit as text
%        str2double(get(hObject,'String')) returns contents of ymin_crop_edit as a double


% --- Executes during object creation, after setting all properties.
function ymin_crop_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ymin_crop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function width_crop_edit_Callback(hObject, eventdata, handles)
% hObject    handle to width_crop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width_crop_edit as text
%        str2double(get(hObject,'String')) returns contents of width_crop_edit as a double


% --- Executes during object creation, after setting all properties.
function width_crop_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width_crop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function height_crop_edit_Callback(hObject, eventdata, handles)
% hObject    handle to height_crop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of height_crop_edit as text
%        str2double(get(hObject,'String')) returns contents of height_crop_edit as a double


% --- Executes during object creation, after setting all properties.
function height_crop_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to height_crop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in crop_chck_box.
function crop_chck_box_Callback(hObject, eventdata, handles)
% hObject    handle to crop_chck_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of crop_chck_box

if get(hObject,'Value')
    handles.crop_batch = 1;
    xmin = 32.775132;
    ymin = 28.542328;
    width = 129.100529;
    height = 369.312169;
    set(handles.xmin_crop_edit, 'String', num2str(xmin));
    set(handles.ymin_crop_edit, 'String', num2str(ymin));
    set(handles.width_crop_edit, 'String', num2str(width));
    set(handles.height_crop_edit, 'String', num2str(height));
    % last : 36.478836, 52.087302, 122.751323x122.751323
    % fibril stef : 112.140212, 8.436508, 55.026455x390.476190
    % fibril wide : 32.775132, 28.542328, 129.100529x369.312169
else % unchecked
    handles.crop_batch =0;
    set(handles.xmin_crop_edit, 'String', num2str(0));
    set(handles.ymin_crop_edit, 'String', num2str(0));
    set(handles.width_crop_edit, 'String', num2str(0));
    set(handles.height_crop_edit, 'String', num2str(0));
end



% --- Executes on button press in crop_phase.
function crop_phase_Callback(hObject, eventdata, handles)
% hObject    handle to crop_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.phase_cell)
    phase1_origin = handles.phase_cell{end};
else % no ph
    fprintf(2,'no phase\n');  return
end
usr = 0;

str1 = get(handles.plot_phase_popup, 'String');
img00 = handles.amp_cell{end};
if length(handles.err_cell) >=1
    err00 = handles.err_cell{end};
else
    err00 = img00*0;
end

while usr~=1 % Variable pour le contentement de l'usager
    
    set(handles.plot_phase_popup, 'String', str1)
    set(handles.plot_phase_popup, 'Value', length(str1));%length(str1)); % reset to init value
    %     handles.x_current = handles.x;
    %     handles.y_current = handles.y;
    plot_phase_popup_Callback(hObject, eventdata, handles);
    axes(handles.plot_phase); %#ok<LAXES>
    %     imagesc(handles.phase1); colormap(hsv); colorbar;
    
    title(handles.plot_phase, 'Select a ROI with the mouse (right click to make it square)')
    rect = getrect(handles.plot_phase); % [xmin ymin width height];
    
    xmin = (rect(1)-min(handles.x))*(length(handles.x)-1)/(max(handles.x)-min(handles.x))+1;
    ymin = (rect(2)-min(handles.y))*(length(handles.y)-1)/(max(handles.y)-min(handles.y))+1;
    width = (rect(3)-min(handles.x))*(length(handles.x)-1)/(max(handles.x)-min(handles.x))+1;
    height = (rect(4)-min(handles.y))*(length(handles.y)-1)/(max(handles.y)-min(handles.y))+1;
    %     xmin = rect(1);
    %     ymin = (rect(2));
    %     width = (rect(3));
    %     height = (rect(4));
    if (round(ymin) > size(phase1_origin, 1) || round(ymin + height - 1) < 0)  ...
            || (round(xmin) > size(phase1_origin, 2) || round(xmin + width - 1) < 0)
        clear rect; disp('rect out of bounds !');
    else    
        ymin00 = ymin; xmin00 = xmin;
        ymin = max(ymin, 1); xmin = max(xmin, 1); 
        width =  min(size(phase1_origin, 2)-xmin +1, width) + min(ymin00, 1)-1 ; 
        height =  min(size(phase1_origin, 1)-ymin +1, height) + min(xmin00, 1)-1 ;  
        handles.phase_cell(end+1) = {phase1_origin(  round(ymin): round(ymin + height - 1) ,round(xmin): round(xmin + width - 1))};
        
        img = img00(round(ymin): round(ymin + height - 1), round(xmin): round(xmin + width - 1), :);
        
        handles.amp_cell(end+1) = {img};
        
        if (size(err00,1)>0 && size(err00,2)>0)
            handles.err_cell{end+1} = err00(max(round(ymin),0): min(round(ymin + height - 1),size(err00,1)), max(round(xmin),0): min(round(xmin + width - 1),size(err00,2)), :);
        end
        handles.x_cell{end+1} = linspace(0, size(handles.phase_cell{end},2)*handles.resx, size(handles.phase_cell{end},2));
        handles.y_cell{end+1} = linspace(0, size(handles.phase_cell{end},1)*handles.resx, size(handles.phase_cell{end},1));
        %     figure;  imagesc(handles.phase1); colormap(hsv); colorbar;
        
        str11 = str1;
        str11{end+1} = sprintf('phase map %d', length(handles.phase_cell)); %#ok<AGROW>
        set(handles.plot_phase_popup, 'String', str11)
        
        set(handles.plot_phase_popup, 'Value', length(str11)); % set to last value
        plot_phase_popup_Callback(hObject, eventdata, handles);
        
        usr = menu('Région OK ?','Oui','Non');
        
        if usr ~=1
            handles.phase_cell(end) = []; % erase rejected phase map
            handles.amp_cell(end) = [];
            str11(end) = [];
            handles.x_cell(end) = []; handles.y_cell(end) = [];
            set(handles.plot_phase_popup, 'String', str11)
            
            set(handles.plot_phase_popup, 'Value', length(str11)); % set to last value
        else
            str2 = get(handles.plot_diff_menu, 'String');
            str2{end+1} = sprintf('Interf. contrast %d', length(handles.amp_cell)); %#ok<AGROW>
            str2{end+1} = sprintf('Err. map %d', length(handles.err_cell)); %#ok<AGROW>
            set(handles.plot_diff_menu, 'String', str2)     
        end
    end
end


% hh = pl_phasemap_ISHG(handles.undocked_fig1, handles.screensize, handles.fact, handles.left_offset_fig, ...
% handles.top_offset_fig, h, handles.phase_test, handles.x, handles.y, ...
%     handles.xTitle_dflt, handles.yTitle_dflt, handles.Titre1, handles.cmap_brbgsat, handles.phi_mat_default, ...
% handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz);

guidata(hObject, handles);


% --- Executes on button press in med_contr_chck.
function med_contr_chck_Callback(hObject, eventdata, handles)
% hObject    handle to med_contr_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of med_contr_chck

% --- Executes on button press in repair_ph_button.
function repair_ph_button_Callback(hObject, eventdata, handles)
% hObject    handle to repair_ph_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject, 'Value')-1
    case 1 % erase
        % % del
        switch questdlg('Sure to erase phases?','Erase all phases (keep 1st one)','Yes','No','No')
            case 'No'
                disp('user canceled erase');
            case 'Yes'
                handles.phase_cell = handles.phase_cell(1);
                handles.amp_cell = handles.amp_cell(1);
                handles.err_cell = handles.err_cell(1);
                handles.x_cell = handles.x_cell(1);
                handles.y_cell = handles.y_cell(1);
                str0={'View plot phase ...', 'Phase map'};
                set(handles.plot_phase_popup, 'Value', 2);
                set(handles.plot_phase_popup, 'String', str0);
                str1={'View plot ...', 'Interf. ctr'};
                set(handles.plot_diff_menu, 'Value', 2);
                set(handles.plot_diff_menu, 'String', str1);
                guidata(hObject, handles);
        end
    case 2 % repair
        switch questdlg('Sure to set phase00 to phase cell?','Homogen. phases (keep phase_cell)','Yes','No','No')
            case 'Yes'
                handles.phase_cell00=handles.phase_cell;
                handles.amp_cell00 = handles.amp_cell;
        end
end
set(hObject, 'Value',1);
guidata(hObject, handles);


% --- Executes on button press in reset_range_chck.
function reset_range_chck_Callback(hObject, eventdata, handles)
% hObject    handle to reset_range_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reset_range_chck

switch get(hObject,'Value') % 1 if reset range
    
    case 1 % reset to 0
        %         handles.range_sat_old = handles.range_sat;
        %         handles.range_sat = 0;
        
        handles.old_min_contr_hist = str2double(get(handles.min_hist_edt, 'String'));
        handles.old_max_contr_hist = str2double(get(handles.max_hist_edt, 'String'));
        
        set(handles.max_hist_edt, 'String', num2str(handles.max_amp00));
        set(handles.min_hist_edt, 'String', '0');
        
        
    case 0 % use last range
        %         handles.range_sat = handles.range_sat_old;
        
        set(handles.max_hist_edt, 'String', handles.old_max_contr_hist);
        set(handles.min_hist_edt, 'String', handles.old_min_contr_hist );
        
end

guidata(hObject, handles);


% --- Executes on selection change in plot_phase_popup.
function plot_phase_popup_Callback(hObject, eventdata, handles)
% hObject    handle to plot_phase_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plot_phase_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plot_phase_popup

ii = get(handles.plot_phase_popup,'Value')-1; % not hObject

if ii <= 0
    ii = get(handles.plot_phase_popup,'Value') - 1; % don't use hObject here ! it does not works with an explicit call to the callback
    if ii == 0
        ii = 1;
    end
end

colorbar(handles.plot_princ,'off')
cla(handles.plot_princ, 'reset');
set(handles.plot_princ, 'Visible', 'off');
set(handles.plot_phase, 'Visible', 'on');

axes(handles.plot_phase);

ii_safe = min(ii, length(handles.phase_cell));
[handles.h_phase, handles.hfig_out] = pl_phasemap_ISHG(handles.undocked_phase, handles.screensize, handles.fact, handles.left_offset_fig, ...
    handles.top_offset_fig, handles.h_phase, handles.phase_cell{ii_safe}, handles.x_cell{ii_safe}, handles.y_cell{ii_safe}, ...
    handles.xTitle_dflt, handles.yTitle_dflt, handles.Titre1, handles.cmap_default1, handles.phi_mat_default, ...
    handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz);

% cmap_p_menu

if ii > 2 % 1 =non sat
    if find(handles.phase_cell{ii_safe}==1.05) % sat. value
        set(handles.cmap_p_menu, 'Value', 3);
        set(handles.h_phase,'CLim',[-1 1.05]);
        colormap(handles.h_phase, handles.cmap_hsvsat);
    end
end

handles.int_x = 1./str2double(handles.nbinsX_edt.String); % int_x is inverse of the nbins

handles.int_y = str2double(handles.nbinsY_edt.String);%% Histogramme 3D

[ handles.data_tot, handles.ylabel_hist3, handles.title_hist3, handles.nbinsx, handles.h_hist ] = ...
    hist_3D_ISHG(handles.phase_cell{ii}, handles.amp_cell{ii}, handles.int_x, handles.int_y, ...
    handles.title_hist1, handles.yaxis_hist1, handles.Counts, handles.phi_mat_default,...
    handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz, handles.undocked_fig2, ...
    handles.h_hist, handles.screensize, handles.fact, handles.left_offset_fig, handles.top_offset_fig, handles.offset_pi2, handles.range_sat );
% Histogram 3D function

handles.max_amp = round(max(max(handles.data_tot(:,2))));
set(handles.max_hist_edt, 'String', num2str(handles.max_amp));
min_amp = round(min(min(handles.data_tot(:,2))));
set(handles.min_hist_edt, 'String', num2str(min_amp));
% % colorbar(handles.plot_phase)
% % % else % with sat
% % handles.h_phase = pl_phasemap_ISHG(handles.undocked_phase, handles.screensize, handles.fact, handles.left_offset_fig, ...
% %     handles.top_offset_fig, handles.h_phase, handles.phase_test, handles.x, handles.y, ...
% %     handles.xTitle_dflt, handles.yTitle_dflt, handles.Titre1, handles.cmap_sat_dflt, handles.phi_mat_default, ...
% %     handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz);
% %
% % handles.ph_hist_current = handles.phase_test_hist;
% % handles.second_hist_current = handles.second_test_hist;


if handles.undocked_fig1
    handles.h_phase_out = handles.h_phase;
end

if handles.flag_save_fig_out
    d= uigetdir;
    
    if d == 0
        %     d=pwd;
        % end
        disp('Canceled saving !')
    else
%         %     if strcmp(get(gcf, 'Name'), 'I_SHG_GUI')
%         %         h_save =
%         %     else
        if isa(handles.fname, 'cell') % many files
            fname = handles.fname{1};
        else
            fname = handles.fname;
        end
        
        handles.hfigsave_graph(1) = handles.hfig_out;
        savefig(handles.hfigsave_graph, fullfile(d, sprintf('%s_ph.fig', fname(1:end-4))) ,'compact');
        
        set(handles.undck_plotphase, 'Value', 0);
        
        handles.undocked_phase = 0;
        handles.h_phase =  handles.plot_phase;

    end
end

handles.flag_save_fig_out = 0;
set_histrange_vis(handles);
% pause(0.1);

guidata(hObject, handles);

function set_histrange_vis(handles)

set(handles.hist_range_button, 'Visible', 'on');
set(handles.min_hist_edt, 'Visible', 'on');
set(handles.max_hist_edt, 'Visible', 'on');
set(handles.ampmin_txt, 'Visible', 'on');
set(handles.ampmax_txt, 'Visible', 'on');

% --- Executes during object creation, after setting all properties.
function plot_phase_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_phase_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nb_img_in_stack_edt_Callback(hObject, eventdata, handles)
% hObject    handle to nb_img_in_stack_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nb_img_in_stack_edt as text
%        str2double(get(hObject,'String')) returns contents of nb_img_in_stack_edt as a double


% --- Executes during object creation, after setting all properties.
function nb_img_in_stack_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nb_img_in_stack_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ph_color_scale_strt_Callback(hObject, eventdata, handles)
% hObject    handle to ph_color_scale_strt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ph_color_scale_strt as text
%        str2double(get(hObject,'String')) returns contents of ph_color_scale_strt as a double

strt_ph_scale = str2double(get(hObject,'String'));

if (strt_ph_scale == -1 && str2double(get(handles.ph_color_scale_end,'String')) == 1)
    colormap(handles.h_phase, hsv);
else
    
    colormap(handles.h_phase, jet);
    
end

if (strt_ph_scale < str2double(get(handles.ph_color_scale_end,'String')))
    
    set(handles.h_phase, 'CLim', [strt_ph_scale, str2double(get(handles.ph_color_scale_end,'String'))]);
    colorbar(handles.h_phase);
end

% --- Executes during object creation, after setting all properties.
function ph_color_scale_strt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ph_color_scale_strt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ph_color_scale_end_Callback(hObject, eventdata, handles)
% hObject    handle to ph_color_scale_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ph_color_scale_end as text
%        str2double(get(hObject,'String')) returns contents of ph_color_scale_end as a double

end_ph_scale = str2double(get(hObject,'String'));

if (end_ph_scale == 1 && str2double(get(handles.ph_color_scale_strt,'String')) == -1)
    colormap(handles.h_phase, hsv);
else
    
    colormap(handles.h_phase, jet);
    
end

if (str2double(get(handles.ph_color_scale_strt,'String')) < end_ph_scale)
    set(handles.h_phase, 'CLim', [str2double(get(handles.ph_color_scale_strt,'String')), end_ph_scale]);
    colorbar(handles.h_phase);
end

% --- Executes during object creation, after setting all properties.
function ph_color_scale_end_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ph_color_scale_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function first_phase_edt_Callback(hObject, eventdata, handles)
% hObject    handle to first_phase_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of first_phase_edt as text
%        str2double(get(hObject,'String')) returns contents of first_phase_edt as a double

set(handles.max_ps_indic_edt, 'String', num2str(str2double(handles.first_phase_edt.String) + 360/2*max(2, str2double(handles.slice_per_step_edt.String)) - str2double(handles.step_ps_edt.String)));

% --- Executes during object creation, after setting all properties.
function first_phase_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to first_phase_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in point_ROI_batch_chck.
function point_ROI_batch_chck_Callback(hObject, eventdata, handles)
% hObject    handle to point_ROI_batch_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of point_ROI_batch_chck

% switch get(handles.point_ROI_batch_chck,'Value')
%     case 1 % point centers
%         handles.point_center_batch_plotting = 1;
%     case 0 % don`t point centers
%         handles.point_center_batch_plotting = 0;
% end


function  handles =  expmiji_phase_util(handles, expboth, str_ph, fname, predef_foldr_put_mijph, predef_filenm_put_mijph)

[foldr_put, title_put] = exp_plot_util( handles, fname, handles.h_phase, handles.plot_phase_popup, predef_foldr_put_mijph, predef_filenm_put_mijph, str_ph, 1);  % 1 for ph_flag
    %         exp_plot_util(handles, fname, hplot, hpopup, callback, foldr_put, str, ph_flag)
if expboth % export also icontr
    exp_plot_util( handles, fname, handles.plot_princ, handles.plot_diff_menu, foldr_put, [title_put, '_ictr'],'ictr', 0);  % 0 for ph_flag
end

% --- Executes on selection change in action_map_menu.
function action_map_menu_Callback(hObject, eventdata, handles)
% hObject    handle to action_map_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns action_map_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from action_map_menu
if ~get(handles.batch_chck, 'Value')% no batch
    handles.predef_foldr_put_mijph=0;
    handles.predef_filenm_put_mijph=0;
end
val = get(handles.action_map_menu,'Value');
if val~=5 % does not need phase
    num_plot = min(length(handles.phase_cell), get(handles.plot_phase_popup,'Value')-1);
    if ~(num_plot> 0 && length(handles.phase_cell) >= num_plot && ~isempty(handles.phase_cell)) % not enough phase stored ??
        im=get(handles.plot_phase, 'Children'); 
        if ~isempty(im)
            handles.phase_cell{max(num_plot,1)}= im.CData;
        else
            msgbox('no phase in plot (and in memory) !');
            global mij_obj %#ok<TLEV>
            if ~mij_obj
                [mij_obj, ~, ~] = miji_save_img(0, handles.cmap_blkredblkgrn, mij_obj, '', '', '', '', '', 1); % 1 for ph_flag
            end
            return;
        end
    end
    phase11 = handles.phase_cell{max(num_plot,1)};
    if ~(isfield(handles, 'amp_cell') && num_plot <= length(handles.amp_cell) && ~isempty(handles.amp_cell))
        im=get(handles.plot_princ, 'Children'); 
        if ~isempty(im); handles.amp_cell{max(num_plot,1)}= im.CData; end
    end
    if ~isempty(handles.amp_cell); amp = handles.amp_cell{max(num_plot,1)}; end
end

if isa(handles.fname, 'cell') % many files
    fname = handles.fname{1};
else
    fname = handles.fname;
end

fname = sprintf('%s_ctr×%d' , fname(1:end-5), 2^(handles.double_contrast+1)); % % fname end-1 to avoid naming mix if several files

switch val
    
    case 2 % save mat
        
        d=uigetdir;
        
% %         phase = handles.phase_cell{end}; % for saving
%         % handles.phase_calib = handles.phase1;
%         %
%         % handles.saved_ref = 1;
        
        save(fullfile(d,sprintf('%s_phase.mat', fname)), 'phase11', 'amp');
        fprintf('\nFile %s saved in %s\n', sprintf('%s_phase', fname), d);
        
    case 3 % save fig
        
        set(handles.undk_p1_chck, 'Value', 1);
        
        if get(handles.undk_p1_chck, 'Value')
        
            try
                get(handles.h_princ,'Children'); % error if deleted
                %  %               disp('ok')
                hf_contr = get(handles.h_princ, 'Parent');
                if strcmp(hf_contr.Name, 'I_SHG_GUI')
                    error('Gui selected')
                end
                
            catch % if it has been deleted

                hf_contr = figure('outerposition',...
                    [min(handles.screensize(3)*(1-handles.fact), handles.left_offset_fig) min(handles.screensize(4)*(1-handles.fact), handles.top_offset_fig) ...
                    handles.screensize(3)*handles.fact handles.screensize(4)*handles.fact]);
                % [left bottom width height]
                handles.h_princ = axes; % create axes in current figure
            end
        end
        
% %         handles.cmap_princ_curr = parula;
% %         num_cell = (ii - size(handles.contr00, 3))/2;
        draw_plots_ISHG( 0, 0, handles.amp_cell{end}, handles.x_cell{end}, handles.y_cell{end}, handles.h_princ, ...
            handles.xTitle_dflt, handles.yTitle_dflt, handles.Titre2, parula , 0, '', 0, 1, ...
            handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz );
        
        handles.hfigsave_graph(2) = hf_contr;
        set(handles.undk_p1_chck, 'Value', 0);
        callback_undk_p1 = get(handles.undk_p1_chck,'Callback');  hgfeval({callback_undk_p1, hObject,eventdata});
        
        set(handles.undck_plotphase, 'Value', 1);
%         str1 = get(handles.plot_phase_popup, 'String');
%         set(handles.plot_phase_popup, 'Value', length(str1));
        
        handles.undocked_phase = 1;
        handles.h_phase =  handles.h_phase_out;
        % guidata(hObject, handles);
        handles.flag_save_fig_out = 1;
        % undck_plotphase_Callback(hObject, eventdata, handles);
        % the plt_popup callback will save the fig
        plot_phase_popup_Callback(hObject, eventdata, handles);
        
    case 4 % export phase to miji
        expboth = get(handles.expboth2miji_chck, 'Value');
        str_ph = 'phmp';
        if (strcmp(handles.cmap_p_menu.String{handles.cmap_p_menu.Value}, 'B/R/B/G of fit') && isfield(handles, 'cmap_brbg_fit') && size(handles.cmap_brbg_fit,2) == 3 ) 
            str_ph = [str_ph, 'brbg'];
            if ~handles.expboth
                set(handles.expboth2miji_chck, 'Value', 0); % not for BRBG
                expboth0=expboth;  expboth = 0;
            end
        end
        if get(handles.curr_fldr_save_chck, 'Value') % checked to take current folder
            pred_fldr = handles.folder_name;
            pred_fn = sprintf('%s_ctr×%d' , 'phmp', 2^(handles.double_contrast+1));
        else
            pred_fldr = 0; pred_fn =0;
        end
        if handles.exp_phase_before_ictr
           handles =  expmiji_phase_util(handles, expboth, str_ph, fname, pred_fldr, pred_fn); % export with exp_plot_util the phase, and the ictr if asked
           if expboth % export also icontr
               hpopup= handles.plot_phase_popup;
               set(hpopup, 'Value', length(handles.plot_phase_popup.String));
                % %     callback(hpopup, 0, handles);
               hgfeval({get(hpopup,'Callback'), hpopup ,0});
           end
        else
            if expboth % export also icontr
                [foldr_put, title_put] = exp_plot_util( handles, fname, handles.plot_princ, handles.plot_diff_menu, 0, 0,'ictr', 0);  % 0 for ph_flag
                exp_plot_util( handles, fname, handles.h_phase, handles.plot_phase_popup, foldr_put, [title_put, '_phmp'], str_ph, 1);  % 1 for ph_flag
            else
                exp_plot_util( handles, fname, handles.h_phase, handles.plot_phase_popup, 0,0, str_ph, 1);  % 1 for ph_flag
            end
        end
        if (exist('expboth0', 'var') && ~handles.expboth)
            set(handles.expboth2miji_chck, 'Value', expboth0); % not for BRBG
        end
    case 5 % export pl-left to miji
        exp_plot_util( handles, fname, handles.plot_princ, handles.plot_diff_menu, 0,0, 'ictr', 0);  % 0 for ph_flag

    case 6 % plot line profile
        
        axes(handles.h_phase);
        
        title(handles.h_phase, 'Select your 2 ends of segment !'); a = ginput( 2);
        x1 = a(1); y1 = a(3); x2 = a(2); y2 = a(4);
        
        im = get(handles.h_phase, 'Children');
        
        figure(126);
        set(126, 'outerposition',...
            [min(handles.screensize(3)*(1-handles.fact), handles.left_offset_fig) min(handles.screensize(4)*(1-handles.fact), handles.top_offset_fig) ...
            handles.screensize(3)*handles.fact handles.screensize(4)*handles.fact]);
        if ~isfield(handles, 'x'); handles.x = handles.x_cell{end}; end
        if ~isfield(handles, 'y'); handles.y = handles.y_cell{end}; end
        subplot(1, 2, 1); imagesc(handles.x, handles.y,im.CData*180); title('Phase in [-180°, 180°]'); % handles.x, handles.y,
        axis image; colormap jet; hc = colorbar; title(hc, 'phi(°)');hold on;
        plot([x1,x2], [y1,y2], 'w+-', 'LineWidth', 2);
        set(gca, 'CLim',  [str2double(get(handles.ph_color_scale_strt,'String'))*180, str2double(get(handles.ph_color_scale_end,'String'))*180]);
        stpX = handles.x(2) - handles.x(1);
        stpY = handles.y(2) - handles.y(1);
        if (stpX ~= 1)
            x11 = length(handles.x)/max(handles.x)* x1;
            x22 = length(handles.x)/max(handles.x)* x2;
        else
            x11 = x1; x22=x2;
        end
        if (stpY ~= 1)
            y11 = y1*length(handles.y)/max(handles.y); 
            y22 = length(handles.y)/max(handles.y)*y2;
        else
            y11 = y1; y22=y2;
        end
        [theProfile] = improfile(im.CData*180, [x11, x22], [y11, y22]);
        subplot(1,2,2); plot(1:length(theProfile),theProfile, 'r+-');
        xlabel('Distance (PX)'); ylabel('Phase contrast (°)'); title('Profile');
        
    case 7 % plot surf phase
        
        ans1 = get(handles.h_phase, 'Children'); aa = ans1.CData;
        figure(110); 
        set(110, 'outerposition',...
            [min(handles.screensize(3)*(1-handles.fact), handles.left_offset_fig) min(handles.screensize(4)*(1-handles.fact), handles.top_offset_fig) ...
            handles.screensize(3)*handles.fact handles.screensize(4)*handles.fact]);ha = axes; surf(aa); daspect([1, 1, (max(max(aa)) - min(min(aa)))/max(size(aa, 1), size(aa, 2))]);
        set(ha, 'CLim',  [str2double(get(handles.ph_color_scale_strt,'String')), str2double(get(handles.ph_color_scale_end,'String'))]);
    
    case 8 % Correct phase tilt
        
%         ii = get(handles.plot_phase_popup,'Value')-1;
% 
%         phase = handles.phase_cell{ii};

% %         ch_detilt = menu('Choice', 'Surface fit (best)', 'Coplanar points', 'Load .mat');
        f=choosedialog;
        [unwrap_flag, ch_detilt] = f.dialog_tiltcorr();
        
        if strcmp(get(handles.plot_phase, 'Visible'), 'on') % phase current
            data11 = phase11;
             ampflag = 0;
            if unwrap_flag
                if unwrap_flag == 1 % else use previous
                    handles = f.unwrap_prepare(handles); % will unwrap
                else % use previous
                     disp('previous unwrap chosen, so it will use handles.im_unwrapped !');
                end
                data11 = handles.im_unwrapped;
            end
        else % amp current
            cc=get(handles.plot_princ, 'Children');data11 = cc.CData; 
% %             data11 = amp;
            ampflag = 1;
        end
        if ch_detilt
            switch ch_detilt
                case 1 % surface fit
                    x= [0 0 0]; y= [0 0 0]; coplanar_points = 0;

                case 2 % coplanar points
                    prompt = {'x1 (Put 0 to select them on plot, Cancel if you want to surface fit the tilt)', 'x2', 'x3', 'y1', 'y2', 'y3'};
                    dlg_title = 'Coplanar points : ';
                    num_lines = 1;
                    % Valeurs par défaut
                    def = {'0', '0', '0', '0', '0', '0'};
                    answer = inputdlg(prompt,dlg_title,num_lines,def);

                    if ~isempty(answer)
                        x(1) = str2double(answer{1});
                        x(2) = str2double(answer{2});
                        x(3) = str2double(answer{3});
                        y(1) = str2double(answer{4});
                        y(2) = str2double(answer{5});
                        y(3) = str2double(answer{6});
                        coplanar_points = 1;
                    else
                        x= [0 0 0]; y= [0 0 0]; coplanar_points = 0;
                    end
            end

            if ch_detilt ~= 3 % % not load, surface
                h = handles.h_phase;
                try, cd(handles.folder_name); end %#ok<TRYNC,NOCOM>
                [IM, nx, ny]=subplane_mp(data11,x',y', coplanar_points, h, unwrap_flag, handles.surf_fit_tilt_auto_bacth, ampflag);

                if coplanar_points
                    fprintf('New x = %f\n', nx);
                    fprintf('New y = %f\n', ny);
                end
            else % load
                [FileName,PathName, FILTERINDEX] = uigetfile({sprintf('*.%s', ext); sprintf('*.%s', ext_var); sprintf('*.%s', ext_var2); '*.*'}, 'select your file of data', fullfile(dirname, sprintf('data.%s', ext)),'MultiSelect','off');
                if ~FILTERINDEX
                    fprintf(2, 'Wrong file !\n');
                end
                matfit = load(fullfile(PathName, FileName));
                if isa( matfit, 'struct') % mat
                    str = fieldnames(matfit);  
                    matfit = getfield(matfit, str{1}); %#ok<GFLD>
                end
                IM = data11 - matfit;
            end
            
            if unwrap_flag % reput in -1 1
                IM = wrapToPi(IM*pi)/pi;
            end
            handles.phase_cell{end+1} = IM;
%             handles.amp_cell
            handles = add_data_in_list(hObject, handles, 1, 1, 1, 0, [1, Inf], [1, Inf]); % doxy, do_amp, do_err, off
            guidata(hObject, handles);
            % menu('New selection ?', 'Choose the tilted plane on img', 'Load tilt )
            add_plot_ph(hObject, handles, 'Phase map de-tilted');

        end

    case 9 % Correct phase by ref, old 10
        [handles, phase_before_corr] = corr_ph_ref_util(hObject, handles);
        h_f=figure(42);
        set(h_f,'Color', [1 1 1], 'outerposition',...
            [min(handles.screensize(3)*(1-handles.fact), handles.left_offset_fig) min(handles.screensize(4)*(1-handles.fact), handles.top_offset_fig) ...
            handles.screensize(3)*handles.fact handles.screensize(4)*handles.fact]);

        subplot(2, 2,1);
        imagesc(phase_before_corr); colormap(hsv); colorbar;caxis([-1,1]);  axis('image');
        title('Before correction')
        subplot(2, 2,2);
        imagesc(handles.phase_calib); colormap(hsv); colorbar; caxis([-1,1]); axis('image');
        title('Calib')
        subplot(2, 2,3);
        imagesc(handles.phase_cell{end}); colormap(hsv); colorbar;caxis([-1,1]);  axis('image');
        title('After correction')

        guidata(hObject, handles);
        
     case 10 % correct oscillations, old 13
            handles.phase_cell{end+1} = osc_corr_movemean(phase11);
            handles = add_data_in_list(hObject, handles, 1, 1, 1, 0, [1, Inf], [1, Inf]);
            guidata(hObject, handles);
            add_plot_ph(hObject, handles, 'phase corr oscillations');
            
     case 11 % unwrap phase, old 9
            f=choosedialog;
            handles = f.unwrap_prepare(handles);  % % unwrap_2D_MP is inside that
            guidata(hObject, handles); % safety !!!
            disp('unwrap saved in handles.im_unwrapped !');
            
            handles.phase_cell00{end+1} = handles.im_unwrapped;
            handles.phase_cell{end+1} = handles.im_unwrapped;
            handles = add_data_in_list(hObject, handles, 1, 1, 1, 0, [1, Inf], [1, Inf]);
            guidata(hObject, handles);
            add_plot_ph(hObject, handles, sprintf('phase map %d', length(handles.phase_cell)));
            colormap(handles.cmap_cubehelix); caxis('auto');
     
    case 12 % % wrap phase
        try
            im = get(handles.h_princ,'Children'); % error if deleted
            hf_contr = get(handles.h_princ, 'Parent');
            if strcmp(hf_contr.Name, 'I_SHG_GUI'); error('Gui selected'); end
            ph=im.CData;
        catch % if it has been deleted
            ph=handles.phase_cell{end};
        end
        ph1=wrapToPi(ph*pi)/pi;
        handles.phase_cell{end+1} = ph1;
        handles = add_data_in_list(hObject, handles, 1, 1, 1, 0, [1, Inf], [1, Inf]);
        guidata(hObject, handles);
        add_plot_ph(hObject, handles, 'phase wrapped');
    
    case 13 % avg curr map, old 11
         if (get(handles.expboth2miji_chck, 'Value') && ~strcmp(get(handles.h_phase, 'Visible'), 'on'))  % avg both phase and ictr
             plot_phase_popup_Callback(handles.plot_phase_popup, eventdata, handles);
         end
         handles = avg_currmap_util(hObject, eventdata, handles,0);% 0 for ask
         
         if get(handles.expboth2miji_chck, 'Value')  % avg both phase and ictr.
             ind= length(handles.plot_diff_menu.String)+0; str=get(handles.plot_diff_menu, 'String');
             if  ~strcmp(str{ind}(1:9), 'Interf. c')
                 ind= ind-1; % last is surely err
             end
            set(handles.plot_diff_menu, 'Value', ind);
            plot_diff_menu_Callback(handles.plot_diff_menu, eventdata, handles);
            avg_currmap_util(hObject, eventdata, handles,1); % % 1 for skip
         end
         
    case 14 % % crop current map
        
        [im, ph, vy, vx] = plot_crop_util(handles, [], 0, 0);
        if ph % phase
            handles.phase_cell{end+1} = im;
            do_amp=1;
        else
            handles.amp_cell{end+1} = im;
            handles.phase_cell{end+1} = handles.phase_cell{end}(vy, vx);
            do_amp=0;
        end
        handles = add_data_in_list(hObject, handles, 1, do_amp, 1, 0, vy, vx);
        guidata(hObject, handles);
        add_plot_ph(hObject, handles, sprintf('Phase map %d', length(handles.phase_cell)));
        guidata(hObject, handles);
        
    case 15 % % reg shift between X and Y lines
        if strcmp(get(handles.plot_phase, 'Visible'), 'on') % phase current
            hh = handles.plot_phase;
            phflag = 1;
        else
            hh = handles.plot_princ;
            phflag = 0;
        end
        cc=get(hh, 'Children');im0 = cc.CData;
        lim1 = []; off_shift = 0;
        disp('put break around line 4040 in func to change off_shift final if you want !');
        if ~handles.dipimage
            try, findshift; %#ok<NOCOM>
            catch ME
            end
            % % get DIPIMAGE at http://www.diplib.org/download
            if strcmp(ME.message, 'Required parameter missing') % ok
               disp('');
            else; run('C:\Program Files\DIPimage 2.9\dipstart.m'); 
            end
            handles.dipimage = 1;
            guidata(hObject, handles);
        end
        [shiftv, im] = reg_shift_advanced_func(im0, lim1, off_shift, []);
        fprintf('shiftv %.1f\n', shiftv);
        yy ='Yes save the imgs'; nn= 'No ';
        choice = questdlg('Is the disp shift correct?', ...
        'save?', ...
        yy,nn,yy);
        tt = 'single';
        switch choice 
            case yy
                addpath(fullfile(handles.genpath1, 'codes Matlab\Variety'));
                ff=fullfile(handles.folder_name, handles.fname); ff=[ff(1:end-4), 'STACKREG.tif'];
                if (~phflag && isfield(handles,'img_3D') && handles.num_images > 0 && get(handles.plot_diff_menu, 'Value')==2) % img 3D
                    data = handles.img_3D(:,round(shiftv)+1:end,:)*0;
                    tt = handles.type_im;
                    for ii =1:size(data,3)
                        [~, data(:,:,ii)] = reg_shift_advanced_func(handles.img_3D(:,:,ii), 0, 0, shiftv);
                    end
                else
                    data = im;
                end
                saveastiff(cast(data,tt), ff); 
%             case nn
        end
    case 16 % correct lines
        as=inputdlg({'# of separated lines(1,...)','thickness up','thickness down', 'dir(X=1,Y=2)'}, 'Lines to avg', 1,{'[]','[]','[]', '1'});
        % prompt,dlg_title,num_lines,defAns
        
        lines=str2num(as{1}); %#ok<ST2NM>
        if isempty(lines); return; end
        thicknesses_top=str2num(as{2}); %#ok<ST2NM>
        if length(thicknesses_top) < length(lines); thicknesses_top = [thicknesses_top, ones(1,length(lines)-length(thicknesses_top))]; end
        thicknesses_btm=str2num(as{3}); %#ok<ST2NM>
        if length(thicknesses_btm) < length(lines); thicknesses_btm = [thicknesses_btm, ones(1,length(lines)-length(thicknesses_btm))]; end
        dir=str2double(as{4}); if dir >1; disp('transpose your mat for column corr, lol'); end
        for k=1:length(lines)
%             sumline_top=0;
%             for ii=0:thicknesses_top(k)-1
%                 sumline_top = sumline_top+phase11(lines(k)-ii, :);
%             end
%             sumline_btm=0;
%             for ii=0:thicknesses_btm(k)-1
%                 sumline_btm = sumline_btm+phase11(lines(k)+ii, :);
%             end
%             sumline =sum(,1);
            phase11(max(1,lines(k)-(thicknesses_top(k)-1)):min(size(phase11,1),lines(k)+(thicknesses_btm(k)-1)), :) = ...
                (phase11(max(1,lines(k)-(thicknesses_top(k)-1)-1), :)+phase11(min(size(phase11,1),lines(k)+(thicknesses_btm(k)-1)+1),:))/2;
        end
        
        handles.phase_cell{end+1} = phase11;
        handles = add_data_in_list(hObject, handles, 1, 1, 1, 0, [1, Inf], [1, Inf]); % doxy, do_amp, do_err, off, vy, vx)
        guidata(hObject, handles);
        add_plot_ph(hObject, handles, 'phase corrlines');
        
    case 17 % real imaginary parts
        figure; ax=subplot(2,2,1); imagesc(phase11); colormap(ax,hsv); title(ax, 'Phase'); axis image; colorbar;
        subplot(2,2,2); imagesc(amp); colormap parula; title( '|\chi^{(2)}|'); axis image; colorbar;
        subplot(2,2,3); imagesc(amp.*cos(phase11*pi)); colormap parula;title( 'Re[\chi^{(2)}]'); axis image; colorbar;
        subplot(2,2,4); imagesc(amp.*sin(phase11*pi)); colormap parula;title( 'Im[\chi^{(2)}]'); axis image; colorbar;
        
    case 18 % treat further
        nm1 = 'treat';
        put_peaks = '1'; do_mijiph = '1'; do_clrwheel = '1'; saved_ref = '1'; 
        fit_batch = '1';
        fit_clrwheelBRBG  = '1';
        range_interf_keep = '1';
        prefixdf = 'raw';
        sel_files = '0';
        use_prev_str ='0';
        plot_hist_off = 0;
         % % '1', '1', '1', '0', '0', '0', '1','raw', '0'
        
        [handles, ~,put_peaks ,do_mijiph,do_clrwheel, sel_files,num_fit, prefix] = ...
prompt_batch_util(handles, '1',nm1 ,put_peaks ,do_mijiph,do_clrwheel,saved_ref ,fit_batch ,fit_clrwheelBRBG,range_interf_keep,prefixdf ,sel_files, use_prev_str);
        handles.nm_clrwheel = sprintf('clrwheel%s.fig',prefix);
        handles.nm_clrwheelBRBG = sprintf('clrwheel%s_corrref.fig',prefix);

        batch_furthertreat_core(hObject, handles, prefix, do_mijiph, put_peaks,do_clrwheel, sel_files,num_fit, plot_hist_off);
end

set(hObject,'Value', 1)

function [foldr_put, title_put] = exp_plot_util(handles, fname, hplot, hpopup, foldr_put, title_put, str, ph_flag)
global mij_obj
cont = 1;
try
    cc = get(hplot, 'Children');
    curr_map = cc.CData;
catch
% %             set(handles.plot_phase_popup, 'Value', length(handles.plot_phase_popup.String)-1)
    ind= length(hpopup.String)-ph_flag+0; hpopupstr=get(hpopup, 'String');
     if (~ph_flag && ~strcmp(hpopupstr{ind}(1:9), 'Interf. c'))
         ind= ind-1; % last is surely err
     end
    set(hpopup, 'Value', ind+ph_flag);
% %     callback(hpopup, 0, handles);
    hgfeval({get(hpopup,'Callback'), hpopup ,0});
    try
        cc = get(hplot, 'Children'); curr_map = cc.CData;
        fprintf(2, '%s was not current plot, I chose the last one !\n', str);
    catch
        cont = 0;
        fprintf(2, '\n %s is not current plot !\n', str);
    end
end
if ph_flag  % phase
    val = get(handles.cmap_p_menu, 'Value');
    if val==3; cmap = 'brbg';
    elseif val==5; cmap = 'hsv sat';
    elseif val==7; cmap = 'parula';
    elseif val==9; cmap = 'Grays';
    else % %strcmp(handles.cmap_p_menu.String{val}, 'B/R/B/G of fit'))
        cmap = handles.cmap_default1; 
    end
else
    val = get(handles.cmap1_menu, 'Value');
    if val==2; cmap = 'brbg';
    elseif val==3; cmap = 'hsv sat';
    elseif val==4; cmap = 'parula';
    elseif val==6; cmap = 'Grays';
    else % %strcmp(handles.cmap_p_menu.String{val}, 'B/R/B/G of fit'))
        cmap = colormap(hplot); 
    end
end
if cont
    [mij_obj, foldr_put, title_put] = miji_save_img(curr_map, cmap, mij_obj, handles.folder_name, fname, foldr_put, title_put, str, ph_flag); % 1 for ph_flag
% %    [mij_obj, foldr_put, title_put] = miji_save_img(curr_map, cmap, mij_obj, folder_name, fname, foldr_put, title_put, str, ph_flag) 
end
 

function handles = add_data_in_list(hObject, handles, doxy, do_amp, do_err, off, vy, vx)
if doxy
    while length(handles.x_cell) < length(handles.phase_cell)
        handles.x_cell{end+1}= handles.x_cell{end-off}(max(1,vx(1)):min(length(handles.x_cell{end-off}), vx(end))); 
        handles.y_cell{end+1}= handles.y_cell{end-off}(max(1,vy(1)):min(length(handles.y_cell{end-off}), vy(end)));
    end
end
if do_amp
    mat=handles.amp_cell{end-off};
    vyy = max(1,vy(1)):min(length(handles.y_cell{end-off}), vy(end));
    vxx = max(1,vx(1)):min(length(handles.x_cell{end-off}), vx(end));
    if length(vyy) > length(1:size(mat,1)); vyy=1:size(mat,1);end
    if length(vxx) > length(1:size(mat,2)); vxx=1:size(mat,2);end
    handles.amp_cell{end+1}= handles.amp_cell{end-off}(vyy, vxx);
end
if (do_err && ~isempty( handles.err_cell))
    ind = max(length(handles.err_cell)-off, 1);
    if (size(handles.err_cell{ind}, 1) < length(handles.y_cell{end-off}) || size(handles.err_cell{ind}, 2) < length(handles.x_cell{end-off}))
        handles.err_cell{end+1} = [];
    else
         handles.err_cell{end+1} = handles.err_cell{ind}(max(1,vy(1)):min(length(handles.y_cell{end-off}), vy(end)), max(1,vx(1)):min(length(handles.x_cell{end-off}), vx(end)));
    end
elseif isempty( handles.err_cell)
    handles.err_cell = {handles.phase_cell{1}*0};
end

    
function add_plot_ph(hObject, handles, name)

str=get(handles.plot_phase_popup, 'String');
str{end+1} = name;
set(handles.plot_phase_popup, 'String', str);
set(handles.plot_phase_popup, 'Value', length(str)); 
% % hgfeval({get(handles.plot_phase_popup,'Callback'), hObject ,0});
plot_phase_popup_Callback(hObject, 0, handles); 

% --- Executes during object creation, after setting all properties.
function action_map_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to action_map_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [handles, phase_before_corr] = corr_ph_ref_util(hObject, handles)

if ~strcmp(get(handles.plot_phase, 'Visible'), 'on'); plot_phase_popup_Callback(handles.plot_phase_popup, 0, handles); end
calib.phase = handles.phase_calib;
calib.amp = handles.amp_calib;
calib.do_realign = handles.do_realign_calib;
diff_ctr = length(handles.plot_diff_menu.String)- handles.plot_diff_menu.Value ;
a= get(handles.plot_phase, 'Children'); phase1 =  a.CData;

func_hdl = load_stack_plot_ISHG;
[calib, handles.amp_cell, handles.phase_cell, ...
handles.x, handles.y, handles.x_cell, handles.y_cell, phase_before_corr, ~] = correct_phase_ref_subtract( handles.amp_cell, ...
handles.phase_cell, handles.resx, handles.x_cell, handles.y_cell, ~get(handles.batch_chck, 'Value'), calib, phase1, diff_ctr, func_hdl.scaling_img, 0); % 0 for ask

handles.phase_calib = calib.phase;
handles.amp_calib = calib.amp;
handles.do_realign_calib = calib.do_realign;
handles.phase_cell00(end+1) = handles.phase_cell(end); % handles.phase_cell is done inside func
if ~isa(handles.amp_cell00, 'cell'); handles.amp_cell00 = {handles.amp_cell00}; end
handles.amp_cell00(end+1) = handles.amp_cell(end);

aa=get(handles.plot_diff_menu, 'String');
aa{end+1} = sprintf('Interf. contrast %d', length(handles.amp_cell));
aa{end+1} = sprintf('Err. map %d', length(handles.err_cell));
set(handles.plot_diff_menu, 'String', aa);

handles.range_sat = 0;

guidata(hObject, handles);

add_plot_ph(hObject, handles,  sprintf('phase map %d', length(handles.phase_cell)));

function handles = avg_currmap_util(hObject, eventdata, handles, skip_prompt)

 try % phase ?
    cc = get(handles.h_phase, 'Children');
    curr_map = cc.CData; ch = 1;
 catch
     try % ctr ?
        cc = get(handles.plot_princ, 'Children');
        curr_map = cc.CData; ch = 2;
     catch
         ch = 0;
     end
 end

 if ch > 0
    if (isa(handles.fname, 'cell') && get(handles.batch_chck, 'Value')) % batch files
        skip_asking_coeff_median = 1;
    else
        skip_asking_coeff_median = 0;
    end
     func_hdl = load_stack_plot_ISHG;
     [curr_map,  handles.mean_phase_meth] = func_hdl.avg_nearest_ut(curr_map, skip_asking_coeff_median, handles.coef_centre, handles.coef_cross, handles.coef_diag, handles.mean_phase_meth, handles.sigma_med, 1, skip_prompt);
     if ~isempty(curr_map) 
         switch ch
             case 1 % phase
                handles.phase_cell00{end+1} = curr_map;
                handles.phase_cell{end+1} = curr_map;
% %         handles.amp_cell00{end+1} = handles.amp_cell{end};
% %                   handles.amp_cell(end+1) = handles.amp_cell(end);
                handles =  add_data_in_list(hObject, handles, 1, 1, 1, 0, [1, Inf], [1, Inf]);
                 guidata(hObject, handles);
                add_plot_ph(hObject, handles, sprintf('phase map %d', length(handles.phase_cell)));
                
             case 2 % left
                aa=get(handles.plot_diff_menu, 'String'); val =get(handles.plot_diff_menu, 'Value');
                if ((length(aa{val})>=9 && strcmp(aa{val}(1:9), 'Interf. c')) || (val ==1 && (length(aa{end})>=9 && strcmp(aa{end}(1:9), 'Interf. c') )))
                    aa{end+1} = sprintf('Interf. contrast %d', length(handles.amp_cell));
                    if ~isa(handles.amp_cell00, 'cell'); handles.amp_cell00 = {handles.amp_cell00}; end
                    handles.amp_cell00(end+1) = {curr_map};
                    handles.amp_cell{end+1} = curr_map;
                elseif strcmp(aa{val}(1:6), 'Rel. E')  
                    handles.err_cell{end+1} = curr_map;
                    aa{end+1} = sprintf('Err. map %d', length(handles.err_cell));
                else
                    aa{end+1} = 'mapavg';
%                             figure; 
                    imagesc(handles.plot_princ, curr_map); colormap gray; axis image;colorbar;
                end
                 handles.phase_cell(end+1) = handles.phase_cell(end);
                set(handles.plot_diff_menu, 'String', aa);
                set(handles.plot_diff_menu, 'Value', length(aa));
                handles =  add_data_in_list(hObject, handles, 1, 0, 1, 0, [1, Inf], [1, Inf]);
                 guidata(hObject, handles);
                plot_diff_menu_Callback(handles.plot_diff_menu, eventdata, handles);
         end
     end
 end

% --- Executes on button press in force_incr_order_chck.
function force_incr_order_chck_Callback(hObject, eventdata, handles)
% hObject    handle to force_incr_order_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of force_incr_order_chck


% --- Executes on button press in lut_contrast_button.
function lut_contrast_button_Callback(hObject, eventdata, handles)
% hObject    handle to lut_contrast_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

state = get(handles.plot_princ, 'Visible');

switch state
    case 'on' % left plot
        imcontrast( handles.plot_princ);
    case 'off' % phase
        imcontrast( handles.plot_phase);
end

% --- Executes on button press in load_phase_button.
function load_phase_button_Callback(hObject, ~, handles)
% hObject    handle to load_phase_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.batch_chck, 'Value')
    if (~isempty(handles.phase_calib) && (sum(handles.phase_calib(isnan(handles.phase_calib))) ~= 0)); saved_ref = '-1'; else;  saved_ref = '1'; end
   handles = batch_util(hObject, handles, '0', 'treat', '1', '1', '1', saved_ref, '0', '1', '1','loaded', '0'); % standard_notload, nm1 ,put_peaks ,do_mijiph,do_clrwheel,saved_ref ,fit_batch ,fit_clrwheelBRBG,range_interf_keep,prefix ,sel_files
% % standard_notload, nm1 ,put_peaks ,do_mijiph,do_clrwheel,saved_ref ,fit_batch ,fit_clrwheelBRBG,range_interf_keep,prefix ,sel_files
else
   [handles,~] = load_phase_util(hObject, handles, 0); % 0 for ask load
end

guidata(hObject, handles);

function [handles,FILTERINDEX] = load_phase_util(hObject, handles, choice_load)

[FILTERINDEX, ~, ~, handles] = load_phase_func(handles, 'data', choice_load, 'tif');
if length(handles.phase_cell) == length(handles.x_cell); doxy=0; else; doxy=1; end
handles = add_data_in_list(hObject, handles, doxy, 0, 1, 0, [1, Inf], [1, Inf]);
handles.phase_cell00{end+1} = handles.phase_cell{end};
handles.amp_cell00=handles.amp_cell{end};

handles.x = 1:handles.resx:size(handles.phase_cell{1}, 2)*handles.resx;
handles.y =1:handles.resy:size(handles.phase_cell{1}, 1)*handles.resy;
handles.x_cell = {handles.x};
handles.y_cell = {handles.y};

guidata(hObject, handles);
if FILTERINDEX
    if ~isfield(handles, 'err_cell')
        handles.err_cell={};
    end
% %     [] = varargout(2:end);
    if ~isfield(handles, 'Titre1') % no load img used
        handles.langue = get(handles.language_menu, 'Value');
        [ handles.Titre1, handles.Titre2, handles.Titre3, handles.Titre4, handles.Titre5, handles.Titre6, handles.Counts, handles.Yaxis1, handles.Yaxis2, handles.Legendehisto ] = string_axis_ISHG( handles.langue );
        handles.title_hist1 = handles.Titre6;
        handles.yaxis_hist1 = handles.Yaxis2;
        %             [handles.xTitle_dflt, handles.yTitle_dflt, handles.phi_mat_default, handles.cmap_brbgsat, ...
        %                 handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz, ...
        %                 ~, ~, ~, ~] = load_axis_param;
    end

    set(handles.action_map_menu, 'Visible', 'on');
    set(handles.expboth2miji_chck, 'Visible', 'on');
    set(handles.analyze_panel, 'Visible', 'on');
    set( handles.hist_button, 'Enable', 'on');
    set(handles.two_d_hist_popup, 'Enable', 'on');
    str_l = {'Phase map', 'Interf. ctr'};
    pop_l = {handles.plot_phase_popup, handles.plot_diff_menu};
    for ii = 1:2
        str=get(pop_l{ii}, 'String');
        if iscell(str)
            str1 = str{1};
        else
            str1 = str;
        end
        str = {str1, sprintf('%s loaded', str_l{ii})};
        set(pop_l{ii}, 'String', str);
        set(pop_l{ii}, 'Value', length(str));
    end
    handles.Titre4_modified = handles.Titre4;
    plot_phase_popup_Callback(handles.plot_phase_popup, 0, handles);
end


% --- Executes on button press in allcontrol_show_chck.
function allcontrol_show_chck_Callback(hObject, ~, handles)
% hObject    handle to allcontrol_show_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allcontrol_show_chck

state = get(hObject,'Value');
switch state
    case 1 % checked
        set(handles.interf_pannel, 'Visible', 'on');
        set(handles.analyze_panel, 'Visible', 'on');
        set(handles.hist_range_button, 'Visible', 'on');
        set(handles.min_hist_edt, 'Visible', 'on');
        set(handles.max_hist_edt, 'Visible', 'on');
        set(handles.ampmin_txt, 'Visible', 'on');
        set(handles.ampmax_txt, 'Visible', 'on');
        set(handles.edt_pearson, 'Visible', 'on');
        set(handles.text_pearson, 'Visible', 'on');
        set(handles.fit_hist_button, 'Visible', 'on');
        set(handles.baseline_hist_edt, 'Visible', 'on'); set(handles.baseline_hist_lbl, 'Visible', 'on');
        set(handles.action_map_menu, 'Visible', 'on');
        set(handles.expboth2miji_chck, 'Visible', 'on');
        set(handles.plot_diff_menu, 'Visible', 'on');
        
    case 0 % unchecked
        set(handles.interf_pannel, 'Visible', 'off');
        set(handles.analyze_panel, 'Visible', 'off');
        set(handles.hist_range_button, 'Visible', 'off');
        set(handles.min_hist_edt, 'Visible', 'off');
        set(handles.max_hist_edt, 'Visible', 'off');
        set(handles.ampmin_txt, 'Visible', 'off');
        set(handles.ampmax_txt, 'Visible', 'off');
        set(handles.edt_pearson, 'Visible', 'off');
        set(handles.text_pearson, 'Visible', 'off');
        set(handles.fit_hist_button, 'Visible', 'off');
        set(handles.baseline_hist_edt, 'Visible', 'off'); set(handles.baseline_hist_lbl, 'Visible', 'off');
        set(handles.action_map_menu, 'Visible', 'off');
        set(handles.expboth2miji_chck, 'Visible', 'off');
        set(handles.plot_diff_menu, 'Visible', 'off');
end


% --- Executes on button press in sel_only_few_fr_per_ps.
function sel_only_few_fr_per_ps_Callback(hObject, eventdata, handles)
% hObject    handle to sel_only_few_fr_per_ps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sel_only_few_fr_per_ps


% --- Executes on button press in batch_chck.
function batch_chck_Callback(hObject, eventdata, handles)
% hObject    handle to batch_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of batch_chck


% --- Executes on button press in ind_files_chck.
function ind_files_chck_Callback(hObject, eventdata, handles)
% hObject    handle to ind_files_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ind_files_chck



function it_max_algo_edt_Callback(hObject, eventdata, handles)
% hObject    handle to it_max_algo_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of it_max_algo_edt as text
%        str2double(get(hObject,'String')) returns contents of it_max_algo_edt as a double

setappdata(0 , 'k_max', str2double(get(handles.it_max_algo_edt,'String'))); % if you put hObject you cannot call it from outside
setappdata(0 , 'param_algo_changed', 1); 

% --- Executes during object creation, after setting all properties.
function it_max_algo_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to it_max_algo_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epsilon_ph_th_edt_Callback(hObject, eventdata, handles)
% hObject    handle to epsilon_ph_th_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epsilon_ph_th_edt as text
%        str2double(get(hObject,'String')) returns contents of epsilon_ph_th_edt as a double

setappdata(0 , 'epsilon_ph_th', str2double(get(handles.epsilon_ph_th_edt,'String')));
setappdata(0 , 'param_algo_changed', 1); 

% --- Executes during object creation, after setting all properties.
function epsilon_ph_th_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsilon_ph_th_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function epsilon_tilt_th_edt_Callback(hObject, eventdata, handles)
% hObject    handle to epsilon_tilt_th_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of epsilon_tilt_th_edt as text
%        str2double(get(hObject,'String')) returns contents of epsilon_tilt_th_edt as a double

setappdata(0 , 'epsilon_tilt_th', str2double(get(handles.epsilon_tilt_th_edt,'String')));
setappdata(0 , 'param_algo_changed', 1); 

% --- Executes during object creation, after setting all properties.
function epsilon_tilt_th_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epsilon_tilt_th_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in reuse_curr_ph_chck.
function reuse_curr_ph_chck_Callback(hObject, eventdata, handles)
% hObject    handle to reuse_curr_ph_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reuse_curr_ph_chck


% --- Executes on button press in invphshft_order_chk.
function invphshft_order_chk_Callback(hObject, eventdata, handles)
% hObject    handle to invphshft_order_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of invphshft_order_chk


% --- Executes on button press in ramps_phshft_flag.
function ramps_phshft_flag_Callback(hObject, eventdata, handles)
% hObject    handle to ramps_phshft_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ramps_phshft_flag


% --- Executes on selection change in two_d_hist_popup.
function two_d_hist_popup_Callback(hObject, eventdata, handles)
% hObject    handle to two_d_hist_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns two_d_hist_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from two_d_hist_popup

ind = get(hObject,'Value')-1;
if ind > 0 % not dflt
    set(hObject,'Value', 1); % dflt
end
switch ind
    case 1 % % 2D hist
        %% hist 2D 
        handles.int_x = 1./str2double(handles.nbinsX_edt.String); % int_x is inverse of the nbins
        num_cell = get(handles.plot_phase_popup, 'Value')-1;
        [handles.hhist1, handles.hist2_xdata, handles.hist2_ydata] = hist2D_ISHG(  handles.fact,  handles.left_offset_fig, ...
        handles.top_offset_fig, handles.phase_cell{num_cell},  handles.int_x,  handles.phi_mat_default,  handles.Counts,  handles.Titre4_modified, ...
        handles.axes_font_size,  handles.xaxis_sz,  handles.yaxis_sz,  handles.title_sz,  handles.clrbr_tl_sz, handles.screensize, handles.offset_pi2, handles.hhist1 );
        set(handles.fit_hist_button, 'Visible', 'on');
        set(handles.baseline_hist_edt, 'Visible', 'on'); set(handles.baseline_hist_lbl, 'Visible', 'on');
        set(handles.edt_pearson, 'Visible', 'on');
        set(handles.text_pearson, 'Visible', 'on');
        set_histrange_vis(handles);
        
    case 2 % % polar hist wheel
       clrwheel_twod_hist_util(handles);
    otherwise % % 3 or 4
        contents = cellstr(get(hObject,'String'));
        name = contents{ind+1};
        switch name
            case handles.str_ratiof_sigma %  'ratio f sigma, shg';
                if ~isfield(handles, 'data_ratio_f')
                    addpath(fullfile(handles.genpath1,'codes Matlab\ratio f'));
    % %                 pp=pwd;
    % %                 cd(fullfile(genpath,'codes Matlab\ratio f'));
                    load('sig_f.mat','data');handles.data_ratio_f = data;
                    load('f_spline.mat','f_spline');handles.f_spline = f_spline;
                    load('f2_spline.mat','f2_spline');handles.f2_spline = f2_spline;
    % %                 cd(pp); % back
                end
                shggraph = handles.shg && length(handles.img_shg)>1;
                calcul_f_func(202, handles.f_spline, handles.f2_spline, handles.data_ratio_f, shggraph, ...
                str2double(handles.stringb1), str2double(handles.stringc1), [], ...
                str2double(handles.stringb2), str2double(handles.stringc2), [], handles.img_shg, handles.img_shg_max, handles.phase_cell{end});
            case handles.str_hist_ctr_shg % 'Corr. hist ctr/shg';
                min1=str2double(handles.min_hist_edt.String);
                max1=str2double(handles.max_hist_edt.String);
                figure;ha=axes;hist3(ha, [handles.img_shg(handles.img_shg>min1 & handles.img_shg<max1),handles.amp_cell{end}(handles.img_shg>min1 & handles.img_shg<max1)], ...
                [str2double(handles.nbinsX_edt.String), str2double(handles.nbinsY_edt.String)]);
                set(gcf, 'renderer', 'zbuffer');set(get(ha, 'child'), 'FaceColor', 'interp', 'CDataMode', 'auto');
                view(ha, [0,90]); colormap hot; colorbar% equivalent to 2D
                xlabel('SHG [a.u]'); ylabel('intef. ctr');

        end
end

 guidata(hObject, handles);


function handles = clrwheel_twod_hist_util(handles)

num_cell = get(handles.plot_phase_popup, 'Value')-1;
if (num_cell > length(handles.phase_cell) || isempty(handles.phase_cell{num_cell}))
    a=get(handles.h_phase, 'Children'); 
    handles.phase_cell{num_cell} =a.CData;
end
while (num_cell>=1 && max(max(handles.phase_cell{num_cell}))>1)
    num_cell = num_cell -1;
end
num_cell = max(1, num_cell);
ph_hist = handles.phase_cell{num_cell}; 
nb_plr = str2double(get(handles.nb_bins_clrwh_edt, 'String')); %500;% nb_plr = round(1*2*str2double(handles.nbinsX_edt.String));
norm_hist = str2double(get(handles.norm_fact_hist, 'String')); %0.7;
% %  if norm_hist <=0: norm_hist = 0;
if norm_hist > 0
    norm_hist = min(norm_hist, 1/norm_hist);
    func_hdl = fit_hist2d_funcs; div_fact = 6; % the width of the window used to remove the first peak, in order to select well the 2nd one
    % see 'div_fact' in the following for better understanding
    diff_peak_max = 1.5; % frac of pi, maximum ratio between the 2 phase peaks accepted
    [hist2_ydata,hist2_xdata]  = histcounts(ph_hist, -1:2/nb_plr:1); 
    [~, ~, hist2_ydata]=func_hdl.norm_peaks_hist(hist2_ydata, div_fact, diff_peak_max, norm_hist);
    hist2_xdata = hist2_xdata*pi;
else; ph_hist=ph_hist(ph_hist<1&ph_hist>-1)*pi; norm_hist = 0; 
end
if (strcmp(handles.cmap_p_menu.String{handles.cmap_p_menu.Value}, 'B/R/B/G of fit') && isfield(handles, 'cmap_brbg_fit') && size(handles.cmap_brbg_fit,2) == 3 ) 
    addpath(fullfile(handles.genpath1,'codes Matlab\2017-02-22.PSHGAnalysis\Sous-fonctions'));
    fcclr = 'w'; alpha_hist = 0.5; 
    if mod(length(handles.cmap_brbg_fit), 2); off = 1; else; off =0; end
    ColorbarCircle; hold on; cc  = handles.cmap_brbg_fit(1:end-off, :); handles.hclrwh=hf;
    cc1=[cc(round((length(cc)-1)*3/4)+1:end-1, :); cc(1:round((length(cc)-1)*3/4), :)];
    colormap(cc1);
    of_fit = 1; fac1 = 1.2;
else
    addpath(fullfile(handles.genpath1, 'codes Matlab\2017-02-22.PSHGAnalysis'));
    N = 1024; alpha_hist = 0.3;  scale = 1;
    [handles.hclrwh, ax] = colorWheelGenerator(1, 1, N, 0, 2, 600, scale);% % background_white, rotate_wheel, N, secondorder_clrwheel, outerRadfact, innerRadfact)
    hold on; off = 55; ft_sz = 48;
    %ff=get(ax,'parent'); %ww=ff.Position(3);htt = ff.Position(4);
    a=text(ax, N/2,1,'\pi/2', 'Interpreter', 'tex','HorizontalAlignment', 'center', 'FontSize',ft_sz); set(a, 'Position', [N/2, (1-off)+N/2*(scale-1), 0]);
    a=text(ax, 1, N/2,'-\pi', 'Interpreter', 'tex','HorizontalAlignment', 'center', 'FontSize',ft_sz);set(a , 'Position', [1-off+N/2*(scale-1), N/2, 0]);
    a=text(ax, N, N/2,'0', 'Interpreter', 'tex','HorizontalAlignment', 'center', 'FontSize',ft_sz); set(a, 'Position', [N+off+N/2*(1-scale), N/2, 0]);
    a=text(ax, N/2, N,'-\pi/2', 'Interpreter', 'tex','HorizontalAlignment', 'center', 'FontSize',ft_sz); set(a, 'Position', [N/2, N+off+N/2*(1-scale), 0]);
    fcclr = 'k';% [100 100 100]/256;%'k';   
    of_fit = 0;
end
hhh=polaraxes;
% hhh=polaraxes('position', [0.0,0.0,1,1]);%[0.25,0.25,2,2]); hhh.ActivePositionProperty = 'position';
if norm_hist
    hp = polarhistogram( 'BinEdges',hist2_xdata,'BinCounts',hist2_ydata);
else % normal
    hp = polarhistogram(ph_hist, nb_plr, 'BinLimits',[-pi,pi]); 
end
% half white
%hp = polarhistogram(ph_hist(ph_hist<0&ph_hist>=-pi) , nb_plr, 'BinLimits',[-pi,pi]);
%hold on;hp2 = polarhistogram(ph_hist(ph_hist>=0&ph_hist<=pi) , nb_plr, 'BinLimits',[-pi,pi]);
set(hp, 'FaceColor', fcclr, 'EdgeColor',...
 fcclr, 'FaceAlpha', alpha_hist, 'EdgeAlpha', alpha_hist);  %set(hhh,  'Visible', 'off');
% half white
 %set(hp2, 'FaceColor', 'w', 'EdgeColor','w', 'FaceAlpha', alpha_hist, 'EdgeAlpha', alpha_hist);  %set(hhh,  'Visible', 'off');
if of_fit
    hold on; h1=get(hhh,'Children'); max1=max(h1.BinCounts);
%             VV = fliplr(handles.fit2dstruct{end});   % % 1=right, 2=left % % handles.fit2dstruct{1}+ handles.fit2dstruct{2}
    VV = handles.fit2dstruct{end};
%             VV=[VV(round(length(VV)/4)+1:round(3*length(VV)/4)), VV(round(3*length(VV)/4)+1:length(VV)), VV(1:round(length(VV)/4))];
    polarplot( -pi:((2*pi)/(length(VV) -1)):pi, VV/max(VV)*max1, 'y', 'LineWidth', 6); 
    set(hhh,  'RLim', get(hhh,  'RLim')*fac1);
end
set(hhh,  'Visible', 'off');


% --- Executes during object creation, after setting all properties.
function two_d_hist_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to two_d_hist_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_ps_indic_edt_Callback(hObject, eventdata, handles)
% hObject    handle to max_ps_indic_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_ps_indic_edt as text
%        str2double(get(hObject,'String')) returns contents of max_ps_indic_edt as a double


% --- Executes during object creation, after setting all properties.
function max_ps_indic_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_ps_indic_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in double_contr_chck.
function double_contr_chck_Callback(hObject, eventdata, handles)
% hObject    handle to double_contr_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clickType = get(ancestor(hObject, 'figure'), 'SelectionType'); % normal
double_contr_util(hObject, handles, clickType);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over double_contr_chck.
function double_contr_chck_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to double_contr_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clickType = get(ancestor(hObject, 'figure'), 'SelectionType');
% % disp(clickType)
double_contr_util(hObject, handles, clickType);

function double_contr_util(hObject, handles, clickType)
% % home-made

double_contrast = str2double(hObject.String(end-1:end));
if isnan(double_contrast)
    double_contrast = str2double(hObject.String(end));
    n = 1;
else
    n=2;
end

if handles.contr_radio.Value == 1
    minctr = 2;
    if sum(strcmp(clickType, {'alt', 'open'}))
        % %'right click action goes here!'
        double_contrast = max(floor(double_contrast/2), minctr);
        
    else
        % normal click
        double_contrast = max(floor(double_contrast*2), minctr);
        
    end
else % % raw
    double_contrast= 1;
end


 set(hObject, 'String', [hObject.String(1:end-n), num2str(double_contrast)]);
 handles.double_contrast = max(0, floor(log(double_contrast)/log(2)-1));
%  disp( handles.double_contrast)
 
 guidata(hObject, handles);


% --- Executes on button press in norm_interframes_chk.
function norm_interframes_chk_Callback(hObject, eventdata, handles)
% hObject    handle to norm_interframes_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of norm_interframes_chk


% --- Executes on button press in norm_interframes_ctr_chk.
function norm_interframes_ctr_chk_Callback(hObject, eventdata, handles)
% hObject    handle to norm_interframes_ctr_chk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of norm_interframes_ctr_chk



function val_maxsat_edt_Callback(hObject, eventdata, handles)
% hObject    handle to val_maxsat_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of val_maxsat_edt as text
%        str2double(get(hObject,'String')) returns contents of val_maxsat_edt as a double


% --- Executes during object creation, after setting all properties.
function val_maxsat_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to val_maxsat_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function val_min_edt_Callback(hObject, eventdata, handles)
% hObject    handle to val_min_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of val_min_edt as text
%        str2double(get(hObject,'String')) returns contents of val_min_edt as a double


% --- Executes during object creation, after setting all properties.
function val_min_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to val_min_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in avg_on_ph_chck.
function avg_on_ph_chck_Callback(hObject, eventdata, handles)
% hObject    handle to avg_on_ph_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of avg_on_ph_chck


% --- Executes on button press in expboth2miji_chck.
function expboth2miji_chck_Callback(hObject, eventdata, handles)
% hObject    handle to expboth2miji_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of expboth2miji_chck



function scaling_Y_edt_Callback(hObject, eventdata, handles)
% hObject    handle to scaling_Y_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scaling_Y_edt as text
%        str2double(get(hObject,'String')) returns contents of scaling_Y_edt as a double


% --- Executes during object creation, after setting all properties.
function scaling_Y_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scaling_Y_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scaling_X_edt_Callback(hObject, eventdata, handles)
% hObject    handle to scaling_X_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scaling_X_edt as text
%        str2double(get(hObject,'String')) returns contents of scaling_X_edt as a double


% --- Executes during object creation, after setting all properties.
function scaling_X_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scaling_X_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function norm_fact_hist_Callback(hObject, eventdata, handles)
% hObject    handle to norm_fact_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of norm_fact_hist as text
%        str2double(get(hObject,'String')) returns contents of norm_fact_hist as a double
val1 = str2double(get(hObject,'String'));
if (isnan(val1) || val1 <0 || val1>1)
    set(hObject,'String',num2str(0));
end

% --- Executes during object creation, after setting all properties.
function norm_fact_hist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to norm_fact_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function baseline_hist_edt_Callback(hObject, eventdata, handles)
% hObject    handle to baseline_hist_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of baseline_hist_edt as text
%        str2double(get(hObject,'String')) returns contents of baseline_hist_edt as a double


% --- Executes during object creation, after setting all properties.
function baseline_hist_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baseline_hist_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in saving_popup.
function saving_popup_Callback(hObject, eventdata, handles)
% hObject    handle to saving_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns saving_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from saving_popup

switch get(hObject,'Value')-1
    case 1 % save diff ctr
        in_gui = 0; % 1 does not work

        save_results_stack_ISHG2( handles.plot_diff_menu, handles.plot_phase, handles.contr, in_gui, handles.undk_p1_chck, handles.screensize, handles.fact, handles.contrast, ...
            handles.cmap_redgreen, handles.xTitle_dflt, handles.yTitle_dflt, handles.left_offset_fig, handles.top_offset_fig, handles.axes_font_size, handles.xaxis_sz, handles.yaxis_sz, handles.title_sz, handles.clrbr_tl_sz,...
            handles.x, handles.y);
        % function to save the results image of inferferometric contrast in tiff images
        guidata(hObject, handles);
        return
    case 2 % save stck shg mod.
        strspec = '';
        img_3D = handles.img_3D;
    case 3 % save stck shg mod. reordered
        if length(handles.order_frame_vect) <= 1
            fprintf(2, 'do contr calc before !!\n'); return;
        end
        strspec = '_reord_increasing';
        N = numel(handles.order_frame_vect);
        order_frame = reshape(handles.order_frame_vect', 1,N);
        order_frame =unique(order_frame,'stable');

        if sum(sort(order_frame)==order_frame) == N  % % already increasing
            img_3D = handles.img_3D;
        else % % not already increasing
            img_3D = handles.img_3D(:,:,order_frame);
        end
end
if isa(handles.fname, 'cell') % many files
    fname = handles.fname{1};
else; fname = handles.fname;
end

nb_avg = str2double(get(handles.avg_img_edt, 'String'));
dir_file = uigetdir(handles.folder_name, sprintf('Choose the folder to put your 3D stack with an AVG of %d', nb_avg));

outputFileFull = fullfile(dir_file, sprintf('%s_%dAVG%s.tif', fname(1:end-4), nb_avg, strspec)); 
addpath(fullfile(handles.genpath1,'codes Matlab\Variety'));
saveastiff(cast(img_3D, handles.type_im), outputFileFull); % is a variety file function
% % if > 4GB, try options.big = true;

% %works, but not for single uint32
% for k = 1:size(handles.img_3D, 3)
%     imwrite(uint16(handles.img_3D(:, :, k)), outputFileFull, 'WriteMode', 'append',  'Compression','none');
% end


% --- Executes during object creation, after setting all properties.
function saving_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saving_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in disp_interm_plots_chck.
function disp_interm_plots_chck_Callback(hObject, eventdata, handles)
% hObject    handle to disp_interm_plots_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of disp_interm_plots_chck


% --- Executes on button press in curr_fldr_save_chck.
function curr_fldr_save_chck_Callback(hObject, eventdata, handles)
% hObject    handle to curr_fldr_save_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of curr_fldr_save_chck


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over load_stck_button.
function load_stck_button_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to load_stck_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clickType = get(ancestor(hObject, 'figure'), 'SelectionType'); % normal
handles = ld_stck3d_util(hObject, handles,clickType);

guidata(hObject, handles);



function nb_bins_clrwh_edt_Callback(hObject, eventdata, handles)
% hObject    handle to nb_bins_clrwh_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nb_bins_clrwh_edt as text
%        str2double(get(hObject,'String')) returns contents of nb_bins_clrwh_edt as a double


% --- Executes during object creation, after setting all properties.
function nb_bins_clrwh_edt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nb_bins_clrwh_edt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
