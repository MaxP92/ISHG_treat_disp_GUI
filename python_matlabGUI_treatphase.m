function python_matlabGUI_treatphase(incr_ordr, nb_slice_per_step, ctr_mult, foldername)
%
% 
guifig=findobj(allchild(groot), 'flat', 'type', 'figure','Name', 'I_SHG_GUI');

slice_per_step_edt=findall(guifig, 'Tag','slice_per_step_edt');
slice_per_step_edt.String=num2str(nb_slice_per_step);
force_incr_order_chck=findall(guifig, 'Tag','force_incr_order_chck');
force_incr_order_chck.Value=incr_ordr; % num2str(
double_contr_chck=findall(guifig, 'Tag','double_contr_chck');
double_contr_chck.String(end)=num2str(ctr_mult);
curr_fldr_save_chck = findall(guifig, 'Tag','curr_fldr_save_chck');
set(curr_fldr_save_chck,'Value',1);  % take previous folder
expboth2miji_chck= findall(guifig, 'Tag','expboth2miji_chck');
set(expboth2miji_chck,'Value',1);  % save also ictr
action_map_menu = findall(guifig, 'Tag','action_map_menu');
set(action_map_menu, 'Value', 4); % save tiff miji

histmode_menu=findall(guifig, 'Tag', 'histmode_menu');
set(histmode_menu, 'Value', 3);
disp_interm_plots_chck=findall(guifig, 'Tag', 'disp_interm_plots_chck');
set(disp_interm_plots_chck, 'Value', 0);

choose_file_tiff=findall(guifig, 'Tag','choose_file_tiff');
set(choose_file_tiff,'String', ['&', foldername]);
choose_file_tiff.Callback(choose_file_tiff,[]); % load with string of button
histmode_menu=findall(guifig, 'Tag','histmode_menu');
histmode_menu.Callback(histmode_menu,[]); % ctr mode

load_stck_button=findall(guifig, 'Tag','load_stck_button');
load_stck_button.Callback(load_stck_button,[]); % load with string of button
int_contrast_button=findall(guifig, 'Tag','int_contrast_button');
int_contrast_button.Callback(int_contrast_button,[]); % calc contr
relative_phase_button=findall(guifig, 'Tag','relative_phase_button');
relative_phase_button.Callback(relative_phase_button,[]); % calc phase + ictr + err + hist
action_map_menu.Callback(action_map_menu,[]);

end

