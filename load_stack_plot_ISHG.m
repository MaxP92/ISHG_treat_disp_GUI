function func_hdl = load_stack_plot_ISHG
func_hdl.avg_nearest_ut = @avg_nearest_util;
func_hdl.scaling_img = @scaling_img;
func_hdl.load_stack_plot_ISHG_func = @load_stack_plot_ISHG_func;
end

function [ img_3D, type_im, num_images,img_shg, ...
     h, slice_per_ps_step, start_phase, diff_phase,...
     x_phase, nonequal_spacing, x_center, y_center, disp_vect, contr_mode] =...
    load_stack_plot_ISHG_func( fname, shg, undocked_fig, h, init_phase, first_stp, snd_step, slice_per_ps_step, ...
    contr_mode, avg_img, near_avg_ini, norm_ini, skip_asking_coeff_median, point_center_batch_plotting, nImage,mImage, ...
    num_images, force_incr_order, axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, screensize, fact, ...
    left_offset_fig, top_offset_fig, select_only_few_frame_per_ps, inv_order, normalize_interframes, scaling)
% [ img_3D, num_images, xTitle_dflt, yTitle_dflt, phi_mat_default, ...
% img_shg, screensize, fact, left_offset_fig, top_offset_fig, ...
% axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, ...
% cmap_redgreen, cmap_blkredblkgrn, cmap_brbgsat, cmap_hsvsat, h, slice_per_ps_step, start_phase, diff_phase, x_phase, nonequal_spacing] =...
% load_stack_plot_ISHG( fname, shg, undocked_fig, h, init_phase, second_phase, third_phase, fourth_phase, snd_step, contr_mode, avg_img, near_avg_ini, norm_ini )
%
% 2015-10 edited by Maxime PINSARD
%
%  To load the 3D image stack and plot the first image

% read the image stack

perturbations = 0;
norm_add_then_mult = 0;  % 1 = two steps
modX = 1; modY = 1; % % modulate along X and/or Y
eta = 0.5; % raio mult vs addition of perturb
freq_sinY = 10; freq_sinX = 5;
coef_centred=1; coef_crossd=1; coef_diagd=1; mean_methd=0; sigma_medd=0; real_sigmad=1; % dflt

if near_avg_ini
    near_avg_ini = 2; % !!
    % performed before if 1, after if 2!
end

% warning('OFF', 'all'); % because of Tiff warning
diff_phase = first_stp;
if isa(fname, 'cell') ; fname0 = fname{1}; else; fname0 =fname; end
InfoImage=imfinfo(fname0); bt = InfoImage.BitDepth;

if bt > 16
    fmt = sprintf('int%d', bt);
else
    fmt = sprintf('uint%d', bt);
end

img_3D=zeros(nImage, mImage, num_images, fmt);

if isa(fname, 'cell') % many files
    for i=1:num_images
        img_3D(:,:,i)=imread(fname{i});
    end
else
    warning('OFF', 'all'); % because of Tiff warning
    TifLink = Tiff(fname, 'r');
    for i=1:num_images
        TifLink.setDirectory(i);
        img_3D(:,:,i)=TifLink.read();
        % %    img_3D(:,:,i)= img_3D(:,:,i) + 30000; % !!!
        % %    fprintf(2, '1000 added : load_stack l.30\n');
    end
    TifLink.close();
    warning('ON', 'all');  % because of Tiff warning
end
type_im = class(img_3D);
img_3D = single(img_3D);
% real_sigma = 0;

%% scaling

if (~isempty(scaling ) && (scaling(1)>1 || scaling(2)>1))
% %     rg = (scaling(1)>1)+(scaling(1)>1);
    img_3D00 = img_3D;
    img_3D = zeros(size(img_3D00)./[scaling(2),scaling(1),  1]);
    for ii=1:size(img_3D,3)
        im_temp = img_3D00(:,:,ii);
        im_temp = scaling_img(im_temp, scaling);
        img_3D(:,:,ii) = im_temp;
    end
    clear img_3D00 im_temp;
% %     a1=reshape(a, 2,4*6/2 );
% % a11=median(a1, 1);
% % step = 1.0;aaa=cell2mat(arrayfun(@(i,j) mean(a(:,i:j),2),1:step:length(a)-step, 1+step:step:length(a),'un',0)); %aaa=aaa(:,1:2:end);
end

%% Avg whole image
if avg_img > 1
    [img_3D, num_images] = avg_whole_img_ISHG (avg_img, num_images, img_3D, select_only_few_frame_per_ps, fname, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, h);
    % function for the average on the whole images
end

start_phase = init_phase;
nonequal_spacing = 0; % default 
contr_chosen = 0; % default 

%% Avg neighbors 
if near_avg_ini == 1 % % Avg before (optional) norm. inter-frames
    [im1, ~] = avg_nearest_util(img_3D, skip_asking_coeff_median, coef_centred, coef_crossd, coef_diagd, mean_methd, sigma_medd, real_sigmad, 0);
    if ~isempty(im1); img_3D = im1; end % % empty if user canceled
end

%% normalization values in range (16 bits or other)

if norm_ini
    A = min(min(min(img_3D))); B = max(max(max(img_3D)));
    for i=1:size(img_3D, 3)
        newMax = 2^15/2; newMin = -2^15/2;
        % %         A = mean(mean(im00(:,:,i))); B = sqrt(var(var(im00(:,:,i))));
        % %         A=(max(max(im00(:,:,i))) + min(min(im00(:,:,i))))/2; B=(max(max(im00(:,:,i))) - min(min(im00(:,:,i))))/2;
        
        img_3D(:,:,i) = (img_3D(:,:,i)-A).*(newMax - newMin)/(B-A) + newMin;
        
        %         A = 250; B = A/2; % sigmoid case
        %         im00(:,:,i) = (newMax - newMin).*1./(1+exp(-((im00(:,:,i)-B)/A))) + newMin;
    end
    
end


if (perturbations && (modX || modY))
    
    fprintf(2, 'pertubations added !!\n');
    img_3D = perturbations_spatial(img_3D, eta, modX, modY, freq_sinY, freq_sinX);
end

if (undocked_fig || point_center_batch_plotting)
    
    try
        get(h,'Children'); % error if deleted
        if point_center_batch_plotting
            err;
        end
    catch % if it has been deleted
        figure('outerposition',...
            [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
            screensize(3)*fact screensize(4)*fact]);
        % [left bottom width height]
        
        h = axes; % create axes in current figure
    end
    axes_font_size = 16;
    xaxis_sz = 16;
    yaxis_sz = 16;
    title_sz=20;
    clrbr_tl_sz = 16;
end

if (h == 0 || ~isvalid(h))
    figure;h = axes; % create axes in current figure    
end

draw_plots_ISHG( 0, 0, img_3D(:,:,1), 1:size(img_3D(:,:,1), 2), 1:size(img_3D(:,:,1), 1), h, ...
    'X (pixels)', 'Y (pixels)', 'X image of the stack', 'gray', 0, 'SHG int. (a. u.)', 0, 1, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );
%  To plot the figures of the ISHG program : hist or image figure
% see the help of this function for entry vars

% Si l'image SHG de la zone d'ou provient le stack est demandée, on l'ouvre ici

if shg==1
    fname2 = uigetfile('.tif','sel. SHG file','shg.tif');
    img_shg = single(imread(fname2));
else
    img_shg = 0; % init for nothing
end


%% normalize interframes, to be insensitive to constant modulation on every interferograms
if normalize_interframes
    img_3D = normalize_interframes_func(img_3D, near_avg_ini, norm_add_then_mult, 0);
end

%% Avg neighbors 

if near_avg_ini == 2 % % Avg after (optional) norm. inter-frames
    [im1, ~] = avg_nearest_util(img_3D, skip_asking_coeff_median, coef_centred, coef_crossd, coef_diagd, mean_methd, sigma_medd, real_sigmad, 0);
    if ~isempty(im1); img_3D = im1; end % % empty if user canceled
end

replot = 1;
if replot
    draw_plots_ISHG( 0, 0, img_3D(:,:,1), 1:size(img_3D(:,:,1), 2), 1:size(img_3D(:,:,1), 1), h, ...
    'X (pixels)', 'Y (pixels)', 'X image of the stack', 'gray', 0, 'SHG int. (a. u.)', 0, 1, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );
end

%% point center for batch plotting

if point_center_batch_plotting
    axes(h); title(h, 'Click on you pattern center'); [x_center,y_center] = ginput(1); 
else
    x_center = 0; y_center=0;
end

%% vector of phase shift

if (isnan(slice_per_ps_step) || snd_step~=round(snd_step) || first_stp~=round(first_stp)) % not integer steps, so this number is not specified
    if first_stp ~= snd_step
        nonequal_spacing = 1; 
    else
    	start_phase = init_phase; diff_phase = first_stp;
        if contr_mode == 1
            slice_per_ps_step = round(num_images*first_stp/180); % useless ??
            if mod(diff_phase, 180) % contr impossible !
                msgbox('You chose contrast mode, but no 180deg spacing could be found !');
                contr_mode = 0;
            end
        else
            slice_per_ps_step = 1;
        end
    end
end

disp_vect = 1:num_images; disp_vect00 = disp_vect; 

if force_incr_order % force recognizing increasing order, so the logic is different : slice_per_ps_step is considered not in a row (pas à la suite)
    switch slice_per_ps_step
        case 1 % just increasing
            contr_mode = 0; % no contr mode, go on
% %         case 2 % 0, 180 but not in a row
% %              % nothing, go on
        case 0 
            % nothing, go on
        otherwise % % case 3 % 0, 180, 360 but not in a row for instance
            if contr_mode
                x_phase = start_phase:diff_phase:start_phase+(round(slice_per_ps_step*num_images/(slice_per_ps_step-1))-1)*diff_phase;
                if (mod(num_images,2) && mod(num_images, 5)) % % odd
                    num_images = num_images -1;
                    img_3D=img_3D(:, :, 1:end-1);
                end
                img_3D00 = img_3D;
                try
                    switch num_images==360/diff_phase*slice_per_ps_step/2
                        case 1 % normal
                            for ii = 1:slice_per_ps_step:num_images
                                for k=0:(slice_per_ps_step-1)
                                    if ~mod(num_images, slice_per_ps_step)% normal
                                        ind = (ii+(slice_per_ps_step-1))/slice_per_ps_step+k*num_images/slice_per_ps_step;
                                    else; ind = ii + k*180/diff_phase;
                                    end
                                    img_3D(:,:,ii+k) = img_3D00(:,:,ind);
                                    disp_vect(ii+k) = disp_vect00(ind);
                                end
                            end
                        case 0 % anormal
                            disp_vect = 1:360/diff_phase*slice_per_ps_step/2;
                            img_3D= zeros(size(img_3D00,1),size(img_3D00,2),360/diff_phase*slice_per_ps_step/2);
                            new_steperps = floor(2*num_images/(360/diff_phase));
                            indices = 1:num_images;ct_newarray=1;
                            ct=1;ii=1;
                            while ii <= num_images
                                for k=0:min((slice_per_ps_step-1), new_steperps-1) 
                                    ind = ct + k*180/diff_phase;

                                    if (ind>num_images || ct_newarray>360/diff_phase*slice_per_ps_step/2)
                                        if (isempty(indices(indices == k+ii))||(k<=1)); ii=num_images;break; else; continue; end
                                    end
                                    img_3D(:,:,ct_newarray) = img_3D00(:,:,ind);
                                    disp_vect(ct_newarray) = disp_vect00(ind);
                                    indices(indices == k+ii) = [];
        %                             disp([ind, ct_newarray])  % (ct-1)*slice_per_ps_step+k
                                    ct_newarray=ct_newarray+1;
                                end
                                ct=ct+1; ii=ii+1;
        %                         disp('--')
                            end
                            
                            img_3D=img_3D(:,:,1:ct_newarray-1-mod(ct_newarray-1,2));
                            disp_vect=disp_vect(1:ct_newarray-1-mod(ct_newarray-1,2));
                            
                    end
                catch ME
                    fprintf('indexes are %d %d ind %d\n\n', ii,k, (ii+(slice_per_ps_step-1))/slice_per_ps_step+k*num_images/slice_per_ps_step);
                    rethrow(ME);
                end
%                 slice_per_ps_step = -1 ; % pass the following

            end
    end
end

if slice_per_ps_step == 0 % special order

    diff_phase = 0;

%     slice_per_ps_step = -1;
    nonequal_spacing = 1;
    prompt = {'Contrast mode (= 1, 0 for raw)', 'Total number of frames if raw, nb of contrast frame if contrast'};
    dlg_title = 'Number of frames';
    deflt = {'1', '3'};
    num_lines = 1;
    answer = inputdlg(prompt,dlg_title,num_lines,deflt);
    contr_chosen = str2double(answer{1});
    nb_frame = str2double(answer{2});

    if contr_chosen % contrast with custom phases
%         prompt = {'Number of CONTRAST frames'};
%         dlg_title = 'Nb';
%         deflt = {'3'};
%         num_lines = 1;
%         answer = inputdlg(prompt,dlg_title,num_lines,deflt);
%         nb_frame = str2double(answer{1});

        x_phase = zeros(1, nb_frame);
        slice_per_ps_step = zeros(2, nb_frame); % rubbish var, sorry !!
        dfltvect = {{'2', '1', '0'}, {'3', '2', '180'}, {'5', '4', '15'}};

        for compt = 1:nb_frame
            prompt = {'Frame # ', ' - Frame #', 'Phase shift (degree)'};
            dlg_title = sprintf('Soustraction of frame numb. %d (indicate number in stack)', compt);
            deflt = dfltvect{compt};
            num_lines = 1;
            answer = inputdlg(prompt,dlg_title,num_lines,deflt);
            slice_per_ps_step(1,compt) = str2double(answer{1});  % rubbish var, sorry !!
            slice_per_ps_step(2, compt) = str2double(answer{2});  % rubbish var, sorry !!
            x_phase(compt) = str2double(answer{3});

        end
        [x_phase , I] = sort(x_phase); % increasing order
        slice_per_ps_step = slice_per_ps_step(:, I);

    else
        slice_per_ps_step = 1; 
    end
end

if slice_per_ps_step > 0 % non-special order

    if slice_per_ps_step == 1 % increasing order

        if first_stp == snd_step % equal spacing

            fprintf(2, 'Warning, regarding the first phases you gave, we are considering that the spacing in phase shift is regular over all the stack\n');
            diff_phase = snd_step;
            if contr_mode == 1 % contrast chosen
                x_phase = start_phase:diff_phase:start_phase+(round(num_images/2)-1)*diff_phase;
            else % raw
                x_phase = start_phase:diff_phase:start_phase+(num_images-1)*diff_phase;
            end
            
        else % irregular spacing
            diff_phase = 180/round(num_images/2); % to fit with line 105 of contrast_by_subtract
            nonequal_spacing = 1;

        end

        slice_per_ps_step = 1; % increasing order
        msg = 'You chose contrast mode : this assumes you have 180-shifted (only) images in your stack, but all is in increasing order. Otherwise, reload stack with "Use Raw" option.';

        if contr_mode == 1 % contrast chosen
            msgbox(msg);
        end

    elseif slice_per_ps_step > 1 % not strictly increasing 


            if first_stp ~= snd_step % unequal
                nonequal_spacing = 1;
                diff_phase = 0;
            end


        if ~nonequal_spacing % equal
            if contr_mode == 1 % contrast chosen
%                 if slice_per_ps_step == 2 % 0, 180, phi 
%                 x_phase = start_phase:diff_phase:start_phase+(round(num_images*(slice_per_ps_step-1)/(slice_per_ps_step))-1)*diff_phase;
                  x_phase0 = start_phase:diff_phase:start_phase+(round(num_images/slice_per_ps_step)-1)*diff_phase;
                  x_phase = x_phase0; % init

                  if slice_per_ps_step > 2
                      for k=1:slice_per_ps_step-2
                          x_phase = [x_phase (x_phase0+k*180)]; %#ok<AGROW> % DON'T USE (end+1)

                      end
                  end
%                 elseif slice_per_ps_step == 3 % 0, 180, 360, phi
%                     x_phase = start_phase:diff_phase:start_phase+(round(2*num_images/3)-1)*diff_phase;
%                 elseif slice_per_ps_step == 4 % 0, 180, 360, 540, phi 
%                     x_phase = start_phase:diff_phase:start_phase+(round(3*num_images/3)-1)*diff_phase;
%                 end
            else % raw
                x_phase = start_phase:diff_phase:start_phase+(num_images-1)*diff_phase;
            end
        end

    else % strange order

        slice_per_ps_step = -1;

        msgbox('Order of frame not increasing in phase shift and not in contrast mode : re-order your stack manually in increasing order or phi, phi+180 or phi, phi+180, phi+360');
        return;
    end

    %% special order with raw - to be verified

    if (slice_per_ps_step == 0 && contr_chosen )  % special order

        nb_tot = nb_frame;
        x_phase = zeros(1, nb_tot);
        nb_tot0=nb_tot;
        nb_tot = nb_tot0/ceil(nb_tot0/10);

        for compt = 1:ceil(nb_tot0/10) % because displaying more than 10 lines is impossible

            title_prompt = cell(1, nb_tot);
            deflt_vect = cell(1, nb_tot);
            for k = 1:nb_tot
                title_prompt{k} = sprintf('Phase of frame # %d', k+nb_tot*(compt-1));
                deflt_vect{k} = sprintf('%d', (k+nb_tot*(compt-1)-1)*30);
            end

            prompt = title_prompt;
            dlg_title = 'Phase of frames';
            deflt = deflt_vect;
            num_lines = 1;
            answer = inputdlg(prompt,dlg_title,num_lines,deflt);

            for k = 1:nb_tot
                x_phase(k+nb_tot*(compt-1)) = str2double(answer{k});

                if contr_mode == 1
                    if slice_per_ps_step == 3 % 0, 180, 360, phi ...
                        x_phase(k+nb_tot*(compt-1)+nb_tot) = str2double(answer{k})+180;
                    end

                else % raw
                    if (slice_per_ps_step == 3 || slice_per_ps_step == 2) % 0, 180, phi ...
                        x_phase(k+nb_tot*(compt-1)+nb_tot) = str2double(answer{k})+180;
                        if slice_per_ps_step == 3 % 0, 180, 360, phi ...
                            x_phase(k+nb_tot*(compt-1)+2*nb_tot) = str2double(answer{k})+360;
                        end
                    end
                end
            end

        end

    end
end

if nonequal_spacing

    if contr_mode == 1 % contrast chosen

        uiwait(msgbox('You chose contrast mode and increasing order : assuming that the second part of the stack is 180-shifted images. If no, stop and re-order your stack correctly.'));
        if slice_per_ps_step == 3 % 0, 180, 360, phi ...
            nb_tot = round(num_images/3);
            msgbox('Indicate only the principal phases (not the phi + 180, phi + 360)');
        elseif slice_per_ps_step == 2 % 0, 180, phi ...
            nb_tot = round(num_images/2);
            msgbox('Indicate only the principal phases (not the phi + 180)');
        elseif slice_per_ps_step == 1 % increasing order
            nb_tot = round(num_images/2);
        end
    else % raw
        nb_tot = num_images;
    end
    x_phase = zeros(1, nb_tot);
    nb_tot0=nb_tot;
    nb_tot = round(nb_tot0/ceil(nb_tot0/10));

    for compt = 1:ceil(nb_tot0/10) % because displaying more than 10 lines is impossible

        title_prompt = cell(1, nb_tot);
        deflt_vect = cell(1, nb_tot);
        for k = 1:nb_tot
            title_prompt{k} = sprintf('Phase of frame # %d', k+nb_tot*(compt-1));
            deflt_vect{k} = sprintf('%d', (k+nb_tot*(compt-1)-1)*30);
        end

        prompt = title_prompt;
        dlg_title = 'Phase of frames';
        deflt = deflt_vect;
        num_lines = 1;
        answer = inputdlg(prompt,dlg_title,num_lines,deflt);

        for k = 1:nb_tot
            x_phase(k+nb_tot*(compt-1)) = str2double(answer{k});

            if contr_mode == 1
                if slice_per_ps_step == 3 % 0, 180, 360, phi ...
                    x_phase(k+nb_tot*(compt-1)+nb_tot) = str2double(answer{k})+180;
                end

            else % raw
                if (slice_per_ps_step == 3 || slice_per_ps_step == 2) % 0, 180, phi ...
                    x_phase(k+nb_tot*(compt-1)+nb_tot) = str2double(answer{k})+180;
                    if slice_per_ps_step == 3 % 0, 180, 360, phi ...
                        x_phase(k+nb_tot*(compt-1)+2*nb_tot) = str2double(answer{k})+360;
                    end
                end
            end
        end

    end

%     else % manual order specifying

    %         prompt = 'Number of frames';
    %         dlg_title = 'Number of frames';
    %         deflt = '5';
    %         num_lines = 1;
    %         answer = inputdlg(prompt,dlg_title,num_lines,deflt);
    %
    %         nb_tot = str2double(answer{1});
    %
    %         prompt = {'Phases'};
    %         dlg_title = 'Phase of frames';
    %         deflt = deflt_vect;
    %         num_lines = 1;
    %         answer = inputdlg(prompt,dlg_title,num_lines,deflt);
    %
    %         for k = 1:nb_tot
    %             x_phase(k+nb_tot*(compt-1)) = str2double(answer{k});
    %             % to be implemented !!!
    %
    %         end

end

if exist('new_steperps', 'var')
    if new_steperps > 2
        x_phase = x_phase(1:(ct_newarray-1-mod(ct_newarray-1,2))*new_steperps/2);
    else
        x_phase = x_phase(1:(ct_newarray-1-mod(ct_newarray-1,2))/new_steperps);
        if new_steperps <2
            contr_mode = 0;
            msgbox('You chose contrast mode, but not enough 180deg spacing could be found !');
        end
    end
    fprintf(2, 'warning: missing frames, I chose a lower slice per steps:%d\n', new_steperps);
    contr_mode=[contr_mode,new_steperps];
end

if inv_order
% %     x_phase = x_phase(end:-1:1);
    img_3D = img_3D(:,:,end:-1:1);
end

if (length(disp_vect00) ~= length(disp_vect) || sum(disp_vect ~= disp_vect00)<length(disp_vect))
    disp('disp_vect00'); disp(disp_vect); 
end
disp('x_phase00'); disp(x_phase)

end

function [img_3D, mean_meth] = avg_nearest_util(img_3D, skip_asking_coeff_median, coef_centred, coef_crossd, coef_diagd, mean_methd, sigma_medd, real_sigmad, skip_prompt)

def = { num2str(mean_methd), num2str(sigma_medd), num2str(real_sigmad), '[1, 1, 1]', '[3, 3]', '0'};
if skip_prompt
    answer2 = def;
else
    prompt2 = {'Method for nearest neighbors : 0 for median, 1 for mean ', ...
        'If median, method to calculate it : 0 for no weight, -1 for weighted LONG, radius>0 for Joint-Histogram Weighted Median Filter exponential',...
        'If Joint-Histogram Weighted Median Filter, standard deviation of the Gaussian kernel (medianWMF 3, 3 is nice usually)', ...
        'if weighted:[Coeff centre, Coeff cross, Coeff diag.]', ...
           'If median, [sz along X, sz along Y]', 'Circular avg (1 for mult. of pi, 2 for already in rad, 3 for already in °), (for phase, if goes around -pi or pi)'};
    dlg_title2 = 'Method';
    num_lines2 = 1;  answer2 = inputdlg(prompt2,dlg_title2,num_lines2,def);
% Les réponses en caractères sont converties en chiffres qui sont enregistrés dans des variables.
    if isempty(answer2)
        warning('user canceled')
        img_3D = []; return;
    end
end
mean_meth = str2double(cell2mat(answer2(1)));
sigma_med = str2double(cell2mat(answer2(2)));
real_sigma = str2double(cell2mat(answer2(3)));
answer3 = answer2{4};
answer_med_sz= answer2{5};
circ_avg= str2double(answer2{6});

if skip_prompt
    skip_asking_coeff_median = 1;
    coef_centre = coef_centred;
    coef_cross = coef_crossd;
    coef_diag = coef_diagd;
else
    if (sigma_med == -1 || mean_meth == 1) % median weighted
% %         prompt3 = {'Coeff centre', 'Coeff cross', 'Coeff diag.'};
% %         dlg_title3 = 'Coeff weight';
% %         num_lines3 = 1;  def = {num2str(coef_centred), num2str(coef_crossd), num2str(coef_diagd)};
% %      answer3 = inputdlg(prompt3,dlg_title3,num_lines3,def);
% %         % Les réponses en caractères sont converties en chiffres qui sont enregistrés dans des variables.
        if isempty(answer3)
    % %         answer3 = def;
            warning('user canceled')
            img_3D = []; return;
        end
%         coef_centre = str2double(cell2mat(answer3(1)));
%         coef_cross = str2double(cell2mat(answer3(2)));
%         coef_diag = str2double(cell2mat(answer3(3)));
        a=num2cell(str2num(answer3)); [coef_centre, coef_cross, coef_diag]=deal(a{:}); %#ok<ST2NM>
%         [coef_centre, coef_cross, coef_diag]
    else
        coef_centre=0; coef_cross=0; coef_diag=0;
    end
end
img_3D = average_nearest( img_3D, coef_centre, coef_cross, coef_diag , mean_meth, sigma_med, real_sigma, skip_asking_coeff_median, answer_med_sz, circ_avg);

end

function im_temp = scaling_img(im_temp, scaling)
% binning X and/or Y
    for kk =1:2
        if scaling(kk) <= 1
            if kk == 1; continue; % go to Y
            else; break;
            end
        end
        if kk == 2; im_temp = im_temp'; end % % y
        im_temp = squeeze(mean(reshape(im_temp,size(im_temp,1),scaling(kk),size(im_temp,2)/scaling(kk)),2));
        if kk == 2; im_temp = im_temp'; end % % y
    end
end


