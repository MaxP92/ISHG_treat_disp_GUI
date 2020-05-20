function [ x, y, contr, rect, img_shg, x_phase, img_3D, mask_val_sat,disp_vect] = contrast_by_subtract_img_ISHG( img_3D, nb_images, filt, crop, resx, resy, ...
    shg, img_shg, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig, slice_per_ps_step, contr_mode, diff_phase, ...
    x_phase, use_range, avg_median, double_contrast, skip_asking_coeff_median, force_incr_order, disp_vect00, normalize_interframes, sat_vals)
% [ x, y, contr, rect, img_shg, x_phase  ] = contrast_by_subtract_img_ISHG( img_3D, slice_per_ps_step_images, filt, crop, resx, resy, ...
% shg, img_shg, undocked_fig, h, screensize, fact, left_offset_fig,
% top_offset_fig, slice_per_ps_step, contr_mode, diff_phase, x_phase, use_range )
%
% 2015.10 edited by Maxime PINSARD
%
%   To calculate contrast by subtraction of images
% it returns 'contr', a 3D array with the data
% contr is nbp1 x nbp2 x 12 ; img3D can be nbp1 x nbp2 x 18 or 36 (15deg or
% 30deg of diff_phase)

rect = 0; % init for nothing
method_crop = 2; % with rect

if (~contr_mode &&  force_incr_order)% mode raw, and increasing order
    slice_per_ps_step = 1;
end

if size(slice_per_ps_step, 1) > 1 % custom order
    diff_var = slice_per_ps_step;
    slice_per_ps_step = -1;
    
    contr = zeros(size(img_3D, 1), size(img_3D, 2), size(diff_var, 2));
    for k=1:size(diff_var, 2)
        contr(:, :, k) = img_3D(:, :, diff_var(1, k)) - img_3D(:, :, diff_var(2, k));
    end
    
end

if crop
    [ x, y, img_3D, rect ] = crop_img_ISHG(img_3D, method_crop, resx, resy, undocked_fig, h, ...
        screensize, fact, left_offset_fig, top_offset_fig );
else
    % if no crop
    x = linspace(0,size(img_3D,2)*resx,size(img_3D,2));
    y = linspace(0,size(img_3D,1)*resy,size(img_3D,1));
end

if (slice_per_ps_step > 1) % =3 --> classic order : 0, 180, 360, phi, phi+180 ...
    % = 2 --> order : 0, 180, phi, phi+180 ...
    % 4 : 0, 180, 360, 540, phi, phi+180 ..
    
%     slice_per_ps_step = -slice_per_ps_step + 4;

if size(img_3D,3) < 360/diff_phase*slice_per_ps_step/2 % % frames missing
    fprintf(2, 'missing frames ! \n');
    pr= cell(1,slice_per_ps_step);def = pr;
%     pr{1} = sprintf('phi(°)');
    for ii=1:slice_per_ps_step
        pr{ii} = sprintf('phi(°)+%d', (ii-1)*180); %,'phi(°)'
        def{ii} = num2str(1:floor((size(img_3D,3)/slice_per_ps_step))); %sprintf('phi(°)+%d', (ii-1)*180);
    end
    answer = inputdlg(pr,'Frames in pshft files',1,def);
    vspec= cell(1,slice_per_ps_step);
    for ii=1:slice_per_ps_step
        vspec{ii} = str2num(answer{ii}); %#ok<ST2NM>
    end
    if length(cell2mat(vspec)) > size(img_3D,3)
        fprintf(2, 'There are less frames than indicated ! \n');
    end
    for ii=1:slice_per_ps_step
        vspec{ii} = vspec{ii}(1:min(length(vspec{ii}), floor(size(img_3D,3)/slice_per_ps_step))); 
    end
%     if slice_per_ps_step>2 ;fact = 2/slice_per_ps_step;
%     else;fact = 1/slice_per_ps_step;
%     end
    x_phase = 0:diff_phase:(360-diff_phase);
    x_phase = x_phase(vspec{1});%1:length(cell2mat(v))*fact);
    img_3D = img_3D(:,:,1:length(cell2mat(vspec)) );
end

    V=cell(1, slice_per_ps_step);
    disp_vect = cell(1, slice_per_ps_step);
    for k = 1:slice_per_ps_step
        V{k} = img_3D(:, :, k:slice_per_ps_step:end); % 1st third : 0, 15, 30 ...
        disp_vect{k}  = disp_vect00(k:slice_per_ps_step:size(img_3D,3));
     % 2nd third : 180, 195, 210 ...
    end
    vect_A1 = V{1}; 
    vect_A2 = V{2};
    
    for k = 1:slice_per_ps_step-1
        if k < length(V)-1
            vect_A1 = cat(3, vect_A1, V{k+1});
            disp_vect{1} = cat(2, disp_vect{1}, disp_vect{k+1});
            disp_vect{2} = cat(2, disp_vect{2}, disp_vect{k+2});
            vect_A2 = cat(3, vect_A2, V{k+2});
        end
    end
%     if slice_per_ps_step == 1 % =1 --> classic order : 0, 180, 360, phi, phi+180 ...
%         V3 = img_3D(:, :, 3:slice_per_ps_step:end); % 3rd third  : 360, 375, 390 ...
%         
%     else
%         
%     end
    
    if filt == 1
        for k = 1:(slice_per_ps_step-1)*nb_images/slice_per_ps_step
            [ vect_A1(:, :, k), vect_A2(:, :, k) ] = filtering_img_stef ( vect_A1(:, :, k), vect_A2(:, :, k) );
            % function to do the filtering
        end  
    end
    
    if contr_mode % mode contrast
        % control where there are 180 deg phase shifted frame
        if (size(vect_A2, 3) ~= size(vect_A1, 3))
            vect_A1 = vect_A1(:,:,2:end);
            disp_vect{1} = disp_vect{1}(2:end);
            fprintf(2, 'Seems to miss one frame : I`m working without the first 2 ones\n');
        end
        contr = vect_A2 - vect_A1;
        disp_vect = cat(1, disp_vect{1} ,  disp_vect{2});
        
    else % mode raw
        
        contr = cat(3, vect_A1, vect_A2(:,:, nb_images/slice_per_ps_step + 1:end));
        disp_vect = cat(2, disp_vect{1}, disp_vect{2}(nb_images/slice_per_ps_step + 1:end));
    end
    
elseif slice_per_ps_step == 1 % increasing order
    
    if (~contr_mode && mod(size(img_3D, 3), 2)) % raw and odd number
        vect1 = 1:floor(size(img_3D, 3)/2);  
    else
        vect1 = 1:1:floor(size(img_3D, 3)/2);
    end
    
    vect2 = floor(size(img_3D, 3)/2)+1:size(img_3D, 3);
    
    if (size(img_3D, 3) - 180/diff_phase == 1+180/diff_phase)
        vect2 = vect2(2:end);
    end
    
    disp_vect{2} = vect2; 
    disp_vect{1} = vect1;
    
    vect_A1 = img_3D(:,:, vect1);
    vect_A2 = img_3D(:,:, vect2);
    
    if filt == 1
        for k = 1:max(size(vect_A1, 3), size(vect_A2, 3))
            ind1 = k; ind2 = k;
            if size(vect_A1, 3) > size(vect_A2, 3)
                ind2 = k-1;
            elseif size(vect_A2, 3) > size(vect_A1, 3)
                ind1 = k-1;
            end
            [ vect_A1(:, :, ind1), vect_A2(:, :, ind2) ] = filtering_img_stef ( vect_A1(:, :, ind1), vect_A2(:, :, ind2) );
            % function to do the filtering
        end
    end
    
    if contr_mode % mode contrast
        % control where there are 180 deg phase shifted frame
        contr = vect_A2 - vect_A1;
        disp_vect = cat(1,disp_vect{1} ,  disp_vect{2});
    else % mode raw
        disp_vect = cat(2, disp_vect{1}, disp_vect{2});
        contr=cat(3, vect_A1, vect_A2);
%         contr(:,:,1+180/diff_phase:end) = vect_A2;
    end
end

if (shg==1 && length(rect) > 1)
    xmin = rect(1); ymin = rect(2); width = rect(3); height = rect(4);
    img_shg = img_shg(round(ymin): round(ymin + height - 1), ...
        round(xmin): round(xmin + width - 1), :);
end

if (double_contrast && slice_per_ps_step > 1 && contr_mode)
    
    fprintf(2, 'WARNING contrast doubled : be sure that the order was 0,180,360, phi ... (in contrast_by_subtract)\n');
    
    for it_doubl_ctr = 1:double_contrast % % double_contrast = 1 for double, 2 for quadruple (0 = normal contrast)
        disp_vect0 = disp_vect;
        VV=cell(1, slice_per_ps_step - it_doubl_ctr);
        disp_vect = cell(1, slice_per_ps_step - it_doubl_ctr);
        for k = 1:(slice_per_ps_step - it_doubl_ctr)
            vectk = (k-1)*(size(img_3D, 3)/slice_per_ps_step)+1:(k)*(size(img_3D, 3)/slice_per_ps_step);
            disp_vect{k} = disp_vect0(:, vectk);
            VV{k} = contr(:, :, vectk); % 1st third : 0, 15, 30 ...
         % 2nd third : 180, 195, 210 ...
        end
        vvect_A1 = VV{1}; % init
        vvect_A2 = VV{2}; % init

        for k = 1:slice_per_ps_step-it_doubl_ctr-2
%             if k < length(VV)-it_doubl_ctr
            vvect_A1 = cat(3, vvect_A1, VV{k+1});
            disp_vect{1}= cat(2, disp_vect{1}, disp_vect{k+1}); 
            vvect_A2 = cat(3, vvect_A2, VV{k+2});
            disp_vect{2}= cat(2, disp_vect{2}, disp_vect{k+2}); 
%             end
        end

        contr = vvect_A1 - vvect_A2;
        disp_vect = cat(1,disp_vect{1} ,  disp_vect{2});

    %     contr =  contr(:,:,1:size(contr,3)/2)-contr(:,:,size(contr,3)/2+1:end);
        % not the good meth
    % %     for k = 1: size(contr, 3)
    % %         contr(:,:,k) = contr(:,:,k) - median(median(contr(:,:,k),1),2);
    % %     end
        x_phase = x_phase(1:length(x_phase)*(slice_per_ps_step - it_doubl_ctr-1)/(slice_per_ps_step - it_doubl_ctr));

    end
end

% % [min_val, max_val] = sat_vals;
mask_val_sat = sum(single( abs(contr) >= sat_vals(1) &  abs(contr) <= sat_vals(2)), 3); maskmask = mask_val_sat < 3;
mask_val_sat(maskmask) = NaN; mask_val_sat(~maskmask) = 1;% 3 values minimum to recover the phase

if avg_median % avg nearest neighbors : median method
    coef_centre=1; coef_cross=1; coef_diag=1; mean_meth=0; sigma_med=0; real_sigma=1; % dflt
    func_hdl = load_stack_plot_ISHG;
    [im1, ~] = func_hdl.avg_nearest_ut(contr, skip_asking_coeff_median, coef_centre, coef_cross, coef_diag, mean_meth, sigma_med, real_sigma, 0);
    if ~isempty(im1); contr = im1; end % % empty if user canceled
end

if normalize_interframes
    norm_add_then_mult = 0;  % 1 = two steps
    contr = normalize_interframes_func(contr, avg_median, norm_add_then_mult, 1); %near_avg_ini, norm_add_then_mult, 0);
end

if use_range ~= 0 % use only certain frames
    contr = contr(:,:, use_range(1):use_range(2):use_range(3));
    x_phase = x_phase(use_range(1):use_range(2):use_range(3));
    disp_vect = disp_vect(:, use_range(1):use_range(2):use_range(3));
end

if isa(disp_vect, 'cell')
    for k=1:length(disp_vect)
        disp('disp_vect'); disp(disp_vect{k}); 
    end
else
    disp('disp_vect'); disp(disp_vect);
    if  (slice_per_ps_step == 1 || ~contr_mode || force_incr_order) % increasing
        disp((disp_vect-1)*diff_phase);
    else
%         disp_vect = disp_vect(disp_vect
        if exist('vspec', 'var')
            for ii =1:length(vspec)
                vspec{ii} = (slice_per_ps_step*vspec{ii}-slice_per_ps_step+ii);
            end
            
            kk=sort(cell2mat(vspec));
        else
            kk=1:nb_images;%0:diff_phase:(diff_phase*nb_images);
        end
        vect1 = diff_phase*floor((kk-1)/slice_per_ps_step)+mod(kk-1,slice_per_ps_step)*180;
        disp(vect1(disp_vect));
    end
end
% % whos disp_vect
disp('x_phase'); disp(x_phase)

end

