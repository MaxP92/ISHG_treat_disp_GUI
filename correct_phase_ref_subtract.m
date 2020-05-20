function [calib, amp_cell, phase_cell, x, y, x_cell, y_cell, phase_before_corr, folder_name] = correct_phase_ref_subtract( amp_cell, ...
phase_cell, resx, x_cell, y_cell, not_batch_chck, calib, phase_before_corr, diff_ctr, scaling_img, choice_load)
% [handles.phase_calib, handles.amp_cell, handles.phase_cell, handles.x,
% handles.y, handles.x_cell, handles.y_cell] =
% correct_phase_ref_subtract(handles.fname, amp_cell, phase_cell, handles.resx
%
% 2017.10.06 by Maxime PINSARD
% 
% correct the phase in I-SHG galvos by subtracting a ref phase
FILTERINDEX = 1; folder_name = 0; 
ext_dflt = 'fig';
if not_batch_chck %~isa(fname00, 'cell')
    disp('Select your .mat file of quartz (calib) phase !')
% %     [fname, folder_name, FILTERINDEX] = uigetfile('*.mat','Select your .mat file of quartz (calib) phase !', 'calib.mat');
    [FILTERINDEX, folder_name, do_realign, calib] = load_phase_func([], 'calib', ext_dflt, choice_load); % 0 for ask
    if ~FILTERINDEX % no file loaded
        if ~isempty(calib.phase) % param of function
            fprintf(2, '\n File not chosen ! Will use the same calib file as before. \n');
            do_realign = calib.do_realign;
        else
            error('File not chosen ! I end.');
        end
%     else % a file loaded
% % %         if isa(fname, 'char') % a correct file has been chosen
% % %             if strcmp(fname(end-2:end), 'mat')
% % %                 calib=load(fullfile(folder_name, fname));
% %                 m2=struct2cell(m);
% %                 phase_calib=m2{1};
%         if ~choice_load
%             do_realign = 2 - menu('Realign ph. map with ref. by 2D cross-correlation ?', 'Yes', 'No');
%         else; do_realign = 0;
%         end
% % %             end
% % %         end
    end
else
    do_realign = calib.do_realign;
end

if isfield(calib, 'phase')
    phase_calib = calib.phase;
else
    error('no phase mat in loaded !')
end
if isfield(calib, 'amp')
    amp_calib = calib.amp;
else
    warning('no amp mat in loaded !')
    amp_calib  = ones(size(phase_calib));
end

if not_batch_chck ~= -1 % otherwise just do beginning
    amp_corr = amp_cell{max(1,length(amp_cell)-diff_ctr)};

    if do_realign
        [phase_before_corr, phase_calib, amp_corr] = realign_pattern(phase_before_corr, phase_calib, amp_corr);
    end

    x = linspace(0,size(phase_before_corr,2)*resx,size(phase_before_corr,2));
    y = linspace(0,size(phase_before_corr,1)*resx,size(phase_before_corr,1));
    y_cell{end+1} = y; x_cell{end+1} = x;

    phase_before_corr(phase_before_corr == 1.05) = 1.05e6; % because otherwise it gets corrected also

    % subtracting the 2 phase plots

    if sum(size(phase_before_corr)) ~= sum(size(phase_calib))
        if (size(phase_before_corr, 1) == size(phase_calib, 2) && size(phase_before_corr, 2) == size(phase_calib, 1))
            phase_calib = phase_calib';
            amp_calib = amp_calib';
            disp('I transposed the img !');
        else
            if (size(phase_calib,1) >= size(phase_before_corr,1) && size(phase_calib,2) >= size(phase_before_corr,2))
                as1 = inputdlg({'calib size larger: crop at 1:upper left corner, 2:upper right, 3:lower left, 4:lower left (else return)', ...
                'do scaling on X (<=1 no)', 'do scaling on Y (<=1 no)'}, ...
                '!! calib sz larger X&Y', 1,{'1', '1','1'});
                if isempty(as1); crop_at = 0; 
                else; crop_at = str2double(as1{1});
                end
                if (crop_at >= 1 && crop_at <= 4)
                    do_scX = str2double(as1{2});
                    do_scY  = str2double(as1{3});
                    if (do_scX || do_scY)
                        scaling = [do_scX , do_scY];
                        phase_calib = scaling_img(phase_calib, scaling);
                        amp_calib =  scaling_img(amp_calib, scaling);
                    end
                    %inputdlg(prompt,dlg_title,num_lines,defAns)
                    if sum(size(phase_before_corr)) ~= sum(size(phase_calib))
                        if (crop_at == 1 || crop_at == 2); vy = 1:size(phase_before_corr,1);
                        elseif (crop_at == 3 || crop_at == 4); vy = size(phase_calib,1)-size(phase_before_corr,1)+1:size(phase_calib,1);
                        end
                        if (crop_at == 3 || crop_at == 4); vx = 1:size(phase_before_corr,2); 
                        elseif (crop_at == 1 || crop_at == 2); vx = size(phase_calib,2)-size(phase_before_corr,2)+1:size(phase_calib,2);
                        end
                        phase_calib = phase_calib(vy, vx);
                        amp_calib = amp_calib(vy, vx);
                        disp('I croped the larger phase calib !!'); % took the upper left corner for croping
                    end
                else; return;
                end
            else % calib smaller
                 as = inputdlg({'correct size diff by: 1:complete with zeros, 2: extend fit, <=0: return'; ...
                'if extend, mat_fit is : 1:Linear 2D, 11: Linear 1D on X, 12: Linear 1D on Y, 2:Parabola 2D, 21: Parabola 1D on X, 22: Parabola 1D on Y'}, ...
                    ['!! calib sz smaller: d', num2str(size(phase_before_corr)), ' VS c', num2str(size(phase_calib))],...
                    1, {'0'; '11'});
                    if ~isempty(as)
                        ch= str2double(as{1});
                    else; ch= 0;
                    end
                    if ~ch; return;  end
% %             inputdlg(prompt,dlg_title,num_lines,defAns)
                switch ch
                    case 2 % extend fit
    %                     ch1 = menu('Nature of mat_fit', 'Linear 2D', 'Extend 2D fit', 'Return');
                          matfitD= str2double(as{2});
                          addpath('C:\Users\pc\Documents\These\codes Matlab\Variety');
                          [phase_calib, ~, ~]=subplane_mp(phase_calib, size(phase_before_corr,2), size(phase_before_corr,1), 0, 0, 0, -matfitD, 0); 
                          % Z, x, y, coplanar_points, h, was_unwrap, surf_fit_tilt_auto_bacth, ampflag
                          if size(amp_calib) == size(phase_calib); amp_calib = phase_calib;
                          else;  [amp_calib, ~, ~]=subplane_mp(amp_calib, size(phase_before_corr,2), size(phase_before_corr,1), 0, 0, 0, -matfitD, 0); 
                          end
                    case 1 % zeros
                        ph00 = zeros(max(size(phase_before_corr), size(phase_calib))); ph01 = ph00;
                        ph00(1:size(phase_before_corr,1), 1:size(phase_before_corr,2)) = phase_before_corr; 
                        phase_before_corr = ph00;
                        ph00(1:size(amp_corr,1), 1:size(amp_corr,2)) = amp_corr; amp_corr = ph00;
                        ph01(1:size(phase_calib,1), 1:size(phase_calib,2)) = phase_calib; phase_calib = ph01;
                        ph01(1:size(amp_calib,1), 1:size(amp_calib,2)) =amp_calib ; amp_calib =ph01;
                end
            end
        end
    end
    
    if (max(phase_calib(:)) > 1 || min(phase_calib(:)) < -1)
        wrap_fitmat=menu('phasecalib is out of -pi,pi: wrap it ?', 'No', 'Yes',...
 'just recast by max division (if you are sure the phase range does not normally exceed -pi,pi')-1; 
        switch wrap_fitmat
            case 1 % wrap
                phase_calib= wrapToPi(phase_calib*pi)/pi ;
            case 2 % re-cast
                M = max(abs(phase_calib(:)));
                if M > 1 % maybe a interf ctr ref loaded
                    phase_calib = phase_calib/M;
                end
        end
        figure(55); imagesc(phase_calib); colorbar; axis image;
    end
    
    phase_corr = phase_before_corr - phase_calib;    
    phase_corr(phase_corr> 1.05e5) = 1.05e6;
    
    m=min(min(amp_calib));
    if m <= 0; amp_calib = amp_calib-m+1; end % avoid division by zero
    m2 = mean2(amp_calib(~isnan(amp_calib)));
   
%     amp_calib(amp_calib < 1) = m2; % normalize to pos values
    amp_corr = m2*amp_corr./amp_calib;
    amp_corr(amp_corr>100*mean2(amp_corr(~isnan(amp_corr)))) = NaN; %  avoid hotspots!!!

    % wrapping the phase (between -1, 1)
    phase_corr(phase_corr<-1 & phase_corr~= 1.05e6)=phase_corr(phase_corr<-1 & phase_corr~= 1.05e6)+2;
    phase_corr(phase_corr>1 & phase_corr~= 1.05e6)=phase_corr(phase_corr>1 & phase_corr~= 1.05e6)-2;

    phase_corr(phase_corr == 1.05e6) = 1.05; % because otherwise it gets corrected also

    phase_cell{end+1} = phase_corr;

    % the modified amp can be not good for hist3D
    % % amp_cell{end+1} = amp_cell{end};
    % % amp_cell{end-1} = amp_corr;
    amp_cell{end+1} = amp_corr;
else; x=0; y=0;
end
calib.phase = phase_calib;
calib.amp = amp_calib;
calib.do_realign = do_realign;
calib.FILTERINDEX = FILTERINDEX;

end