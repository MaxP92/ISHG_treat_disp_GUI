function [FILTERINDEX, folder_name, realign_cross_correl, handles] = load_phase_func(handles, tit, ext_dflt, not_ask)
% [FILTERINDEX, handles] = load_phase_func(handles)
% % % Max Pinsard 2018.10.01
%  is used by load_phase and correct phase by ref

FILTERINDEX = 0; realign_cross_correl = 0;
% % varargout = {};

if ~not_ask
% not_ask = menu(sprintf('Load phase map %s...', tit),...
%     'From .mat file', ...
%     'From .fig file', 'From .tif file', '!!end please!!');
    not_ask = inputdlg({sprintf('Load phase map %s (tif, fig or mat)', tit),'Realign ph. map with ref. by 2D cross-correlation ?', ...
     'put - sign to data (1 for ph only, 2 for amp only, 3 for both)'},'Input',1, {ext_dflt, '0', '0'});
% else % not ask
end
if ~isempty(not_ask)
    if isa(not_ask, 'cell')
        extend = not_ask{1};
        realign_cross_correl = str2double(not_ask{2});
        minus_sign = str2double(not_ask{3});
    else; extend = not_ask;realign_cross_correl =0; minus_sign = 0; % dflt
    end
    % Handle response
%     switch not_ask
%         case 1%'From .mat file'
%             extend = '.mat';
%         case 2%'From .fig file'
%             extend = '.fig';
%         case 3%'From .tif file'
%             extend = '.tif';
%         case  4%  'end please'
%             error('user ends');
%     end
    msg_win_ld = sprintf('Select your ph .%s file %s!', extend, tit);
    name_dflt = sprintf('phase%s.%s', tit, extend);
    ext = sprintf('*.%s', extend);
    
    [fname, folder_name, FILTERINDEX] = uigetfile(ext, msg_win_ld, name_dflt);
    
    if ~FILTERINDEX
        folder_name=0; return;
    else % good
        extread = fname(end-2:end);
        if ~strcmp(extread, extend); extend = extread; end % user changed his mind
        cd(folder_name);
        if ~isempty(handles) % GUI used
            handles.fname = fname; handles.folder_name = folder_name;
        end
        ff=fullfile(folder_name, fname);
        disp(tit); disp(ff);
        % Handle response
        switch extend
            case 'tif'
                phase1 = imread(ff);
                if (~isempty(handles) && (~isfield(handles, 'x_cell') || isempty(handles.x_cell)))
                    handles.x_cell = {1:size(phase1, 2)};
                    handles.y_cell = {1:size(phase1, 1)};
                end
                [fna, fr_name, FX] = uigetfile(ext, 'Select your INTERF.CTR .tif file !', 'ctr.tif');
                if FX
                    amp = imread(fullfile(fr_name, fna));
                    if ~isempty(handles) % GUI used
% %                         handles.amp_cell = {amp};
                        set(handles.plot_diff_menu, 'Visible', 'on');
                        
                    end
                else % no amp loaded
                    amp = ones(size(phase1));
%                     if ~isempty(handles) % GUI used
%                         handles.amp_cell = {amp};
%                     end
                    fprintf(2, 'interf. ctr not found !! using ones.\n');
                end
                
            case 'mat'
                ld= load(ff);
                if ~isempty(handles) % GUI used
                    if isfield(ld, 'phase'); phase1=ld.phase;
                    else; phase1=struct2cell(ld);phase1=phase1{1};
                    end
                    if isfield(ld, 'amp');amp = ld.amp; else; amp=phase1*0+1; end
                    
                    if (~isfield(handles,'x_cell') ||  isempty(handles.x_cell))
                        handles.x_cell = {1:size(phase1, 2)};
                        handles.y_cell = {1:size(phase1, 1)};
                    end
                end
            case 'fig'
                figs = openfig(ff);
                children_ph=get(findall(figs(1),'type','axes'), 'Children'); % interf. contr
                
                phase1=children_ph.CData;
                xx = children_ph.XData;
                yy = children_ph.YData;
                if ~isempty(handles) % GUI used
                    handles.x_cell = {xx};
                    handles.y_cell = {yy};
                end
                if length(figs) > 1
                    children_contr=findall(figs(2),'type','axes'); % interf. contr
                    amp = children_contr.CData;
                else
                    ld_cntr = questdlg('Load .fig interf. contrast (close window if no)?', 'Load contr',...
                        'Yes',...
                        'use same as phase', 'use same as phase with sign -', 'use same as phase');
                    if ~isempty(handles); handles.contr00 = []; end
                    if strcmp(ld_cntr, 'Yes')
                        [fname, folder_name, FILTERINDEX] = uigetfile(ext, 'contr.fig !', 'contr.fig');
                        uiopen(fullfile(folder_name, fname), 1);
                        children=get(gca, 'Children');
                        amp = children.CData;
                    elseif strcmp(ld_cntr, 'use same as phase')
                        amp = phase1;
                    elseif strcmp(ld_cntr, 'use same as phase with sign -')
                        amp = -phase1; minus_sign = min(1, minus_sign); amp = amp-min(amp(:)) +1;
                    elseif (length(ld_cntr) < 1 || ~FILTERINDEX)
                        amp = phase1*0 + 1;
                        amp(phase1 >1) = 0;
                    end
                    
                end
        end
        
        if (minus_sign == 1 || minus_sign == 3)
            phase1 = -phase1;
        end
        if (minus_sign == 2 || minus_sign == 3)
            amp = -amp;
        end
        if ~isempty(handles) % GUI used
            handles.amp_cell = {amp};
            if ~isfield(handles, 'contr00')
                handles.contr00 = {ones(size(phase1))};
            end
        end
        
        if isempty(handles) % NO GUI used
            if ~exist('ld', 'var')
                ld.phase = phase1;  ld.amp = amp;
            end
            handles = ld;
        else % GUI
            handles.phase_cell = {phase1};
        end
    end
    
    % %     varargout = {fname, folder_name}
end

