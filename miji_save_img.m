function [mij_obj, foldr_put, title_put] = miji_save_img(curr_map, cmap, mij_obj, folder_name, fname, foldr_put, title_put, str, ph_flag) 
%
%
% 2017.10.16 Maxime PINSARD
% 
% needs MIJI (or ImageJ-Matlab - modify script accordingly)
path_fiji = 'C:\Users\pc\Documents\fiji-win64\Fiji.app';

if isa(cmap, 'char')%length(cmap) <= 1
% %     switch cmap
    disp(['using predefined cmap ', cmap]);
    cmap_char = cmap;
else % cmap special
    %% refine LUT

    [nb1, ~] = size(cmap);

    if nb1 < 256 % usually 64

        xx = (1:nb1)';
        xx_refined = linspace(1, nb1, 256)';

        reds = interp1(xx, cmap(:,1), xx_refined);
        greens = interp1(xx, cmap(:,2), xx_refined);
        blues = interp1(xx, cmap(:,3), xx_refined);

        cmap_new = [reds, greens, blues];

    elseif nb1 > 256

    %     cmap_new = zeros(256, 3);
    %     cmap = [zeros(floor((256-mod(size(cmap, 1), 256))/2), 3); cmap; zeros(ceil((256-mod(size(cmap, 1), 256))/2), 3)];

        xx = (1:nb1)';
        xx_refined = linspace(1, nb1, nb1 + 256-mod(size(cmap, 1), 256))';

        reds = interp1(xx, cmap(:,1), xx_refined);
        greens = interp1(xx, cmap(:,2), xx_refined);
        blues = interp1(xx, cmap(:,3), xx_refined);

        cmap2 = [reds, greens, blues];
        cmap_new = downsample(cmap2, size(cmap2, 1)/256);
    %     for i =1:3
    %        
    %         cmap_new(:, i) = mean(reshape(cmap2(:, i), 256, []),2);
    %         cmap_new(:, i) = cmap_new(:, i)/max(cmap_new(:, i) );
    %     end

    else % == 256

        cmap_new = cmap;
    end

    cmap_new(1,:) = cmap(1,:);
    cmap_new(end,:) = cmap(end,:);

    %% LUT save

    cell_cmap = [(0:size(cmap_new, 1)-1)', round(cmap_new*255)];
    path_lut_fiji = fullfile(path_fiji, '\luts');

    fileID = fopen(fullfile(path_lut_fiji, 'temp_map.lut'),'w');
    fprintf(fileID,'%s %s %s %s\n', 'Index Red Green Blue');
    fclose(fileID);

    fileID2 = fopen(fullfile(path_lut_fiji, 'temp_map.lut'),'a');
    formatSpec = '\n%d %d %d %d';

    [nrows,~] = size(cell_cmap);
    for row = 1:nrows
        fprintf(fileID2, formatSpec, cell_cmap(row, :));
    end

    fclose(fileID2);
    cmap_char = 'temp_map';
end
%% miji
new_mij = 0; % default
try
    MIJ.run('Close All');
catch
    new_mij = 1;
end
if ~isa(mij_obj,  'MIJ')
    new_mij = 1;
end

if new_mij
    disp('No instance of Miji running')
    try
        MIJ.exit;
    catch
        disp('No instance of Miji to close')
    end

    if isa(mij_obj,  'MIJ')
        clear mij_obj
    end

    path_ml = 'C:\Program Files\MATLAB\R2017b\java\jar';
    path_mij =fullfile(path_ml, 'mij.jar');
    javaaddpath(path_mij);
    addpath(fullfile(path_fiji, 'scripts'));
    disp('Wait while Miji is opening ...');
    warning('off', 'MATLAB:mir_warning_unrecognized_pragma');
    mij_obj = Miji(true);
    % % methods(MIJ) to get all functions !!
end

switch ph_flag 
    case 0 % interf. contr.
        name = sprintf('%s_%s', fullfile(folder_name, fname), str);
    case 1 % phase
        name = sprintf('%s_%s', fullfile(folder_name, fname), str);
end

MIJ.createImage(name, curr_map, true);
switch ph_flag
    case 1 % phase
        ij.IJ.runMacro('setMinAndMax(-1, 1)');
end

MIJ.run(cmap_char);

try % if the user canceled
    if (~isa(foldr_put,'char') || ~isa(title_put,'char') ) % ask fldr
        MIJ.run('Save'); % cannot put a path or it will show no dlg %, sprintf('save=[%s_.tif]', fullfile(folder_name, fname(1:end-4))));
        ij.IJ.runMacro('getDirectory("image")');  %saveAs("Tiff", "C:/Users/pc/Desktop/0.tif");
        s=strsplit(char(MIJ.getLog), '\n'); s=s(~cellfun(@isempty,s)); foldr_put=s{end};
        title_put = char(MIJ.getCurrentTitle); title_put = title_put(1:max(length(title_put)-9, 1)); %title_put
    else
        MIJ.run('Save', sprintf('path=[%s]',fullfile(foldr_put, sprintf('%s.tif',title_put))));
    end
        %     errr
    
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:Java:GenericException')
        rethrow(ME)
    end
end

disp(fullfile( foldr_put, title_put))