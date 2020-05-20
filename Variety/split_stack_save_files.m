% 201.3.9 Maxime PINSARD
%
% split a stack in TIFF into different substacks, with a recursive name at each files

batch = 11;
bits16 = 1;

for bbtch = 1:batch
    clc
    clearvars
    close all

    %% load your file

    [fname, folder_name, FILTERINDEX] = uigetfile('*.tif','Choose your stack or several single images !','MultiSelect','on');
    if ~FILTERINDEX
        error('File not chosen ! Program ends.');
    end
    cd(folder_name);

    fname1 = fname;

    if iscell(fname) % several files
        msgbox('You chose several files, so it will consider that each file contains one image and will concatenate them !')
        nb_files = length(fname);
        one_stack = 0;
        fname1=fname{1};
    %     fprintf('Your files are being read ... be patient\n')

    else % one big stack

        one_stack = 1;
        nb_files = 1;
    end

    % ch00 = menu('Choose mode', 'SHG', 'iSHG (stacks with phase-shifts)', 'several files with stack to average')-1; % 0 for SHG, 1 for iSHG
    % 
    % crop_batch = 2-menu('Do you want to crop (see code below for parameters?)', 'Yes', 'No');

    warning('OFF', 'all');
    InfoImage=imfinfo(fname1); % disp a warning 'Division by zero when processing YResolution'
    warning('ON', 'all');
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    num_images=length(InfoImage);
    img_3D=zeros(nImage,mImage,max(num_images, nb_files),'uint16');

    if ~one_stack % several files
        for ii=1:nb_files
            img_3D(:,:,ii) = imread(fname{ii},1);
        end

    else % one big stack

        warning('OFF', 'all'); % because of Tiff warning
        TifLink = Tiff(fname1, 'r');
        for i=1:num_images
            TifLink.setDirectory(i);
            img_3D(:,:,i) = TifLink.read();
        end
        TifLink.close();
        warning('ON', 'all');
    end

    max_im = max(max(max(img_3D)));
    min_im = min(min(min(img_3D)));
    img_3D00 = img_3D;

    %% define param

    prompt = {sprintf('Number of sub-stacks (must be a divider of %d', size(img_3D, 3)), 'Suffix to add to file name (beginning)', 'Step', 'unit'};
    dlg_title = 'Parameters';
    num_lines = 1;
    % Valeurs par défaut
    def = {num2str(size(img_3D, 3)/2), 'Z=', '0.3', 'um'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    nb_substacks = str2double(answer{1});
    suffix_name = answer{2};
    step_suffix = str2double(answer{3});
    unit_suffix = answer{4};

    % outputFileName = [fname '_croped.tif'];
    % 
    % sizee = 250;

    %% if avg

    if nb_substacks == 1
        img_3D = median(img_3D00, 3);
    end


    %% split into several stacks and save
    minN=min(min(min(img_3D)));
    maxN=max(max(max(img_3D)));
    
    switch bits16
        case 1
            img_3D_towrite = uint16((img_3D-minN)*(max_im - min_im)/(maxN-minN)+ min_im);
        case 0
            img_3D_towrite = img_3D;
    end

    for ind_substack = 1:nb_substacks

        outputFileName = sprintf('%s_%s%.3f%s.tif', fname1(1:end-4), suffix_name, step_suffix*(ind_substack-1), unit_suffix ); 

        for k = (ind_substack-1)*size(img_3D, 3)/nb_substacks+1:ind_substack*size(img_3D, 3)/nb_substacks

            imwrite(img_3D_towrite(:, :, k), outputFileName, 'WriteMode', 'append',  'Compression','none');
        end
    end

    if ~one_stack % several files
       size_max = 100;

        for ii = 1:floor(size(img_3D, 3)/size_max) % because matlab won't save stacks with too many slices, you have to concatenate them after in ImageJ
            for k = 1 + (ii-1)*size_max:(ii)*size_max
                imwrite(img_3D_towrite(:, :, k), sprintf('%s_ALL_%d.tif', fname1(1:end-4), ii), 'WriteMode', 'append',  'Compression','none');
            end
        end

        % remaining indices
        for k = (ii)*size_max + 1:size(img_3D, 3)
            imwrite(img_3D_towrite(:, :, k), sprintf('%s_ALL_%d.tif', fname1(1:end-4), ii+1), 'WriteMode', 'append',  'Compression','none');
        end

    end
end