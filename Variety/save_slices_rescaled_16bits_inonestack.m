% 201.3.9 Maxime PINSARD
%
% re-save a stack, with its slices rescaled independently (relative intensity info
% is lost !)

clc
clearvars
close all

%% load your file

[fname, folder_name, FILTERINDEX] = uigetfile('*.tif','sel. files ! !','MultiSelect','on');
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
img_3D=zeros(nImage,mImage,max(num_images, nb_files));

if ~one_stack % several files
    for ii=1:nb_files
        img_3D(:,:,ii) = double(imread(fname{ii},1));
    end
    
else % one big stack
    warning('OFF', 'all'); % because of Tiff warning
    TifLink = Tiff(fname1, 'r');
    for i=1:num_images
        TifLink.setDirectory(i);
        img_3D(:,:,i) = double(TifLink.read());
    end
    TifLink.close();
    warning('ON', 'all');
end

img_3D00 = img_3D;

%% re-scaling

for i=1:num_images 
    min1_fwd_new = min(min(img_3D(:,:,i))); max1_fwd_new = max(max(img_3D(:,:,i)));
    img_3D(:,:,i) = (img_3D(:,:,i) - min1_fwd_new)*(2^16-1)/(max1_fwd_new - min1_fwd_new);
    imwrite(uint16(img_3D(:, :, i)), sprintf('%s_intensity_info_lost.tif', fname1(1:end-4)), 'WriteMode', 'append',  'Compression','none');
end




