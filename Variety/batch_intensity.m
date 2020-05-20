
%% batch intensity

clear all
clc

% path00 = 'D:\Documents\';
path00 = 'C:\Users\pc\Documents';

addpath(fullfile(path00, 'These\codes Matlab\Codes_I-SHG\MP'));

[fname, folder_name, FILTERINDEX] = uigetfile('*.tif','Select your stack of images !','MultiSelect','on');
if ~FILTERINDEX
    error('File not chosen ! Program ends.');
end
cd(folder_name);

ch00 = menu('Choose mode', 'SHG', 'iSHG (stacks with phase-shifts)', 'several files with stack to average', 'SHG with ind. images')-1; % 0 for SHG, 1 for iSHG

crop_batch = menu('Do you want to crop (see code below for parameters?)', 'No', 'Yes', 'Point center every frame')-1 ;

if crop_batch == 2
    point_center_batch_plotting = 1;
end


roi_same = 1;

switch ch00
    case 0 % SHG
        
        %         tableau{k,1} = fname{k};
        
        warning('OFF', 'all');
        InfoImage=imfinfo(fname); % disp a warning 'Division by zero when processing YResolution'
        warning('ON', 'all');
        mImage=InfoImage(1).Width;
        nImage=InfoImage(1).Height;
        num_images=length(InfoImage);
        img_3D=zeros(nImage,mImage,num_images,'uint16');
        
        warning('OFF', 'all'); % because of Tiff warning
        TifLink = Tiff(fname, 'r');
        for i=1:num_images
            TifLink.setDirectory(i);
            img_3D(:,:,i)=TifLink.read();
        end
        TifLink.close();
        warning('ON', 'all');
        
        if crop_batch
            
            % xmin = 58.701058; ymin = 50.500000; width = 82.539683; height = 82.539683;
            if (~roi_same || ~exist('rect', 'var'))
                figure(101); hfwd = imagesc(img_3D(:,:,1)); colormap(gray); colorbar;
                suptitle( 'Select an ROI with the mouse (right click to make it square)')
                rect = getrect(101);
            end
            
            % [cropped_img1 rect] = imcrop(tmp1); % crop by
            % rectangle instead of free hand
            img_3D00 = img_3D;
            xmin = rect(1); ymin= rect(2); width= rect(3); height= rect(4);
            % img_3D = imcrop(img_3D,rect);
            % bwd_crop = imcrop(bwd,rect);
            % parameters to be changed accordingly
            img_3D = img_3D(round(ymin): round(ymin + height - 1), round(xmin): round(xmin + width - 1), :);
        end
        
        im = double(img_3D);
        
        tableau = squeeze(mean(mean(im)));
        
        tableau_sqrt = sqrt(tableau);
        
    case 1 % iSHG
        
        for k=1:length(fname)
            
            tableau{k,1} = fname{k};
            
            warning('OFF', 'all');
            InfoImage=imfinfo(fname{k}); % disp a warning 'Division by zero when processing YResolution'
            warning('ON', 'all');
            mImage=InfoImage(1).Width;
            nImage=InfoImage(1).Height;
            num_images=length(InfoImage);
            img_3D=zeros(nImage,mImage,num_images,'uint16');
            
            warning('OFF', 'all'); % because of Tiff warning
            TifLink = Tiff(fname{k}, 'r');
            for i=1:num_images
                TifLink.setDirectory(i);
                img_3D(:,:,i)=TifLink.read();
            end
            TifLink.close();
            warning('ON', 'all');
            
            if crop_batch
                height = 12; width = 20;
                % xmin = 58.701058; ymin = 50.500000; width = 82.539683; height = 82.539683;
                if (~roi_same || ~exist('rect', 'var'))
                    figure(101); hfwd = imagesc(img_3D(:,:,1)); colormap(gray) ; colorbar; axis image;
                    if point_center_batch_plotting
                        title('Click on you pattern center'); [x_center,y_center] = ginput(1);
                    else
                        suptitle( 'Select an ROI with the mouse (right click to make it square)')
                        rect = getrect(hfwd);
                    end
                end
                
                % [cropped_img1 rect] = imcrop(tmp1); % crop by
                % rectangle instead of free hand
                img_3D00 = img_3D;
                if point_center_batch_plotting
                    img_3D = img_3D(max(1, y_center - round(height/2)): min(size(img_3D, 1), y_center + round(height/2)), ...
                        max(1, x_center - round(width/2)): min(size(img_3D, 2), x_center + round(width/2)), :);
                else
                    img_3D = imcrop(img_3D,rect);
                end
                % bwd_crop = imcrop(bwd,rect);
                % parameters to be changed accordingly
                %img_3D = img_3D(round(ymin): round(ymin + height - 1), round(xmin): round(xmin + width - 1), :);
            end
            
            im = double(img_3D);
            
            tableau{k,2} = mean(mean(mean(im)));
            
        end
        
    case 3 % SHG with images to treat 1 by 1
        
        %         tableau{k,1} = fname{k};
        
        warning('OFF', 'all');
        InfoImage=imfinfo(fname); % disp a warning 'Division by zero when processing YResolution'
        warning('ON', 'all');
        mImage=InfoImage(1).Width;
        nImage=InfoImage(1).Height;
        num_images=length(InfoImage);
        img_3D00=zeros(nImage,mImage,num_images,'uint16');
        
        %         pas_content = 1;
        
        warning('OFF', 'all'); % because of Tiff warning
        TifLink = Tiff(fname, 'r');
        for i=1:num_images
            TifLink.setDirectory(i);
            img_3D00(:,:,i)=TifLink.read();
        end
        TifLink.close();
        warning('ON', 'all');
        %             img_3D00(:,:,i) = img_3D(:,:,i);
        for i=1:num_images
            if crop_batch
                pas_content = 1;
                while pas_content
                    % xmin = 58.701058; ymin = 50.500000; width = 82.539683; height = 82.539683;
                    
                        figure(101); hfwd = imagesc(img_3D00(:,:,i)); colormap(gray); colorbar;
                        suptitle( 'Select an ROI with the mouse (right click to make it square)')
                    if (~roi_same || ~exist('rect', 'var'))
                        rect = getrect(101);
                        xmin = rect(1); ymin= rect(2); width= rect(3); height= rect(4);
                    end
                    
                    if point_center_batch_plotting
                        title('Click on you pattern center'); [x_center,y_center] = ginput(1);
                        %             else
                        %                 suptitle( 'Select an ROI with the mouse (right click to make it square)')
                        %                 rect = getrect(hfwd);
                        
                    end
                    
                    
                    % [cropped_img1 rect] = imcrop(tmp1); % crop by
                    % rectangle instead of free hand
                    %                     img_3D00 = img_3D;
                    if point_center_batch_plotting
                        img_3D = img_3D00(max(1, y_center - round(height/2)): min(size(img_3D00, 1), y_center + round(height/2)), ...
                            max(1, x_center - round(width/2)): min(size(img_3D00, 2), x_center + round(width/2)), i);
                    else
                        img_3D = imcrop(img_3D00(:,:,i),rect);
                    end
                    
                    try
                        close(102)
                    end
                    
                    figure(102); imagesc(img_3D); axis image; colormap gray;
                    
                    pas_content = menu('Content ?', 'Y', 'Center not good', 'Size (+center) not good')-1;
                    
                    if pas_content == 2
                        clear rect
                    end
                    
                end
                
            end
            
            im = double(img_3D);
            
            tableau(i) = squeeze(mean(mean(im)));
            
            tableau_sqrt(i) = sqrt(tableau(i));
            % [cropped_img1 rect] = imcrop(tmp1); % crop by
            % rectangle instead of free hand
            %             img_3D00 = img_3D;
            %             xmin = rect(1); ymin= rect(2); width= rect(3); height= rect(4);
            %             % img_3D = imcrop(img_3D,rect);
            %             % bwd_crop = imcrop(bwd,rect);
            %             % parameters to be changed accordingly
            %             img_3D = img_3D(round(ymin): round(ymin + height - 1), round(xmin): round(xmin + width - 1), :);
        end
        
        
        
        
    case 2 % several files with stack to average
        
        for k=1:length(fname)
            
            tableau{k,1} = fname{k};
            
            warning('OFF', 'all');
            InfoImage=imfinfo(fname{k}); % disp a warning 'Division by zero when processing YResolution'
            warning('ON', 'all');
            mImage=InfoImage(1).Width;
            nImage=InfoImage(1).Height;
            num_images=length(InfoImage);
            img_3D=zeros(nImage,mImage,num_images,'uint16');
            
            warning('OFF', 'all'); % because of Tiff warning
            TifLink = Tiff(fname{k}, 'r');
            for i=1:num_images
                TifLink.setDirectory(i);
                img_3D(:,:,i)=TifLink.read();
            end
            TifLink.close();
            warning('ON', 'all');
            
            
            if crop_batch
                
                xmin = 58.701058; ymin = 50.500000; width = 82.539683; height = 82.539683;
                % parameters to be changed accordingly
                img_3D(:,:,i) = img_3D(round(ymin): round(ymin + height - 1), round(xmin): round(xmin + width - 1), :);
            end
            
            im = median(double(img_3D), 3);
            
            tableau{k,2} = mean(mean(im));
            tableau{k,3} = median(median(im));
            nbins = length(im(:));
            [N,edges] = histcounts(im(:),nbins);
            [~, II] = max(N);
            tableau{k,4} = edges(II);
            
            gain_specF_320 = exp(coef.p1*log(320) + coef.p2);
            gain_specF_380 = exp(coef.p1*log(380) + coef.p2);
        end
end

t=tableau(1:2:end); % fwd
t2=tableau(2:2:end);

return;

%% crop a stack that moves

[fname, folder_name, FILTERINDEX] = uigetfile('*.tif','For pointing center in a stack !','MultiSelect','on');
if ~FILTERINDEX
    error('File not chosen ! Program ends.');
end
cd(folder_name);

% ch00 = menu('Choose mode', 'SHG', 'iSHG (stacks with phase-shifts)', 'several files with stack to average')-1; % 0 for SHG, 1 for iSHG
%
% crop_batch = 2-menu('Do you want to crop (see code below for parameters?)', 'Yes', 'No');

roi_same = 0;
img_3D_new = [];

warning('OFF', 'all');
InfoImage=imfinfo(fname); % disp a warning 'Division by zero when processing YResolution'
warning('ON', 'all');
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
num_images=length(InfoImage);
img_3D=zeros(nImage,mImage,num_images,'uint16');

warning('OFF', 'all'); % because of Tiff warning
TifLink = Tiff(fname, 'r');
for i=1:num_images
    TifLink.setDirectory(i);
    img_3D(:,:,i) = TifLink.read();
end
TifLink.close();
warning('ON', 'all');

img_3D00 = img_3D;

outputFileName = [fname '_croped.tif'];

sizee = 250;

for k=1:num_images
    
    
    if (~roi_same || ~exist('rect', 'var'))
        hfwd =figure(101);  imagesc(img_3D00(:,:,k)); colormap(gray) ; colorbar;
        suptitle( 'Select CENTER of pattern')
        point = ginput(1);
        rect = [point(1)- round(sizee/2), point(2)-round(sizee/2), sizee, sizee];
    end
    
    
    img_3D_new = cat(3, img_3D_new, imcrop(img_3D00(:,:, k),rect));
    
    imwrite(img_3D_new(:, :, k), outputFileName, 'WriteMode', 'append',  'Compression','none');
end

%% Assemblate two stacks of different gain

clearvars

volt_PMT = [300, 400, 500, 600, 800, 1000];
gain_PMT = [1000, 7000, 30000, 110000, 800000, 4000000];

coef = fit(log(volt_PMT)', log(gain_PMT)','poly1');

prompt = {            'Gain stack # 1','Gain stack # 2'};
dlg_title = 'Gain';
num_lines = 1;
% Valeurs par défaut
def = {'500', '500'};
answer = inputdlg(prompt,dlg_title,num_lines,def);

volt_stck1 = str2double(cell2mat(answer(1)));
volt_stck2 = str2double(cell2mat(answer(2)));
% % %         gain_specF_calib = exp(coef.p1*log(volt_specF_calib) + coef.p2);
% % %         gain_specB_calib = exp(coef.p1*log(volt_specB_calib) + coef.p2);
gain_stck1 = exp(coef.p1*log(volt_stck1) + coef.p2);
gain_stck2 = exp(coef.p1*log(volt_stck2) + coef.p2);

%

for ii = 1:2
    
    [fname{ii}, folder_name, FILTERINDEX] = uigetfile('*.tif', 'Select your stack of images !', 'Multiselect', 'on');
    if ~FILTERINDEX
        error('File not chosen ! Program ends.');
    end
    cd(folder_name);
    
    warning('OFF', 'all');
    InfoImage=imfinfo(fname{ii}); % disp a warning 'Division by zero when processing YResolution'
    warning('ON', 'all');
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    num_images=length(InfoImage);
    %     img_3D=zeros(nImage,mImage,num_images,'uint16');
    
    warning('OFF', 'all'); % because of Tiff warning
    TifLink = Tiff(fname{ii}, 'r');
    ct = 1;
    for i=2:2:num_images
        TifLink.setDirectory(i);
        if ii == 2
            img_3D_2(:,:,ct) = gain_stck1/gain_stck2*double(TifLink.read());
        else
            img_3D(:,:,ct) = double(TifLink.read());
        end
        ct = ct+1;
    end
    TifLink.close();
    warning('ON', 'all');
    
end

im_tot = cat(3, img_3D, img_3D_2);

outputFileName = sprintf('%s_all.tif', fname{1});

for i=1:size(im_tot, 3)
    
    im_tot(:,:,i) = (im_tot(:,:,i)-min(min(min(im_tot))))*(2^16-1)/(max(max(max(im_tot))) - min(min(min(im_tot))));
    
    imwrite(uint16(im_tot(:, :, i)), outputFileName, 'WriteMode', 'append',  'Compression','none');
    
end

%
