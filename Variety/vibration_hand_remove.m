% 2018.9.17, Maxime PINSARD

cd('C:\Users\pc\Documents\These\Biomedical\Muscle - ishg\2018-8 larvae\Z=-40um');
[fname, folder_name, FILTERINDEX] = uigetfile('*.tif','Select your image !','i.tif');
if ~FILTERINDEX
    error('File not chosen ! Program ends.');
end

img00 = imread(fullfile(folder_name, fname)); % uint16
img_3D = img00; mm=mean(mean(img00));

img_3D(11:12, :) = [ones(2, 8)*mm, img_3D(11:12, 1:end-8)];
img_3D(15:16, :) = [ones(2, 3)*mm, img_3D(15:16, 1:end-3)]; % to the right
img_3D(14, :) = [img_3D(14, 9:end), ones(1, 8)*mm]; % to the left


img00 = img_3D;

%%
imwrite(img_3D, sprintf('%s_re-set.tif', fname(1:end-4)), 'WriteMode', 'append',  'Compression','none');