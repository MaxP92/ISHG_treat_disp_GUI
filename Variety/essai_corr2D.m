% essai corr 2D
% Maxime PINSARD
% 17/04/2016

close all
clear all
clc

nb = 50;
taille_thick = 2;
taille = 20;
template = 0.2*ones(nb);
template(ceil(nb/2)-taille_thick+1:ceil(nb/2)+taille_thick-1,ceil(nb/2)-taille:ceil(nb/2)+taille) = 0.6;
template(ceil(nb/2)-taille:ceil(nb/2)+taille,ceil(nb/2)-taille_thick+1:ceil(nb/2)+taille_thick-1) = 0.6;
offsettemplate = 0.2*ones(nb);
offset = [10 12];
offsettemplate(min(ceil(nb/2)+offset(1), size(offsettemplate,1)),min(ceil(nb/2)-taille+offset(2), size(offsettemplate,2)):min(ceil(nb/2)+taille+offset(2), size(offsettemplate,2))) = 0.6;
offsettemplate(min(ceil(nb/2)-taille+offset(1), size(offsettemplate,1)):min(ceil(nb/2)+taille+offset(1), size(offsettemplate,1)), min(ceil(nb/2)+offset(2), size(offsettemplate,2))) = 0.6;
% offsetTemplate((1:size(template,1))+offset(1), ...
%     (1:size(template,2))+offset(2)) = template;

% offsetTemplate = offsetTemplate();

imagesc(offsettemplate)
colormap gray
figure
imagesc(template)
colormap gray
% axis equal

cc = xcorr2(offsettemplate,template);
[max_cc, imax] = max(abs(cc(:)));
[ypeak, xpeak] = ind2sub(size(cc),imax(1));
corr_offset = [(ypeak-size(template,1)) (xpeak-size(template,2))]