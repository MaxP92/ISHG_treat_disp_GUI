
close all;

I = imread('imgs/image1.png');
figure(1), imshow(I);

tic;
res = jointWMF(I,I,10,25.5,256,256,1,'exp');
toc;

figure(2), imshow(res);