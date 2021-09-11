%% NN Activation Comparison
clc;clear;close all;
addpath('/home/arun/Documents/MATLAB/ImageDB/');
%%
imD=imread('DA_PP.jpg');



net = squeezenet;
im = imread('face.jpg');
imgSize = size(im);
im1=imresize3(imD,imgSize);

analyzeNetwork(net)

% laname='relu_conv1';
laname='conv1';
act11 = activations(net,im,laname);
act21 = activations(net,im1,laname);

sz = size(act11);
act1 = reshape(act11,[sz(1) sz(2) 1 sz(3)]);
act2 = reshape(act21,[sz(1) sz(2) 1 sz(3)]);
I1 = imtile(mat2gray(act1),'GridSize',[8 8]);
I2 = imtile(mat2gray(act2),'GridSize',[8 8]);
%%
figure(11);
imshow(I1)
figure(12);
imshow(I2)
figure(13);
imshow(imabsdiff(I1,I2));
%%
diffmatrix=zeros(sz);
for fn=1:sz(3)
    diffmatrix(:,:,fn)=imabsdiff(act1(:,:,fn),act2(:,:,fn));
end
diffm=mean(diffmatrix,3);
figure;imshow(diffm,[]);
%%
% figure(1);
% % fn=15;
% for fn=1:sz(3)
% subplot(131),
% imshow(act1(:,:,fn),[]);
% subplot(132),
% imshow(act2(:,:,fn),[]);
% subplot(133),
% imshow(imabsdiff(act1(:,:,fn),act2(:,:,fn)),[]);
% title(num2str(fn));
% pause(0.25);
% end