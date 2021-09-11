clc;clear;close all;
%%
nhooper=6;
nhood=8;

se = strel('disk',nhood);

se1 = strel('diamond',nhood+round(nhood/nhooper));
test1=se1.Neighborhood(4:end-3,4:end-3);
test2=imresize(test1,[2*nhood-1, 2*nhood-1]);
%test2=test2(2:end-1,2:end-1);
%%
figure(1),
subplot(131),
imshow(test1);
subplot(132),
imshow(test2);
subplot(133),
imshowpair(se.Neighborhood,test2,'diff');