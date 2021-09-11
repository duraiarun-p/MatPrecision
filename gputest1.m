clc;clear;close all;
%%
load mri;
%%
D=gpuArray(D);
tic
parfor i=1:siz(3)
%     i;
%     i1=min(i+1,siz(3));
    ind=i:min(i+1,siz(3));
    Dt(:,:,i)=D(:,:,2)-D(:,:,1);
%     Do=
end
time=toc
Dt1=gather(Dt);
% figure(1),imshow(D(:,:,));

