clc;clear;close all;
%%
% A = randi(m,n);
% windowL = 25;
% mbSize = 3;
% for k=1:v:size(A,2)
%       if (k+windowL+-1) > size(A,2) || (k+v+windowL-1) > size(A,2)
%           break;
%       end
%       dataA = A(:,k:k+windowL-1);
%       fprintf('a) %d:%d ',k,k+windowL-1);
%       % do something
%       dataB = A(:,k+v+windowL:k+v+windowL-1);
%       % do something
%       fprintf('b) %d:%d\n',k+windowL,k+v+windowL-1);
% end

% for ind = 1:mbSize:size(A,2)
%     Window = A(ind:ind+windowL-1);
% end
m=5;n=5;%image
A = eye(n,m);B = eye(n,m);mbSize = 3;%i=2;j=2;

% refBlkVer = i + m1;
% refBlkHor = j + n1;
% Aw=A(refBlkVer:refBlkVer+mbSize-1,refBlkHor:refBlkHor+mbSize-1);


firstindex_i=1;lastindex_i=m;
firstindex_j=1;lastindex_j=n;
m1=-1;n1=-1;
% AT=A;
% if m1<0 || n1<0
%     AT=fliplr(A);
%     m1=abs(n1);
%     n1=abs(m1);
% end

bwi=max(firstindex_i,firstindex_i+m1);
bwie=min(lastindex_i,lastindex_i+m1);
bwj=max(firstindex_j,firstindex_j+n1);
bwje=min(lastindex_j,lastindex_j+n1);



awi=max(firstindex_i,firstindex_i-m1);
awie=min(lastindex_i,lastindex_i-m1);
awj=max(firstindex_j,firstindex_j-n1);
awje=min(lastindex_j,lastindex_j-n1);

Aw1=A(awi:awie,awj:awje);%static-current
Bw1=B(bwi:bwie,bwj:bwje);%moving-refblk
% delje=n-awje;
% delie=m-awie;
% deli=firstindex_i-awi;
% delj=firstindex_j-awj;
% bwie=lastindex_i-delie;
% bwje=lastindex_j-delje;
% bwi=firstindex_i-deli;
% bwj=firstindex_j-delj;


% Aw1=fliplr(Aw1);