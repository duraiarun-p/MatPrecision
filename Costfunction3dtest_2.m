clc;clear;close all;
%%
rsiz=5;    
Rm=randi(2,rsiz);
Rm1=randi(2,rsiz);
%     Rmt=omegdiif_DrDave(Rm,rsiz,rsiz);
%     Rmt1=omegdeiff(Rm,rsiz,rsiz);
%     Rmt2=omegadd(Rm,rsiz,rsiz);
%     Rmt3=omegadd1(Rm,Rm1,rsiz,rsiz);
%     Rmt4=omegadd2(Rm,Rm1,rsiz,rsiz);
r1=Rm1;r=Rm;Bl1=zeros(rsiz);Bl2=zeros(rsiz);
up=r1(1:end-1,1:end)-r(2:end,1:end);
Bl1(1:end-1,1:end)=up;
% %Down
down=r1(2:end,1:end)-r(1:end-1,1:end);
Bl1(2:end,1:end)=Bl1(2:end,1:end)-down;
% %Right
right=r1(1:end,2:end)-r(1:end,1:end-1);
Bl1(1:end,2:end)=Bl1(1:end,2:end)-right;
% %Left
left=r1(1:end,1:end-1)-r(1:end,2:end);
Bl1(1:end,1:end-1)=Bl1(1:end,1:end-1)-left;
% %leftup45
leftup45=r1(2:end,2:end)-r(1:end-1,1:end-1);
Bl1(2:end,2:end)=Bl1(2:end,2:end)-leftup45;
% %rightup45
rightup45=r(2:end,1:end-1)-r(1:end-1,2:end);
Bl1(2:end,1:end-1)=Bl1(2:end,1:end-1)-rightup45;
% %leftdown45
leftdown45=r1(1:end-1,1:end-1)-r(2:end,2:end);
Bl1(1:end-1,1:end-1)=Bl1(1:end-1,1:end-1)-leftdown45;
% %rightdown45
rightdown45=r1(1:end-1,2:end)-r(2:end,1:end-1);
Bl1(1:end-1,2:end)=Bl1(1:end-1,2:end)-rightdown45;


up=myadd(r1(1:end-1,1:end),r(2:end,1:end));
Bl2(1:end-1,1:end)=up;
% %Down
down=myadd(r1(2:end,1:end),r(1:end-1,1:end));
Bl2(2:end,1:end)=myadd(Bl2(2:end,1:end),down);
% %Right
right=myadd(r1(1:end,2:end),r(1:end,1:end-1));
Bl2(1:end,2:end)=myadd(Bl2(1:end,2:end),right);
% %Left
left=myadd(r1(1:end,1:end-1),r(1:end,2:end));
Bl2(1:end,1:end-1)=myadd(Bl2(1:end,1:end-1),left);
% %leftup45
leftup45=myadd(r1(2:end,2:end),r(1:end-1,1:end-1));
Bl2(2:end,2:end)=myadd(Bl2(2:end,2:end),leftup45);
% %rightup45
rightup45=myadd(r(2:end,1:end-1),r(1:end-1,2:end));
Bl2(2:end,1:end-1)=myadd(Bl2(2:end,1:end-1),rightup45);
% %leftdown45
leftdown45=myadd(r1(1:end-1,1:end-1),r(2:end,2:end));
Bl2(1:end-1,1:end-1)=myadd(Bl2(1:end-1,1:end-1),leftdown45);
% %rightdown45
rightdown45=myadd(r1(1:end-1,2:end),r(2:end,1:end-1));
Bl2(1:end-1,2:end)=myadd(Bl2(1:end-1,2:end),rightdown45);

%%
function vo=myadd(v1,v2)
vo=v1-v2;
end
% function BlocOd=omegdiif_DrDave(omeg,LenX,LenY)
% 
% BlocOd=zeros(LenX,LenY);
% for dispX=-1:1
%     for dispY=-1:1
%         if (dispX==0) && (dispY==0)
%             continue
%         end
%         BlocOd(max(1-dispY,1):(end-max(dispY,0)),max(1-dispX,1):(end-max(dispX,0))) = ...
%             BlocOd(max(1-dispY,1):(end-max(dispY,0)), max(1-dispX,1):(end-max(dispX,0))) - ...
%             abs(omeg(max(1-dispY,1):(end-max(dispY,0)), max(1-dispX,1):(end-max(dispX,0))) - ...
%             omeg(max(1+dispY,1):(end+min(dispY,0)), max(1+dispX,1):(end+min(dispX,0))));
%     end
% end
% 
% end
% 
% function Bl1=omegdeiff(r,M,N)
% 
% % [M,N]=size(r);
% Bl1=zeros(M,N);
% 
% Bl1(1:end-1,1:end)=r(1:end-1,1:end)-r(2:end,1:end);
% % % Bl(2:end,1:end)=Bl(2:end,1:end)+r(1:end-1,1:end);
% % %Right
% Bl1(1:end,2:end)=Bl1(1:end,2:end)-r(1:end,1:end-1);
% % %Left
% Bl1(1:end,1:end-1)=Bl1(1:end,1:end-1)-r(1:end,2:end);
% % Diagnonal 4 directions
% Bl1(2:end,2:end)=Bl1(2:end,2:end)-r(1:end-1,1:end-1);
% Bl1(2:end,1:end-1)=Bl1(2:end,1:end-1)-r(1:end-1,2:end);
% Bl1(1:end-1,1:end-1)=Bl1(1:end-1,1:end-1)-r(2:end,2:end);
% Bl1(1:end-1,2:end)=Bl1(1:end-1,2:end)-r(2:end,1:end-1);
% % %Down
% Bl1(2:end,1:end)=Bl1(2:end,1:end)-r(1:end-1,1:end);
% Bl1(end,1:end)=Bl1(end,1:end)-r(end,1:end);
% 
% Bl1=abs(Bl1);
% 
% end
% function Bl1=omegadd(r,M,N)
% 
% % [M,N]=size(r);
% Bl1=zeros(M,N);
% 
% Bl1(1:end-1,1:end)=r(1:end-1,1:end)+r(2:end,1:end);
% % % Bl(2:end,1:end)=Bl(2:end,1:end)+r(1:end-1,1:end);
% % %Right
% Bl1(1:end,2:end)=Bl1(1:end,2:end)+r(1:end,1:end-1);
% % %Left
% Bl1(1:end,1:end-1)=Bl1(1:end,1:end-1)+r(1:end,2:end);
% % Diagnonal 4 directions
% Bl1(2:end,2:end)=Bl1(2:end,2:end)+r(1:end-1,1:end-1);
% Bl1(2:end,1:end-1)=Bl1(2:end,1:end-1)+r(1:end-1,2:end);
% Bl1(1:end-1,1:end-1)=Bl1(1:end-1,1:end-1)+r(2:end,2:end);
% Bl1(1:end-1,2:end)=Bl1(1:end-1,2:end)+r(2:end,1:end-1);
% % %Down
% Bl1(2:end,1:end)=Bl1(2:end,1:end)+r(1:end-1,1:end);
% Bl1(end,1:end)=Bl1(end,1:end)+r(end,1:end);
% 
% Bl1=abs(Bl1);
% 
% end
% function Bl1=omegadd1(r1,r,M,N)
% 
% % [M,N]=size(r);
% % Bl1=zeros(M,N);
% Bl1=r1;
% Bl1(1:end-1,1:end)=Bl1(1:end-1,1:end)+r(2:end,1:end);
% % % Bl(2:end,1:end)=Bl(2:end,1:end)+r(1:end-1,1:end);
% % %Right
% Bl1(1:end,2:end)=Bl1(1:end,2:end)+r(1:end,1:end-1);
% % %Left
% Bl1(1:end,1:end-1)=Bl1(1:end,1:end-1)+r(1:end,2:end);
% % Diagnonal 4 directions
% Bl1(2:end,2:end)=Bl1(2:end,2:end)+r(1:end-1,1:end-1);
% Bl1(2:end,1:end-1)=Bl1(2:end,1:end-1)+r(1:end-1,2:end);
% Bl1(1:end-1,1:end-1)=Bl1(1:end-1,1:end-1)+r(2:end,2:end);
% Bl1(1:end-1,2:end)=Bl1(1:end-1,2:end)+r(2:end,1:end-1);
% % %Down
% Bl1(2:end,1:end)=Bl1(2:end,1:end)+r(1:end-1,1:end);
% Bl1(end,1:end)=Bl1(end,1:end)+r(end,1:end);
% 
% Bl1=abs(Bl1);
% 
% end
% function Bl1=omegadd2(r1,r,M,N)
% 
% % [M,N]=size(r);
% Bl1=zeros(M,N);
% % %Up
% Bl1(1:end-1,1:end)=Bl1(1:end-1,1:end)+r1(1:end-1,1:end)+r(2:end,1:end);
% % %Right
% Bl1(1:end,2:end)=Bl1(1:end,2:end)+r1(1:end,2:end)+r(1:end,1:end-1);
% % %Left
% Bl1(1:end,1:end-1)=Bl1(1:end,1:end-1)+r1(1:end,1:end-1)+r(1:end,2:end);
% % Diagnonal 4 directions
% % Bl1(2:end,2:end)=Bl1(2:end,2:end)+r1(2:end,2:end)+r(1:end-1,1:end-1);
% % Bl1(2:end,1:end-1)=Bl1(2:end,1:end-1)+r(2:end,1:end-1)+r(1:end-1,2:end);
% % Bl1(1:end-1,1:end-1)=Bl1(1:end-1,1:end-1)+r1(1:end-1,1:end-1)+r(2:end,2:end);
% % Bl1(1:end-1,2:end)=Bl1(1:end-1,2:end)+r1(1:end-1,2:end)+r(2:end,1:end-1);
% % %Down
% Bl1(2:end,1:end)=Bl1(2:end,1:end)+r1(2:end,1:end)+r(1:end-1,1:end);
% Bl1(end,1:end)=Bl1(end,1:end)+r1(end,1:end)+r(end,1:end);
% 
% Bl1=abs(Bl1);
% 
% end