figure(1);
plot(smooth(V(:,1)));hold on;
plot(smooth(V(:,2)));
plot(smooth(V(:,3)));
plot(smooth(V(:,4)));
plot(smooth(V(:,5)));
plot(smooth(V(:,6)));
hold off;
legend({'AniM','BLM2','b-splM','Ani','BLM1','b-spl'});