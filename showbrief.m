figure(1);
plot(dadi(:,1));
ylim([-50,50]);
ylabel('voltage(V)');
xlabel('Time(0.01ms)');
title('voltage');

figure(2);
plot(dadi(:,2));
ylim([-400,400]);
ylabel('current(A)')
xlabel('Time(0.01ms)');
title('current');


% figure(3);
% plot(gm(:,1));
% title('gm volt');
% 
% figure(4);
% plot(gm(:,2));
% title('gm curr');