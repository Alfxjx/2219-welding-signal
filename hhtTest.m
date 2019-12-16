fs = 100000;
t_v = (0:length(imf_volt)-1)/fs;
t_c = (0:length(imf_curr)-1)/fs;

% figure(7);
% hht(imf_volt,fs);
% ylim([0 5000]);
% xlabel('Time(s)');
% title('volt hht');

figure(8);
hht(imf_curr,fs);
ylim([0 6000]);
xlabel('Time(s)');
title('curr hht');

figure(9);
[hs,f,t] = hht(imf_curr);
mesh(seconds(t),f,hs,'EdgeColor','none','FaceColor','interp');
xlabel('Time')
ylabel('Frequency')
zlabel('Instantaneous Energy')
