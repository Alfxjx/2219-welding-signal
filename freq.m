% 计算一下数据的频率

dadi = importdata('1.txt',',',0);
gm = importdata('gm.txt',',',0);

start = 350000;
step = 50000;

dadi_volt = dadi(start:(start+step),1);
dadi_curr = dadi(start:(start+step),2);

gm_volt = gm(start:(start+step),1);
gm_curr = gm(start:(start+step),2);

% 采样频率100kHz
f = 100*1000;

fy1 = fft(dadi_volt,512);
fy2 = fft(dadi_curr,512);
fy3 = fft(gm_volt,512);
fy4 = fft(gm_curr,512);

w = f*(0:256);
figure(1);
plot(w,fy1(1:257));
figure(2);
plot(w,fy2(1:257));
figure(3);
plot(w,fy3(1:257));
figure(4);
plot(w,fy4(1:257));

