%给定频率分别为10Hz和35Hz的两个正弦信号相叠加的复合信号
%采样频率fs=2048Hz的信号
%表达式如下：y=5sin(2*pi*10t)+5*sin(2*pi*35t)
function fftfenxi
clear;clc;
N=2048;
%fft默认计算的信号是从0开始的
t=linspace(1,2,N);deta=t(2)-t(1);fs=1/deta;
x=5*sin(2*pi*10*t)+5*sin(2*pi*35*t);
N1=256;N2=512;w1=0.2*2*pi;w2=0.3*2*pi;w3=0.4*2*pi;
%
x=(t>=-200&t<=-200+N1*deta).*sin(w1*t)+(t>-200+N1*deta&t<=-200+N2*deta).*sin(w2*t)+(t>-200+N2*deta&t<=200).*sin(w3*t);
y = x;
m=0:N-1;
f=1./(N*deta)*m;%可以查看课本就是这样定义横坐标频率范围的

%%%%%%%%%%%%%%%%%%%%%%%%%画出傅里叶频谱图%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%下面计算的Y就是x(t)的傅里叶变换数值
%Y=exp(i*4*pi*f).*fft(y)%将计算出来的频谱乘以exp(i*4*pi*f)得到频移后[-2,2]之间的频谱值
Y=fft(y);
z=sqrt(Y.*conj(Y));
plot(f(1:100),z(1:100));
title('幅频曲线')

xiangwei=angle(Y);
figure(2)
plot(f,xiangwei)
title('相频曲线')

figure(3)
plot(t,y,'r')
%axis([-2,2,0,1.2])
title('原始信号')

%%%%%%%%%%%%%%%%%%%%用Hilbert变换直接求该信号的瞬时频率 %%%%%%%%%%%%%%%%%%%%%
clear;clc;clf;
%假设待分析的函数是z=t^3
N=2048;
%fft默认计算的信号是从0开始的
t=linspace(1,2,N);deta=t(2)-t(1);fs=1/deta;
x=5*sin(2*pi*10*t)+5*sin(2*pi*35*t);
z=x;
hx=hilbert(z);
xr=real(hx);xi=imag(hx);

%计算瞬时振幅
sz=sqrt(xr.^2+xi.^2);
%计算瞬时相位
sx=angle(hx);
%计算瞬时频率
dt=diff(t);
dx=diff(sx);
sp=dx./dt;
plot(t(1:N-1),sp)
title('瞬时频率')

%%%%%%%%%小结：傅里叶变换不能得到瞬时频率，即不能得到某个时刻的频率值。
%Hilbert变换是求取瞬时频率的方法，但如果只用Hilbert变换求出来的瞬时频率也不准确。（出现负频，实际上负频没有意义！）

%%%%%%%%%%%%%%%%%%%%%%%用HHT求取信号的时频谱与边际谱%%%%%%%%%%%%%%%%% 
function HHT
clear;clc;clf;
N=2048;
%fft默认计算的信号是从0开始的
t=linspace(1,2,N);deta=t(2)-t(1);fs=1/deta;
x=5*sin(2*pi*10*t)+5*sin(2*pi*35*t);
z=x;
c=emd(z);

%%%%%%%%%%%%%%%%%%画出每个IMF分量及最后一个剩余分量residual的图形%%%%%%%%%%%%%%%%%
%画出imf矩阵的第一列元素，即原始数据
subplot(m+1,1,1)
plot(t,z)
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['signal','Amplitude'])

%画出全部的IMF分量的图形
for i=1:m-1
subplot(m+1,1,i+1);
set(gcf,'color','w')
plot(t,c(i,:),'k')
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['imf',int2str(i)])
end

%画出最后一个剩余分量residual的波形
subplot(m+1,1,m+1);
set(gcf,'color','w')
plot(t,c(m,:),'k')
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['r',int2str(m-1)])

%%%%%%%%%%%%%%%%%%%%%%%画出每个IMF分量及剩余分量residual的幅频曲线%%%%%%%%%%%%%%%%%%%%%%
%画出第一列原始数据的辐频曲线
figure;
subplot(m+1,1,1)
set(gcf,'color','w')
[f,z]=fftfenxi(t,z);
plot(f,z,'k')
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['initial signal',int2str(m-1),'Amplitude'])

%画出每个IMF分量的辐频曲线
for i=1:m-1
subplot(m+1,1,i+1);
set(gcf,'color','w')
[f,z]=fftfenxi(t,c(i,:));
plot(f,z,'k')
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['imf',int2str(i),'Amplitude'])
end

%画出最后一个剩余分量residual的辐频曲线
subplot(m+1,1,m+1);
set(gcf,'color','w')
[f,z]=fftfenxi(t,c(m,:));
plot(f,z,'k')
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['r',int2str(m-1),'Amplitude'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%计算HHT时频谱和边际谱%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A,fa,tt]=hhspectrum(c);
[E,tt1]=toimage(A,fa,tt,length(tt));

figure(3)
disp_hhs(E,tt1) %二维图显示HHT时频谱，E是求得的HHT谱
pause

E=flipud(E);

for k=1:size(E,1)
bjp(k)=sum(E(k,:))*1/fs; 
end

f=(1:N-2)/N*(fs/2);

figure(5)
plot(f,bjp);
xlabel('频率 / Hz');
ylabel('信号幅值');
title('信号边际谱')%要求边际谱必须先对信号进行EMD分解



function [A,f,tt] = hhspectrum(x,t,l,aff)
error(nargchk(1,4,nargin));

if nargin < 2
t=1:size(x,2);
end

if nargin < 3
l=1;
end

if nargin < 4
aff = 0;
end

if min(size(x)) == 1
if size(x,2) == 1
x = x';

if nargin < 2
t = 1:size(x,2);
end
end
Nmodes = 1;
else
Nmodes = size(x,1);
end

lt=length(t);
tt=t((l+1):(lt-l));


for i=1:Nmodes
an(i,:)=hilbert(x(i,:)')';
f(i,:)=instfreq(an(i,:)',tt,l)';
A=abs(an(:,l+1:end-l));

if aff
disprog(i,Nmodes,max(Nmodes,100))
end

end





function disp_hhs(im,t,inf)
% DISP_HHS(im,t,inf)
% displays in a new figure the spectrum contained in matrix "im"
% (amplitudes in log).
% inputs : - im : image matrix (e.g., output of "toimage")
% - t (optional) : time instants (e.g., output of "toimage") 
% - inf (optional) : -dynamic range in dB (wrt max)
% default : inf = -20

% utilisation : disp_hhs(im) ; disp_hhs(im,t) ; disp_hhs(im,inf) 
% disp_hhs(im,t,inf)

figure
colormap(bone)
colormap(1-colormap);

if nargin==1
inf=-20;
t = 1:size(im,2);
end

if nargin == 2
if length(t) == 1
inf = t;
t = 1:size(im,2);
else
inf = -20;
end
end

if inf >= 0
error('inf doit etre < 0')
end

M=max(max(im));

im = log10(im/M+1e-300);
inf=inf/10;
imagesc(t,fliplr((1:size(im,1))/(2*size(im,1))),im,[inf,0]);
set(gca,'YDir','normal')
xlabel(['time'])
ylabel(['normalized frequency'])
title('Hilbert-Huang spectrum')

function [f,z]=fftfenxi(t,y)
L=length(t);N=2^nextpow2(L);
%fft默认计算的信号是从0开始的
t=linspace(t(1),t(L),N);deta=t(2)-t(1);
m=0:N-1;
f=1./(N*deta)*m;
%下面计算的Y就是x(t)的傅里叶变换数值
%Y=exp(i*4*pi*f).*fft(y)%将计算出来的频谱乘以exp(i*4*pi*f)得到频移后[-2,2]之间的频谱值
Y=fft(y);
z=sqrt(Y.*conj(Y));