%����Ƶ�ʷֱ�Ϊ10Hz��35Hz�����������ź�����ӵĸ����ź�
%����Ƶ��fs=2048Hz���ź�
%���ʽ���£�y=5sin(2*pi*10t)+5*sin(2*pi*35t)
function fftfenxi
clear;clc;
N=2048;
%fftĬ�ϼ�����ź��Ǵ�0��ʼ��
t=linspace(1,2,N);deta=t(2)-t(1);fs=1/deta;
x=5*sin(2*pi*10*t)+5*sin(2*pi*35*t);
N1=256;N2=512;w1=0.2*2*pi;w2=0.3*2*pi;w3=0.4*2*pi;
%
x=(t>=-200&t<=-200+N1*deta).*sin(w1*t)+(t>-200+N1*deta&t<=-200+N2*deta).*sin(w2*t)+(t>-200+N2*deta&t<=200).*sin(w3*t);
y = x;
m=0:N-1;
f=1./(N*deta)*m;%���Բ鿴�α������������������Ƶ�ʷ�Χ��

%%%%%%%%%%%%%%%%%%%%%%%%%��������ҶƵ��ͼ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��������Y����x(t)�ĸ���Ҷ�任��ֵ
%Y=exp(i*4*pi*f).*fft(y)%�����������Ƶ�׳���exp(i*4*pi*f)�õ�Ƶ�ƺ�[-2,2]֮���Ƶ��ֵ
Y=fft(y);
z=sqrt(Y.*conj(Y));
plot(f(1:100),z(1:100));
title('��Ƶ����')

xiangwei=angle(Y);
figure(2)
plot(f,xiangwei)
title('��Ƶ����')

figure(3)
plot(t,y,'r')
%axis([-2,2,0,1.2])
title('ԭʼ�ź�')

%%%%%%%%%%%%%%%%%%%%��Hilbert�任ֱ������źŵ�˲ʱƵ�� %%%%%%%%%%%%%%%%%%%%%
clear;clc;clf;
%����������ĺ�����z=t^3
N=2048;
%fftĬ�ϼ�����ź��Ǵ�0��ʼ��
t=linspace(1,2,N);deta=t(2)-t(1);fs=1/deta;
x=5*sin(2*pi*10*t)+5*sin(2*pi*35*t);
z=x;
hx=hilbert(z);
xr=real(hx);xi=imag(hx);

%����˲ʱ���
sz=sqrt(xr.^2+xi.^2);
%����˲ʱ��λ
sx=angle(hx);
%����˲ʱƵ��
dt=diff(t);
dx=diff(sx);
sp=dx./dt;
plot(t(1:N-1),sp)
title('˲ʱƵ��')

%%%%%%%%%С�᣺����Ҷ�任���ܵõ�˲ʱƵ�ʣ������ܵõ�ĳ��ʱ�̵�Ƶ��ֵ��
%Hilbert�任����ȡ˲ʱƵ�ʵķ����������ֻ��Hilbert�任�������˲ʱƵ��Ҳ��׼ȷ�������ָ�Ƶ��ʵ���ϸ�Ƶû�����壡��

%%%%%%%%%%%%%%%%%%%%%%%��HHT��ȡ�źŵ�ʱƵ����߼���%%%%%%%%%%%%%%%%% 
function HHT
clear;clc;clf;
N=2048;
%fftĬ�ϼ�����ź��Ǵ�0��ʼ��
t=linspace(1,2,N);deta=t(2)-t(1);fs=1/deta;
x=5*sin(2*pi*10*t)+5*sin(2*pi*35*t);
z=x;
c=emd(z);

%%%%%%%%%%%%%%%%%%����ÿ��IMF���������һ��ʣ�����residual��ͼ��%%%%%%%%%%%%%%%%%
%����imf����ĵ�һ��Ԫ�أ���ԭʼ����
subplot(m+1,1,1)
plot(t,z)
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['signal','Amplitude'])

%����ȫ����IMF������ͼ��
for i=1:m-1
subplot(m+1,1,i+1);
set(gcf,'color','w')
plot(t,c(i,:),'k')
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['imf',int2str(i)])
end

%�������һ��ʣ�����residual�Ĳ���
subplot(m+1,1,m+1);
set(gcf,'color','w')
plot(t,c(m,:),'k')
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['r',int2str(m-1)])

%%%%%%%%%%%%%%%%%%%%%%%����ÿ��IMF������ʣ�����residual�ķ�Ƶ����%%%%%%%%%%%%%%%%%%%%%%
%������һ��ԭʼ���ݵķ�Ƶ����
figure;
subplot(m+1,1,1)
set(gcf,'color','w')
[f,z]=fftfenxi(t,z);
plot(f,z,'k')
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['initial signal',int2str(m-1),'Amplitude'])

%����ÿ��IMF�����ķ�Ƶ����
for i=1:m-1
subplot(m+1,1,i+1);
set(gcf,'color','w')
[f,z]=fftfenxi(t,c(i,:));
plot(f,z,'k')
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['imf',int2str(i),'Amplitude'])
end

%�������һ��ʣ�����residual�ķ�Ƶ����
subplot(m+1,1,m+1);
set(gcf,'color','w')
[f,z]=fftfenxi(t,c(m,:));
plot(f,z,'k')
set(gca,'fontname','times New Roman')
set(gca,'fontsize',14.0)
ylabel(['r',int2str(m-1),'Amplitude'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����HHTʱƵ�׺ͱ߼���%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A,fa,tt]=hhspectrum(c);
[E,tt1]=toimage(A,fa,tt,length(tt));

figure(3)
disp_hhs(E,tt1) %��άͼ��ʾHHTʱƵ�ף�E����õ�HHT��
pause

E=flipud(E);

for k=1:size(E,1)
bjp(k)=sum(E(k,:))*1/fs; 
end

f=(1:N-2)/N*(fs/2);

figure(5)
plot(f,bjp);
xlabel('Ƶ�� / Hz');
ylabel('�źŷ�ֵ');
title('�źű߼���')%Ҫ��߼��ױ����ȶ��źŽ���EMD�ֽ�



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
%fftĬ�ϼ�����ź��Ǵ�0��ʼ��
t=linspace(t(1),t(L),N);deta=t(2)-t(1);
m=0:N-1;
f=1./(N*deta)*m;
%��������Y����x(t)�ĸ���Ҷ�任��ֵ
%Y=exp(i*4*pi*f).*fft(y)%�����������Ƶ�׳���exp(i*4*pi*f)�õ�Ƶ�ƺ�[-2,2]֮���Ƶ��ֵ
Y=fft(y);
z=sqrt(Y.*conj(Y));