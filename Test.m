X=input('please input X : X= '); %输入元胞数
Y=input('please input Y : Y= ');
T0=ones(X,Y)*923;   %初始温度               %由925改来
Tcool=-30;   % 冷却速率-3K/s
s0=zeros(X,Y);  % 标注元胞状态，0液
U=zeros(X,Y);  %表示晶向
T=zeros(X,Y);
Q=48;  %晶向总数
Tend=880;
Teut=922;
Ii=imagesc(U);
ti=title('t=0','Fontsize',14);
k=1;
k0=0.33;
ddx=3e-6;
c0=2.8;
dt=0.005;         %时间由0.05改成0.001
Dl=3e-9;
Ta=0.1;
Tn=0.5;
nmax=5.5e10;

Nnuc=round(nmax*2*(X+Y-1)*ddx^2);  %计算四壁中能够形核的元胞数目
pi=Nnuc/(X+Y-1)/2;
mat1=rand(X,2); %给第一列和最后一列赋随机值
% 让它们和形核概率比较，通过大小确定能否形核
mat2=rand(Y,2);
for i=1:X
    if mat1(i,1)<pi
        s0(i,1)=2;   %s0=2 用来标记型壁能够形核的元胞
    elseif mat1(i,2)<pi
        s0(i,Y)=2;
    end
end
for i=1:Y
    if mat2(i,1)<pi
        s0(1,i)=2;
    elseif mat2(i,2)<pi
        s0(X,i)=2;
    end
end
        
Ta2=0.1;
Tn2=5;
nmax2=3.5e8;

Nnuc2=round(nmax2*(X-2)*(Y-2)*ddx^2);  %计算内部能够形核的元胞数目
pi2=Nnuc2/(X-2)/(Y-2);
mat=rand(X-2,Y-2);
for i=1:X-2
    for j=1:Y-2
        if mat(i,j)<pi2
            s0(i+1,j+1)=3;  %s0=3用来标记内部能够形核的元胞
        end
    end
end
% 给预定的形核元胞一个符合高斯分布形核模型的过冷度
for i=1:X
    for j=1:Y
        if s0(i,j)==2
            T(i,j)=normrnd(Tn,Ta);
        elseif s0(i,j)==3
            T(i,j)=normrnd(Tn2,Ta2);   %返回一个符合正态分布的数值
        end
    end
end


T1=zeros(X,Y);
T2=zeros(X,Y);
u=s0;
U1=U; % 标记邻居的晶粒方向
%dx=[-1 0 0 1];
%dy=[0 -1 1 0];

dx=[-1 0 1 0];
dy=[0 -1 0 1];


fs=zeros(X,Y);  %元胞初始状态
fs0=zeros(X,Y);
Cl=ones(X,Y)*c0;  
Cs=zeros(X,Y);
Cl0=zeros(X,Y);  
Cs0=zeros(X,Y);

E=0.3;
st0=0;
GibbT=2e-7;
ML=-3;%待修正

fsx=zeros(X,Y);
fsy=zeros(X,Y);
fsxx=zeros(X,Y);
fsyy=zeros(X,Y);
fsxy=zeros(X,Y);
Dfs=zeros(X,Y);
st1=zeros(X,Y);
clx=zeros(X,Y);
cly=zeros(X,Y);
coe=zeros(X,Y);
cb=0;
%fst1=zeros(n);


while k<100
    u=s0;
for i=1:X
    for j=1:Y
        T1(i,j)=T0(i,j)+Tcool*dt;
       if Tend<T1(i,j)<Teut & (s0(i,j)==2 | s0(i,j)==3)
            T2(i,j)=Teut-T1(i,j);
            if T2(i,j)>T(i,j) & U(i,j)==0
                U(i,j)=round(rand(1)*Q);
                
                u(i,j)=4;
                fs0(i,j)=1;
                
                Cl0(i,j)=0;
                Cs0(i,j)=k0*c0;
                
                XX=i+dx;
                YY=j+dy;
                % 周期边界条件
                XX=mod(XX-1,X)+1;
                YY=mod(YY-1,Y)+1;                
                for m=1:4
                    if s0(XX(m),YY(m))==0
                    u(XX(m),YY(m))=1;
                    U1(XX(m),YY(m))=U(i,j);
                    end
                end
            end 
       elseif  Tend<T1(i,j)<Teut  & s0(i,j)==1
           U(i,j)=U1(i,j);   % 让该邻居形核并再标记它的邻居
          if i==1&j==1
              
           fsx(1,1)=(fs(1,2)-1)/ddx/2;
           fsy(1,1)=(fs(2,1)-1)/ddx/2;
           fsxx(1,1)=(fs(1,2)+1-2*fs(1,1))/ddx^2;
           fsyy(1,1)=(1+fs(2,1)-2*fs(1,1))/ddx^2;
           fsxy(1,1)=(1-fs(2,2))/4/ddx^2;
          elseif i==1&j==Y
           
           fsx(1,Y)=(1-fs(1,Y-1))/ddx/2;
           fsy(1,Y)=(1-fs(2,Y))/ddx/2;
           fsxx(1,Y)=(1+fs(1,Y-1)-2*fs(1,Y))/ddx^2;
           fsyy(1,Y)=(1+fs(2,Y)-2*fs(1,Y))/ddx^2;
           fsxy(1,Y)=(fs(2,Y-1)-1)/4/ddx^2;
          elseif i==X&j==1
           
           fsx(X,1)=(fs(X,2)-1)/ddx/2;
           fsy(X,1)=(fs(X-1,1)-1)/ddx/2;
           fsxx(X,1)=(fs(X,2)+1-2*fs(X,1))/ddx^2;
           fsyy(X,1)=(1+fs(X-1,1)-2*fs(X,1))/ddx^2;
           fsxy(X,1)=(fs(X-1,2)-1)/4/ddx^2;
          elseif i==X&j==Y
           
           fsx(X,Y)=(1-fs(X,Y-1))/ddx/2;
           fsy(X,Y)=(fs(X-1,Y)-1)/ddx/2;
           fsxx(X,Y)=(fs(X,Y-1)+1-2*fs(X,Y))/ddx^2;
           fsyy(X,Y)=(1+fs(X-1,Y)-2*fs(X,Y))/ddx^2;
           fsxy(X,Y)=(1-fs(X-1,Y-1))/4/ddx^2;
          elseif i==1&j~=1&j~=Y
           
           fsx(1,2:Y-1)=(fs(1,3:Y)-fs(1,1:Y-2))/ddx/2;
           fsy(1,2:Y-1)=(1-fs(2,2:Y-1))/ddx/2;
           fsxx(1,2:Y-1)=(fs(1,3:Y)+fs(1,1:Y-2)-2*fs(1,2:Y-1))/ddx^2;
           fsyy(1,2:Y-1)=(1+fs(2,2:Y-1)-2*fs(1,2:Y-1))/ddx^2;
           fsxy(1,2:Y-1)=(fs(2,3:Y)-fs(2,1:Y-2))/4/ddx^2;
          elseif i==X&j~=1&j~=Y
           
           fsx(X,2:Y-1)=(fs(X,3:Y)-fs(X,1:Y-2))/ddx/2;
           fsy(X,2:Y-1)=(fs(X-1,2:Y-1)-1)/ddx/2;
           fsxx(X,2:Y-1)=(fs(X,3:Y)+fs(X,1:Y-2)-2*fs(X,2:Y-1))/ddx^2;
           fsyy(X,2:Y-1)=(1+fs(X-1,2:Y-1)-2*fs(X,2:Y-1))/ddx^2;
           fsxy(X,2:Y-1)=(fs(X-1,3:Y)-fs(X-1,1:Y-2))/4/ddx^2;
           
          elseif j==1&i~=1&i~=X
           
           fsx(2:X-1,1)=(fs(2:X-1,2)-1)/ddx/2;
           fsy(2:X-1,1)=(fs(1:X-2,1)-fs(3:X,1))/ddx/2;
           fsxx(2:X-1,1)=(fs(2:X-1,2)+1-2*fs(2:X-1,1))/ddx^2;
           fsyy(2:X-1,1)=(fs(2:X-1,1)+fs(3:X,1)-2*fs(2:X-1,1))/ddx^2;
           fsxy(2:X-1,1)=(fs(1:X-2,2)-fs(3:X,2))/4/ddx^2;
          elseif j==Y&i~=1&i~=X
           
           
           fsx(2:X-1,Y)=(1-fs(2:X-1,Y-1))/ddx/2;
           fsy(2:X-1,Y)=(fs(1:X-2,Y)-fs(3:X,Y))/ddx/2;
           fsxx(2:X-1,Y)=(fs(2:X-1,Y-1)+1-2*fs(2:X-1,Y))/ddx^2;
           fsyy(2:X-1,Y)=(fs(1:X-2,Y)+fs(3:X,Y)-2*fs(2:X-1,Y))/ddx^2;
           fsxy(2:X-1,Y)=(fs(3:X,Y-1)-fs(1:X-2,Y-1))/4/ddx^2;
          else
           
          % for ii=2:n-1
          %     for jj=2:n-1
          %         fsx(ii,jj)=(fs(ii,jj+1)-fs(ii,jj-1))/dx/2;
          %         fsy(ii,jj)=(fs(ii-1,jj)-fs(ii+1,jj))/dx/2;
          %         fsxx(ii,jj)=(fs(ii,jj+1)+fs(ii,jj-1)-2*fs(ii,jj))/dx^2;
          %         fsyy(ii,jj)=(fs(ii-1,jj)+fs(ii+1,jj)-2*fs(ii,jj))/dx^2;
          %         fsxy(ii,jj)=(fs(ii-1,jj+1)-fs(ii-1,jj-1)-fs(ii+1,jj+1)+fs(ii+1,jj-1))/4/dx^2;
          %    end
          % end   
          fsx(i,j)=(fs(i,j+1)-fs(i,j-1))/ddx/2;                                            %去掉了绝对值
          fsy(i,j)=(fs(i-1,j)-fs(i+1,j))/ddx/2;
          fsxx(i,j)=(fs(i,j+1)+fs(i,j-1)-2*fs(i,j))/ddx^2;
          fsyy(i,j)=(fs(i-1,j)+fs(i+1,j)-2*fs(i,j))/ddx^2;
          fsxy(i,j)=(fs(i-1,j+1)-fs(i-1,j-1)-fs(i+1,j+1)+fs(i+1,j-1))/4/ddx^2;
          end
           
          %界面斜率
          Kx=(2*fsx(i,j)*fsy(i,j)*fsxy(i,j)-fsxx(i,j)*fsy(i,j)^2-fsyy(i,j)*fsx(i,j)^2)/((fsx(i,j)^2+fsy(i,j)^2)^1.5); 
          if fsx(i,j) <fsy(i,j)
              st1(i,j)=atan(fsx(i,j) /fsy(i,j));
          else st1(i,j)=atan(fsy(i,j) /fsx(i,j));
          end
          T2(i,j)=Teut-T1(i,j);                                                            %增加了T2的计算
          fst1=1-15*E*cos(4*(st1(i,j)-st0));
          C6=Cl(i,j)-(T2(i,j)-GibbT*Kx*fst1)/ML;                                                %把cl(i,j)改成了C0
          
            %计算固相率的变化
            XX=i+dx;
            YY=j+dy;
            XX=mod(XX-1,X)+1;
            YY=mod(YY-1,Y)+1; 
            Cb=0;
            for i0=1:4
                if s0(XX(i0),YY(i0))==0
                    Cb=C6-Cl(XX(i0),YY(i0))+Cb;
                end
            end
            DC=Dl*dt*Cb/ddx^2;
            Dfs(i,j)=DC/C6/(1-k0);
            %考虑界面凝固率的变化
             if  Dfs(i,j)>=(1-fs(i,j))  %界面完全凝固
                Dfs(i,j)=1-fs(i,j);
                                                                              %删去了u(i,j)=4；
                Cl0(i,j)=0;
                Cs0(i,j)=(Cs(i,j)*fs(i,j)+k0*Cl(i,j)*Dfs(i,j))/(fs(i,j)+Dfs(i,j));
                                                                              %删去了fs0(i,j)=1; 
                cb=Cl(i,j)*(1-k0)*Dfs(i,j);  
            elseif  Dfs(i,j)<(1-fs(i,j))  %液相还有剩余
                CL=(Cl(i,j)*(1-fs(i,j))-k0*Cl(i,j)*Dfs(i,j))/(1-fs(i,j)-Dfs(i,j));
                if  CL<=C6
                    Cl0(i,j)=CL;
                    Cs0(i,j)=(Cs(i,j)*fs(i,j)+k0*Cl(i,j)*Dfs(i,j))/(fs(i,j)+Dfs(i,j));
                                                                              %删去了fs0(i,j)=fs(i,j)+Dfs(i,j); 
                    cb=0;
                elseif  CL>C6
                     Cl0(i,j)=C6;
                     Cs0(i,j)=(Cs(i,j)*fs(i,j)+k0*Cl(i,j)*Dfs(i,j))/(fs(i,j)+Dfs(i,j));
                                                                              %删去了 fs0(i,j)=fs(i,j)+Dfs(i,j);
                     cb=(CL-C6)*(1-fs(i,j)-Dfs(i,j)); %界面凝固排出的溶质
                end
             end
             
            % 把排除的溶质扩散到周围的液相中
            XX=i+dx;
            YY=j+dy;
            XX=mod(XX-1,X)+1;
            YY=mod(YY-1,Y)+1;                                                                     %新添加的部分
     
            m=0;
            for im0=1:4
                if s0(XX(im0),YY(im0))==0
                    m=m+1;
                end
            end            
            for in0=1:4
                if s0(XX(in0),YY(in0))==0
                    Cl0(XX(in0),YY(in0))= Cl(XX(in0),YY(in0))+cb/m;
                end
            end
            %结束界面胞的生长
        end
        Cl=Cl0;
        Cs=Cs0;
        fs=fs0;
        s0=u;
    end
end
for i=2:X-1
    for j=2:Y-1
    %扩散模块 Diffuse
         if s0(i,j)==0
            m=0;
            XX=i+dx;
            YY=j+dy;
            XX=mod(XX-1,X)+1;
            YY=mod(YY-1,Y)+1; 
            for i0=1:4
                if s0(XX(i0),YY(i0))==0
                    coe(XX(i0),YY(i0))=1;
                else coe(XX(i0),YY(i0))=0;
                end
                m=m+coe(XX(i0),YY(i0));
            end
          Cl0(i,j)=Cl(i,j)+Dl*dt*(Cl(i-1,j)*coe(i-1,j)+Cl(i+1,j)*coe(i+1,j)+Cl(i,j+1)*coe(i,j+1)+Cl(i,j-1)*coe(i,j-1)-m*Cl(i,j))/ddx^2;
           %转变模块 transition module
         elseif s0(i,j)==1
            fs0(i,j)=fs(i,j)+Dfs(i,j);                                 %加上了%+Dfs(i,j)
            if fs(i,j)>=1
                fs0(i,j)=1;
                u(i,j)=4;
                XX=i+dx;
                YY=j+dy;
                XX=mod(XX-1,X)+1;
                YY=mod(YY-1,Y)+1; 
                for i0=1:4
                    if s0(XX(i0),YY(i0))==0
                        u(XX(i0),YY(i0))=1;
                        U1(XX(i0),YY(i0))=U(i,j);                                        %有错误，把m改成i0
                    end
                end
            end   
         end
     end
end
set(Ii,'CData',U);
set(ti,'string',['t= ',num2str(k)]);
pause(0.1);
k=k+1;
T0=T1;
Cl=Cl0;
Cs=Cs0;
fs=fs0;                                                                                   %把fs0赋给fs
s0=u;
end    
            
            