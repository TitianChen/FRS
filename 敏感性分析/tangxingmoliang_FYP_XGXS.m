clear
clc
N=2048;  %%FFT步数
pi=3.14159265;
fiv=fopen('rg160.txt','r');
fic=fopen('respons.txt','r');
fiw=fopen('sa.txt','w');
%%%%读入地震波
aag=fscanf(fiv,'%f %f',[2,inf]);
aag=aag';
tt=aag(:,1);  %时间
ag=aag(:,2);   %地震输入加速度
AAG=fft(ag,N);   %傅里叶变换，得到地震动输入各谐波成分的加速度峰值  ?????????
fs=1/(tt(2)-tt(1)); %采样频率
DDW=2.0*pi*fs/(N-1);
for ii=1:N/2
	TW(ii)=DDW*(ii-1);    %地震动输入各谐波成分频率值
end
for ii=N/2+1:N
    TW(ii)=-DDW*(N-ii)-1;
end
%%%%反应谱控制频率
freq=fscanf(fic,'%f',[1,inf]);  %w1-w106
nih=length(freq);
ccc=0.02; %%%反应谱阻尼值   
np = 500;%%%%%%
fredom=2;  %自由度数
vvs = 2163.2;%剪切波速置为2163.2m/s;
%%%%%结构材料参数
k0_k0=[4.1717E+10 -4.1717E+10;-4.1717E+10 4.1717E+10]; %%结构刚度阵
EEII = linspace(0.5,2,np);

kr0=[0 0;0 4.1717E+10]; %%固定地基，将1点刚度置0
m=[8.43E+06 0;0 5.37e7];
[vv0,dd0]=eig(kr0,m);
w_grid=abs(sqrt(dd0))/2/3.1415926;  %%%固定地基条件下结构的自振频率，自然频率，应等于4.436Hz
dens=2650;  %地基密度;
prxy=0.31;  %地基泊松比;  
E=zeros(fredom,1); 
E(:,1)=1;    %%%位移变换向量

%%开始大循环
for IVS=1:np
 vs=2163.2;
 m=[8.43E+06 0;0 5.37e7];
 ggdd=dens*vs*vs;  %%地基剪切模量
 k_g=32*(1-prxy)*19.8*ggdd/(7-8*prxy); %%弹簧阻尼器等效弹簧刚度
 c_g=0.576*19.8*k_g*sqrt(dens/ggdd);  %%弹簧阻尼器等效阻尼系数
 
 k0 = k0_k0*(EEII(IVS));
 
 k=k0;
 k(1,1)=k(1,1)+k_g;  %%整体刚度阵
 [vv,dd]=eig(k,m);   %%%结构地基体系条件下结构的自振频率，圆频率
 www0=abs(sqrt(dd)); 

 a1=2*0.05/(www0(1)+www0(2));
 a0=2*0.05/(www0(1)+www0(2))*(www0(1)*www0(2));
 c0=a0*m+a1*k;   %%结构瑞雷阻尼阵
 c=c0;
 c(1,1)=c0(1,1)+c_g;%%整体阻尼阵
 
 for ii=1:fredom
    for j=fredom:-1:1                  
       vv(j,ii)=vv(j,ii)/vv(1,ii);   %振型归一化
    end
    w0=www0(ii,ii);  
    ww00(IVS,ii)=w0/2/pi;
    mmode=vv(:,ii);
    M=mmode'*m*mmode;
    K=mmode'*k*mmode;
    C=mmode'*c*mmode;
    %damp=C/2/w0/M;
    damp=0.05;
    r=mmode'*m*E;
    p=r/M;   %%a(t)+2*cc*w0*v(t)+w0*w0*u(t)=-p*ag
    
    for kk=1:N
     w=TW(kk);    %谐波分量
     Hw(kk,ii)=p*w*w/(w0*w0-w*w+2*i*ccc*w0*w); %传递函数（相对）
     A0(kk)=Hw(kk,ii)*AAG(kk);
    end    
    A0=A0';
    U0(:,ii)=A0;    %%%广义坐标
 end
  for ii=1:N
   X(ii,1)=vv(1,1)*U0(ii,1)+vv(1,2)*U0(ii,2);   %%%结构频域响应相对加速度值
   X(ii,2)=vv(2,1)*U0(ii,1)+vv(2,2)*U0(ii,2);
  end 
  AAAG(:,1)=X(:,1)+AAG;  %%%绝对加速度响应
  AAAG(:,2)=X(:,2)+AAG;
  %%求对应的反应谱值
  for ii=1:fredom
   for ih=1:nih
    Freq=2*pi*freq(ih); 
    for kk=1:N
     w=TW(kk);    %谐波分量
     HHw(kk)=(Freq*Freq+i*2*ccc*Freq*w)/(Freq*Freq-w*w+2*i*ccc*Freq*w);   %传递函数（绝对）
     AAAW(kk)=HHw(kk)*AAAG(kk,ii);
    end
    AAT=ifft(AAAW);
    AATR=abs(real(AAT));
    AATI=imag(AAT);
    AMAX=AATR(1); 
    for kk=2:N
     if AATR(kk)>AMAX
       AMAX=AATR(kk);
     end
    end
    Sa(ih,ii)=AMAX;
   end  
  end
  Sa_1(:,IVS)=Sa(:,1)/9.81;
  Sa_2(:,IVS)=Sa(:,2)/9.81;
end


b=k0_k0(1,1)*EEII;

for ii=1:81
a(ii)=freq(ii);  %w
z(ii,:)=Sa_2(ii,:);
end
%a:频率;
%z:反应谱
%plot(a,z);
[xx,yy]=meshgrid(b,a);
surf(xx,yy,z);
hold on
shading interp

   
   

