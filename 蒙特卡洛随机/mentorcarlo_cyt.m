clear
clc
N=2048;  %%FFT����
pi=3.14159265;
fiv=fopen('rg160.txt','r');
fic=fopen('respons.txt','r');
fiw=fopen('sa.txt','w');
%%%%�������
aag=fscanf(fiv,'%f %f',[2,inf]);
aag=aag';
tt=aag(:,1);  %ʱ��
ag=aag(:,2);   %����������ٶ�
AAG=fft(ag,N);   %����Ҷ�任���õ����������г���ɷֵļ��ٶȷ�ֵ  ?????????
fs=1/(tt(2)-tt(1)); %����Ƶ��
DDW=2.0*pi*fs/(N-1);
for ii=1:N/2
	TW(ii)=DDW*(ii-1);    %���������г���ɷ�Ƶ��ֵ
end
for ii=N/2+1:N
    TW(ii)=-DDW*(N-ii)-1;
end

%%%%��Ӧ�׿���Ƶ��
freq=fscanf(fic,'%f',[1,inf]);  %w1-w106
nih=length(freq);
ccc=0.02; %%%��Ӧ������ֵ   

fredom=2;  %���ɶ���
%%%��ȷ���Ա���Ϊ�ػ����в���
%{
np=200;%��Ҫ�ģ���
vs_ave=1;           %1
vs_gam=0.2/1.282;    %��һ����ֵ 90%�ı�֤��:1.282 %95%��֤��;:1.645;
vvs=normrnd(vs_ave,vs_gam,[1 np]);%��ֵΪvs_ave ��׼��Ϊvx_gam
vs0=2163.2;%���в���ȷ��ֵ;
vvs=vs0*vvs;
for j=1:np      %%���򣬼��в�����С����
 for k=j+1:np
  if vvs(k)<vvs(j)
     v=vvs(k);
     vvs(k)=vvs(j);
     vvs(j)=v;
  end
 end
end
vvs(np+1)=vs0;
np = np+1;
%}

np=20;      
vvs=linspace(500,2500,np);%���в���;
%np = 1;
%vvs = 2163.2;%���в�����Ϊ2163.2m/s;


%%%%%�ṹ���ϲ���
k0=[4.1717E+10 -4.1717E+10;-4.1717E+10 4.1717E+10]; %%�ṹ�ն���
kr0=[0 0;0 4.1717E+10]; %%�̶��ػ�����1��ն���0
m=[8.43E+06 0;0 5.37e7];
[vv0,dd0]=eig(kr0,m);
w_grid=abs(sqrt(dd0))/2/3.1415926;  %%%�̶��ػ������½ṹ������Ƶ�ʣ���ȻƵ�ʣ�Ӧ����4.436Hz
dens=2650;  %�ػ��ܶ�;
prxy=0.31;  %�ػ����ɱ�;  
E=zeros(fredom,1); 
E(:,1)=1;    %%%λ�Ʊ任����


%%��ʼ��ѭ��
for IVS=1:np
 vs=vvs(IVS);
 ggdd=dens*vs*vs;  %%�ػ�����ģ��
 k_g=32*(1-prxy)*19.8*ggdd/(7-8*prxy); %%������������Ч���ɸն�
 c_g=0.576*19.8*k_g*sqrt(dens/ggdd);  %%������������Ч����ϵ��

 k=k0;
 k(1,1)=k(1,1)+k_g;  %%����ն���
 [vv,dd]=eig(k,m);   %%%�ṹ�ػ���ϵ�����½ṹ������Ƶ�ʣ�ԲƵ��
 www0=abs(sqrt(dd)); 

 a1=2*0.05/(www0(1)+www0(2));
 a0=2*0.05/(www0(1)+www0(2))*(www0(1)*www0(2));
 c0=a0*m+a1*k;   %%�ṹ����������
 c=c0;
 c(1,1)=c0(1,1)+c_g;%%����������
 
 for ii=1:fredom
    for j=fredom:-1:1                  
       vv(j,ii)=vv(j,ii)/vv(1,ii);   %���͹�һ��
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
     w=TW(kk);    %г������
     Hw(kk,ii)=p*w*w/(w0*w0-w*w+2*i*ccc*w0*w); %���ݺ�������ԣ�
     A0(kk)=Hw(kk,ii)*AAG(kk);
    end    
    A0=A0';
    U0(:,ii)=A0;    %%%��������
  end
   
  for ii=1:N
   X(ii,1)=vv(1,1)*U0(ii,1)+vv(1,2)*U0(ii,2);   %%%�ṹƵ����Ӧ��Լ��ٶ�ֵ
   X(ii,2)=vv(2,1)*U0(ii,1)+vv(2,2)*U0(ii,2);
  end 
  AAAG(:,1)=X(:,1)+AAG;  %%%���Լ��ٶ���Ӧ
  AAAG(:,2)=X(:,2)+AAG;
 
  %%���Ӧ�ķ�Ӧ��ֵ
  for ii=1:fredom
   for ih=1:nih
    Freq=2*pi*freq(ih); 
    for kk=1:N
     w=TW(kk);    %г������
     HHw(kk)=(Freq*Freq+i*2*ccc*Freq*w)/(Freq*Freq-w*w+2*i*ccc*Freq*w);   %���ݺ��������ԣ�
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
 
%%���
fprintf(fiw,'%s  ','Ƶ��\����');
for IVS=1:np
 fprintf(fiw,'%10.2f   ',vvs(IVS));
end
 fprintf(fiw,'\n');

for ii=1:106
 fprintf(fiw,'%8.3f',freq(ii));
 for IVS=1:np
  fprintf(fiw,'%13.4f',Sa_2(ii,IVS))/9.81;
 end  
  fprintf(fiw,'\n');
end

%��άͼ
b=vvs;  %vs
%b=b';
%a=freq;  %w
%a=a';
%z=Sa_2;
for ii=1:81
a(ii)=freq(ii);  %w
z(ii,:)=Sa_2(ii,:);
end

%a:Ƶ��;
%z:��Ӧ��
%plot(a,z);


[xx,yy]=meshgrid(b,a);
surf(xx,yy,z);
hold on
shading interp


% 
% b=vvs;  %vs
% %b=b';
% a=freq;  %w
% %a=a';
% z2=Sa_2;
% [xx,yy]=meshgrid(b,a);
% surf(xx,yy,z2);
% hold on
% shading interp
% 
% 
% 
figure
 b=vvs;  %vs
 %b=b';
 nn=[19 29 33 34 35 36 37 38 41 43 46 50];
 mm=length(nn);
 SSA2=zeros(mm,np);
 k=0;
 for ii=1:mm
  nnn=nn(ii);
  k=k+1;
  %a1(k)=freq(nnn);  %w
  %SSA2(k,:)=Sa_2(nnn,:);
  a1=freq(nnn); 
  SSA2=Sa_2(nnn,:);
  z1=a1*ones(1,length(b));
  plot3(z1,b,SSA2)
  hold on
 end


%plot(vvs,ww00(:,1)/2/pi)
%plot(vvs,ww00(:,2)/2/pi)




   
   
   

