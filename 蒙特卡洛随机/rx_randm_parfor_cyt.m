%��������ͳ�� ȷ����  �Ŵ�1.5��2���ĳɹ�
%����������������� RXģ�� ASCE4-98�ػ�����ģ�� ˲̬��Ӧ����
%������ʼ��
clear
clc

nwav=2800;         %����    ****
aahh0=0.1*9.81;    %��Ƶ���ֵ  ****

vs_ave=800.2;        %ȷ���� ���в���   2163.2m/s  ****

skce=zeros(12,12);

eexi=zeros(5,1);  %����¥���յ�����ֵ �̶�Ϊ5��ֵ
eexi(1)=0.005;     %****
eexi(2)=0.02;      %****
eexi(3)=0.05;      %****
eexi(4)=0.07;      %****
eexi(5)=0.1;       %****

%������Ԫ��ϵ  AA/SAX ���Ǽ������
%R,1,1204,115400,115400,39.16,39.16,0,   
%     AA   IX     IY
%RMORE,0,1000000,1.111,1.111,0,0,
%          IXX     SkAX  SkAY
% faix=12*E*IX*SkAY/G/A/l/l   ��uyӰ��
% faiy=12*E*IY*SkAX/G/A/l/l   ��uxӰ��
% E=G*2*(1+u)
% 1-1 1-5 1-7 1-11  7-1 7-5 7-7 7-11 E.AA*IY/(1+faiy)
% 2-2 2-4 2-8 2-10  8-2 8-4 8-8 8-10 E.AA*IX/(1+faix)
% 3-3 3-9 9-3 9-9   E.AA
% 4-2 4-8 10-2 10-8 E.AA*IX/(1+faix)  4-4 10-10 E.AA.IX*(4+faix)/(1+faix)
% 4-10 10-4  E.AA.IX*(2-faix)/(1+faix)
% 5-1 5-7 11-1 11-7 E.AA*IY/(1+faiy)  5-5 11-11 E.AA.IY*(4+faiy)/(1+faiy)
% 5-11 11-5  E.AA.IY*(2-faiy)/(1+faiy)
% 6-6 6-12 12-6 12-12 E.IXX

%������� ���ٶ�ʱ�� ĳһ��ˮƽX����
fiv=fopen('rg160x.txt','r');%ˮƽ�����;
kcm=fscanf(fiv,'%d %f',[2,1]);
nstep=kcm(1);   %���𲨵ĵ��
if nwav>nstep
    nwav=nstep; %nwav����ʱ���
end
deltk=kcm(2);   %ʱ�䲽��
ecm1=fscanf(fiv,'%f %f',[2,nstep]);
ecm1=ecm1';
fclose(fiv);
   %���ڵ�����ٶȷ�ֵ
   amx=0;
   for i=1:nstep
      aii=ecm1(i,1)-(i-1)*deltk;
      if aii>1e-3
          disp('���� aii>1e-3');          
      end
      if abs(ecm1(i,2)) > amx
          amx=abs(ecm1(i,2));
      end
   end
amx=aahh0/amx;       %������𲨷�ֵ��Ϊaahh0 ��λm/s2
for i=1:nstep
    ecm1(i,2)=ecm1(i,2)*amx;
end
%-- Ϊ���ͳ�������
yaa=zeros(nwav,1);
exi=0.02;
for i=1:nwav
    yaa(i)=ecm1(i,2);
end
saatmp=calresp(nwav,deltk,exi,yaa);

%�������е������������...
fir=fopen('distb_all.txt','r');
  %fprintf(fiw,'%8d %8d %8d %8d\n',np,nnp,npl,numnd);
kcm=fscanf(fir,'%d %d %d %d',[4,1]);
nprnd=kcm(1);  %��ʵ���������
nnpoornd=kcm(2); %��10�ı�����ĸ���
nnpllrnd=kcm(3); %��10������
numrnd=kcm(4);   %�������������
vsrd=zeros(nprnd,numrnd);
for ikk=1:numrnd
    disp(['ikk=',num2str(ikk)]);
    %fprintf(fiw,'%5d %13.5e %13.5e %13.5e %13.5e %5d\n',ikk,radmp(ikk,np+1),radmp(ikk,np+2),radmp(ikk,np+3),radmp(ikk,np+4),radmp(ikk,np+5));
    kkcm=fscanf(fir,'%d %f %f %f %f %d',[6,1]);
    ikk1=kkcm(1);  %�ڼ����������
    if ikk1 ~= ikk
        disp('���� ikk1 ~= ikk');
        return;
    end
    for iil=1:nnpllrnd
        na=(iil-1)*10;
        kcm0=fscanf(fir,'%f %f %f %f %f %f %f %f %f %f',[10,1]);  %10��1��
        for i=1:10
            if na+i > nprnd
                break
            end
            vsrd(na+i,ikk)=kcm0(i);  %��¼
            if kcm0(i) <=0
                disp('���� kcm0(i) <=0');
                return;
            end
        end
    end %iil
end %ikk
fclose(fir);

%��������
fid=fopen('shuju_all.txt','r');
tcm=fscanf(fid,'%d %d',[2,1]);
np=tcm(1);  %�ڵ���
ne=tcm(2);  %��Ԫ��
tcm1=fscanf(fid,'%d %f %f %f',[4,np]);
tcm1=tcm1';  %�ڵ�����
tcm2=fscanf(fid,'%d %d',[2,ne]);
tcm2=tcm2';  %��Ԫ�ڵ��
tcmpp=fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11,ne]);
tcmpp=tcmpp';  %��Ԫ ȷ���Բ��� ia 1,ib 2,e 3,u 4,kexi 5,area 6,ix 7,iy 8,skx 9,sky 10, siipp 11
tcm3=fscanf(fid,'%d %f %f %f %f %f %f',[7,12]);
tcm3=tcm3';  %�ʵ����
skkele=zeros(12,12,ne);
ngg=np*6;  %�����ɶ�
for i=1:np
    i1=tcm1(i,1);
    i2=tcm3(i,1);
    if i1~=i
        disp('����1');
    end 
    if i2~=i
        disp('����2');
    end     
end
fclose(fid);

%-- ����һ�¸�����Ԫ�ĳ���
tcml=zeros(ne,1);
for i=1:ne
    i1=tcm2(i,1);
    i2=tcm2(i,2);
    tcml(i)=abs(tcm1(i2,4)-tcm1(i1,4));
end

iclc=1;   %---------------------------- ��ѭ�� ��ѭ��
nprnd5=nprnd+5;
allrsp=zeros(1500,nprnd5);   %5������
arfhh=zeros(10,nprnd5);
tic
%tStart=tic;
%matlabpool close;
% matlabpool(4);      %parallel using 4 threads

for iclc=1:nprnd5
    
    tmpi0=zeros(1500,1);
    tmpi1=zeros(10,1);
    
    %disp(iclc);

% %     tlaps=toc(tStart);
%         fidc=fopen(['cout/fil' num2str(iclc) ],'w');  %shanchu  ****
% %        fidc=fopen(['filcc.txt'],'wt');
% %        %fprintf(fidc,'%8d\n',numct);
% %        fprintf(fidc,'%13.4e\n',tlaps);
%        fclose(fidc);
%     %length(dir('cout'));

 
    %-- ����
    ngg=np*6;  %�����ɶ�
    skk=zeros(ngg,ngg);
    smm=zeros(ngg,ngg);
    
    %-- ������� vsrd=zeros(nprnd,numrnd);
   %���в��� 1  vs_rd
   %�ṹ���ϵ���ģ��  2  e_rd
   %�ṹ�������� 3 damp_rd
   %�ʵ������(���ܶ�) 4-15  12���� ת������(ͬ��) mass_rd
   %����Ԫ�ļ������ 11���� 16-26 sa_rd
   %����Ԫ�Ĺ��Ծ� 11���� 27-37 ������һ�� ii_rd   
    mass_rd=zeros(12,1);
    sa_rd=zeros(11,1);
    ii_rd=zeros(11,1);
   if iclc<=nprnd   %�������
      vs_rd=vsrd(iclc,1);
      e_rd=vsrd(iclc,2);
      damp_rd=vsrd(iclc,3);
      for i=1:12
         mass_rd(i)=vsrd(iclc,3+i);  %�൱���ܶ�
      end
      for i=1:11
         sa_rd(i)=vsrd(iclc,15+i);
         ii_rd(i)=vsrd(iclc,26+i);
      end
   else
      %ȷ��+1, �ػ�ģ��1/1.5 +2, 1.5 +3, 1/2 +4, 2 +5
      vs_rd=1;
      e_rd=1;
      damp_rd=1;
      for i=1:12
         mass_rd(i)=1;  %�൱���ܶ�
      end
      for i=1:11
         sa_rd(i)=1;
         ii_rd(i)=1;
      end

      if iclc==nprnd+2  %�ػ�����ģ������1.5��  1/1.5
          vs_rd=1/sqrt(1.5);
      end
      if iclc==nprnd+3  %�ػ�����ģ������1.5��
          vs_rd=sqrt(1.5);
      end
      if iclc==nprnd+4  %�ػ�����ģ������2��   1/2
          vs_rd=sqrt(0.5);
      end
      if iclc==nprnd+5  %�ػ�����ģ������2��
          vs_rd=sqrt(2.0);
      end       
   end
    
    %-- ��Ԫie������  tcm3(12,7) ������������
    for ii=1:np
        i1=6*(tcm3(ii,1)-1); %������ɶ�
        for i=1:6
            smm(i1+i,i1+i)=smm(i1+i,i1+i)+tcm3(ii,i+1)*mass_rd(ii);
        end  %i   
    end %ii

    %-- ��Ԫie�ն��� skkele(12,12,ne)
    for ie=1:ne
        
%         if ie==ne
%             ie=ne;
%         end
        
        i1=6*(tcm2(ie,1)-1);
        i2=6*(tcm2(ie,2)-1);
        
        %��Ԫ ȷ���Բ��� ia 1,  ib 2,  e 3,  u 4,  kexi 5,  area 6,  ix 7,  iy 8,  skx 9,  sky 10 
% faix=12*E*IX*SkAY/G/A/l/l   ��uyӰ��
% faiy=12*E*IY*SkAX/G/A/l/l   ��uxӰ��
% E=G*2*(1+u)
% 1-1 1-5 1-7 1-11  7-1 7-5 7-7 7-11 E.AA*IY/(1+faiy)
% 2-2 2-4 2-8 2-10  8-2 8-4 8-8 8-10 E.AA*IX/(1+faix)
% 3-3 3-9 9-3 9-9   E.AA
% 4-2 4-8 10-2 10-8 E.AA*IX/(1+faix)  4-4 10-10 E.AA.IX*(4+faix)/(1+faix)
% 4-10 10-4  E.AA.IX*(2-faix)/(1+faix)
% 5-1 5-7 11-1 11-7 E.AA*IY/(1+faiy)  5-5 11-11 E.AA.IY*(4+faiy)/(1+faiy)
% 5-11 11-5  E.AA.IY*(2-faiy)/(1+faiy)
% 6-6 6-12 12-6 12-12 E.IXX
        esl=tcml(ie);        %��Ԫ����
        ee=tcmpp(ie,3)*e_rd;        %����ģ��
        emu=tcmpp(ie,4);        %���ɱ�
        ekexi=tcmpp(ie,5)*damp_rd;     %�����
        eara=tcmpp(ie,6);     %��� 
        eix=tcmpp(ie,7)*ii_rd(ie);       %ix
        eiy=tcmpp(ie,8)*ii_rd(ie);       %iy
        eskx=tcmpp(ie,9)/sa_rd(ie);      %skx
        esky=tcmpp(ie,10)/sa_rd(ie);     %skx
        eiip=tcmpp(ie,11)*ii_rd(ie);         %ip ����ֱz���ת������
        
        if ie==1
            ekeximin=ekexi;
        else 
            if ekeximin>ekexi
                ekeximin=ekexi;
            end
        end
                
        faiy=24*(1+emu)*eiy*eskx/eara/esl/esl;  %��ux����
        faix=24*(1+emu)*eix*esky/eara/esl/esl; 
        
        scle=zeros(12,12);  %���ڴ洢��Ԫ�ն���
        scle(1,1)=12*ee*eiy/(1+faiy)/esl/esl/esl;        %ux1
        scle(1,5)=6*ee*eiy/(1+faiy)/esl/esl;
        scle(1,7)=-12*ee*eiy/(1+faiy)/esl/esl/esl; 
        scle(1,11)=6*ee*eiy/(1+faiy)/esl/esl;
        
        scle(7,7)=12*ee*eiy/(1+faiy)/esl/esl/esl;  %ux2
        scle(7,11)=-6*ee*eiy/(1+faiy)/esl/esl;
        
        scle(2,2)=12*ee*eix/(1+faix)/esl/esl/esl;    %uy1
        scle(2,4)=-6*ee*eix/(1+faix)/esl/esl; 
        scle(2,8)=-12*ee*eix/(1+faix)/esl/esl/esl; 
        scle(2,10)=-6*ee*eix/(1+faix)/esl/esl;
        
        scle(8,8)=12*ee*eix/(1+faix)/esl/esl/esl;  %uy2
        scle(8,10)=6*ee*eix/(1+faix)/esl/esl;
        
        scle(3,3)=ee*eara/esl;     %uz1
        scle(3,9)=-ee*eara/esl; 
  
        scle(9,9)=ee*eara/esl;     %uz2
        
        scle(4,4)=(4+faix)*ee*eix/(1+faix)/esl;   %ct x1
        scle(4,8)=6*ee*eix/(1+faix)/esl/esl; 
        scle(4,10)=(2-faix)*ee*eix/(1+faix)/esl;
        
        scle(10,10)=(4+faix)*ee*eix/(1+faix)/esl;   %ct x2
         
        scle(5,5)=(4+faiy)*ee*eiy/(1+faiy)/esl;    %ct y1
        scle(5,7)=-6*ee*eiy/(1+faiy)/esl/esl;
        scle(5,11)=(2-faiy)*ee*eiy/(1+faiy)/esl; 
        
        scle(11,11)=(4+faiy)*ee*eiy/(1+faiy)/esl;  %ct y2
        
        scle(6,6)=ee/2/(1+emu)*eiip/esl;   %ct z1
        scle(6,12)=-ee/2/(1+emu)*eiip/esl; 
        
        scle(12,12)=ee/2/(1+emu)*eiip/esl;  %ct z2
        
        for i=1:12
            for k=1:i-1
                scle(i,k)=scle(k,i);
            end
        end
        
        for i=1:6
            for k=1:6
                %i1 i1  i k
                skk(i1+i,i1+k)=skk(i1+i,i1+k)+scle(i,k);
                %i2 i2  6+i 6+k
                skk(i2+i,i2+k)=skk(i2+i,i2+k)+scle(6+i,6+k);
                %i1 i2  i 6+k
                skk(i1+i,i2+k)=skk(i1+i,i2+k)+scle(i,6+k);
                %i2 i1  6+i k
                skk(i2+i,i1+k)=skk(i2+i,i1+k)+scle(6+i,k);
            end
        end
        
 
    end


    
    %-- ����̶��ػ���һ������Ƶ��
    %Լ������ �޳�1-6�к��� 1�ŵ�����ǰ��
    nggn=6*(np-1);
    skkn=zeros(nggn,nggn);
    smmn=zeros(nggn,nggn);
    for i=1:nggn
    for k=1:nggn
        skkn(i,k)=skk(6+i,6+k);
        smmn(i,k)=smm(6+i,6+k);
    end
    end

    %�����������ֵ+��������
    % [V,D]=eig(A,B)��
    %  AV=BVD ����V���������� (B-1/D*A)V=0   V��������{},{},{}... �Ѿ���֤
    %       (skkn-w2*smmn)V=0
    [vv1,dd1]=eig(smmn,skkn);
    fhz=zeros(nggn,1);
    for i=1:nggn
      fhz(i)=sqrt(1/dd1(i,i))/2/3.1415926;  %Hz ֻ��˳�򲻶� �̶��ػ�
    end

    %��������ֵ����������������(ֻ���м���!) �����ʵ����ɶ����Ӧ
    %��ͬƵ��ֵ �������ڵ��x\y\z����ֵ˳������ ��fhz()�е�hz��С��������
    for ii=1:nggn
    for jj=ii+1:nggn
        if fhz(ii)>fhz(jj)
            cztp=fhz(ii);
            fhz(ii)=fhz(jj);
            fhz(jj)=cztp;
            
            for k=1:nggn
                cztp=vv1(k,ii);
                vv1(k,ii)=vv1(k,jj);
                vv1(k,jj)=cztp;
            end
        end
    end
    end

    %��֤������ �̶��ػ� 1����10��
    cfh1=zeros(nggn,1);
    cfh2=zeros(nggn,1);
    for i=1:nggn
      cfh1(i)=vv1(i,1);
      cfh2(i)=vv1(i,10);
    end
    aaxm=cfh1'*skkn*cfh2;     %��֤��
    aaxm00=abs(aaxm);
    f1=fhz(1)                 %1��Ƶ�� hz
    
    jxx=0;
    jzz=0;  %�������Ťת���Ծ�
    for i=1:np
      ia=(i-1)*6;
      jxx=jxx+smm(ia+4,ia+4);
      jxx=jxx+smm(ia+1,ia+1)*(tcm1(i,4)-tcm1(1,4))^2;
      jzz=jzz+smm(ia+6,ia+6);
    end
    jyy=jxx;

   %=== ʩ�ӵ����������ػ�ģ�� ===

   %a.1  �ػ����ϵĹ̶����� ****
   prxy=0.31;
   dens=2650;
   ggdd=dens*vs_rd*vs_ave*vs_rd*vs_ave;  %��������в��ټ������ģ��
   rr0=19.8;

   %a.2 �ػ����ϵ��µĵ���������������
   kkx=32*(1-prxy)*rr0*ggdd/(7-8*prxy);
   kkz=4*rr0*ggdd/(1-prxy);
   kkf=8*rr0^3*ggdd/(1-prxy)/3;
   kkr=16*rr0^3*ggdd/3;

   betaf=3*(1-prxy)*jxx/(8*dens*rr0^5);
   ckx=0.576*rr0*kkx*sqrt(dens/ggdd);
   ckz=0.85*rr0*kkz*sqrt(dens/ggdd);
   ckf=0.3*rr0*kkf/(1+betaf)*sqrt(dens/ggdd);
   ckr=sqrt(kkr*jzz)/(1+2*jzz/(dens*rr0^5));

   kky=kkx;
   kkfx=kkf;
   kkfy=kkf;
   cky=ckx;
   ckfx=ckf;
   ckfy=ckf;
   
   %a.3 ���ӵ�ģ�͵ĸն�����
   nggn=6*np;
   skka=skk;
   smma=smm;
   skka(1,1)=skka(1,1)+kkx;   %x
   skka(2,2)=skka(2,2)+kky;   %y
   skka(3,3)=skka(3,3)+kkz;   %z
   skka(4,4)=skka(4,4)+kkfx;  %ҡ��fx
   skka(5,5)=skka(5,5)+kkfy;  %ҡ��fy
   skka(6,6)=skka(6,6)+kkr;   %Ťתrz

   [vv2,dd2]=eig(smma,skka);
   fhza=zeros(nggn,1);
   for i=1:nggn
      fhza(i)=sqrt(1/dd2(i,i))/2/3.1415926;  %���Եػ� Hz ֻ��˳�򲻶� �̶��ػ�
      %fhza(i)=abs(fhza(i));
   end
   
   for ii=1:nggn
      for jj=ii+1:nggn
        if fhza(ii)>fhza(jj)
            cztp=fhza(ii);
            fhza(ii)=fhza(jj);
            fhza(jj)=cztp;
            
            for k=1:nggn
                cztp=vv2(k,ii);
                vv2(k,ii)=vv2(k,jj);
                vv2(k,jj)=cztp;
            end
         end
       end
   end
   
    %��֤������ ���Եػ� 1����10��
    cfha=zeros(nggn,1);
    cfhb=zeros(nggn,1);
    for i=1:nggn
      cfha(i)=vv2(i,1);
      cfhb(i)=vv2(i,10);
    end
    aaxm=cfha'*skka*cfhb;     %��֤�� ****
    aaxm11=abs(aaxm);
    f1c=fhza(1)   
   
    tmpi1(1)=real(fhz(1));
    tmpi1(2)=imag(fhz(1));
    tmpi1(3)=real(fhz(20));
    tmpi1(4)=imag(fhz(20));
    tmpi1(5)=real(fhza(1));
    tmpi1(6)=imag(fhza(1));
    tmpi1(7)=real(fhza(20));
    tmpi1(8)=imag(fhza(20));   
    tmpi1(9)=aaxm00;
    tmpi1(10)=aaxm11;    
    
    arfhh(:,iclc)=tmpi1;
      
   %a.4 ȡ1��20�׵�Ƶ�ʼ���������[C]
   %  ��=2*��1*��2*��/(��1+��2) 
   %  ��=2*��/(��1+��2)
   ww1=2*3.1415926*abs(fhza(1));
   ww2=2*3.1415926*abs(fhza(20));
   alfa=2*ww1*ww2*ekeximin/(ww1+ww2);
   beta=2*ekeximin/(ww1+ww2);           %Ӧ�ð���ʵ�Ĳ������� ������õ�Ԫ�������� ????

   scca=alfa*smma+beta*skk;   %������ (�ն��󲿷ֲ�Ҫ���ǵػ������򲿷���)
   scca(1,1)=scca(1,1)+ckx;   %x
   scca(2,2)=scca(2,2)+cky;   %y
   scca(3,3)=scca(3,3)+ckz;   %z
   scca(4,4)=scca(4,4)+ckfx;  %ҡ��fx
   scca(5,5)=scca(5,5)+ckfy;  %ҡ��fy
   scca(6,6)=scca(6,6)+ckr;   %Ťתrz

   %a.5 Newmarkϵ��
   ala=0.5;
   blb=0.25;
   a0k=1/(blb*deltk*deltk);
   a1k=ala/(blb*deltk);
   a2k=1/(blb*deltk);
   a3k=1/(2*blb)-1;
   a4k=ala/blb-1;
   a5k=deltk/2*(ala/blb-2);
   a6k=deltk*(1-ala);
   a7k=ala*deltk;

   %a.6 �����Ķ������� ����ʱ�̷���
   sktp=skka+a0k*smma+a1k*scca; 
   skin=inv(sktp);

   xa=zeros(nwav,1);   %ʱ��
   ya2=zeros(nwav,1);    %���ٶ���Ӧ 2�ŵ�
   ya8=zeros(nwav,1);    %���ٶ���Ӧ 8�ŵ�
   ya12=zeros(nwav,1);   %���ٶ���Ӧ 12�ŵ�
   ut0=zeros(nggn,1);
   ut1=zeros(nggn,1);
   vt0=zeros(nggn,1);
   vt1=zeros(nggn,1);
   at0=zeros(nggn,1);
   at1=zeros(nggn,1);

   iddk=zeros(nggn,1);
   for i=1:np
      iddk(6*i-6+1)=1;   %ָ��x��
   end

   for istp=1:nwav
      xa(istp)=(istp-1)*deltk;
    
      %disp(['istp=',num2str(istp)]);
    
      w1=a0k*ut0+a2k*vt0+a3k*at0;
      w2=a1k*ut0+a4k*vt0+a5k*at0;
      xxaagg=ecm1(istp,2);    
      rri=-smma*iddk*xxaagg+smma*w1+scca*w2;         %****
    
      ut1=skin*rri;
    
      at1=a0k*(ut1-ut0)-a2k*vt0-a3k*at0;
      vt1=vt0+a6k*at0+a7k*at1;
    
      ya2(istp)=at1(6*2-6+1)+xxaagg;      %2�ŵ��x�� ���Լ��ٶ���Ӧ   ****
      ya8(istp)=at1(6*8-6+1)+xxaagg;      %8�ŵ��x�� ���Լ��ٶ���Ӧ   ****
      ya12(istp)=at1(6*12-6+1)+xxaagg;    %12�ŵ��x�� ���Լ��ٶ���Ӧ  ****
    
      ut0=ut1;
      vt0=vt1;
      at0=at1;    
   end

   %a.7 ����ya2��ya8��ya12�ڵ�¥����

   for iexi=1:5
       exi=eexi(iexi);
       ssaa2=calresp(nwav,deltk,exi,ya2);  %����ssaa(nih,3) ��Ӧ�� 1Ƶ�� 2���Լ��ٶȷ�Ӧ��(m/s2)  ssaa2(1,3)Ϊnihֵ
       ssaa8=calresp(nwav,deltk,exi,ya8);
       ssaa12=calresp(nwav,deltk,exi,ya12);
       nih=ssaa2(1,3);
       iaa=(iexi-1)*300;
       for i=1:nih
          tmpi0(iaa+i)=ssaa2(i,2); 
          tmpi0(iaa+100+i)=ssaa8(i,2); 
          tmpi0(iaa+200+i)=ssaa12(i,2); 
       end
       allrsp(:,iclc)=tmpi0;
   end %iexi
 
end %iclc
% matlabpool close;   %close all the matlabpools
toc

ficcc=fopen('fhz1_20.txt','wt');
for iclc=1:nprnd+5
    f11a=arfhh(1,iclc);
    f11b=arfhh(2,iclc);
    f22a=arfhh(3,iclc);
    f22b=arfhh(4,iclc);
    f33a=arfhh(5,iclc);
    f33b=arfhh(6,iclc);
    f44a=arfhh(7,iclc);
    f44b=arfhh(8,iclc);
    aaxm00=arfhh(9,iclc);
    aaxm11=arfhh(10,iclc);
    fprintf(ficcc,' %7d | %13.4e %13.4e | %13.4e %13.4e ||  %13.4e %13.4e | %13.4e %13.4e ||  %13.4e %13.4e \n',iclc,f11a,f11b,f22a,f22b,f33a,f33b,f44a,f44b,aaxm00,aaxm11);
end
fclose(ficcc);

   fibb=fopen('rspp_all3s.txt','wt');   %�ṹ2��8��12��ľ��Լ��ٶ���Ӧ
   nih=saatmp(1,3);
   osp=zeros(10,1);
   fprintf(fibb,'%5d%8d\n',nih,nprnd);
   for i=1:nih
       fprintf(fibb,'%13.5f\n',saatmp(i,1));
   end
   for iexi=1:5
       iaa=(iexi-1)*300;
       fprintf(fibb,'%13.7f\n',eexi(iexi));
       for iclc=1:nprnd+5
           %300��Ԫ��
           for il=1:30
               na=(il-1)*10;
               for i=1:10
                   osp(i)=allrsp(iaa+na+i,iclc);
               end
               fprintf(fibb,'%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n', ...
                  osp(1),osp(2),osp(3),osp(4),osp(5),osp(6),osp(7),osp(8),osp(9),osp(10));
           end
       end %iclc
   end %iexi
   fclose(fibb);
   
   nprnd
  







