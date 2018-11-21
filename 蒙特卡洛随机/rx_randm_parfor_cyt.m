%包含概率统计 确定性  放大1.5和2倍的成果
%所有随机参数条件下 RX模型 ASCE4-98地基动力模型 瞬态响应分析
%参数初始化
clear
clc

nwav=2800;         %地震波    ****
aahh0=0.1*9.81;    %设计地震动值  ****

vs_ave=800.2;        %确定性 剪切波速   2163.2m/s  ****

skce=zeros(12,12);

eexi=zeros(5,1);  %计算楼层普的阻尼值 固定为5个值
eexi(1)=0.005;     %****
eexi(2)=0.02;      %****
eexi(3)=0.05;      %****
eexi(4)=0.07;      %****
eexi(5)=0.1;       %****

%和梁单元关系  AA/SAX 这是剪切面积
%R,1,1204,115400,115400,39.16,39.16,0,   
%     AA   IX     IY
%RMORE,0,1000000,1.111,1.111,0,0,
%          IXX     SkAX  SkAY
% faix=12*E*IX*SkAY/G/A/l/l   对uy影响
% faiy=12*E*IY*SkAX/G/A/l/l   对ux影响
% E=G*2*(1+u)
% 1-1 1-5 1-7 1-11  7-1 7-5 7-7 7-11 E.AA*IY/(1+faiy)
% 2-2 2-4 2-8 2-10  8-2 8-4 8-8 8-10 E.AA*IX/(1+faix)
% 3-3 3-9 9-3 9-9   E.AA
% 4-2 4-8 10-2 10-8 E.AA*IX/(1+faix)  4-4 10-10 E.AA.IX*(4+faix)/(1+faix)
% 4-10 10-4  E.AA.IX*(2-faix)/(1+faix)
% 5-1 5-7 11-1 11-7 E.AA*IY/(1+faiy)  5-5 11-11 E.AA.IY*(4+faiy)/(1+faiy)
% 5-11 11-5  E.AA.IY*(2-faiy)/(1+faiy)
% 6-6 6-12 12-6 12-12 E.IXX

%读入地震波 加速度时程 某一个水平X方向
fiv=fopen('rg160x.txt','r');%水平向地震波;
kcm=fscanf(fiv,'%d %f',[2,1]);
nstep=kcm(1);   %地震波的点号
if nwav>nstep
    nwav=nstep; %nwav计算时间段
end
deltk=kcm(2);   %时间步长
ecm1=fscanf(fiv,'%f %f',[2,nstep]);
ecm1=ecm1';
fclose(fiv);
   %调节地震加速度峰值
   amx=0;
   for i=1:nstep
      aii=ecm1(i,1)-(i-1)*deltk;
      if aii>1e-3
          disp('错误 aii>1e-3');          
      end
      if abs(ecm1(i,2)) > amx
          amx=abs(ecm1(i,2));
      end
   end
amx=aahh0/amx;       %输入地震波峰值改为aahh0 单位m/s2
for i=1:nstep
    ecm1(i,2)=ecm1(i,2)*amx;
end
%-- 为最后统计输出用
yaa=zeros(nwav,1);
exi=0.02;
for i=1:nwav
    yaa(i)=ecm1(i,2);
end
saatmp=calresp(nwav,deltk,exi,yaa);

%读入所有的随机样本数据...
fir=fopen('distb_all.txt','r');
  %fprintf(fiw,'%8d %8d %8d %8d\n',np,nnp,npl,numnd);
kcm=fscanf(fir,'%d %d %d %d',[4,1]);
nprnd=kcm(1);  %真实随机数个数
nnpoornd=kcm(2); %凑10的倍数后的个数
nnpllrnd=kcm(3); %凑10的行数
numrnd=kcm(4);   %共几组随机数字
vsrd=zeros(nprnd,numrnd);
for ikk=1:numrnd
    disp(['ikk=',num2str(ikk)]);
    %fprintf(fiw,'%5d %13.5e %13.5e %13.5e %13.5e %5d\n',ikk,radmp(ikk,np+1),radmp(ikk,np+2),radmp(ikk,np+3),radmp(ikk,np+4),radmp(ikk,np+5));
    kkcm=fscanf(fir,'%d %f %f %f %f %d',[6,1]);
    ikk1=kkcm(1);  %第几列随机数字
    if ikk1 ~= ikk
        disp('错误 ikk1 ~= ikk');
        return;
    end
    for iil=1:nnpllrnd
        na=(iil-1)*10;
        kcm0=fscanf(fir,'%f %f %f %f %f %f %f %f %f %f',[10,1]);  %10个1行
        for i=1:10
            if na+i > nprnd
                break
            end
            vsrd(na+i,ikk)=kcm0(i);  %记录
            if kcm0(i) <=0
                disp('错误 kcm0(i) <=0');
                return;
            end
        end
    end %iil
end %ikk
fclose(fir);

%读入数据
fid=fopen('shuju_all.txt','r');
tcm=fscanf(fid,'%d %d',[2,1]);
np=tcm(1);  %节点数
ne=tcm(2);  %单元数
tcm1=fscanf(fid,'%d %f %f %f',[4,np]);
tcm1=tcm1';  %节点坐标
tcm2=fscanf(fid,'%d %d',[2,ne]);
tcm2=tcm2';  %单元节点号
tcmpp=fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11,ne]);
tcmpp=tcmpp';  %单元 确定性参数 ia 1,ib 2,e 3,u 4,kexi 5,area 6,ix 7,iy 8,skx 9,sky 10, siipp 11
tcm3=fscanf(fid,'%d %f %f %f %f %f %f',[7,12]);
tcm3=tcm3';  %质点参数
skkele=zeros(12,12,ne);
ngg=np*6;  %总自由度
for i=1:np
    i1=tcm1(i,1);
    i2=tcm3(i,1);
    if i1~=i
        disp('错误1');
    end 
    if i2~=i
        disp('错误2');
    end     
end
fclose(fid);

%-- 计算一下各个单元的长度
tcml=zeros(ne,1);
for i=1:ne
    i1=tcm2(i,1);
    i2=tcm2(i,2);
    tcml(i)=abs(tcm1(i2,4)-tcm1(i1,4));
end

iclc=1;   %---------------------------- 大循环 大循环
nprnd5=nprnd+5;
allrsp=zeros(1500,nprnd5);   %5种阻尼
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

 
    %-- 参数
    ngg=np*6;  %总自由度
    skk=zeros(ngg,ngg);
    smm=zeros(ngg,ngg);
    
    %-- 随机数字 vsrd=zeros(nprnd,numrnd);
   %剪切波速 1  vs_rd
   %结构材料弹性模量  2  e_rd
   %结构材料阻尼 3 damp_rd
   %质点的质量(或密度) 4-15  12个点 转动惯量(同样) mass_rd
   %梁单元的剪切面积 11个梁 16-26 sa_rd
   %梁单元的惯性矩 11个梁 27-37 各方向一致 ii_rd   
    mass_rd=zeros(12,1);
    sa_rd=zeros(11,1);
    ii_rd=zeros(11,1);
   if iclc<=nprnd   %随机因子
      vs_rd=vsrd(iclc,1);
      e_rd=vsrd(iclc,2);
      damp_rd=vsrd(iclc,3);
      for i=1:12
         mass_rd(i)=vsrd(iclc,3+i);  %相当于密度
      end
      for i=1:11
         sa_rd(i)=vsrd(iclc,15+i);
         ii_rd(i)=vsrd(iclc,26+i);
      end
   else
      %确定+1, 地基模量1/1.5 +2, 1.5 +3, 1/2 +4, 2 +5
      vs_rd=1;
      e_rd=1;
      damp_rd=1;
      for i=1:12
         mass_rd(i)=1;  %相当于密度
      end
      for i=1:11
         sa_rd(i)=1;
         ii_rd(i)=1;
      end

      if iclc==nprnd+2  %地基剪切模量降低1.5倍  1/1.5
          vs_rd=1/sqrt(1.5);
      end
      if iclc==nprnd+3  %地基剪切模量增大1.5倍
          vs_rd=sqrt(1.5);
      end
      if iclc==nprnd+4  %地基剪切模量降低2倍   1/2
          vs_rd=sqrt(0.5);
      end
      if iclc==nprnd+5  %地基剪切模量增大2倍
          vs_rd=sqrt(2.0);
      end       
   end
    
    %-- 单元ie质量阵  tcm3(12,7) 集成总质量阵
    for ii=1:np
        i1=6*(tcm3(ii,1)-1); %起点自由度
        for i=1:6
            smm(i1+i,i1+i)=smm(i1+i,i1+i)+tcm3(ii,i+1)*mass_rd(ii);
        end  %i   
    end %ii

    %-- 单元ie刚度阵 skkele(12,12,ne)
    for ie=1:ne
        
%         if ie==ne
%             ie=ne;
%         end
        
        i1=6*(tcm2(ie,1)-1);
        i2=6*(tcm2(ie,2)-1);
        
        %单元 确定性参数 ia 1,  ib 2,  e 3,  u 4,  kexi 5,  area 6,  ix 7,  iy 8,  skx 9,  sky 10 
% faix=12*E*IX*SkAY/G/A/l/l   对uy影响
% faiy=12*E*IY*SkAX/G/A/l/l   对ux影响
% E=G*2*(1+u)
% 1-1 1-5 1-7 1-11  7-1 7-5 7-7 7-11 E.AA*IY/(1+faiy)
% 2-2 2-4 2-8 2-10  8-2 8-4 8-8 8-10 E.AA*IX/(1+faix)
% 3-3 3-9 9-3 9-9   E.AA
% 4-2 4-8 10-2 10-8 E.AA*IX/(1+faix)  4-4 10-10 E.AA.IX*(4+faix)/(1+faix)
% 4-10 10-4  E.AA.IX*(2-faix)/(1+faix)
% 5-1 5-7 11-1 11-7 E.AA*IY/(1+faiy)  5-5 11-11 E.AA.IY*(4+faiy)/(1+faiy)
% 5-11 11-5  E.AA.IY*(2-faiy)/(1+faiy)
% 6-6 6-12 12-6 12-12 E.IXX
        esl=tcml(ie);        %单元长度
        ee=tcmpp(ie,3)*e_rd;        %弹性模量
        emu=tcmpp(ie,4);        %泊松比
        ekexi=tcmpp(ie,5)*damp_rd;     %阻尼比
        eara=tcmpp(ie,6);     %面积 
        eix=tcmpp(ie,7)*ii_rd(ie);       %ix
        eiy=tcmpp(ie,8)*ii_rd(ie);       %iy
        eskx=tcmpp(ie,9)/sa_rd(ie);      %skx
        esky=tcmpp(ie,10)/sa_rd(ie);     %skx
        eiip=tcmpp(ie,11)*ii_rd(ie);         %ip 饶竖直z轴的转动惯量
        
        if ie==1
            ekeximin=ekexi;
        else 
            if ekeximin>ekexi
                ekeximin=ekexi;
            end
        end
                
        faiy=24*(1+emu)*eiy*eskx/eara/esl/esl;  %对ux贡献
        faix=24*(1+emu)*eix*esky/eara/esl/esl; 
        
        scle=zeros(12,12);  %用于存储单元刚度阵
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


    
    %-- 计算固定地基的一阶自振频率
    %约束处理 剔除1-6行和列 1号点在最前面
    nggn=6*(np-1);
    skkn=zeros(nggn,nggn);
    smmn=zeros(nggn,nggn);
    for i=1:nggn
    for k=1:nggn
        skkn(i,k)=skk(6+i,6+k);
        smmn(i,k)=smm(6+i,6+k);
    end
    end

    %计算广义特征值+特征向量
    % [V,D]=eig(A,B)：
    %  AV=BVD 其中V是特征向量 (B-1/D*A)V=0   V是列向量{},{},{}... 已经验证
    %       (skkn-w2*smmn)V=0
    [vv1,dd1]=eig(smmn,skkn);
    fhz=zeros(nggn,1);
    for i=1:nggn
      fhz(i)=sqrt(1/dd1(i,i))/2/3.1415926;  %Hz 只是顺序不对 固定地基
    end

    %做个特征值和特征向量的排序(只换列即可!) 行与质点自由度相对应
    %相同频率值 按顶部节点的x\y\z振型值顺序排列 按fhz()中的hz由小到大排列
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

    %验证正交性 固定地基 1阶与10阶
    cfh1=zeros(nggn,1);
    cfh2=zeros(nggn,1);
    for i=1:nggn
      cfh1(i)=vv1(i,1);
      cfh2(i)=vv1(i,10);
    end
    aaxm=cfh1'*skkn*cfh2;     %验证用
    aaxm00=abs(aaxm);
    f1=fhz(1)                 %1阶频率 hz
    
    jxx=0;
    jzz=0;  %绕竖轴的扭转惯性矩
    for i=1:np
      ia=(i-1)*6;
      jxx=jxx+smm(ia+4,ia+4);
      jxx=jxx+smm(ia+1,ia+1)*(tcm1(i,4)-tcm1(1,4))^2;
      jzz=jzz+smm(ia+6,ia+6);
    end
    jyy=jxx;

   %=== 施加弹簧阻尼器地基模型 ===

   %a.1  地基材料的固定参数 ****
   prxy=0.31;
   dens=2650;
   ggdd=dens*vs_rd*vs_ave*vs_rd*vs_ave;  %由随机剪切波速计算剪切模量
   rr0=19.8;

   %a.2 地基材料导致的弹簧与阻尼器参数
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
   
   %a.3 叠加到模型的刚度阵中
   nggn=6*np;
   skka=skk;
   smma=smm;
   skka(1,1)=skka(1,1)+kkx;   %x
   skka(2,2)=skka(2,2)+kky;   %y
   skka(3,3)=skka(3,3)+kkz;   %z
   skka(4,4)=skka(4,4)+kkfx;  %摇摆fx
   skka(5,5)=skka(5,5)+kkfy;  %摇摆fy
   skka(6,6)=skka(6,6)+kkr;   %扭转rz

   [vv2,dd2]=eig(smma,skka);
   fhza=zeros(nggn,1);
   for i=1:nggn
      fhza(i)=sqrt(1/dd2(i,i))/2/3.1415926;  %弹性地基 Hz 只是顺序不对 固定地基
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
   
    %验证正交性 弹性地基 1阶与10阶
    cfha=zeros(nggn,1);
    cfhb=zeros(nggn,1);
    for i=1:nggn
      cfha(i)=vv2(i,1);
      cfhb(i)=vv2(i,10);
    end
    aaxm=cfha'*skka*cfhb;     %验证用 ****
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
      
   %a.4 取1和20阶的频率计算阻尼阵[C]
   %  α=2*ω1*ω2*ξ/(ω1+ω2) 
   %  β=2*ξ/(ω1+ω2)
   ww1=2*3.1415926*abs(fhza(1));
   ww2=2*3.1415926*abs(fhza(20));
   alfa=2*ww1*ww2*ekeximin/(ww1+ww2);
   beta=2*ekeximin/(ww1+ww2);           %应该按真实的材料阻尼 来核算该单元的阻尼阵 ????

   scca=alfa*smma+beta*skk;   %阻尼阵 (刚度阵部分不要考虑地基无限域部分了)
   scca(1,1)=scca(1,1)+ckx;   %x
   scca(2,2)=scca(2,2)+cky;   %y
   scca(3,3)=scca(3,3)+ckz;   %z
   scca(4,4)=scca(4,4)+ckfx;  %摇摆fx
   scca(5,5)=scca(5,5)+ckfy;  %摇摆fy
   scca(6,6)=scca(6,6)+ckr;   %扭转rz

   %a.5 Newmark系数
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

   %a.6 完整的动力计算 动力时程分析
   sktp=skka+a0k*smma+a1k*scca; 
   skin=inv(sktp);

   xa=zeros(nwav,1);   %时间
   ya2=zeros(nwav,1);    %加速度响应 2号点
   ya8=zeros(nwav,1);    %加速度响应 8号点
   ya12=zeros(nwav,1);   %加速度响应 12号点
   ut0=zeros(nggn,1);
   ut1=zeros(nggn,1);
   vt0=zeros(nggn,1);
   vt1=zeros(nggn,1);
   at0=zeros(nggn,1);
   at1=zeros(nggn,1);

   iddk=zeros(nggn,1);
   for i=1:np
      iddk(6*i-6+1)=1;   %指定x向
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
    
      ya2(istp)=at1(6*2-6+1)+xxaagg;      %2号点的x向 绝对加速度响应   ****
      ya8(istp)=at1(6*8-6+1)+xxaagg;      %8号点的x向 绝对加速度响应   ****
      ya12(istp)=at1(6*12-6+1)+xxaagg;    %12号点的x向 绝对加速度响应  ****
    
      ut0=ut1;
      vt0=vt1;
      at0=at1;    
   end

   %a.7 计算ya2、ya8、ya12节点楼层谱

   for iexi=1:5
       exi=eexi(iexi);
       ssaa2=calresp(nwav,deltk,exi,ya2);  %返回ssaa(nih,3) 反应谱 1频率 2绝对加速度反应谱(m/s2)  ssaa2(1,3)为nih值
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

   fibb=fopen('rspp_all3s.txt','wt');   %结构2、8、12点的绝对加速度响应
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
           %300个元素
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
  







