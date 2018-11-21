%楼层反应谱的后处理;
fibb=fopen('rspp_all3sx.txt','r');   %结构2、8、12点的随机样本 绝对加速度响应
kcm=fscanf(fibb,'%d %d',[2,1]);
nih=kcm(1);    %频率点
nprnd=kcm(2);  %样本总数
ssaa=zeros(nih,4);      %1频率
osp=zeros(10,1);        %10个元素
spall=zeros(450,1);     %1个样本
sta2=zeros(nih,nprnd+15);
sta8=zeros(nih,nprnd+15);
sta12=zeros(nih,nprnd+15);
spssc=zeros(nprnd,1);
for i=1:nih
    fii=fscanf(fibb,'%f',[1,1]);
    ssaa(i,1)=fii;
end

for iexi=1:5  %5种阻尼比
    
    ekexi=fscanf(fibb,'%f',[1,1]);   %阻尼比 
    
    ia=0;
    for iclc=1:nprnd+5
        
        ia=ia+1;
        if(ia==500)
            ia=0;
            disp(['iclc=',num2str(iclc)]);
        end
    
        %-- 450个元素
        for il=1:45
            na=(il-1)*10;
            osp=fscanf(fibb,'%f %f %f %f %f %f %f %f %f %f',[10,1]);
            for i=1:10
                spall(na+i)=osp(i);
            end
        end
        %-- 分配到几个 2 8 12
        for i=1:nih
            sta2(i,iclc)=spall(i);
            sta8(i,iclc)=spall(150+i);
            sta12(i,iclc)=spall(300+i);
        end
    end %iclc
    
    % +1确定; +2 1/1.5 +3 1.5;  +4  1/2  +5  2;
    %计算3个楼层的 随机楼层谱的均值 上下限制
    xa=zeros(nih,1);
    yaq=zeros(nih,1);
    yave=zeros(nih,1);
    yamax=zeros(nih,1);
    yamin=zeros(nih,1);
    ya1=zeros(nih,1);   %90%
    ya2=zeros(nih,1);   %90%数数
    ya3=zeros(nih,1);   %1.5倍 包络
    ya4=zeros(nih,1);   %2倍   包络
    ya5=zeros(nih,1);   %1/1.5
    ya6=zeros(nih,1);   %1.5
    ya7=zeros(nih,1);   %1/2
    ya8=zeros(nih,1);   %2 
    ya9=zeros(nih,1);    %95%
    ya10=zeros(nih,1);   %95%数数
    ifl=3;
    for ifl=1:3
    
       disp(['ifl=',num2str(ifl)]);
    
       %-- 读入
       if ifl==1
          stss=sta2;   %2点
       else
          if ifl==2
            stss=sta8;   %8点
          else
            if ifl==3
                stss=sta12;    %12点
            end
          end
       end
    
    
       for ih=1:nih
          %-- 平均值 ＋6
          ave=0;
          for i=1:nprnd
            ave=ave+stss(ih,i);
          end
          ave=ave/nprnd;            %平均值
          stss(ih,nprnd+6)=ave;
                
          %-- 标准差 ＋7
          asum=0;
          for i=1:nprnd
            asum=asum+(stss(ih,i)-ave)*(stss(ih,i)-ave);
          end
          asegm=sqrt(asum/nprnd);    %样本空间总体  [标准差]
          stss(ih,nprnd+7)=asegm; 
       
          %-- 上下限制 +8下  +9上
          amin=stss(ih,1);
          amax=stss(ih,1);
          for i=2:nprnd
            if amin>stss(ih,i)
                amin=stss(ih,i);
            end
            if amax<stss(ih,i)
                amax=stss(ih,i);
            end
          end
          stss(ih,nprnd+8)=amin;
          stss(ih,nprnd+9)=amax;
       end %ih
   
       %-- 绘制图
       for ih=1:nih
          xa(ih)=ssaa(ih,1);
          yaq(ih)=stss(ih,nprnd+1)/9.81;     %1 确定
          yave(ih)=stss(ih,nprnd+6)/9.81;    %7 平均值
          yamin(ih)=stss(ih,nprnd+8)/9.81;   %8 下
          yamax(ih)=stss(ih,nprnd+9)/9.81;   %9 上
       end
    
       %-- 3种计算90%保证率的结果 1.27倍 90%  1.645倍 95%   2.58倍 99%   3倍 99.99%
       %-- 标准差+1.27  +10
       idkrnd=round(0.9*nprnd);            %****
       idkrnd_95=round(0.95*nprnd);
       for ih=1:nih
         stss(ih,nprnd+10)=stss(ih,nprnd+6)+1.282*stss(ih,nprnd+7);    %**** 90%保证率
         stss(ih,nprnd+11)=stss(ih,nprnd+6)-1.282*stss(ih,nprnd+7);    %**** 10%保证率
         stss(ih,nprnd+12)=stss(ih,nprnd+6)+1.645*stss(ih,nprnd+7);    %**** 95%保证率
       end
       %-- 排序
%        for ih=1:nih
%           %ih
%           for i=1:nprnd
%             spssc(i)=stss(ih,i);
%           end
%           ckcc=sort(spssc);
%           for i=2:nprnd
%             if ckcc(i)<ckcc(i-1)
%                 disp('错误 ckcc(i)<ckcc(i-1)');
%             end
%           end
%           stss(ih,nprnd+11)=ckcc(idkrnd);        %**** 90%保证率 数数...
%           stss(ih,nprnd+13)=ckcc(idkrnd_95);        %**** 90%保证率 数数... 
%        end     

       for ih=1:nih
         ya1(ih)=stss(ih,nprnd+10)/9.81;   %+10 几倍标准差  90% 
         ya2(ih)=stss(ih,nprnd+11)/9.81;   %+11 数数        10%
         ya9(ih)=stss(ih,nprnd+12)/9.81;   %+10 几倍标准差   95%
         ya10(ih)=stss(ih,nprnd+13)/9.81;   %+11 数数
       end
%     ya5=zeros(nih,1);   %1/1.5
%     ya6=zeros(nih,1);   %1.5
%     ya7=zeros(nih,1);   %1/2
%     ya8=zeros(nih,1);   %2       
       for ih=1:nih
           ya3(ih)=yaq(ih);
           ya4(ih)=yaq(ih);
           if ya3(ih)<stss(ih,nprnd+2)/9.81
               ya3(ih)=stss(ih,nprnd+2)/9.81;  %1/1.5
           end
           ya5(ih)=stss(ih,nprnd+2)/9.81;      %1/1.5
           if ya3(ih)<stss(ih,nprnd+3)/9.81
               ya3(ih)=stss(ih,nprnd+3)/9.81;  %1.5
           end  
           ya6(ih)=stss(ih,nprnd+3)/9.81;      %1.5 
           if ya4(ih)<stss(ih,nprnd+4)/9.81
               ya4(ih)=stss(ih,nprnd+4)/9.81;  %1/2
           end
           ya7(ih)=stss(ih,nprnd+4)/9.81;      %1/2 
           if ya4(ih)<stss(ih,nprnd+5)/9.81
               ya4(ih)=stss(ih,nprnd+5)/9.81;  %2
           end  
           ya8(ih)=stss(ih,nprnd+5)/9.81;      %2 
       end
    
       %plot(xa,yaq,'k',xa,yamin,'r-',xa,yamax,'r',xa,yaave,'r')   %绘制出匹配来 红色(目标)
       %-- 输出到文件
       %if iexi==1   %第一个阻尼比
       if ifl==1
               switch iexi
                   case 1
                       fioo=fopen('cmp_2nd_1.txt','wt');       %2
                   case 2
                       fioo=fopen('cmp_2nd_2.txt','wt');       %2
                   case 3
                       fioo=fopen('cmp_2nd_3.txt','wt');       %2
                   case 4
                       fioo=fopen('cmp_2nd_4.txt','wt');       %2
                   case 5
                       fioo=fopen('cmp_2nd_5.txt','wt');       %2
               end
       else
              if ifl==2
                  switch iexi
                      case 1
                          fioo=fopen('cmp_8nd_1.txt','wt');       %8
                      case 2
                          fioo=fopen('cmp_8nd_2.txt','wt');       %8
                      case 3
                          fioo=fopen('cmp_8nd_3.txt','wt');       %8
                      case 4
                          fioo=fopen('cmp_8nd_4.txt','wt');       %8
                      case 5
                          fioo=fopen('cmp_8nd_5.txt','wt');       %8
                  end
              else
                  if ifl==3
                      switch iexi
                          case 1
                              fioo=fopen('cmp_12nd_1.txt','wt');       %12
                          case 2
                              fioo=fopen('cmp_12nd_2.txt','wt');       %12
                          case 3
                              fioo=fopen('cmp_12nd_3.txt','wt');       %12
                          case 4
                              fioo=fopen('cmp_12nd_4.txt','wt');       %12
                          case 5
                              fioo=fopen('cmp_12nd_5.txt','wt');       %12
                      end
                  end
              end
       end
       %end  %iexi==1
       
       %xa(ih),yaq(ih),yave(ih),yamin(ih),yamax(ih),ya1(ih),ya2(ih),ya3(ih),ya4(ih)
       %频率   确定性   平均值   最小值     最大值    90%保证率 90%数数  1.5包络  2包络 
       %,ya5(ih),ya6(ih),ya7(ih),ya8(ih)
       %,1/1.5     1.5     1/2    2   
       fprintf(fioo,'%13.5f \n',ekexi);
       for ih=1:nih
        fprintf(fioo,'%13.5f %13.5e %13.4e %13.4e %13.4e %13.4e %13.4e %13.4e %13.4e %13.4e %13.4e |  %13.4e %13.4e   |  %13.4e %13.4e \n',xa(ih),yaq(ih),yave(ih),yamin(ih),yamax(ih),ya1(ih),ya2(ih),ya9(ih),ya10(ih),ya3(ih),ya4(ih),ya5(ih),ya6(ih),ya7(ih),ya8(ih));
       end
       fclose(fioo);
       
    
       xxkk=zeros(41,1);
       xxoo=zeros(40,1);
       if ifl==2 &&  iexi==2   %8号点 特殊分析 2%阻尼比  **** 
          for ih=1:nih
            iaaa=0;
            if abs(xa(ih)-2.1) <1e-3  %****
                for i=1:nprnd
                    spssc(i)=stss(ih,i);
                end
                yaave=stss(ih,nprnd+6);
                yagam=stss(ih,nprnd+7);
                yamin=stss(ih,nprnd+8);
                yamax=stss(ih,nprnd+9);
                fopc=fopen('8_21.txt','wt');
                fopcp=fopen('8_21_sa.txt','wt');
                iaaa=1;
            else
                if abs(xa(ih)-3.3)<1e-3   %****
                    for i=1:nprnd
                       spssc(i)=stss(ih,i);
                    end   
                    yaave=stss(ih,nprnd+6);
                    yagam=stss(ih,nprnd+7);
                    yamin=stss(ih,nprnd+8);
                    yamax=stss(ih,nprnd+9);    
                    fopc=fopen('8_33.txt','wt');
                    fopcp=fopen('8_33_sa.txt','wt');
                    iaaa=1;
                else
                  if abs(xa(ih)-4.2)<1e-3   %****
                      for i=1:nprnd
                       spssc(i)=stss(ih,i);
                      end   
                      yaave=stss(ih,nprnd+6);
                      yagam=stss(ih,nprnd+7);
                      yamin=stss(ih,nprnd+8);
                      yamax=stss(ih,nprnd+9);    
                      fopc=fopen('8_42.txt','wt');
                      fopcp=fopen('8_42_sa.txt','wt');
                      iaaa=1;
                  else
                     if abs(xa(ih)-3)<1e-3   %****
                      for i=1:nprnd
                       spssc(i)=stss(ih,i);
                      end   
                      yaave=stss(ih,nprnd+6);
                      yagam=stss(ih,nprnd+7);
                      yamin=stss(ih,nprnd+8);
                      yamax=stss(ih,nprnd+9);    
                      fopc=fopen('8_30.txt','wt');
                      fopcp=fopen('8_30_sa.txt','wt');
                      iaaa=1; 
                     end
                  end
                end
            end
            
            if iaaa==1
               detxx=(yamax-yamin)/40;    %分40块
               xxkk(1)=yamin;
               for ik=1:40
                   xxkk(ik+1)=xxkk(ik)+detxx;
                   xxoo(ik)=(xxkk(ik)+xxkk(ik+1))/2;
               end
               iiccd=zeros(40,1);
               for i=1:nprnd
                   if spssc(i) <yamin*0.999
                       disp('错误 spssc(i) <yamin*0.99');
                   end
                   if spssc(i) >yamax*1.001
                       disp('错误 spssc(i) >yamax*1.001');
                   end                  
                   ist=0;
                   for ik=1:40
                       if spssc(i)>=xxkk(ik) && spssc(i)<=xxkk(ik+1)
                           ist=ik; 
                           break;
                       end
                   end
                   if ist==0
                       %disp('错误 ist==0');
                       if spssc(i)<=xxkk(1)
                           ist=1;
                       else spssc(i)>=xxkk(41)
                           ist=40;
                       end
                   end
                   if ist==0
                       disp('错误 ist==0');
                       stop
                   end
                   iiccd(ist)=iiccd(ist)+1;    %数个数
               end
               
               %-- 输出
               fprintf(fopc,'%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %7d %13.5e\n',yaave/9.81,yagam/9.81,ya1(ih),ya2(ih),yamin/9.81,yamax/9.81,nprnd,detxx/9.81);
               for ih=1:40
                   fprintf(fopc,'%13.5e %7d \n',xxoo(ih)/9.81,iiccd(ih));
               end
               for ih=1:nprnd
                   fprintf(fopcp,'%15.7e \n',spssc(ih)/9.81);
               end
               fclose(fopc);
            end %iaaa=1
            
          end %ih
          
          fiaa=fopen('minganxing_8.txt','wt');   %结构2、8、12点的绝对加速度响应  验证了calresp.m的结果
          for i=1:nih
              fprintf(fiaa,'%13.5f %13.5e %13.5e %13.5e\n',ssaa(i),stss(i,nprnd+1)/9.81,stss(i,nprnd+10)/9.81,stss(i,nprnd+11)/9.81);
          end
          fclose(fiaa);
       end %ifl==2  8号点 特殊分析
       
       if ifl==3&&iexi==2
         ficc=fopen('minganxing_12.txt','wt');   %结构2、8、12点的绝对加速度响应  验证了calresp.m的结果
          for i=1:nih
              fprintf(ficc,'%13.5f %13.5e %13.5e %13.5e\n',ssaa(i),stss(i,nprnd+1)/9.81,stss(i,nprnd+10)/9.81,stss(i,nprnd+11)/9.81);
          end
          fclose(ficc);  
       end
       
        if ifl==3&&iexi==5
         fidd=fopen('minganxing_12_0.05.txt','wt');   %结构2、8、12点的绝对加速度响应  验证了calresp.m的结果
          for i=1:nih
              fprintf(fidd,'%13.5f %13.5e %13.5e %13.5e\n',ssaa(i),stss(i,nprnd+1)/9.81,stss(i,nprnd+10)/9.81,stss(i,nprnd+11)/9.81);
          end
          fclose(fidd);  
        end
        if ifl==2&&iexi==5
         fiee=fopen('minganxing_8_0.05.txt','wt');   %结构2、8、12点的绝对加速度响应  验证了calresp.m的结果
          for i=1:nih
              fprintf(fiee,'%13.5f %13.5e %13.5e %13.5e\n',ssaa(i),stss(i,nprnd+1)/9.81,stss(i,nprnd+10)/9.81,stss(i,nprnd+11)/9.81);
          end
          fclose(fiee);  
       end
    
    end %ifl  1:3
   
end %iexi
fclose(fibb);

% %读入确定性的楼层谱结果
% tmp=zeros(4,1);
% fiaa=fopen('rspp_queding_nod3s.txt','r');   %结构2、8、12点的绝对加速度响应  验证了calresp.m的结果
% for i=1:nih
%     tmp=fscanf(fiaa,'%f %f %f %f',[4,1]);
%     if abs(tmp(1)-ssaa(i,1))>1e-4
%        disp('错误 abs(tmp(1)-ssaa(i,1))>1e-4');
%     end
%     sta2(i,nprnd+1)=tmp(2);
%     sta8(i,nprnd+1)=tmp(3);
%     sta12(i,nprnd+1)=tmp(4);
% end
% fclose(fiaa);
 

ifl


