clear
clc
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
rr_all=zeros(6,nih);
slop=zeros(6,nih);
fopc=fopen('sensibility.txt','w');
fopcc=fopen('sensibility1.txt','w');
rank=zeros(nprnd,6);
Sa8_2=zeros(150,nprnd);
xx=zeros(nprnd,1);
yy=zeros(nprnd,1);
for i=1:nih
    fii=fscanf(fibb,'%f',[1,1]);
    ssaa(i,1)=fii;
end

for iexi=1:2  %5种阻尼比
    
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
end
%  iicc=150*3+150; %%%8号点2%阻尼
   for ii=1:nih
       for jj=1:nprnd
         Sa8_2(jj,ii)=sta8(ii,jj)/9.81;
       end
   end
        
   %敏感性   

   
  
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
   for ii=1:3    %%分析6种随机数
    rank(:,ii)=vsrd(:,ii);
   end
   rank(:,4)=vsrd(:,4);
   rank(:,5)=vsrd(:,16);
   rank(:,6)=vsrd(:,27);
   
  
 for jj=1:nih
    yy=zeros(nprnd,1);
    yy(:,1)=Sa8_2(:,jj);   %%响应结果
    yy=yy';
    ikk=ssaa(jj,1);
        
    ave_yy=0;  
    for kk=1:nprnd
       ave_yy=ave_yy+yy(1,kk);
    end
    ave_yy=ave_yy/nprnd;            %平均值
          
    asum_yy=0;
    for kk=1:nprnd
      asum_yy=asum_yy+(yy(1,kk)-ave_yy)*(yy(1,kk)-ave_yy);
    end
    asegm_yy=sqrt(asum_yy/nprnd); 
          
   for ii=1:6  %%分析6种随机数
%     disp(['随机数',num2str(ii)]);
%     fprintf(fopc,'%10d\n',ii);
      xx=zeros(nprnd,1);
      xx(:,1)=rank(:,ii);     %%随机数
      xx=xx';
      ave_xx=0;
      for kk=1:nprnd
        ave_xx=ave_xx+xx(1,kk);
      end
      ave_xx=ave_xx/nprnd;            %平均值
      
      asum_xx=0;
      for kk=1:nprnd
        asum_xx=asum_xx+(xx(1,kk)-ave_xx)*(xx(1,kk)-ave_xx);
      end
      asegm_xx=sqrt(asum_xx/nprnd); 
          
                  
           asum_cov=0;
           for kk=1:nprnd
            asum_cov=asum_cov+(xx(1,kk)-ave_xx)*(yy(1,kk)-ave_yy);
           end
           rr=asum_cov/(asegm_xx*asegm_yy*nprnd);
           
            rr_all(ii,jj)=rr;
           
%            fprintf(fopc,'%10.5f',rr);
            
          if abs(ssaa(jj,1)-3.15)<1e-3||abs(ssaa(jj,1)-2.9)<1e-3 
           A=polyfit(xx,yy,1);
           zz=polyval(A,xx);
           slop(ii,jj)=A(1,1);
%            figure;
%            plot(xx,yy,'r*',xx,zz,'b');
           plot(xx,yy,'r*','markersize',2);
            hold on;
           plot(xx,zz,'b');
          if ii==1
              xlabel(strcat('随机参数G与确定值的缩放比例'),'fontname','宋体','FontSize',7);
           else if ii==2
                  xlabel(strcat('随机参数E与确定值的缩放比例'),'fontname','宋体','FontSize',7);
               else if ii==3
                       xlabel(strcat('随机参数damp与确定值的缩放比例'),'fontname','宋体','FontSize',7);
                   else if ii==4
                           xlabel(strcat('随机参数mass与确定值的缩放比例'),'fontname','宋体','FontSize',7);
                       else if ii==5
                               xlabel(strcat('随机参数Sa与确定值的缩放比例'),'fontname','宋体','FontSize',7);
                           else 
                               xlabel(strcat('随机参数I与确定值的缩放比例'),'fontname','宋体','FontSize',7);
                           end
                       end
                   end
               end
           end
           ylabel('楼层反应谱值Sa(g)','fontname','宋体','FontSize',7);
%            xlabel('横轴x','fontname','宋体','fontsize',12,'fontweight','bold')
           set(gca,'FontName','Times New Roman','FontSize',4.6);
           set(gcf,'Units','centimeters','Position',[10 10 4.9 3.8])
%              set(1,'visible','off');
            aabb=[ii,ikk];
%             print(gcf,'-djpeg',strcat(num2str(aabb),'.fig'));
            saveas(gcf,strcat(num2str(aabb),'.fig'));
            close gcf
          end
         end
      end
  
   for jj=1:nih
   fprintf(fopc,'%10.3f ',ssaa(jj,1));
  end     
   fprintf(fopc,'\n');
  for ii=1:6
     for jj=1:nih 
        fprintf(fopc,'%10.3f ',rr_all(ii,jj));
     end
     fprintf(fopc,'\n');
  end 
   fprintf(fopc,'\n');
   for ii=1:6
     for jj=1:nih 
        fprintf(fopc,'%10.3f ',slop(ii,jj));
     end
     fprintf(fopc,'\n');
  end  
  
  
  
  
   for jj=1:nih
   fprintf(fopcc,'%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n',ssaa(jj,1),rr_all(1,jj),rr_all(2,jj),rr_all(3,jj),rr_all(4,jj),rr_all(5,jj),rr_all(6,jj));
  end     
        
  fclose(fopc);
  fclose(fopcc);

