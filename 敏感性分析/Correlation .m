%计算相关性;   
fopc=fopen('sensibility.txt','w');
   rank=zeros(nprnd,6);
   Sa8_2=zeros(150,nprnd);
   xx=zeros(nprnd,1);
   yy=zeros(nprnd,1);
   
   for ii=1:3    %%分析6种随机数
    rank(:,ii)=vsrd(:,ii);
   end
   rank(:,4)=vsrd(:,4);
   rank(:,5)=vsrd(:,16);
   rank(:,6)=vsrd(:,27);
   
   iicc=150*3+150; %%%8号点2%阻尼
   for ii=1:150
       for jj=1:nprnd
         Sa8_2(jj,ii)=allrsp(iicc+ii,jj);
       end
   end
        
   %敏感性
   
       
 for jj=1:nih
          
   for ii=1:6  %%分析6种随机数
%       disp(['随机数',num2str(ii)]);
%       fprintf(fopc,'%10d\n',ii);
      xx=zeros(nprnd,1);
      xx(:,1)=rank(:,ii);
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
          
          
        yy=zeros(nprnd,1);
        yy(:,1)=Sa8_2(:,jj);
        yy=yy';
        if abs(saatmpx(jj,1)-2.1)<1e-3 
           ikk=2.1;
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
           
           asum_cov=0;
           for kk=1:nprnd
            asum_cov=asum_cov+(xx(1,kk)-ave_xx)*(yy(1,kk)-ave_yy);
           end
           rr=asum_cov/(asegm_xx*asegm_yy*nprnd);
           
           disp(['fih=',num2str(ikk)]);
           disp(['相关系数',num2str(rr)]);
           fprintf(fopc,'%10.1f %10.5f\n',ikk,rr);
                      
%            A=polyfit(xx,yy,1);
%            zz=polyval(A,xx);
%            plot(xx,yy,'r*',xx,zz,'b');
%            hold on;
%            disp(['线性系数',num2str(A)]);
          
        else
           if abs(saatmpx(jj,1)-3.3)<1e-3 
               ikk=3.3;
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
           
               asum_cov=0;
               for kk=1:nprnd
                asum_cov=asum_cov+(xx(1,kk)-ave_xx)*(yy(1,kk)-ave_yy);
               end
               rr=asum_cov/(asegm_xx*asegm_yy*nprnd);
            
              disp(['fih=',num2str(ikk)]);
              disp(['相关系数',num2str(rr)]);
              fprintf(fopc,'%10.1f %10.5f\n',ikk,rr);
        
              
%            A=polyfit(xx,yy,1);
%            zz=polyval(A,xx);
%            plot(xx,yy,'r*',xx,zz,'b');
%            hold on;
%            disp(['线性系数',num2str(A)]);
           else
               if abs(saatmpx(jj,1)-4.2)<1e-3 
                   ikk=4.2;
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
           
                   asum_cov=0;
                   for kk=1:nprnd
                     asum_cov=asum_cov+(xx(1,kk)-ave_xx)*(yy(1,kk)-ave_yy);
                   end
                   rr=asum_cov/(asegm_xx*asegm_yy*nprnd);
                   
                   disp(['fih=',num2str(ikk)]);
                   disp(['相关系数',num2str(rr)]);
                   fprintf(fopc,'%10.1f %10.5f\n',ikk,rr);
                                     
%                   A=polyfit(xx,yy,1);
%                   zz=polyval(A,xx);
%                   plot(xx,yy,'r*',xx,zz,'b');
%                   hold on;
%                   disp(['线性系数',num2str(A)]);
               end
           end
        end
   end
  end
   
      
   
   
   
   


