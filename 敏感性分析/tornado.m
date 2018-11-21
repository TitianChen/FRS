G = 1.225;%���в���
E = 1.421;%����ģ��
damp = 1.327;%�����
mass = 0.971;%�ڵ�����
pG = 0.55;
pE = 0.6;
pdamp = 0.95;
pmass = 0.55;
x1 = [-G*pG,G*(1-pG),-E*pE,E*(1-pE),-damp*pdamp,damp*(1-pdamp),-mass*pmass,mass*(1-pmass)];
y1 = [1,1,2,2,3,3,4,4];
for i = 1:4
    xx1 = [x1(2*i-1),x1(2*i)];
    yy1 = [y1(2*i-1),y1(2*i)];
    plot(xx1,yy1,'b',xx1,yy1-0.1,'b',xx1,yy1+0.1,'b','Linewidth',6);
    hold on;
end
plot([0,0],[0,5],'k');
set(gca,'yticklabel',{'','','���в���','','����ģ��','','�����','','�ڵ�����'});