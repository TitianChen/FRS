%{
%���в���&��Ӧ�� ���ϵ������
%vvs:���в���;1*np
%a:Ƶ��;
%z:��Ӧ��;
zz_max = [];
for i = 1:np
    z_vvs = z(:,i);
    zz_max = [zz_max max(z_vvs)];
end
plot(vvs,zz_max);
xlabel('�ػ����в���m/s');
ylabel('���Ӧ���ٶ�');
vvs_coef = corrcoef(vvs,zz_max)
%}

%{
%�ڵ�����&��Ӧ�� ���ϵ������
%b:�ڵ�����;1*np
%a:Ƶ��;
%z:��Ӧ��;
np = 500;
zz_max = [];
for i = 1:np
    z_b = z(:,i);
    zz_max = [zz_max max(z_b)];
end
zz_max = zz_max(1:np);
plot(b,zz_max);
xlabel('�ڵ�����');
ylabel('���Ӧ���ٶ�');
m_coef = corrcoef(b,zz_max);
%}

%{
%�����&��Ӧ�� ���ϵ������
%b:�����;1*np
%a:Ƶ��;
%z:��Ӧ��;
zz_max = [];
for i = 1:np
    z_b = z(:,i);
    zz_max = [zz_max max(z_b)];
end
zz_max = zz_max(1:np);
plot(b,zz_max);
xlabel('�����');
ylabel('���Ӧ���ٶ�');
m_coef = corrcoef(b,zz_max)
%}


%����ģ��&��Ӧ�� ���ϵ������
%bb:����ģ��;1*np
%a:Ƶ��;
%z:��Ӧ��;
zz_max = [];
for i = 1:np
    z_b = z(:,i);
    zz_max = [zz_max max(z_b)];
end
bb = 4e10*EEII;
zz_max = zz_max(1:np);
plot(bb,zz_max);
xlabel('����ģ��');
ylabel('���Ӧ���ٶ�');
m_coef = corrcoef(bb,zz_max)

