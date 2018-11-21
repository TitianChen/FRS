%{
%剪切波速&反应谱 相关系数分析
%vvs:剪切波速;1*np
%a:频率;
%z:反应谱;
zz_max = [];
for i = 1:np
    z_vvs = z(:,i);
    zz_max = [zz_max max(z_vvs)];
end
plot(vvs,zz_max);
xlabel('地基剪切波速m/s');
ylabel('最大反应加速度');
vvs_coef = corrcoef(vvs,zz_max)
%}

%{
%节点质量&反应谱 相关系数分析
%b:节点质量;1*np
%a:频率;
%z:反应谱;
np = 500;
zz_max = [];
for i = 1:np
    z_b = z(:,i);
    zz_max = [zz_max max(z_b)];
end
zz_max = zz_max(1:np);
plot(b,zz_max);
xlabel('节点质量');
ylabel('最大反应加速度');
m_coef = corrcoef(b,zz_max);
%}

%{
%阻尼比&反应谱 相关系数分析
%b:阻尼比;1*np
%a:频率;
%z:反应谱;
zz_max = [];
for i = 1:np
    z_b = z(:,i);
    zz_max = [zz_max max(z_b)];
end
zz_max = zz_max(1:np);
plot(b,zz_max);
xlabel('阻尼比');
ylabel('最大反应加速度');
m_coef = corrcoef(b,zz_max)
%}


%弹性模量&反应谱 相关系数分析
%bb:弹性模量;1*np
%a:频率;
%z:反应谱;
zz_max = [];
for i = 1:np
    z_b = z(:,i);
    zz_max = [zz_max max(z_b)];
end
bb = 4e10*EEII;
zz_max = zz_max(1:np);
plot(bb,zz_max);
xlabel('弹性模量');
ylabel('最大反应加速度');
m_coef = corrcoef(bb,zz_max)

