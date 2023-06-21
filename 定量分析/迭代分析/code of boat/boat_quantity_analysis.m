%Step1:船体数学模型建立
%导入船体点集数据
load('X.mat');
load('Y.mat');
load('Z.mat');
%定义船体参数 单位：cm g
%以下参数都基于对于船体实物的测量给出
boat.L = 74; %长
boat.W = 16.5; %宽
boat.HB = boat.W / 2; %半宽
boat.D = 7; %深
density_water = 1;    % g / cm^3
boat.mass = 653.8; %船体外壳总质量
%x为长，y为宽，z为深，做切片（微分）；切片数量为num
num=10;
dy = boat.W/num; 
dz = boat.D/num; 
dx = boat.L/num;
%构建船体外壳曲面模型
mesh.xs = 0:dx:boat.L;
mesh.ys = -boat.HB:dy:boat.HB; 
mesh.zs = 0:dz:boat.D;
[mesh.ygrid,mesh.zgrid] = meshgrid(mesh.ys,mesh.zs);
[Xi,Yi]=meshgrid(mesh.xs,mesh.ys);
Zi=griddata(X,Y,Z, Xi,Yi);%griddata:对二维或三维散点数据插值——增加有效数据量，减少误差
surf(Xi,Yi,Zi);
shading flat;
axis('equal');

%Step2:计算船体重心的三维坐标:
COM_x=0;%x轴上的重心
COM_y=0;%y轴上的重心
COM_z=0;%z轴上的重心
tnt=0;%有效点个数
i=find(isnan(Zi));
Zi(i)=0;%去掉NaN点的影响
for i=1:num
    for j=1:num
            COM_x=COM_x+Xi(i,j);
            COM_y=COM_y+Yi(i,j);        
            COM_z=COM_z+Zi(i,j);
            tnt=tnt+1;
    end
end
%加入除船壳外的结构
%可乐1
m1=356;
boat.finalmass=boat.mass+m1;
x1=34;
y1=0;
z1=9.8;
%可乐2
m2=356;
boat.finalmass=boat.finalmass+m2;
x2=44;
y2=0;
z2=9.5;
%桅杆
m3=116;
boat.finalmass=boat.finalmass+m3;
x3=39;
y3=0;
z3=49.5;
%上层外壳
m4=116.1;
boat.finalmass=boat.finalmass+m4;
x4=0;
y4=39;
z4=12.5;
%尾部甲板
m5=49.4;
boat.finalmass=boat.finalmass+m5;
x5=64.5;
y5=0;
z5=6.2;
%电池
m6=185.1;
boat.finalmass=boat.finalmass+m6;
x6=12;
y6=0;
z6=5.5;
%电机1
m7=190.4;
boat.finalmass=boat.finalmass+m7;
x7=53;
y7=4;
z7=3;
%电机2
m8=191.6;
boat.finalmass=boat.finalmass+m8;
x8=53;
y8=-4;
z8=3;
%电调1
m9=49.4;
boat.finalmass=boat.finalmass+m9;
x9=20.5;
y9=6;
z9=7;
%电调2
m10=49.6;
boat.finalmass=boat.finalmass+m10;
x10=20.5;
y10=-6;
z10=7;
%船桨1
m11=43.5;
boat.finalmass=boat.finalmass+m11;
x11=73.5;
y11=3.8;
z11=1;
%船桨2
m12=40.3;
boat.finalmass=boat.finalmass+m12;
x12=73.5;
y12=-3.8;
z12=1;
%前端甲板
m13=18.2;
boat.finalmass=boat.finalmass+m13;
x13=20.5;
y13=0;
z13=6.3;
%最终的重心三维坐标
density=boat.mass/tnt;
COM_x=(COM_x*density+m1*x1+m2*x2+m3*x3+m4*x4+m5*x5+m6*x6+m7*x7+m8*x8+m9*x9+m10*x10+m11*x11+m12*x12+m13*x13)/boat.finalmass;
COM_y=(COM_y*density+m1*y1+m2*y2+m3*y3+m4*y4+m5*y5+m6*y6+m7*y7+m8*y8+m9*y9+m10*y10+m11*y11+m12*y12+m13*y13)/boat.finalmass;
COM_z=(COM_z*density+m1*z1+m2*z2+m3*z3+m4*z4+m5*z5+m6*z6+m7*z7+m8*z8+m9*z9+m10*z10+m11*z11+m12*z12+m13*z13)/boat.finalmass;

%Step3:计算船体浮心三维坐标
%初始化：COB_x,COB_y,COB_z记录每个切片的浮心（2D）
COB_x=[];
COB_y=[];
COB_z=[];
MMass=[];%每个浮心点的权重（质量）
tot_mass=0;%总质量
Y=[];
%遍历1-180度的所有theta值
for theta=1:1:180
    for i=1 :num   %对每个x积分
        ZZ=[];%记录该X下的Z轴坐标
        for j=1 :num+1 %讨论y从-boat.HB到boat.HB
            t=calculate(i*dx,mesh.ys(1,j));
            ZZ=[ZZ t];
        end
        boat.hull = mesh.zgrid > ZZ; %在所计算的Z的值之上的就是船的截面
        %不同theta下的水线表示
        %确定重力浮力平衡以迭代调整水线
        d =water_line(mesh,theta,boat,dx,dy,dz,num);%寻找函数零点
        y = mesh.ys;
        z = -tand(theta).*y+d;
        if theta>90 && theta<=180
            water = (mesh.zgrid-boat.D/2) < z;
        elseif theta<=90 && theta>=0
            water = (mesh.zgrid-boat.D/2) > z;
        end
        %定义淹没区域
        sub_region = flipud(boat.hull)&water;%矩阵交集,求取排水体积;&位与运算（都是1才得1）
        %计算船体浮心
        COB = centerOfMass2(sub_region,mesh);
        COB_z=[COB_z COB(1,2)];
        COB_y=[COB_y COB(1,1)];
        COB_x=[COB_x i*dx];
        MMass=[MMass matrixSum(sub_region)];
        tot_mass=tot_mass+matrixSum(sub_region);
    end
    %得出三维浮心坐标
    ANS=centerOfMass3(COB_x,COB_y,COB_z,MMass,tot_mass,num);%计算三维浮心
    COB_x=ANS(1,1);
    COB_y=ANS(1,2);
    COB_z=ANS(1,3);

    %Step4:求恢复力臂
    k1 = -tand(theta);
    if ((k1*COB_z+COB_y-COM_y)/k1)-COM_z<0
        flag1=1;
    else
        flag1=-1;
    end
    l = flag1*abs(COB_y+k1*COB_z-k1*COM_z-COM_y)/(k1^2+1)^0.5;%恢复力臂
    Y=[Y,l-4];
end    

%Step5:绘制扶正力臂曲线
x=180:-1:1;
x1=1:1:180;
y1=interp1(x,Y,x1,'spline');
plot(x1,y1);