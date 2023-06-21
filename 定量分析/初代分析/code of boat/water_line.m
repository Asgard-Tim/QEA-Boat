function D_line = water_line(mesh,theta,boat,dx,dy,dz,num)   
     number = 1;
     i=1;
     d = -0.02;
     d_gap = 0.07;
     density_water = 1000; 
while number<100
     y = mesh.ys;
     z = -tand(theta).*y+d;
     if theta>90 && theta<=180
        water = (mesh.zgrid-boat.D/2) < z;
     elseif theta<=90 && theta>=0
        water = (mesh.zgrid-boat.D/2) > z;
     end                      %水线下的水体矩阵

     sub_region = flipud(boat.hull)&water;   %矩阵交集，求取排水体积

    force=sub_region.*density_water*dz*dy*dx;    %计算单个微分方块的力
   
    lift= matrixSum(force);                     %整体浮力
    up_force=lift-(boat.finalmass/num);                %浮力剩余值

    if up_force<0                           
        d=d+d_gap/i;                        
    elseif up_force>0
        d=d-d_gap/i;
    else 
        break;
    end

    if up_force>-0.01&&up_force<0.01              %判断阈值，水线收敛
       break;
    end
    i=i*2;                                      %缩短步长
    number=number+1;
end  
    D_line = d;
    
end
