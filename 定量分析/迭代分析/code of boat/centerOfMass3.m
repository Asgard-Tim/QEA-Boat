function ANS=centerOfMass3(XMass,YMass,ZMass,MMass,m,num)
x=0;
y=0;
z=0;
    for i=1:num
        x=x+XMass(1,i)*MMass(1,i)/m;
        y=y+YMass(1,i)*MMass(1,i)/m;
        z=z+ZMass(1,i)*MMass(1,i)/m;
    end
    ANS=[x,y,z];
end

