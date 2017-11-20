m=3;
n=2; 
c=[10;-7];
z=c*-1;
b=[1;2;3];
%from previous iteration
xb=[2;1;4];
A=[1 -1 0 1 0; 2 -1 0 0 1; 0 1 1 0 0];
N=[1 -0;0 1;0 0];
c=[c;zeros(m,1)];
bas=n+1:m+n; 
nbas=1:n;
B=A(:,[1,2,3]);


while any(z(:)<0)
    %Step2
    [tj,j]=min(z);
    %step3
    dx=inv(B)*N(:,j);
    %step4
    [t,tx]=max(dx./xb);
    %step5
    i=bas(tx);
    %step6
    dz=-transpose(inv(B)*N)*[0;0;1];
    %step7
    s=z(j)/dz(j);
    %step8
    xb=xb-(inv(t).*dx);
    z=z-s.*dz;
    %Xb and Zn
    xb(xb==0)=inv(t);
    z(z==0)=s;  
    disp(xb);
    disp(z);
    %Missing a way to update the basis,A,N etc
end

for y = size(z)
    disp(z(y) +"x" + y )
end

for k = size(b)
    disp(xb(k) +"x" + k )
end
