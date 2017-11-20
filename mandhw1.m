m=3;
n=2; 
c=[4;3];
z=c*-1;
%Initial case xb = b
b=[1;3;5];
A=[1 -1; 2 -1; 0 1];
N=A;
A=[A eye(m)]; 
c=[c;zeros(m,1)];
bas=n+1:m+n; 
nbas=1:n;
B=A(:,bas);

while any(z(:)<0)
    %Step2
    [tj,j]=min(z);
    %step3
    dx=inv(B)*N(:,j);
    %step4
    [t,tx]=max(dx./b);
    %step5
    i=bas(tx);
    %step6
    dz=-1*transpose(N(j,:));
    %step7
    s=z(j)/dz(j);
    %step8
    b=b-(t.*dx);
    z=z-s.*dz;
    %Xb and Zn
    b(b==0)=t;
    z(z==0)=s;  
    disp(b);
    disp(z);
end