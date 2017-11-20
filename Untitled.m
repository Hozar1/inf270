m=3;
n=2; 
c=[4;3];
z=c*-1;
b=[1;3;5];
A=[1 -1; 2 -1; 0 1];
N=A;
A=[A eye(m)]; 
c=[c;zeros(m,1)];
bas=n+1:m+n; 
nbas=1:n;


while any(z(:)<0)
    %Step2
    j=find(min(z));
    %step3
    dx=N(:,j);
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
    b(b==0)=t;
    z(z==0)=s;  
end