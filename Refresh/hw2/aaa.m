m=3;
n=2;
c=[4;3];
z=c*-1;
b=[1;3;5];
A=[1 -1; 2 -1; 0 1];
A=[A eye(m)];
c=[c;zeros(m,1)];

bas=n+1:m+n;
nbas=1:n;
B=A(:,bas);
N=A(:,nbas);
cb=c(nbas);
cx=c(bas);

bas1=n+1:m+n;
nbas1=1:n;
B1=A(:,bas1);
N1=A(:,nbas1);
cb1=c(nbas1);
cx1=c(bas1);


[L,U,P] = lu(A)
%General Simplex
while any(cb(:)>0)
    
    %Step2
    [tj,j]=max(cb);
    ji=j;
    j=nbas(j);
    %step3
    ej=zeros(size(cb));
    ej(ji)=1;
    
    dx=inv(B)*N*ej;
    %step4
    
    
    xb =inv(B) * b;
    [t,tx]=(max(dx./xb));
    t=inv(t);
    %step5 dont think this is correct
    i=bas(tx);
    %step 6
    ei=zeros(size(cx));
    ei(i==bas)=1;
    dz=-1 *transpose(inv(B)*N)*ei;
    %step7
    s=-tj/dz(ji);
    %step8
    xb=xb-t*dx;
    cb=-cb-s*dz;
    %step 9
    tempB = B(:,i==bas);
    B(:,bas==i) = N(:,ji);
    N(:,nbas==j) = tempB;
    bas(bas==i) = j;
    nbas(nbas==ji) = i;
    cb(cb==0)=s;
    cb=-cb
    xb(xb==0)=t
    
end



cx=c(bas);
f=transpose(cx)*xb;
disp(xb);
disp(cb);
disp(f);