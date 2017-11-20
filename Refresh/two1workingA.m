 
c=[6;8;5;9];
z=c*-1;
b=[5;3];
A=[2 1 1 3;1 3 1 2];
[m,mn]=size(b);
[nm,n]=size(A);
A=[A eye(m)]; 
c=[c;zeros(m,1)];
    
bas=n+1:m+n;
nbas=1:n;
B=A(:,bas);
N=A(:,nbas);
cb=c(nbas);
cx=c(bas);


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
    if t<= 0
       break; 
    end
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
    cb=-cb;
    xb(xb==0)=t  ;
    
end
cx=c(bas);
f=transpose(cx)*xb;
disp("Xb");
disp(xb);
disp("Cb");
disp(cb);
disp("f=" +f);