



c=[-2;-1];
z=c*-1;
b=[-1;-2;1];
A=[-1 1;-1 -2;0 1];
Aux=A;
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
%Auxilliary
Aux(:,n+1)=1
Aux=[Aux eye(m)];
cux(m,1)=0;
cux(m+1,1)=-1;
cux=[cux;zeros(m+1,1)];
bux=b;
basux=n+2:m+n+1;
nbasux=1:n+1;
Bux=Aux(:,basux);
Nux=Aux(:,nbasux);
cbux= cux(nbasux);
cxux=cux(basux);
while any(cux(:)<0)
    
    %Step2
    [tj,j]=max(cbux);
    ji=j;
    j=nbasux(j);
    %step3
    ej=zeros(size(cbux));
    ej(ji)=1;
    dx=inv(Bux)*Nux*ej;
    %step4
    xb =inv(Bux) * bux;
    [t,tx]=(max(dx./xb));
    t=inv(t);
    if t<= 0
        disp("Unbounded");
        break;
    end
    if isinf(t)
        disp("Infeasible");
        break;
    end
    %step5 dont think this is correct
    i=basux(tx);
    %step 6
    ei=zeros(size(cxux));
    ei(i==basux)=1;
    dz=-1 *transpose(inv(Bux)*Nux)*ei;
    %step7
    s=-tj/dz(ji);
    %step8
    xb=xb-t*dx;
    cbux=-cbux-s*dz;
    %step 9
    tempB = Bux(:,i==basux);
    Bux(:,basux==i) = Nux(:,ji);
    Nux(:,nbasux==j) = tempB;
    
    basux(basux==i) = j;
    nbasux(nbasux==ji) = i;
    cbux(cbux==0)=s;
    cbux=-cbux;
    xb(xb==0)=t  ;
    
    
end
f=transpose(cxux)*xb;
if f == 0
    disp("Infeasible")
else
    
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
            disp("Unbounded");
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
end