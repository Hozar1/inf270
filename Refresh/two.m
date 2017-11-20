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

bas1=n+1:m+n;
nbas1=1:n;
B1=A(:,bas1);
N1=A(:,nbas1);
cb1=c(nbas1);
cx1=c(bas1);

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
    cb=-cb;
    xb(xb==0)=t;
    
end

%LUFact

[L,U,P] = lu(B);
% aj=P*b;
% dxy=linsolve(L,aj);
% dxx=linsolve(U,dxy);




EtaFile = cell(21,1);
EtaCounter = 1;

Ka=[(nnz(L)+nnz(U)/m),sqrt(m),20];
K=min(Ka);
while any(cb1(:)>0)
    
    
    dxy=linsolve(L,aj);
    dxx=linsolve(U,dxy);
    
    %Step2
    [tj,j]=max(cb1);
    ji=j;
    j=nbas(j);
    %step3
    ej=zeros(size(cb1));
    ej(ji)=1;
    %step4
    [t,tx]=(max(dxx./xb));
    t=inv(t);
    %step5 dont think this is correct
    i=bas(tx);
    %step 6
    ei=zeros(size(cx));
    ei(i==bas)=1;
    %aj == dx   
    aj=N1*ej
    ai=B1*ei;
    
    
    %dZn
    v=linsolve(transpose(B1),ei);
    dZn=-1*transpose(N1)*v;
    s=-tj/dZn(ji);
    cb1=-cb1-s*dZn;
    
    cb1(cb1==0)=s;
    
    
    
    %Update B and update A? update C
    B1=B1 + (aj-ai)*transpose(ei);
    N1=N1 + (ai-aj)*transpose(ej);
    Eta = eye(m) + (dxx -ei)*transpose(ei);
    
    
   % cb1=-cb1
  %  xb1(xb1==0)=t
    EtaFile{EtaCounter,1} = Eta;
    EtaCounter = EtaCounter +1;
end


cx=c(bas);
f=transpose(cx)*xb;
disp(xb);
disp(cb);
disp(f);