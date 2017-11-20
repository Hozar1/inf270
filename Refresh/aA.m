m=3;
n=2; 
c=[4;3];
z=c*-1;
b=[1;3;5];
A=[1 -1; 2 -1; 0 1];
N=A;
A=[A eye(m)]; 
c=[c;zeros(m,1)];
    
Abas=[3 4 5];
Bbas=[1 4 5];
Cbas=[1 2 5];
Dbas=[1 2 3];

nbas=[1 2];
B=A(:,Dbas);

%xb basic variables
xb =inv(B) * b;
cb=c(bas);
disp ("    Xb   Xn");
disp (xb + "  "+  cb);
%xn NBV

