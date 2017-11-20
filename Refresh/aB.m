m=3;
n=2; 
c=[4;3];
z=c*-1;
b=[1;3;5];
A=[1 -1; 2 -1; 0 1];
A=[A eye(m)]; 
c=[c;zeros(m,1)];
    
Abas=[3 4 5];
Bbas=[1 4 5];
Cbas=[1 2 5];
Dbas=[1 2 3];
bas=Cbas;
nbas=[4 3];
B=A(:,bas);
N=A(:,nbas);
%xb basic variables
xb =inv(B) * b;
cb=c(nbas);
cx=c(bas);
f=transpose(cx)*xb;
disp("For the basis ");
g=fprintf(' %d ', bas);
disp(" we get the objective value =" + f);

str="";
disp("    " + str);
str="";
v=1;
while v<size(A,2) -1
    str="w"+v+"=";
    str(end+1) = xb(v);
    w=1;
    while w< size(B,1) +1
        str(end+1) = (-1*B(v,w) + "x"+w);
        w=w+1;
    end
    disp(str);
    v=v+1;
end

f1=transpose(cx)*xb;
f2=transpose(cx) *inv(B) * b
disp(f1 + "=" + f2); 
