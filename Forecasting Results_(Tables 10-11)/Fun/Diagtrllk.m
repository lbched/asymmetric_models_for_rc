function L = Diagtrllk(x0,yt,ytpm,ytn,n,T)
A1=diag((x0(1:n)));
A2=diag((x0((n+1):(2*n))));
B=diag((x0((2*n+1):end)));
St=Diagtrrec(A1,A2,B,yt,ytpm,ytn,n,T);
L=zeros(T,1);
for j=1:T
d=det(St(:,:,j));
a=trace(St(:,:,j)\yt(:,:,j));
L(j)=-0.5*log(d)-0.5*(a);
end
end