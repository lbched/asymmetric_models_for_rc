function L = Diagsymllk(x0,yt,n,T)
A=diag((x0(1:n)));
B=diag((x0((n+1):(2*n))));
St=Diagsymrec(A,B,yt,n,T);
L=zeros(T,1);
for j=1:T
d=det(St(:,:,j));
a=trace(St(:,:,j)\yt(:,:,j));
L(j)=-0.5*log(d)-0.5*(a);
end
end