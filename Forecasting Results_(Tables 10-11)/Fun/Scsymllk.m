function L = Scsymllk(x0,yt,n,T)
A = (x0(1));
B = (x0(2));
St = Scsymrec(A,B,yt,n,T);
L = zeros(T,1);
for j = 1:T
d = det(St(:,:,j));
a = trace(St(:,:,j)\yt(:,:,j));
L(j) = -0.5*log(d)-0.5*(a);
end
end