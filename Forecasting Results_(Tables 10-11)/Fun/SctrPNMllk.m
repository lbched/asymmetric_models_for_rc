function L = SctrPNMllk(x0,yt,ytp,ytn,ytm,n,T)
A1 = (x0(1));
A2 = (x0(2));
A3 = (x0(3));
B = (x0(4));
St = SctrPNMrec(A1,A2,A3,B,yt,ytp,ytn,ytm,n,T);
L = zeros(T,1);
for j = 1:T
d = det(St(:,:,j));
a = trace(St(:,:,j)\yt(:,:,j));
L(j) = -0.5*log(d)-0.5*(a);
end
end