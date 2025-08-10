function L = SctrPNtauMllk(x0,yt,ytp,ytn,ytmp,ytmn,n,T)
A1 = (x0(1));
A2 = (x0(2));
A3 = (x0(3));
A4 = (x0(4));
B = (x0(5));
St = SctrPNtauMrec(A1,A2,A3,A4,B,yt,ytp,ytn,ytmp,ytmn,n,T);
L = zeros(T,1);
for j = 1:T
d = det(St(:,:,j));
a = trace(St(:,:,j)\yt(:,:,j));
L(j) = -0.5*log(d)-0.5*(a);
end
end