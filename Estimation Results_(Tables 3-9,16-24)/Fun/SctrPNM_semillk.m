function L = SctrPNM_semillk(x0,yt,ytp,ytn,ytm,ytp1,ytn1,ytm1,n,T)
A1 = (x0(1));
A2 = (x0(2));
A3 = (x0(3));
A4 = (x0(4));
A5 = (x0(5));
A6 = (x0(6));
B = (x0(7));
St = SctrPNM_semirec(A1,A2,A3,A4,A5,A6,B,yt,ytp,ytn,ytm,ytp1,ytn1,ytm1,n,T);
L = zeros(T,1);
for j = 1:T
d = det(St(:,:,j));
a = trace(St(:,:,j)\yt(:,:,j));
L(j) = -0.5*log(d)-0.5*(a);
end
end