function St = SctrPNM_semirec(A1,A2,A3,A4,A5,A6,B,yt,ytp,ytn,ytm,ytp1,ytn1,ytm1,n,T)
St = zeros(n,n,T);
S1 = mean(yt(:,:,:),3);
P = mean(ytp(:,:,:),3);
N1 = mean(ytn(:,:,:),3);
M = mean(ytm(:,:,:),3);
P1 = mean(ytp1(:,:,:),3);
N2 = mean(ytn1(:,:,:),3);
M1 = mean(ytm1(:,:,:),3);
C = (1-B^2)*S1-(A1^2*P*S1^(-1)+A2^2*N1*S1^(-1)+A3^2*M*S1^(-1)+A4^2*P1*S1^(-1)+A5^2*N2*S1^(-1)+A6^2*M1*S1^(-1))*S1;
St(:,:,1) = S1;
for t = 2:T
St(:,:,t) = C+A1^2*ytp(:,:,t-1)+A2^2*ytn(:,:,t-1)+A3^2*ytm(:,:,t-1)+A4^2*ytp1(:,:,t-1)+A5^2*ytn1(:,:,t-1)+A6^2*ytm1(:,:,t-1)+B^2*St(:,:,t-1);
end
end