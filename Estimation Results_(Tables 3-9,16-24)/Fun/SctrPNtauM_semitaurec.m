function St = SctrPNtauM_semitaurec(A1,A2,A3,A4,A5,A6,A7,A8,B,yt,ytp,ytn,ytmp,ytmn,ytp1,ytn1,ytmp1,ytmn1,n,T)
St = zeros(n,n,T);
S1 = mean(yt(:,:,:),3);
P = mean(ytp(:,:,:),3);
N1 = mean(ytn(:,:,:),3);
MP = mean(ytmp(:,:,:),3);
MN = mean(ytmn(:,:,:),3);
P1 = mean(ytp1(:,:,:),3);
N2 = mean(ytn1(:,:,:),3);
MP1 = mean(ytmp1(:,:,:),3);
MN1 = mean(ytmn1(:,:,:),3);
C = (1-B^2)*S1-(A1^2*P*S1^(-1)+A2^2*N1*S1^(-1)+A3^2*MP*S1^(-1)+A4^2*MN*S1^(-1)+A5^2*P1*S1^(-1)+A6^2*N2*S1^(-1)+A7^2*MP1*S1^(-1)+A8^2*MN1*S1^(-1))*S1;
St(:,:,1) = S1;
for t = 2:T
St(:,:,t) = C+A1^2*ytp(:,:,t-1)+A2^2*ytn(:,:,t-1)+A3^2*ytmp(:,:,t-1)+A4^2*ytmn(:,:,t-1)+A5^2*ytp1(:,:,t-1)+A6^2*ytn1(:,:,t-1)+A7^2*ytmp1(:,:,t-1)+A8^2*ytmn1(:,:,t-1)+B^2*St(:,:,t-1);
end
end