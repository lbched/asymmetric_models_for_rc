function St = SctrPNtauMrec(A1,A2,A3,A4,B,yt,ytp,ytn,ytmp,ytmn,n,T)
St = zeros(n,n,T);
S1 = mean(yt(:,:,:),3);
P = mean(ytp(:,:,:),3);
N1 = mean(ytn(:,:,:),3);
MP = mean(ytmp(:,:,:),3);
MN = mean(ytmn(:,:,:),3);
C = (1-B^2)*S1-(A1^2*P*S1^(-1)+A2^2*N1*S1^(-1)+A3^2*MP*S1^(-1)+A4^2*MN*S1^(-1))*S1;
St(:,:,1) = S1;
for t = 2:T
St(:,:,t) = C+A1^2*ytp(:,:,t-1)+A2^2*ytn(:,:,t-1)+A3^2*ytmp(:,:,t-1)+A4^2*ytmn(:,:,t-1)+B^2*St(:,:,t-1);
end
end