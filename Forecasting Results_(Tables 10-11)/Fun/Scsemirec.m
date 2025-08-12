function St = Scsemirec(A1,A2,A3,B,yt,ytp,ytn,ytm,n,T)
St = zeros(n,n,T);
S1 = mean(yt(:,:,:),3);
P = mean(ytp(:,:,:),3);
N1 = mean(ytn(:,:,:),3);
M = mean(ytm(:,:,:),3);
C = (1-B^2)*S1-(A1^2*P*S1^(-1)+A2^2*N1*S1^(-1)+A3^2*M*S1^(-1))*S1;
St(:,:,1) = S1;
for t = 2:T
St(:,:,t) = C+A1^2*ytp(:,:,t-1)+A2^2*ytn(:,:,t-1)+A3^2*ytm(:,:,t-1)+B^2*St(:,:,t-1);
end
end