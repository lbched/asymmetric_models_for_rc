function St = Scsymrec(A,B,yt,n,T)
St = zeros(n,n,T);
S1 = mean(yt(:,:,:),3);
C = (1-A^2-B^2)*S1;
St(:,:,1) = S1;
for t = 2:T
St(:,:,t) = C+A^2*yt(:,:,t-1)+B^2*St(:,:,t-1);
end
end