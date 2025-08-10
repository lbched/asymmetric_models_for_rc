function St = PLTtrrec(A1,A2,B,yt,ytpm,ytn,n,T)
Qt = zeros(n*(n+1)/2,T);
St = zeros(n,n,T);
Q1 = vech(mean(yt(:,:,:),3),n);
P=kron(mean(ytpm(:,:,:),3)^(1/2)*mean(yt(:,:,:),3)^(-1/2),mean(ytpm(:,:,:),3)^(1/2)*mean(yt(:,:,:),3)^(-1/2));
N1=kron(mean(ytn(:,:,:),3)^(1/2)*mean(yt(:,:,:),3)^(-1/2),mean(ytn(:,:,:),3)^(1/2)*mean(yt(:,:,:),3)^(-1/2));
A1A1=kron(A1,A1);
BB=kron(B,B);
A2A2=kron(A2,A2);
C=(eye(n*(n+1)/2)-(elimination(n)*A1A1*duplication(n))*(elimination(n)*P*duplication(n))-(elimination(n)*A2A2*duplication(n))*(elimination(n)*N1*duplication(n))-(elimination(n)*BB*duplication(n)))*Q1;
Qt(:,1) = Q1;
for t=2:T
Qt(:,t) = C+(elimination(n)*A1A1*duplication(n))*vech(ytpm(:,:,t-1),n)+(elimination(n)*A2A2*duplication(n))*vech(ytn(:,:,t-1),n)+(elimination(n)*BB*duplication(n))*Qt(:,t-1);
end
for j=1:T
St(:,:,j)=buildSymmetric(Qt(:,j),n);
end
end