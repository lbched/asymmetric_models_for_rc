function St = PLTsymrec1(A,B,yt,n,T,Q1)
Qt = zeros(n*(n+1)/2,T);
St = zeros(n,n,T);
AA=kron(A,A);
BB=kron(B,B);
C=(eye(n*(n+1)/2)-(elimination(n)*AA*duplication(n))-(elimination(n)*BB*duplication(n)))*Q1;
Qt(:,1) = Q1;
for t=2:T
Qt(:,t) = C+(elimination(n)*AA*duplication(n))*vech(yt(:,:,t-1),n)+(elimination(n)*BB*duplication(n))*Qt(:,t-1);
end
for j=1:T
St(:,:,j)=buildSymmetric(Qt(:,j),n);
end
end