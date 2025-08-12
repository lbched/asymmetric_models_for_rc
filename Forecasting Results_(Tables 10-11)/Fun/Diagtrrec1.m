function St = Diagtrrec1(A1,A2,B,ytpm,ytn,n,T,Q1,P,N1)
Qt = zeros(n*(n+1)/2,T);
St = zeros(n,n,T);
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