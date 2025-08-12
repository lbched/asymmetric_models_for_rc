function [vechA]=vech(A,k)
% Put A into vech of A
vechA=zeros(k*(k+1)/2,1);
l=1;
for j=1:k
for i=j:k
vechA(l)=A(i,j);
l=l+1;
end
end
end