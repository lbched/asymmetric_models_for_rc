function [L]=elimination(k)
% Create elimination matrix
I=eye(k*(k+1)/2);
D=zeros(k^2,k*(k+1)/2);
for j=1:k*(k+1)/2
for i=1:k
if (i>=j)
D((j-1)*k+i,:)=I((j-1)*(k-j/2)+i,:);
D((i-1)*k+j,:)=I((j-1)*(k-j/2)+i,:);
end
end
end
L=(D'*D)\D';
end