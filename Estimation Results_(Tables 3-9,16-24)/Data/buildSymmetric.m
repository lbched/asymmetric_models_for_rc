function K = buildSymmetric(V,n)
if length(V) ~= (n*(n+1)/2)
error('V is not of a length consistent with n')
end
K = zeros(n,n);
ind = 0;
for cind = 1:n
for rind = cind:n
ind = ind + 1;
K(rind,cind) = V(ind);
if rind > cind
K(cind,rind) = V(ind);
end
end
end
end