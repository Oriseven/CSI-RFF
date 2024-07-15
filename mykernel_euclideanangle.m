function G=mykernel_euclideanangle(U,V)
G=zeros(size(U,1),size(V,1));
for i=1:size(U,1)
    for j=1:size(V,1)
        num =sum(U(i,1:52).*V(j,1:52)) + sum(U(i,53:end).*V(j,53:end));
        den = norm(U(i,:))*norm(V(j,:));
        dis = abs(1 - num / den);
        G(i,j) = exp(-dis^2);
    end
end
end