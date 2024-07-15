function G=mykernel_euclidean(U,V)
G=zeros(size(U,1),size(V,1));
for i=1:size(U,1)
    for j=1:size(V,1)
        G(i,j) = exp(-norm(U(i,:)-V(j,:))^2);
    end
end
end