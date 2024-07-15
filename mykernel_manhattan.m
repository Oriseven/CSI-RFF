function G=mykernel_manhattan(U,V)
G=zeros(size(U,1),size(V,1));
for i=1:size(U,1)
    for j=1:size(V,1)
        G(i,j) = exp(-sum(abs(U(i,:)-V(j,:)))^2);
    end
end
end