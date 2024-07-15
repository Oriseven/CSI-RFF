function G=mykernel_chebyshev(U,V)
G=zeros(size(U,1),size(V,1));
for i=1:size(U,1)
    for j=1:size(V,1)
        G(i,j) = exp(-max(abs(U(i,:)-V(j,:)))^2);
    end
end
end