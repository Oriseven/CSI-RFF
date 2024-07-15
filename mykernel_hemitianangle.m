function G=mykernel_hemitianangle(U,V)
G=zeros(size(U,1),size(V,1));
 
for i=1:size(U,1)
    for j=1:size(V,1)
        x = U(i,1:52)+1i*U(i,53:end);
        y = V(j,1:52)+1i*V(j,53:end);

        dsq=abs(sum(x'.*y.'));   
        dis = 1 - dsq / (norm(x)*norm(y));
        G(i,j) = exp(-dis^2);
    end
end
end