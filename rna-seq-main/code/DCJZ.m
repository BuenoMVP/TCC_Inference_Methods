function B=DCJZ(A)
m=size(A)
B=zeros(m)
for i=1:m
    for j=1:m
        if A(i,j)~=0
            B(i,j)=A(i,j)
            B(j,i)=A(i,j)
        else
        end
    end
end
end