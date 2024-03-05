function ME=ME(W)
W(1,:)=[];
[x,y]=size(W);
 D=pdist(W);
E=squareform(D);
sdist=zeros(x,x);
for i=1:x
    for j=i:x
        if j-i==1
           sdist(i,j)=E(i,j);
        elseif j-i>1 
            sdist(i,j)=sdist(i,j-1)+E(j-1,j);
        else
            sdist(i,j)=0;
        end
    end
end
sd=sdist+sdist';
sdd=sd+diag(ones(1,x));
L=E./sdd;
ME=eigs(L,1)/x;        