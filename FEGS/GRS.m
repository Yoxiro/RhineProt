function g=GRS(seq,P,V)
load M
l_seq=length(seq);
k=size(M,1);
cha='ACDEFGHIKLMNPQRSTVWY';
for j=1:k
    DPC=zeros(20,20);
    c{1}=[0 0 0];
    d=zeros(1,3);
    for i=1:l_seq
        if i==1
            x=seq(i)==M(j,:);
            c{i+1}=c{i}+x*P;
        else
            x=seq(i)==M(j,:);
            if all(x==0)
                d=d*(i-2)/(i-1);
                c{i+1}=c{i}+[0 0 1]+d; 
            elseif all(y==0)
                d=d*(i-2)/(i-1);
                c{i+1}=c{i}+x*P+d;
            else
                d=d*(i-2)/(i-1)+V{y==1,x==1}/(i-1);
                c{i+1}=c{i}+x*P+d;
            end
        end
        y=x;
    end
    g{j}=cell2mat(c');
    clear c
end