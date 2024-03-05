function [AAC,DPC]=SAD(seq,a)
len_seq=length(seq);
len_a=length(a);
for i=1:len_a
    c{i}=seq==a(i);
    AAC(i,1)=sum(c{i})/len_seq; 
end
if len_seq~=1
    for i=1:len_a
        for j=1:len_a
            DPC(i,j)= sum(([c{j}(1,2:len_seq) 0]*2-c{i})==1)/(len_seq-1);  
        end
    end
else
     DPC=zeros(20,20);
end


