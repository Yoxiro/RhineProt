function FV=FEGS(data)
[P,V]=coordinate;
s=fastaread(strcat(num2str(data),'.fasta'));
data={s(:).Sequence};
l=length(data);
parfor i=1:l
    g_p{i}=GRS(data{i},P,V);
    for u=1:158
        EL(i,u)=ME(g_p{i}{u});
    end
end
char='ARNDCQEGHILKMFPSTWYV';
parfor i=1:l
    [AAC,DIC]=SAD(data{i},char);
    FA(i,:)=AAC';
    FD(i,:)=DIC(:)';
end
FV=[EL FA FD];