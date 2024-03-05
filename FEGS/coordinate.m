function [P,V]=coordinate
for i=1:20
    pt{i}=[cos(i*2*pi/20) sin(i*2*pi/20) 1];
end
P=cell2mat(pt');
for i=1:20
    for j=1:20
        V{i,j}=pt{i}+(1/4)*(pt{j}-pt{i});
    end
end