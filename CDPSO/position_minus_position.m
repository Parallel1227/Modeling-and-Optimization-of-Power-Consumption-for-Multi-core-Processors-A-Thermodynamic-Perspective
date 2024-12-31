function change=position_minus_position(best,pop)
for i=1:size(best,1)
    for j=1:size(best,2)
        change(i,j)=find(pop(i,:)==best(i,j));
        temp=pop(i,j);
        pop(i,j)=pop(i,change(i,j));
        pop(i,change(i,j))=temp;
    end
end
end
