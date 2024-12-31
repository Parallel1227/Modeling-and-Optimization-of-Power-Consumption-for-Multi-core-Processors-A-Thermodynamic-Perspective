function pop = position_plus_velocity(pop,v)
for i=1:size(pop,1)
    for j=1:size(pop,2)
        if v(i,j)~=0
            temp=pop(i,j);
            pop(i,j)=pop(i,v(i,j));
            pop(i,v(i,j))=temp;
        end
    end
end
end
