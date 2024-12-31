function change = constant_times_velocity(constant,change)
for i=1:size(change,1)
    for j=1:size(change,2)
        if rand>constant
            change(i,j)=0;
        end
    end
end
end
