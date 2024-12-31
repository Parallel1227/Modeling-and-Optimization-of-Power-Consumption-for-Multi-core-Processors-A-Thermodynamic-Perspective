clc;
close all;
clear all;

n=12;
c1=0.1;                         
c2=0.075;                       
w=1;                            
m=100;                          
pop=zeros(m,n);                 
v=zeros(m,n);                   
gen=1;                          
genmax=100;                     
fitness=zeros(m,1);             
Pbest=zeros(m,n);               
Pbest_fitness=zeros(m,1);       
Gbest=zeros(genmax,n);          
Gbest_fitness=zeros(genmax,1);  
Length_ave=zeros(genmax,1);     
ws=1;                           
we=0.8;                         

for i=1:m
    pop(i,:)=randperm(n);
    v(i,:)=randperm(n);
end

for i=1:m
    fitness(i)=Heat_Trans_2DSTopt(pop(i,:));
end

Pbest_fitness=fitness;
Pbest=pop;
[Gbest_fitness(1),min_index]=min(fitness);
Gbest(1,:)=pop(min_index,:);
Length_ave(1)=mean(fitness);

while gen<genmax
    gen
    gen=gen+1;
    w = ws - (ws-we)*(gen/genmax)^2;

    change1=position_minus_position(Pbest,pop);
    change1=constant_times_velocity(c1,change1);

    change2=position_minus_position(repmat(Gbest(gen-1,:),m,1),pop);
    change2=constant_times_velocity(c2,change2);

    v=constant_times_velocity(w,v);

    for i=1:m
        for j=1:n
            if change1(i,j)~=0
                v(i,j)=change1(i,j);
            end
            if change2(i,j)~=0
                v(i,j)=change2(i,j);
            end
        end
    end

    pop=position_plus_velocity(pop,v);

    fitness=zeros(m,1);
    for i=1:m
        fitness(i)=Heat_Trans_2DSTopt(pop(i,:));
    end

    for i=1:m
        if fitness(i)<Pbest_fitness(i)
            Pbest_fitness(i)=fitness(i);
            Pbest(i,:)=pop(i,:);
        end
    end

    [minvalue,min_index]=min(fitness);
    if minvalue<Gbest_fitness(gen-1)
        Gbest_fitness(gen)=minvalue;
        Gbest(gen,:)=pop(min_index,:);
    else
        Gbest_fitness(gen)=Gbest_fitness(gen-1);
        Gbest(gen,:)=Gbest(gen-1,:);
    end
    Length_ave(gen)=mean(fitness);
end

[Shortest_Length,index] = min(Gbest_fitness);
Shortest_Route = Gbest(index,:);

figure;
plot(1:genmax,Gbest_fitness,'b',1:genmax,Length_ave,'r:')
