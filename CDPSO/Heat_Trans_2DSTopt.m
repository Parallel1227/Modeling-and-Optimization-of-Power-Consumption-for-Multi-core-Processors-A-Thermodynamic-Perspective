
function diss_sum = Heat_Trans_2DSTopt(S_index)
% S_index = [1 3 5 7 2 4 12 10 9 11 8 6];
Ntr = [8.8e+8  1.4e+8  6.1e+8   2.5e+8	1.1e+9   7.2e+8];
t_opr = [1 1 1 1 1 1 1 1 1 1 1 1];
t_opr = t_opr*0.1;
Ntr = repmat(Ntr,1,2);
Ntr = Ntr(S_index);
% t_opr = t_opr(S_index);
Ntr_order = 1:1:12;
Ntr_order = Ntr_order(S_index);
Ntr = reshape(Ntr,3,4);
t_opr = reshape(t_opr,3,4);
Ntr_order = reshape(Ntr_order,3,4);
Rbit = 1;
f0 = 3;
R = 8e+7;
diss = zeros(4,3);

delta_t = 0.1;
t_max = 100;
T_env = 300;
T_ini = 300;
t_1 = zeros(1,6);
t_2 = zeros(1,6);
times = zeros(1,6);
flag = zeros(12*(t_max/delta_t),6);
t_opr = t_opr./delta_t;
T_chip1 = zeros(t_max/delta_t,1);
T_chip2 = zeros(t_max/delta_t,1);
T_chip3 = zeros(t_max/delta_t,1);
T_chip4 = zeros(t_max/delta_t,1);
T_chip5 = zeros(t_max/delta_t,1);
T_chip6 = zeros(t_max/delta_t,1);
T_chip7 = zeros(t_max/delta_t,1);
T_chip8 = zeros(t_max/delta_t,1);
T_chip9 = zeros(t_max/delta_t,1);
T_chip10 = zeros(t_max/delta_t,1);
T_chip11 = zeros(t_max/delta_t,1);
T_chip12 = zeros(t_max/delta_t,1);
t1 = 1;
t2 = 1;
t3 = 1;
t4 = 1;
t5 = 1;
t6 = 1;
t7 = 1;
t8 = 1;
t9 = 1;
t10 = 1;
t11 = 1;
t12 = 1;
t_diss1 = 1;
t_diss2 = 1;
t_diss3 = 1;
t_diss4 = 1;
t_diss5 = 1;
t_diss6 = 1;
t_diss7 = 1;
t_diss8 = 1;
t_diss9 = 1;
t_diss10 = 1;
t_diss11 = 1;
t_diss12 = 1;
diss_chip_sum = zeros(3,2);
l_x1 = 0.005;%核x方面尺寸 m
l_x2 = 0.005;%核x方面间距 m
l_y1 = 0.005;%核y方面尺寸 m
l_y2 = 0.005;%核y方面间距 m
l_c = 0.002;%核厚度 m
l_sub = 0.0001;%衬底厚度 m
l_pcb = 0.002;%PCB厚度 m
lam_chip = 102;%核导热系数 W/(m*K)
lam_air = 0.026;%空气传热系数 W/(m*K)
lam_sub = 101;%衬底导热系数 W/(m*K)
lam_pcb = 0.6;%PCB导热系数 W/(m*K)
alpha_air = 26.3;%空气自然对流换热系数 W/(m*K)
alpha_h = 50;%核与衬底换热系数 W/(m*K)
alpha_sink = 401;%芯片与散热片换热系数 W/(m*K)
m_chip = 0.001;%核质量kg
c_chip = 700;%核/硅 比热容 J/(kg*K)

T = [T_env  T_env   T_env   T_env   T_env   T_env;
     T_env  T_ini   T_ini   T_ini   T_ini   T_env;
     T_env  T_ini   T_ini   T_ini   T_ini   T_env;
     T_env  T_ini   T_ini   T_ini   T_ini   T_env;
     T_env  T_env   T_env   T_env   T_env   T_env];



R_chip_sub = (0.5*l_c + 0.5*l_sub)/(l_x1*l_y1*alpha_h);
R_x = (l_x1+l_x2)/l_y1/l_sub/lam_sub + R_chip_sub*2;
R_y = (l_y1+l_y2)/l_x1/l_sub/lam_sub + R_chip_sub*2;

flag_t = zeros(1,6);
for i = 2:1:4
    for j = 2:1:5
        if mod(Ntr_order(i-1,j-1),6) == 0
            index_order = 6;
        else
            index_order = mod(Ntr_order(i-1,j-1),6);
        end
        if flag_t(index_order) == 0
            t_1(index_order) = t_opr(i-1,j-1);
        else
            t_2(index_order) = t_opr(i-1,j-1);
        end
        flag_t(index_order) = 1;
        times(index_order) = floor(t_max/(t_1(index_order)+t_2(index_order))./delta_t);
    end
end

 for t = 1:1:12*(t_max/delta_t)
     for i = 2:1:4
         for j = 2:1:5
             if mod(Ntr_order(i-1,j-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(i-1,j-1),6);
             end
             Q_x = (T(i+1,j)-T(i,j))/R_x + (T(i-1,j)-T(i,j))/R_x;
             Q_y = (T(i,j+1)-T(i,j))/R_y + (T(i,j-1)-T(i,j))/R_y;
             Q_air = (T_env - T(i,j))*(l_x1*l_y1*alpha_sink);
             Q_pcb = (T_env - T(i,j))*l_x1*l_y1/(1/alpha_h+1/alpha_air+l_sub/lam_sub+l_pcb/lam_pcb);
             if flag(t,index_order) == 0
                 if mod(floor(t/12),(t_1(index_order)+t_2(index_order))) < t_1(index_order)
                     S(i-1,j-1) = Chip_Power(Ntr(i-1,j-1),Rbit,f0,R,T(i,j));
                 else
                     S(i-1,j-1) = 0;
                 end
             else
                 if mod(floor(t/12),(t_1(index_order)+t_2(index_order))) < t_1(index_order)
                     S(i-1,j-1) = 0;
                 else
                     S(i-1,j-1) = Chip_Power(Ntr(i-1,j-1),Rbit,f0,R,T(i,j));
                 end
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                 diss(i-1,j-1) = S(i-1,j-1);
             end
             T(i,j) = T(i,j) + ((Q_x + Q_y + Q_air + Q_pcb + S(i-1,j-1)))/m_chip/c_chip*delta_t;
             flag(t,index_order) = 1;
         end
     end
     switch (mod(t,12))
         case 1
             T_chip1(t1) = T(2,2);
             if mod(Ntr_order(2-1,2-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(2-1,2-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip1(t_diss1) = diss(1,1);
                diss_chip_sum(1,1) = diss_chip_sum(1,1) + (diss_chip1(t_diss1)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss1 = t_diss1 + 1;
             end
             t1 = t1 + 1;
         case 2
             T_chip2(t2) = T(2,3);
             if mod(Ntr_order(2-1,3-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(2-1,3-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip2(t_diss2) = diss(1,2);
                diss_chip_sum(1,2) = diss_chip_sum(1,2) + (diss_chip2(t_diss2)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss2 = t_diss2 + 1;
             end
             t2 = t2 + 1;
         case 3
             T_chip3(t3) = T(2,4);
             if mod(Ntr_order(2-1,4-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(2-1,4-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip3(t_diss3) = diss(1,3);
                diss_chip_sum(1,1) = diss_chip_sum(1,1) + (diss_chip3(t_diss3)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss3 = t_diss3 + 1;
             end
             t3 = t3 + 1;
         case 4
             T_chip4(t4) = T(2,5);
             if mod(Ntr_order(2-1,5-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(2-1,5-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip4(t_diss4) = diss(1,4);
                diss_chip_sum(1,2) = diss_chip_sum(1,2) + (diss_chip4(t_diss4)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss4 = t_diss4 + 1;
             end
             t4 = t4 + 1;
         case 5
             T_chip5(t5) = T(3,2);
             if mod(Ntr_order(3-1,2-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(3-1,2-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip5(t_diss5) = diss(2,1);
                diss_chip_sum(2,1) = diss_chip_sum(2,1) + (diss_chip5(t_diss5)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss5 = t_diss5 + 1;
             end
             t5 = t5 + 1;
         case 6
             T_chip6(t6) = T(3,3);
             if mod(Ntr_order(3-1,3-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(3-1,3-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip6(t_diss6) = diss(2,2);
                diss_chip_sum(2,2) = diss_chip_sum(2,2) + (diss_chip6(t_diss6)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss6 = t_diss6 + 1;
             end
             t6 = t6 + 1;
         case 7
             T_chip7(t7) = T(3,4);
             if mod(Ntr_order(3-1,4-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(3-1,4-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip7(t_diss7) = diss(2,3);
                diss_chip_sum(2,1) = diss_chip_sum(2,1) + (diss_chip7(t_diss7)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss7 = t_diss7 + 1;
             end
             t7 = t7 + 1;
         case 8
             T_chip8(t8) = T(3,5);
             if mod(Ntr_order(3-1,5-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(3-1,5-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip8(t_diss8) = diss(2,4);
                diss_chip_sum(2,2) = diss_chip_sum(2,2) + (diss_chip8(t_diss8)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss8 = t_diss8 + 1;
             end
             t8 = t8 + 1;
         case 9
             T_chip9(t9) = T(4,2);
             if mod(Ntr_order(4-1,2-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(4-1,2-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip9(t_diss9) = diss(3,1);
                diss_chip_sum(3,1) = diss_chip_sum(3,1) + (diss_chip9(t_diss9)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss9 = t_diss9 + 1;
             end
             t9 = t9 + 1;
         case 10
             T_chip10(t10) = T(4,3);
             if mod(Ntr_order(4-1,3-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(4-1,3-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip10(t_diss10) = diss(3,2);
                diss_chip_sum(3,2) = diss_chip_sum(3,2) + (diss_chip10(t_diss10)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss10 = t_diss10 + 1;
             end
             t10 = t10 + 1;
         case 11
             T_chip11(t11) = T(4,4);
             if mod(Ntr_order(4-1,4-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(4-1,4-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip11(t_diss11) = diss(3,3);
                diss_chip_sum(3,1) = diss_chip_sum(3,1) + (diss_chip11(t_diss11)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss11 = t_diss11 + 1;
             end
             t11 = t11 + 1;
         otherwise
             T_chip12(t12) = T(4,5);
             if mod(Ntr_order(4-1,5-1),6) == 0
                 index_order = 6;
             else
                 index_order = mod(Ntr_order(4-1,5-1),6);
             end
             if floor(floor(t/12)/(t_1(index_order)+t_2(index_order))) == times(index_order)-1
                diss_chip12(t_diss12) = diss(3,4);
                diss_chip_sum(3,2) = diss_chip_sum(3,2) + (diss_chip12(t_diss12)*delta_t)/(t_1(index_order)+t_2(index_order))./delta_t;
                t_diss12 = t_diss12 + 1;
             end
             t12 = t12 + 1;
     end
 end
t = 1:1:(t_max/delta_t);
diss_sum = sum(sum(diss_chip_sum));
 
