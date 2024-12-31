function P = Chip_Power(ntr,rbit,f0,R,T)

P_static = 30.99*T - 8406;
P_dynamic = -0.3457*T + 567.2;
P = 2 * ntr * rbit * f0 * R * P_dynamic + ntr * rbit * R * P_static;
P = P * 4.14e-21;