from pylab import *
import numpy as np
import scipy
from scipy.linalg import null_space
from scipy import integrate
from matplotlib import rcParams
import pandas as pd

config = {
    "font.family": 'serif',
    "mathtext.fontset": 'stix',
    "font.serif": ['SimSun'],
}
rcParams.update(config)


# plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签


def Fermi(x):
    return 1.0 / (exp(x) + 1)


def Bose(x):
    if (x > 1e-4):
        return 1.0 / (exp(x) - 1)
    else:
        return 1e5


def Gauss_in0(x, a):
    sigma2 = 0.005
    mu = a
    return 1 / sqrt(2 * np.pi * sigma2) * exp(-1.0 / 2 / sigma2 * (x - mu) * (x - mu))


def NOT_propagation(Vin,Vout,V_D):
    A=zeros((4,4))
    c=zeros(4)

    Gamma_l=0.2
    Gamma_r=0.2
    Gamma=0.2
    Gamma_g=0.2

    E_P=Vin
    E_N=1.5*V_D-Vin
    mu_l=0.0
    mu_r=mu_l-V_D
    kBT=1
    Cg=200.0

    mu_g=0.0-Vout
    k_Nl=Gamma_l*Fermi((E_N-mu_l)/kBT)
    k_lN=Gamma_l*(1.0-Fermi((E_N-mu_l)/kBT))
    k_rP=Gamma_r*(1.0-Fermi((E_P-mu_r)/kBT))
    k_Pr=Gamma_r*Fermi((E_P-mu_r)/kBT)
    k_Ng=Gamma_g*Fermi((E_N-mu_g)/kBT)
    k_gN=Gamma_g*(1.0-Fermi((E_N-mu_g)/kBT))
    k_Pg=Gamma_g*Fermi((E_P-mu_g)/kBT)
    k_gP=Gamma_g*(1.0-Fermi((E_P-mu_g)/kBT))
    if(E_N>E_P):
        k_NP=Gamma*Bose((E_N-E_P)/kBT)
        k_PN=Gamma*(1+Bose((E_N-E_P)/kBT))
    else:
        k_PN=Gamma*Bose((E_P-E_N)/kBT)
        k_NP=Gamma*(1+Bose((E_P-E_N)/kBT))

    A[1][0]=k_Pr+k_Pg
    A[2][0]=k_Nl+k_Ng
    A[0][0]=-A[1][0]-A[2][0]
    A[0][1]=k_rP+k_gP
    A[2][1]=k_NP
    A[3][1]=k_Nl+k_Ng
    A[1][1]=-A[0][1]-A[2][1]-A[3][1]
    A[0][2]=k_lN+k_gN
    A[1][2]=k_PN
    A[3][2]=k_Pr+k_Pg
    A[2][2]=-A[0][2]-A[1][2]-A[3][2]
    A[1][3]=k_lN+k_gN
    A[2][3]=k_rP+k_gP
    A[3][3]=-A[1][3]-A[2][3]

    p=null_space(A)
    sum=p[0][0]+p[1][0]+p[2][0]+p[3][0]
    p[0][0]=p[0][0]/sum
    p[1][0]=p[1][0]/sum
    p[2][0]=p[2][0]/sum
    p[3][0]=p[3][0]/sum
    p_N=p[2][0]+p[3][0]
    p_P=p[1][0]+p[3][0]
    J1=k_gN*p_N-k_Ng*(1-p_N)
    J2=k_gP*p_P-k_Pg*(1-p_P)
    J3=k_Nl*(1.0-p_N)-k_lN*p_N
    J4=k_rP*p_P-k_Pr*(1.0-p_P)

    Vout-=1.0*(J1+J2)*tint/Cg
    dissipation=1.0*J3*tint*(mu_l-mu_g)-1.0*J4*tint*(mu_r-mu_g)

    return Vout,dissipation

Vout_NOT1 = 0.0
dissipation_NOT1 = 0.0

Vout_NOT2 = 0.0
dissipation_NOT2 = 0.0


V_D = 5.0
tint = 1000
T = 1001000

# Vout=0.0


Ntot = int(T / tint)
output1 = np.zeros(Ntot)
diss1 = np.zeros(Ntot)
output2 = np.zeros(Ntot)
diss2 = np.zeros(Ntot)
time = np.zeros(Ntot)
error1 = np.zeros(Ntot)
error2 = np.zeros(Ntot)

for i in range(Ntot):
    Vout_NOT1, dissipation_NOT1 = NOT_propagation(0.0, Vout_NOT1, V_D)
    output1[i] = Vout_NOT1
    diss1[i] = diss1[i - 1] + dissipation_NOT1
    error1[i] = 1 - integrate.quad(Gauss_in0, 0.90 * V_D, np.inf, args=(output1[i]))[0]

    Vout_NOT2, dissipation_NOT2 = NOT_propagation(5.0, Vout_NOT2, V_D)
    output2[i] = Vout_NOT2
    diss2[i] = diss2[i - 1] + dissipation_NOT2
    error2[i] = 1 - integrate.quad(Gauss_in0, -1 * np.inf, 0.10 * V_D, args=(output2[i]))[0]

    time[i] = i * tint


# x = np.linspace(0,100)
ntime = np.zeros(21)
ndiss1 = np.zeros(21)
ndiss2 = np.zeros(21)
noutput1 = np.zeros(21)
noutput2 = np.zeros(21)
nerror1=np.zeros(21)
nerror2=np.zeros(21)

for i in range(21):
    ntime[i] = time[i * 50]
    ndiss1[i] = diss1[i * 50]
    ndiss2[i] = diss2[i * 50]
    noutput1[i] = output1[i * 50]
    noutput2[i] = output2[i * 50]
    nerror1[i] = error1[i * 50]
    nerror2[i] = error2[i * 50]



Vout_NOT1 = 5.0
dissipation_NOT1 = 0.0

Vout_NOT2 = 5.0
dissipation_NOT2 = 0.0


V_D = 5.0
tint = 1000
T = 1001000

# Vout=0.0


Ntot = int(T / tint)
output1 = np.zeros(Ntot)
diss1 = np.zeros(Ntot)
output2 = np.zeros(Ntot)
diss2 = np.zeros(Ntot)
time = np.zeros(Ntot)
error1 = np.zeros(Ntot)
error2 = np.zeros(Ntot)

for i in range(Ntot):
    Vout_NOT1, dissipation_NOT1 = NOT_propagation(0.0, Vout_NOT1, V_D)
    output1[i] = Vout_NOT1
    diss1[i] = diss1[i - 1] + dissipation_NOT1
    error1[i] = 1 - integrate.quad(Gauss_in0, 0.90 * V_D, np.inf, args=(output1[i]))[0]


    Vout_NOT2, dissipation_NOT2 = NOT_propagation(5.0, Vout_NOT2, V_D)
    output2[i] = Vout_NOT2
    diss2[i] = diss2[i - 1] + dissipation_NOT2
    error2[i] = 1 - integrate.quad(Gauss_in0, -1 * np.inf, 0.10 * V_D, args=(output2[i]))[0]

    time[i] = i * tint

# x = np.linspace(0,100)
ntime = np.zeros(21)
ndiss1 = np.zeros(21)
ndiss2 = np.zeros(21)
noutput1 = np.zeros(21)
noutput2 = np.zeros(21)
nerror1 = np.zeros(21)
nerror2 = np.zeros(21)

for i in range(21):
    ntime[i] = time[i * 50]
    ndiss1[i] = diss1[i * 50]
    ndiss2[i] = diss2[i * 50]
    noutput1[i] = output1[i * 50]
    noutput2[i] = output2[i * 50]
    nerror1[i] = error1[i * 50]
    nerror2[i] = error2[i * 50]

dfdata = {
    '时间':ntime,
    '当前输入为1错误率':nerror1,
    '当前输入为0错误率':nerror2
}

df = pd.DataFrame(dfdata)
df.to_excel('output.xlsx', index=False)

# dfdata = {
#     '时间':ntime,
#     '前一输入为1当前输入为1输出':noutput1,
#     '前一输入为1当前输入为0输出':noutput2,
#     '前一输入为1当前输入为1能耗':ndiss1,
#     '前一输入为1当前输入为0能耗':ndiss2,
# }
#
# df = pd.DataFrame(dfdata)
# df.to_excel('output.xlsx', index=False)




