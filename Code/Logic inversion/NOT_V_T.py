from pylab import *
import numpy as np
import scipy
from scipy.linalg import null_space
from scipy import integrate
from matplotlib import rcParams

config = {
    "font.family": 'Times new roman',
    "font.size": '13',
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


def NOT_propagation(Vin,Vout,V_D,kBT):
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
Vout_NOT2 = 0.0
Vout_NOT3 = 0.0
dissipation_NOT1 = 0.0
dissipation_NOT2 = 0.0
dissipation_NOT3 = 0.0


V_D = 5.0
tint = 1000
T = 2001000

# Vout=0.0


Ntot = int(T / tint)
output1 = np.zeros(Ntot)
diss1 = np.zeros(Ntot)
output2 = np.zeros(Ntot)
diss2 = np.zeros(Ntot)
output3 = np.zeros(Ntot)
diss3 = np.zeros(Ntot)
time = np.zeros(Ntot)
kBT1 = 1.0
kBT2 = 1.33
kBT3 = 1.66

for i in range(Ntot):
    Vout_NOT1, dissipation_NOT1 = NOT_propagation(0.0, Vout_NOT1, V_D, kBT1)
    Vout_NOT2, dissipation_NOT2 = NOT_propagation(0.0, Vout_NOT2, V_D, kBT2)
    Vout_NOT3, dissipation_NOT3 = NOT_propagation(0.0, Vout_NOT3, V_D, kBT3)
    output1[i] = Vout_NOT1
    output2[i] = Vout_NOT2
    output3[i] = Vout_NOT3
    diss1[i] = diss1[i - 1] + dissipation_NOT1
    diss2[i] = diss2[i - 1] + dissipation_NOT2
    diss3[i] = diss3[i - 1] + dissipation_NOT3

    time[i] = i * tint

# x = np.linspace(0,100)
ntime = np.zeros(21)
ndiss1 = np.zeros(21)
noutput1 = np.zeros(21)
ndiss2 = np.zeros(21)
noutput2 = np.zeros(21)
ndiss3 = np.zeros(21)
noutput3 = np.zeros(21)

for i in range(21):
    ntime[i] = time[i * 100]
    ndiss1[i] = diss1[i * 100]
    noutput1[i] = output1[i * 100]
    ndiss2[i] = diss2[i * 100]
    noutput2[i] = output2[i * 100]
    ndiss3[i] = diss3[i * 100]
    noutput3[i] = output3[i * 100]


# plt.rc('font',family = 'Times New Roman')
# plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
plt.figure(1, figsize=(5, 3.4), dpi=300)
plt.tick_params(labelsize=13)
plt.rc('font', size=13,family='Times New Roman')  # 将字体大小更改为12
ax = plt.gca()
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
ax.spines['right'].set_linewidth(1.5)
ax.spines['top'].set_linewidth(1.5)
plt.plot(ntime, noutput1, label='$T_{i, j}=300K$', marker="+", color='#4C68B0')
plt.plot(ntime, noutput2, label='$T_{i, j}=400K$', marker="*", color='#548D8A')
plt.plot(ntime, noutput3, label='$T_{i, j}=500K$', marker="x", color='#D3A518')
plt.hlines(4.75, 0, 2000000, colors="r", linestyles="dashed",label='Judgment threshold')
# plt.plot(ntime, noutput3, label='$\mathrm{X}=0, \mathrm{Y}=1$', marker="*", color='limegreen')
# plt.plot(ntime, noutput4, label='$\mathrm{X}=0, \mathrm{Y}=0$', marker="1", color='slategrey')
# plt.rcParams['font.sans-serif']=['SimHei'] #用来正常显示中文标签
# plt.title('前一时刻输入为11时，XOR门在不同输入下的时间/输出电压曲线')
plt.xlabel("Operation time($β\hbar$)", fontsize=13,fontproperties='Times New Roman')
plt.ylabel("Output voltage($V_{T}$)", fontsize=13,fontproperties='Times New Roman')
plt.legend(loc = 'lower right', fontsize=11)
plt.grid()
plt.xlim(0, 2000000)
plt.ylim(4, 5)
plt.yticks(fontproperties='Times New Roman')
plt.xticks(fontproperties='Times New Roman')
ax = plt.gca()  #获取当前图像的坐标轴信息
xfmt = ScalarFormatter(useMathText=True)
xfmt.set_powerlimits((0, 0))  # Or whatever your limits are . . .
gca().xaxis.set_major_formatter(xfmt)
# plt.text(-200000, 5, '(a)', fontsize=15, fontproperties='Times New Roman')


plt.tight_layout()
plt.savefig('NOT_输出电压_温度.svg', format='svg')


