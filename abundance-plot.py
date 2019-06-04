import matplotlib
import matplotlib.pyplot as plt
import numpy as np

wardle_output = np.loadtxt('output.d', skiprows=49)
mhd = np.loadtxt('data/ten.xq',skiprows=3)
mhd_rho = np.loadtxt('data/rhon.xq',skiprows=3)
mhd_time = mhd[:,0]
mhd_temp = mhd[:,1]
mhd_n = mhd_rho[:,1]/(2*1.6737236e-24)

# f = open('output.d')
# lines = f.readlines()
# title = lines[19]

time = wardle_output[:, 0]
dg_file = np.loadtxt('data/dg.xq',skiprows=3)
dg = dg_file[:,1]
# species = wardle_output[0:, 1:]

HCO = wardle_output[:, 2]
HCOp = wardle_output[:,3]
CO = wardle_output[:, 6]
H2 = wardle_output[:,3]
HCN = wardle_output[:, 4]
HNC = wardle_output[:, 5]
CH3OH = wardle_output[:, 7]

fig,ax=plt.subplots()

#ax.loglog(10**time, 10**HNCO, label=r"$HNCO$",linewidth=1.0)
ax.loglog(10**time, 10**HCO, label=r"$HCO$",linewidth=1.0)
ax.loglog(10**time, 10**HCOp, label=r"$HCO+$",linewidth=1.0)
# ax.loglog(10**time, 10**CO, label=r"$CO$",linewidth=1.0)
#ax.loglog(10**time, 10**HNC, label=r"$HNC$",linewidth=1.0)
# ax.loglog(10**time, 10**CH3OH, label=r"$CH_{3}OH$", linewidth=1.0)
# ax.loglog(10**time, 10**C, label=r"$C$",linewidth=1.0)
# ax.loglog(10**time, 10**O, label=r"$O$",linewidth=1.0)
# ax.loglog(10**time, 10**COp, label=r"$CO+$", linewidth=1.0)
ax.set_xlabel('t (yrs)')
ax.set_ylabel(r'$X_{species}$')
ax2=ax.twinx()
ax2.loglog(mhd_time/(3.154e7),mhd_temp,label="T",linestyle='--',linewidth=0.5)
ax2.set_ylabel('T (K)')
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
legend = ax.legend(h1+h2, l1+l2, loc='best')
# legend.get_frame().set_facecolor('#ffffff')
legend.set_zorder(20)
#ax2.set_xlim([1e2,1e3])
plt.savefig('Tabundances.png', dpi=500, bbox_inches='tight')
plt.close()

# fig,ax=plt.subplots()

# ax.loglog(10**time, 10**HCN, label=r"$HCN$")
# ax.loglog(10**time, 10**HNC, label=r"$HNC$")
# ax.loglog(10**time, 10**CO, label=r"$CS$")
# ax.loglog(10**time, 10**NH3, label=r"$NH_{3}$")
# ax.loglog(10**time, 10**OH, label=r"$OH$")
# ax.loglog(10**time, 10**ELECTR, label=r"$ELECTR$")
# ax.set_xlabel('t (yrs)')
# ax.set_ylabel(r'$X_{species}$')
# ax.legend(loc='best')
# ax2=ax.twinx()
# ax2.loglog(mhd_time/(3.154e7),mhd_n,label="T")
# ax2.set_ylabel(r'n (cm$^{-3}$)')
# ax2.set_xlim([1e1,1e3])
# plt.savefig('nabundances.png', dpi=500,bbox_inches='tight')
# plt.close()
'''
fig,ax=plt.subplots()

ax.loglog(10**time, 10**HCN, label=r"$HCN$",linewidth=1.0)
#ax.loglog(10**time, 10**HNC, label=r"$HNC$",linewidth=1.0)
ax.loglog(10**time, 10**NH3, label=r"$NH_{3}$",linewidth=1.0)
ax.set_xlabel('t (yrs)')
ax.set_ylabel(r'$X_{species}$')
ax2=ax.twinx()
ax2.loglog(mhd_time/(3.154e7),dg,label="d:g",linestyle='--',linewidth=0.5)
ax2.set_ylabel('d:g')
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
legend = ax.legend(h1+h2, l1+l2, loc='best')
# legend.get_frame().set_facecolor('#ffffff')
legend.set_zorder(20)
ax2.set_xlim([1e2,1001])
plt.title('Fixed '+r'$T$'+', varying '+r'$d:g$')
# fig.set_size_inches(8, 6)
plt.savefig('Habundances.jpeg', dpi=500,bbox_inches='tight')
plt.close()


fig,ax=plt.subplots()

ax.loglog(10**time, 10**CO, label=r"$CO$",linewidth=1.0)
ax.loglog(10**time, 10**OH, label=r"$OH$",linewidth=1.0)
ax.loglog(10**time, 10**H2O, label=r"$H_{2}O$",linewidth=1.0)
ax.set_xlabel('t (yrs)')
ax.set_ylabel(r'$X_{species}$')
ax2=ax.twinx()
ax2.loglog(mhd_time/(3.154e7),dg,label="d:g",linestyle='--',linewidth=0.5)
ax2.set_ylabel('d:g')
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
legend = ax.legend(h1+h2, l1+l2, loc='best')
# legend.get_frame().set_facecolor('#ffffff')
legend.set_zorder(20)
ax2.set_xlim([1e2,1001])
plt.title('Fixed '+r'$T$'+', varying '+r'$d:g$')
plt.savefig('Oabundances.png', dpi=500,bbox_inches='tight')
plt.close()
'''

plt.figure(1)
plt.loglog(mhd_time/(3.154e7),mhd_temp)
plt.xlabel('time/yrs')
plt.ylabel('T/K')
plt.savefig('temp.png',dpi=500)
plt.close()
