import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Read in the data files
t_n,ten = np.loadtxt('data/ten.xq',unpack=True) # neutral temp
t,tei = np.loadtxt('data/tei.xq',skiprows=3,unpack=True) # ion temp
t,tee = np.loadtxt('data/tee.xq',skiprows=3,unpack=True) # electron temp
t,rhon = np.loadtxt('data/rhon.xq',unpack=True) # neutral density
t,rhoi = np.loadtxt('data/rhoi.xq',skiprows=3,unpack=True) # ion density
t,rhoe = np.loadtxt('data/rhoe.xq',skiprows=3,unpack=True) # electron density
t,rhog= np.loadtxt('data/rhog.xq',skiprows=3,unpack=True) # grain density
t_dg,dg = np.loadtxt('data/dg.xq', skiprows=3, unpack=True)  # grain density

# Find location to plot from in second subplot
t_lim = np.where(t >= 5e9)[0][0]

# Convert mass density to number density
non=rhon/(2*1.67e-24)
noi=rhoi/(2*1.67e-24)
noe=rhoe/(2*1.67e-24)
nog=rhog/(2*1.67e-24)

# Plot the temperature
fig1=plt.figure(1)
ax1=fig1.add_subplot(211)
ax1.loglog(t_n,ten,label=r"$T_{neutral}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$T$'+' '+r'$(K)$')
plt.legend(loc='best')
plt.title('Temperature variations through the Wardle Instability')
plt.tight_layout()

ax2=fig1.add_subplot(212)
ax2.loglog(t_n[t_lim:],ten[t_lim:],label=r"$T_{neutral}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$T$'+' '+r'$(K)$')
plt.legend(loc='best')
plt.tight_layout()

plt.savefig('plots/T-neutral.png')
plt.close(1)


fig1=plt.figure(2)
ax1=fig1.add_subplot(211)
ax1.loglog(t,tei,label=r"$T_{ion}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$T$'+' '+r'$(K)$')
plt.legend(loc='best')
plt.title('Temperature variations through the Wardle Instability')
plt.tight_layout()

ax2=fig1.add_subplot(212)
ax2.loglog(t[t_lim:],tei[t_lim:],label=r"$T_{ion}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$T$'+' '+r'$(K)$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('plots/T-ion.png')
plt.close(2)


fig1=plt.figure(3)
ax1=fig1.add_subplot(211)
ax1.loglog(t,tee,label=r"$T_{electron}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$T$'+' '+r'$(K)$')
plt.legend(loc='best')
plt.title('Temperature variations through the Wardle Instability')
plt.tight_layout()

ax2=fig1.add_subplot(212)
ax2.loglog(t[t_lim:],tee[t_lim:],label=r"$T_{electron}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$T$'+' '+r'$(K)$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('plots/T-electron.png')
plt.close(3)


# Plot the density
fig1=plt.figure(4)
ax1=fig1.add_subplot(211)
ax1.loglog(t,non,label=r"$n_{neutral}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$n$'+' '+r'$(cm^{-3})$')
plt.legend(loc='best')
plt.title('Density variations through the Wardle Instability')
plt.tight_layout()

ax2=fig1.add_subplot(212)
ax2.loglog(t[t_lim:], non[t_lim:], label = r"$n_{neutral}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$n$'+' '+r'$(cm^{-3})$')

plt.legend(loc='best')
plt.tight_layout()
plt.savefig('plots/n-neutral.png')
plt.close(4)


fig1=plt.figure(5)
ax1=fig1.add_subplot(211)
ax1.loglog(t,rhoi,label=r"$\rho_{ion}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$\rho$'+' '+r'$(g\ cm^{-3})$')
plt.legend(loc='best')
plt.title('Density variations through the Wardle Instability')
plt.tight_layout()

ax2=fig1.add_subplot(212)
ax2.loglog(t[t_lim:],rhoi[t_lim:],label=r"$\rho_{ion}$")
plt.ylabel(r'$\rho$'+' '+r'$(g\ cm^{-3})$')
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.legend(loc='best')

ax3 = ax2.twinx()
ax3.loglog(t[t_lim:],noi[t_lim:])
plt.ylabel(r'$n$'+' '+r'$(cm^{-3})$',rotation=270)
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('plots/rho-ion.png')
plt.close(5)


fig1=plt.figure(6)
ax1=fig1.add_subplot(211)
ax1.loglog(t,rhoe,label=r"$\rho_{electron}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$\rho$'+' '+r'$(g\ cm^{-3})$')
plt.legend(loc='best')
plt.title('Density variations through the Wardle Instability')
plt.tight_layout()

ax2=fig1.add_subplot(212)
ax2.loglog(t[t_lim:],rhoe[t_lim:],label=r"$\rho_{electron}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$\rho$'+' '+r'$(g\ cm^{-3})$')
plt.legend(loc='best')

ax3 = ax2.twinx()
ax3.loglog(t[t_lim:],noe[t_lim:])
plt.ylabel(r'$n$'+' '+r'$(cm^{-3})$',rotation=270)
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('plots/rho-electron.png')
plt.close(6)


fig1=plt.figure(7)
ax1=fig1.add_subplot(211)
ax1.loglog(t,rhog,label=r"$\rho_{grain}$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$\rho$'+' '+r'$(g\ cm^{-3})$')
plt.legend(loc='best')
plt.title('Density variations through the Wardle Instability')
plt.tight_layout()

ax2=fig1.add_subplot(212)
ax2.loglog(t[t_lim:],rhog[t_lim:],label=r"$\rho_{grain}$")
plt.ylabel(r'$\rho$'+' '+r'$(g\ cm^{-3})$')
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.legend(loc='best')

ax3 = ax2.twinx()
ax3.loglog(t[t_lim:],nog[t_lim:])
plt.ylabel(r'$n$'+' '+r'$(cm^{-3})$',rotation=270)
plt.xlabel(r'$t$'+' '+r'$(s)$')

plt.tight_layout()
plt.savefig('plots/rho-grain.png')
plt.close(7)


fig1 = plt.figure(8)
ax1 = fig1.add_subplot(211)
ax1.loglog(t_dg, dg, label=r"$d:g$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$d:g$')
plt.legend(loc='best')
plt.title('d:g variations through the Wardle Instability')
plt.tight_layout()

ax2 = fig1.add_subplot(212)
ax2.loglog(t_dg[t_lim:], dg[t_lim:], label=r"$d:g$")
plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.ylabel(r'$d:g$')

plt.legend(loc='best')
plt.tight_layout()
plt.savefig('plots/dg.png')
plt.close(4)


fig1 = plt.figure(9)
ax1 = fig1.add_subplot(111)
ax1.loglog(t_n/(60*60*24*365), ten, 'b-', label="T")
plt.xlabel(r'$t$'+' '+r'$(yrs)$')
plt.ylabel(r'$T$'+' '+r'$(K)$', color='b')
ax1.tick_params('y', colors='b')
plt.xlim([1e2,1e3])
plt.legend(loc='upper left')

ax2 = ax1.twinx()
ax2.loglog(t_dg/(60*60*24*365), dg, 'r--', label="d:g")
plt.ylabel(r'$d:g$', color='r')
plt.xlim([1e2, 1e3])
ax2.tick_params('y', colors='r')
# plt.xlabel(r'$t$'+' '+r'$(yrs)$')

# ax2.loglog(t[t_lim:], rhog[t_lim:], label=r"$\rho_{grain}$")
# plt.ylabel(r'$\rho$'+' '+r'$(g\ cm^{-3})$')
# plt.xlabel(r'$t$'+' '+r'$(s)$')
plt.legend(loc='upper right')

plt.tight_layout()
plt.savefig('plots/combined.png')
plt.close(7)
