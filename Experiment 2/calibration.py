from mossbauer_common import *
counts = np.loadtxt("calibration_100s.TKA")

plt.stairs(counts, np.arange(len(counts)+1), color='k')
plt.xlabel("Channel")
plt.ylabel("Counts")
plt.xlim(left=0)

plt.annotate("6.3 keV", (350, 170))
plt.annotate("14.4 keV", (900, 90))

#plt.plot([6.3,6.3],[0,100])
#plt.plot([14.4,14.4],[0,100])
plt.annotate("Escape Peak", (4000, 60))
plt.annotate("121.9 keV", (6500, 150))


plt.savefig("calibration_plot.png", dpi=300)
plt.show()