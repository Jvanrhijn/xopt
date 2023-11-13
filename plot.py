import numpy as np
import matplotlib.pyplot as plt
import matplotlib
font = {'family' : 'normal',
        'size'   : 14}
matplotlib.rc('font', **font)




def read_output(path):
    with open(path, "r", encoding="utf-16") as output:
        ls = output.readlines()
        lsplit = map(lambda l: l.split(), ls)
        lsplit = map(lambda l: (int(l[0]), float(l[1])), lsplit)
        ns, lens = list(map(list, zip(*lsplit)))
        return ns, lens


fig, ax = plt.subplots(nrows=1, ncols=2)

for f, m in zip(["xopt", "yopt", "2opt"], ['o', 'D', '*']):
    ns, lens = read_output("out_" + f)
    ns = np.array(ns)
    lens = np.array(lens)

    ax[0].plot(ns, lens, marker=m, color="black", label=f)
    ax[1].plot(ns, lens/(np.sqrt(ns)), marker=m, color="black", label=f)


ax[0].set_ylabel(r"Average tour length")
ax[1].set_ylabel(r"Average tour length / $\sqrt{n}$")

for a in ax:
    a.set_xlabel("n")
    a.legend()
    a.set_ylim(0)

plt.show()

