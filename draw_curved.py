import matplotlib.pyplot as plt
import numpy as np


def get_label(a):
    return str(round(a[1], 3)) + " * x + " + str(round(a[0], 3))


file = open("len_time_graph_1", "r")
l, t = [], []
for line in file:
    l1, t1 = map(float, line.split())
    t.append(t1)
file.close()


file = open("len_curved_1", "r")
for line in file:
    l1 = float(line)
    l.append(l1 * 8)
file.close()

plt.figure(figsize=(12, 7))

# a = np.polyfit(t, l, 1)
# plt.plot(t, l, color="blue", marker=".", markersize=3.0, linestyle="", label=get_label(a))
plt.plot(t, l, color="blue", marker=".", markersize=3.0, linestyle="")

plt.legend()
plt.ylabel('Average length of curved protofilaments in microtubule, nm')
plt.xlabel('Time')
plt.savefig('curved_len(time)_microtubes.png')
plt.show()
