import matplotlib.pyplot as plt
import numpy as np


def get_label(a):
    return str(round(a[0], 3)) + " * x + " + str(round(a[1], 3))


file = open("len_each_time_graph_1", "r")
l1, l2, l3, l4, t = [], [], [], [], []
for line in file:
    _l1, _l2, _l3, _l4, t1 = map(float, line.split())
    l1.append(_l1)
    l2.append(_l2)
    l3.append(_l3)
    l4.append(_l4)
    t.append(t1)
file.close()


a1 = np.polyfit(t, l1, 1)
a2 = np.polyfit(t, l2, 1)
a3 = np.polyfit(t, l3, 1)
a4 = np.polyfit(t, l4, 1)

plt.figure(figsize=(12, 7))
plt.plot(t, l1, color="blue", marker=".", markersize=3.0, linestyle="-", label=get_label(a1))
plt.plot(t, l2, color="green", marker=".", markersize=3.0, linestyle="-", label=get_label(a2))
plt.plot(t, l3, color="red", marker=".", markersize=3.0, linestyle="-", label=get_label(a3))
plt.plot(t, l4, color="orange", marker=".", markersize=3.0, linestyle="-", label=get_label(a4))

plt.legend()
plt.title('Length(time)')
plt.ylabel('Length of each protofilament in microtubule')
plt.xlabel('Time')
plt.savefig('len(time)_microtubes.png')
plt.show()
