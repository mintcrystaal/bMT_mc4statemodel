import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


# считать наклон между каждыми 2 соседними. Если резко падает, то катастрофа. когда конец катастрофы? найти параметры для диффузии. см раздел quantification ...


file = open("len_time_graph_1", "r")
length, t = [], []
for line in file:
    l1, t1 = map(float, line.split())
    length.append(l1 * 8)
    t.append(t1)
file.close()

minpeakprominence = 20
peaks, properties = find_peaks(length, prominence=minpeakprominence)
print(peaks)
print(properties)
