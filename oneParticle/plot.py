import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("output1.txt")
data2 = np.loadtxt("output2.txt")
xs1 = data1[:, 0::2]
ys1 = data1[:, 1::2]
xs2 = data2[:, 0::2]
ys2 = data2[:, 1::2]

N = len(xs1)

for i in range(N):
    plt.plot(xs1[i], ys1[i], color="red")
    plt.plot(xs2[i], ys2[i], color="blue")
plt.show()
