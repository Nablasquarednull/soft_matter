import numpy as np
import matplotlib.pyplot as plt

b = 3.0
N = 100
T = 1
epsilon = 1e-5

x = np.zeros((T,N+1))
y = np.zeros((T,N+1))
z = np.zeros((T,N+1))

filename='simulation_FJC_b=%.1f_N=%d_T=%d.xyz'%(b,N,T)
with open(filename,'r') as f:
    for t in range(T):
        f.readline()
        f.readline()
        for n in range(N+1):
            data = f.readline().split()
            x[t,n] = float(data[1])
            y[t,n] = float(data[2])
            z[t,n] = float(data[3])

# Structure Factor
k = np.arange(0.01, 1, 0.01)
I = np.zeros((len(k), N+1, N+1))

for l in range(len(k)):
    kval = k[l]
    for j in range(N+1):
        R_j = np.array([x[t,j], y[t,j], z[t,j]])
        for i in range(N+1):
            R_i = np.array([x[t,i], y[t,i], z[t,i]])
            diff = np.linalg.norm(R_j - R_i)
            if diff < 1e-12:
                I[l,j,i] = 1.0
            else:
                I[l,j,i] = np.sin(kval * diff) / (kval * diff)

# the structure factor is done, now we plot it and compare with the Guinier approximation

#We first compute the guinier approx
R_g = np.sqrt(N * b**2 / 6)
Guinier = ((N+1)**2) *(1 - ((k *R_g)**2)/3)





S_sum = I.sum(axis=(1,2))
plt.figure()
plt.plot(k, S_sum, 'o-',label = 'simulation')
plt.plot(k,Guinier,label= 'Guinier approx.')
plt.xlim(0.001,0.3)
plt.ylim(bottom=0)
plt.ylim(top=12000)
plt.legend()
plt.plot()
plt.xlabel('k')

plt.ylabel(r'Sum over all i,j of $sin(k|R_i-R_j|)/(k|R_i-R_j|)$')
plt.savefig('Guinier_approx.png')
plt.show()
