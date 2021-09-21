#!/usr/bin/env python

import numpy as np

n_com = 750
n_diffenv = 750
rho_eta = 0.8

if n_diffenv < n_com:
    env = np.linspace(0,1,n_diffenv)
    E = np.zeros(n_com)
    c = 0
    for i in env:
        x = n_com / n_diffenv
        if x.is_integer():
            for j in range(int(x)):
                E[c] = i
                c += 1
        else:
            print("Error: n_com / n_diffenv must give an integer")
            break
else:
    E = np.linspace(0,1,n_com)

window = 0.01

# swap environmental conditions randomly until autocorrelation reaches rho_eta 
trial = True
while trial:
    copy_E = np.array(E)
    for i in range(1000):
        autocorr = np.corrcoef(copy_E[1:],copy_E[:-1])[0,1]
        if autocorr >= rho_eta + window:
            ind1,ind2 = np.random.choice(range(n_com),2,replace=False)
            e1 = copy_E[ind1]
            e2 = copy_E[ind2]
            copy_E[ind1] = e2
            copy_E[ind2] = e1
        elif autocorr < rho_eta:
            print("Error")
            break
        else: # i.e. (rho_eta <= autocorr < rho_eta + window)
            trial = False
            np.savetxt("env_rhoeta{0}_e{1}.txt".format(rho_eta,n_diffenv),copy_E)
