#!/usr/bin/env python3

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

topology = "SmallWorld"
topology_mat = np.loadtxt("../../data/raw/5nets_networks_N/N25/"+topology+"_undir")
fluc="0.99"
g = nx.Graph(topology_mat)
pos=nx.spring_layout(g)
variances = [] # variance of Ks across space at each time
flucs = ["0.5","0.7","0.9","0.99"]
fig2,ax2=plt.subplots(figsize=[12,4])
for fluc in flucs:
	Kfile = np.loadtxt("../../data/interim/K_outputs/Ksp5_Complete_undir_N25e25_hetero0.8_fluc"+fluc+"_std1.0_d0.5_r1")
	time_steps = [0,500,1000,1500,2000]
	n_graphs = len(time_steps)
	c = 0
	fig1,axs=plt.subplots(1,n_graphs,figsize=[16,5])
	for time in time_steps:
		Kvec = Kfile[:,time] # n_com x 1 vector of carrying capacity
		_=nx.set_node_attributes(g,Kvec,"K")
		_=nx.draw(g,node_color=Kvec,pos=pos,cmap=plt.cm.summer,node_size=50,ax=axs[c])
		_=axs[c].set_title("t="+str(time))
		c += 1
	output="../../data/processed/img_"+topology+"_fluc"+fluc+".png"
	_=fig1.savefig(output)

	mins,maxs,means,vars = [],[],[],[]
	x = range(Kfile.shape[1]) # timesteps
	for time in x:
		Kvec = Kfile[:,time]
		mins.append(Kvec.min())
		maxs.append(Kvec.max())
		means.append(Kvec.mean())
		vars.append(Kvec.var())

	_=ax2.plot(x,vars,label=fluc)
	_=ax2.set_ylim([0,1])
_=ax2.legend()
output="../../data/processed/plot_"+topology+"_fluc"+fluc+".png"
_=fig2.savefig(output)
