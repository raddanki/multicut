from collections import defaultdict 
from Graph import Graph
from scipy.optimize import linprog
import numpy as np
import random

def input_graph():

	n = 5
	G = Graph(n)
	
	
	# G.Edges = [[0,1,0,0,1],[1,0,1,1,1],[0,1,0,1,0],[0,1,1,0,1],[1,1,0,1,0]]
	G.Edges = [[0,1,1,1,1],[1,0,1,1,1],[1,1,0,1,1],[1,1,1,0,1],[1,1,1,1,0]]
	G.Wt = [[0,1,1,1,1],[1,0,1,1,1],[1,1,0,1,1],[1,1,1,0,1],[1,1,1,1,0]]
	
	SS = [[0,1],[0,3],[2,3]]


	return G, SS

def random_graph(n=7, k = 3, p=0.5, maxwt=100):


	G = Graph(n)
	E =[[0 for j in range(0,n)] for i in range(0,n)]
	W =[[0 for j in range(0,n)] for i in range(0,n)]
	 
	for i in range(0,n):
		for j in range(i+1,n):
			if random.random() <= p:
				E[i][j] = E[j][i] = 1
				W[i][j] = W[j][i] = random.randint(1,maxwt)

	G.Edges = E
	G.Wt = W

	SS = [[-1,-1] for i in range(0,k)]

	for i in range(0, k):
		
		si = random.randint(0, n-1)
		ti = random.randint(0, n-1)

		while si == ti:
			ti = random.randint(0,n-1) 

		SS[i][0],SS[i][1] = si,ti

	return G, SS

def lp(G, SS):

	A_ub = []
	b_ub = []
	c = []
	Elist = []


	for i in range(0, G.n):
		for j in range(i+1, G.n):
			if G.Edges[i][j] == 1: 
				c.append(G.Wt[i][j])
				Elist.append([i,j])


	for i in range(0,len(SS)):

		G.Paths = []
		G.generate_path(SS[i][0],SS[i][1],[False for t in range(0,G.n)], [])

		for path in G.Paths:

			b_ub.append(-1)

			ll = len(A_ub)
			A_ub.append([0 for i in range(0, len(c))])

			path_edges = []
			for k in range(0,len(path)-1):
				path_edges.append((path[k],path[k+1]))

			count = 0
			for k in range(0,G.n):
				for l in range(k+1,G.n):
					if G.Edges[k][l] > 0:
						if (k,l) in path_edges or (l,k) in path_edges:
							A_ub[ll][count] = -1
						else :
							A_ub[ll][count] = 0
						
						count = count + 1


	# print(len(A_ub))
	res = linprog(c, A_ub=A_ub, b_ub=b_ub,bounds=(0,None),method='interior-point')
	# print(res)
	
	E = [[0 for i in range(0, G.n)] for j in range(0, G.n)]
	for i in range(0, len(c)):
		x1 = int(Elist[i][0])
		x2 = int(Elist[i][1])
		
		E[x1][x2] = E[x2][x1] = res.x[i]

	return E,res
	
def removeEdges(F, E, E_lp):
	for (i,j) in F:
		E[i][j] = 0
		E[j][i] = 0
		E_lp[i][j] = 0
		E_lp[j][i] = 0

	return E, E_lp

def multicut(G, SS):

	cutwt = 0.0
	E_lp, res = lp(G, SS)
	E = G.Edges.copy()
	W = G.Wt.copy()
	r = 0

	for i in range(0, len(SS)):
		if G.connected(SS[i][0],SS[i][1], E) == True:

			r = random.uniform(0,0.5)
			dist = G.shortestpath(SS[i][0],E_lp)
			F = G.cutEdges( E_lp, dist, r)

			cutwt = cutwt + np.sum([W[i][j] for (i,j) in F])
			E, E_lp = removeEdges(F, E, E_lp)

	return cutwt, res.fun

# G, SS = input_graph()
G, SS = random_graph()


cutwt, opt_lb = multicut(G, SS)
print('approx alg : ',cutwt, '\nlower bound for opt : ',opt_lb)





