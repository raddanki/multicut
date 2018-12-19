from collections import defaultdict 
from Graph import Graph
from scipy.optimize import linprog
import numpy as np
import random
import math

"""
Give any custom input graph for testing. 

Output should have
1. Graph type object {Refer Graph.py file for class definition}
2. List of s-t pairs 

"""
def input_graph():

	n = 5
	G = Graph(n)

	G.Edges = [[0,1,1,1,1],[1,0,1,1,1],[1,1,0,1,1],[1,1,1,0,1],[1,1,1,1,0]]
	G.Wt = [[0,1,1,1,1],[1,0,1,1,1],[1,1,0,1,1],[1,1,1,0,1],[1,1,1,1,0]]
	
	SS = [[0,1],[0,3],[2,3]]

	return G, SS

"""
Generate a random graph with each edge being included with probability p.

Weights in the graph are randomly initialized with integers from [1, maxwt].

s-t pairs are randomly chosen among the vertices in the graph.

"""
def random_graph(n=7, k = 3, p=0.5, maxwt=100):


	# random weighted graph
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


	# random s-t pairs for multi-cut
	SS = [[-1,-1] for i in range(0,k)]

	for i in range(0, k):
		
		si = random.randint(0, n-1)
		ti = random.randint(0, n-1)

		while si == ti:
			ti = random.randint(0,n-1) 

		SS[i][0],SS[i][1] = si,ti

	return G, SS

"""
Input:  Graph G and a set of source-sink pairs

Write an LP for multicut by generating all possible paths between every source-sink

Solve LP using lp solver

Output: Fractional values for edges in LP solution 

"""
def lp(G, SS):

	A_ub = []
	b_ub = []
	c = []
	Elist = []

	# instead of creating a $n^2$ length vector, we will only create $|E| = m$ length vector
	for i in range(0, G.n):
		for j in range(i+1, G.n):
			if G.Edges[i][j] == 1: 
				c.append(G.Wt[i][j])
				Elist.append([i,j])


	#generating constraints for lp by considering all paths for every source-sink pair
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

	#solve lp
	res = linprog(c, A_ub=A_ub, b_ub=b_ub,bounds=(0,None),method='interior-point')
	
	#read the fractional values of lp for every edge in G
	E = [[0 for i in range(0, G.n)] for j in range(0, G.n)]
	for i in range(0, len(c)):
		x1 = int(Elist[i][0])
		x2 = int(Elist[i][1])
		
		E[x1][x2] = E[x2][x1] = res.x[i]

	return E,res
	
#remove edges of F from E and E_lp by setting the values to 0
def removeEdges(F, E, E_lp):
	for (i,j) in F:
		E[i][j] = 0
		E[j][i] = 0
		E_lp[i][j] = 0
		E_lp[j][i] = 0

	return E, E_lp

#find radius using ball growing
def get_radius(G, k, E, W, V_opt, dist):
		
	Br = []

	for i in range(0,G.n):
		if dist[i] < 0.5:
			Br.append(i)

	radius = [dist[i] for i in Br]
	radius.sort()

	for r in radius : 

		V_r = float(V_opt)/k
		c_r = 0.0
		for k in range(0, len(Br)) :
			for l in range(k+1, len(Br)):
				
				u = Br[k]
				v = Br[l]

				if E[u][v] > 0 : 

					V_r = V_r + W[u][v]*E[u][v]

		for u in range(0, len(E)):
			for v in range(0, len(E[u])):
				if E[u][v] > 0 and (u in Br and v not in Br):
					V_r = V_r + W[u][v]*(r-dist[u])
					c_r = c_r + W[u][v]
		
		if r > 0 and c_r <= 2*V_r*math.log(k+1) :
			return np.nextafter(r,0) 

	return np.nextafter(0.5,0)


"""
Input : G, SS

Find a fractional solution by solving LP

Round the fractional values to give a valid cut

Output : Approx-Multicut, LP-objective value

"""
def multicut(G, SS):

	cutwt = 0.0
	E_lp, res = lp(G, SS)
	E = G.Edges.copy()
	W = G.Wt.copy()
	r = 0

	#volume of OPT
	V_opt = res.fun

	#multicut approximation algorithm
	for i in range(0, len(SS)):
		if G.connected(SS[i][0],SS[i][1], E) == True:

			dist = G.shortestpath(SS[i][0],E_lp)

			#r = random.uniform(0,0.5)			
			r = get_radius(G, len(SS), E_lp, W, V_opt , dist)
			F = G.cutEdges( E_lp, dist, r)

			cutwt = cutwt + np.sum([W[i][j] for (i,j) in F])
			E, E_lp = removeEdges(F, E, E_lp)

	return cutwt, res.fun

# G, SS = input_graph()
G, SS = random_graph(n=7, k = 3, p=0.5, maxwt=100)

#cutwt : weight of the multicut from the approximation algorithm
#opt_lb : Objective value of LP relaxation obtained from the lp solver. It serves as a lower bound for OPT cut.
cutwt, opt_lb = multicut(G, SS)

print('approx alg : ',cutwt, '\nlower bound for opt : ',opt_lb)





