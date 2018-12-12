import sys
class Graph:


	def __init__(self, n):
		self.n = n
		
		self.Edges = [[0 for j in range(0,n)] for i in range(0,n)]
		self.Wt = [[1.0 for j in range(0,n)] for i in range(0,n)]

		self.lpG  = [[0 for j in range(0,n)] for i in range(0,n)]
		self.lpwt = [[0.0 for j in range(0,n)] for i in range(0,n)]
		
		self.Paths = []
		

	#add an edge of weight u, v into the graph
	def addEdge(self,u,v, w_uv):
		self.Edges[u].append(v)
		self.wt[u,v] = w_uv
		self.wt[v,u] = w_uv

	#add edges into list from adjacency matrix
	def addEdges(self, E):
		for i in range(0, len(E)):
			for j in range(0, len(E[i])):
				if not E[i][j] == 0:
					self.addEdge(i, j, E[i][j])
					self.lpG[i][j] = 1


	#delete an edge from edge list E
	def delEdge(self,u,v, E):
		E[u].remove(v)


	#check if two nodes are connected in the graph with edge list E. We use BFS.
	def connected(self,u,v,E):

		visited = [False for i in range(0, self.n)]
		queue = []

		queue.append(u)
		visited[u] = True

		found = False
		while len(queue) > 0:
			s = queue.pop(0)
			for t in [ i for i in range(0, len(E[s])) if  not E[s][i] == 0] :
				if visited[t] == False:
					visited[t] = True
					if t == v:
						found = True
						return found
					queue.append(t)


		return found

	#generate all the paths between source "s" and a fixed destination "v". All such paths are collected in self.Paths variable.
	def generate_path(self,u,v, visited, path):


		visited[u] = True
		path.append(u)

		if u == v:
			self.Paths.append(path.copy())
		else:
			for t in [ i for i in range(0, len(self.Edges[u])) if  not self.Edges[u][i] == 0]:
				if visited[t] == False:
					self.generate_path(t, v, visited, path)

		path.pop()
		visited[u] = False


	#return min distance vertex which is not yet included in shortest path tree
	def min_vertex(self, spTree, dist):

		dmin = sys.maxsize
		dv = -1

		for i in range(0, self.n):
			if spTree[i] ==False and dist[i] < dmin:
				dmin = dist[i]
				dv = i 

		return dv

	#find all shortest path distances from a vertex u in graph with edges E using dijkstra's algorithm
	def shortestpath(self, u, E):
		dist = [sys.maxsize for i in range(0, self.n)]
		spTree = [False for i in range(0,self.n)]

		dist[u] = 0

		for i in range(0, self.n):
			v = self.min_vertex(spTree, dist)

			spTree[v] = True

			for j in range(0,len(E[v])):
				if E[v][j] > 0 and spTree[j] == False and dist[j] > dist[v] + E[v][j]  :
					dist[j] = dist[v] + E[v][j]



		return dist

	#find cut set in the graph i.e, all edges in E which are at a shortest distance of "r" from a vertex u.
	def cutEdges(self, E, dist, r):
		Br = []
		for i in range(0,self.n):
			if dist[i] < r:
				Br.append(i)

		F = []

		for i in range(0, self.n):
			for j in range(i+1, self.n):
				if E[i][j] > 0 and ((i in Br and j not in Br) or (i not in Br and j in Br)):
					F.append((i,j))

		return F




