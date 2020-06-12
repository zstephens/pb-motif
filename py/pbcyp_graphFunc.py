#
#	various functions for analysis of directed graphs
#


def longestPath_DAG(g,s,w):
	"""
	return the longest (highest scoring) path in a directed acyclic graph,
	starting from a particular node

	g - input graph
	s - starting node
	w - graph weights
	returns: list of nodes in top-scoring path
	"""
	reverse_graph = {v[0]:[] for v in g}
	for v in g:
		for v2 in v[1]:
			reverse_graph[v2].append(v[0])
	reverse_graph = {k:list(set(reverse_graph[k])) for k in reverse_graph.keys()}
	start = False
	dist  = {v[0]:0 for v in g}
	dist[s[0]] = 999999999
	w_mat = {v[0]:[] for v in g}
	for v in g[::-1]:
		if start:
			if len(reverse_graph[v[0]]):
				w_mat[v[0]] = sorted([(dist[u]+w[u][v[0]],u) for u in reverse_graph[v[0]]],reverse=True)
				dist[v[0]]  = w_mat[v[0]][0][0]
		if v[0] == s[0]:
			start = True
	maxV = max(dist.values())
	# traceback
	outPath = []
	if maxV > 0:
		maxVal = [v for v in dist.keys() if dist[v] == maxV][0]
		vvv = maxVal
		outPath.append(vvv)
		while True:
			if len(w_mat[vvv]):
				vvv = w_mat[vvv][0][1]
				outPath.append(vvv)
			else:
				break
	return outPath[::-1]


def exhaustive_DAG(g,w,p):
	"""
	find the top-scoring path in a directed acyclic graph, brute force trying
	every possible starting position

	g - input graph
	w - graph weights
	p - prize for each node
	returns: list paths, sorted by score
	"""
	topPaths = []
	for start in g:
		topPaths.append([p[start[0]],[start[0]]])
		dag_path = longestPath_DAG(g,start,w)
		if len(dag_path) >= 2:
			for i in range(len(dag_path)):
				for j in range(i+1,len(dag_path)+1):
					path    = dag_path[i:j]
					myPrize = p[path[0]]
					for k in range(1,len(path)):
						myPrize +=  w[path[k-1]][path[k]]
					topPaths.append([myPrize,[n for n in path]])
	return sorted(topPaths,reverse=True)[0]


def get_connected_subgraphs(A):
	"""
	return a list of all connected subgraphs given an adjacency matrix

	A - adjacency matrix
	returns: list of list of points, clustered
	"""
	clusters_out = []
	not_visited = {i:True for i in range(len(A))}
	for k in not_visited.keys():
		if not_visited[k]:
			newCluster = dfs_adj(A,k)
			clusters_out.append([n for n in newCluster])
			for n in newCluster:
				not_visited[n] = False
	return clusters_out

#
if __name__ == '__main__':
	test_A = [[0,1,0,0,0],
	          [1,0,0,0,0],
	          [0,0,0,1,1],
	          [0,0,1,0,1],
	          [0,0,0,1,0]]
	print(get_connected_subgraphs(test_A))

