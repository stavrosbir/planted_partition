from math import log, ceil, sqrt
import random as r
# from itertools import chain

def ceiling(a):
	return int(ceil(a))

def sgn(a):
	if a < 0:
		return -1
	elif a == 0:
		return 0
	else:
		return 1

def do_with_probability(p):
	return r.random() < p

def planted_partition_graph(n, a, b):
	graph = [[] for i in xrange(n)]
	for i in xrange(n/2):
		# print_progress(100.0*i/(n/2))
		for j in xrange(i+1, n/2):
			if do_with_probability(a/n):
				graph[i].append(j)
				graph[j].append(i)
			if do_with_probability(a/n):
				ii = i + n/2; jj = j + n/2
				graph[ii].append(jj)
				graph[jj].append(ii)
		for j in xrange(n/2, n):
			if do_with_probability(b/n):
				graph[i].append(j)
				graph[j].append(i)
	# print "\r100%   "
	return graph

def slow_R(n):
	return 2*ceiling(log(log(log(n))))

# def subtract_list(a, b):
# 	return [v for v in a if v not in b]

# def intersection(a, b):
# 	return [v for v in a if v in b]

# def union(a, b):
# 	return list(set(a) | set(b))

# def is_subset(a, b):
# 	for v in a:
# 		if v not in b:
# 			return False
# 	return True

# def remove_vertices(G, vertices, n):
# 	for i in xrange(n):
# 		if i in vertices:
# 			G[i] = []
# 		else:
# 			G[i] = subtract_list(G[i], vertices)

def find_w_star(G, vertices_prime2, vertices_prime, n):
	num = ceiling(sqrt(log(log(n))))
	min_val = n
	for i in vertices_prime2:
		number_of_neighbors_in_prime = sum(1 for v in G[i] if v in vertices_prime)
		val = abs(number_of_neighbors_in_prime - num)
		if val < min_val:
			min_val = val
			min_pos = i
	w_star = min_pos
	S_star = set(G[w_star]) & vertices_prime
	return w_star, S_star

def dijsktra(G, s, n):
	nodes = range(n)
	results = [n]*n
	results[s] = 0

	visited = [False]*n
	visited[s] = True
	queue = [s]

	while queue:
		u = queue.pop(0)
	
		new_weight = results[u]+1

		for v in G[u]:
			if not visited[v]:
				queue.append(v)
				visited[v] = True

				if new_weight < results[v]:
					results[v] = new_weight

	return results

def get_indexes_of_value(l, val):
	return set(i for i in xrange(len(l)) if l[i]==val)

def get_indexes_below_value(l, val):
	return set(i for i in xrange(len(l)) if l[i]<val)

# def matrix_combine4(C, n):
# 	result = []
# 	for l in C:
# 		for i in xrange(n):
# 			result.append(list(chain(l[0][i], l[1][i], l[2][i], l[3][i])))
# 	return result

# Sparse lib
def chain_dict4(a, b, c, d, n):
	e = dict(a)
	for i in b:
		e[n+i] = b[i]
	for i in c:
		e[2*n+i] = c[i]
	for i in d:
		e[3*n+i] = d[i]
	return e

def matrix_combine4_sparse(C, n):
	result = []
	for l in C:
		for i in xrange(n):
			result.append(chain_dict4(l[0][i], l[1][i], l[2][i], l[3][i], n))
	return result

def intersection(l1, l2):
	return (value for value in l1 if value in l2)

def dot_product_sparse(a, b):
	return sum(a[i] * b[i] for i in intersection(a, b))

def matrix_vector_mult_sparse(A, b, n):
	vector = {}
	for i in xrange(n):
			value = dot_product_sparse(A[i], b)
			if value != 0:
				vector[i] = value
	return vector

def matrix_powering_sparse(G, V, S_star, S_u, d, k, n):
	# indicator vector z of S_star
	z = {v : 1 for v in S_star}

	# G induced by V
	G_partial = [[] for i in xrange(n)]
	for i in V:
		G_partial[i] = [j for j in G[i] if j in V]

	# adjajency matrix
	A1, A2 = [], []
	for u in xrange(n):
		row1 = {v : 1-d/n for v in G_partial[u]}
	 	A1.append(row1)
	 	row2 = {v : -d/n for v in xrange(n) if v != u and v not in G[u]}
	 	A2.append(row2)

	# identity matrix
	I = [{i : 1} for i in xrange(n)]

	# diagonal matrix
	D1 = [{i : (-(1-d/n)**2)*len(G_partial[i])} for i in xrange(n)]
	D2 = [{i : (-(1-d/n)**2)*(len(G_partial[i])-1)} for i in xrange(n)]

	# zeros matrix
	Zero = [{} for i in xrange(n)]

	# M goes here
	ID1 = [{i : (-(d/n)**2)*(n-1-len(G_partial[i]))} for i in xrange(n)]
	ID2 = [{i : (-(d/n)**2)*(n-2-len(G_partial[i]))} for i in xrange(n)]

	M = matrix_combine4_sparse([[A1, D2, A2, ID1], [I, Zero, Zero, Zero], [A1, D1, A2, ID2], [Zero, Zero, I, Zero]], n)

	# M_hat goes here
	M_hat = matrix_combine4_sparse([[A1, D1, A2, ID1], [Zero]*4, [Zero]*4, [Zero]*4], n)

	# Q_cal zero goes here
	Q_cal0 = matrix_combine4_sparse([[I, Zero, I, Zero]], n)

	# Finally

	result = matrix_vector_mult_sparse(Q_cal0, z, n)
	for i in xrange(k-1):
		result = matrix_vector_mult_sparse(M, result, n)
	result = matrix_vector_mult_sparse(M_hat, result, n)

	return sum(result[v] for v in intersection(S_u, V) if v in result)

# -----

def random_ksi():
	return 2*r.random()-1

def randspin():
	return 1 if r.random() < 0.5 else -1

def belief_propagation(G, n, a, b):
	# step 0
	R = slow_R(n)
	s = (a-b)/2
	d = (a+b)/2
	delta = (1-d/(s**2))/2
	s_prime = (1-delta)*s
	d_prime = (1-delta)*d
	kappa = 10
	k = ceiling(log(n)-1)

	# step 1
	vertices = set(xrange(n))
	vertices_prime2 = set(r.sample(xrange(n), ceiling(sqrt(n))))
	vertices_prime = vertices - vertices_prime2
	# G_prime = remove_vertices(G, vertices_prime2, n)

	# step 2
	w_star, S_star = find_w_star(G, vertices_prime2, vertices_prime, n)

	# step 3
	distances = {u : dijsktra(G, u, n) for u in vertices_prime}
	S = {u : get_indexes_of_value(distances[u], R) for u in vertices_prime}

	# step 4
	U, V = {}, {}
	for j in xrange(1, ceiling(log(n))):
		num = ceiling(n*delta) - ceiling(sqrt(n))
		U[j] = set(r.sample(vertices_prime, num))
		V[j] = vertices_prime - U[j]

	# step 5
	tau = {}
	for j in xrange(1, ceiling(log(n))):
		for u in vertices_prime:
			tau[(j, u)] = sgn(matrix_powering_sparse(G, V[j], S_star, S[u], d, k, n)+kappa*(s_prime**(k+R+1))*len(S_star)*random_ksi()/d_prime)

	# step 6
	J = {}
	spins = [0]*n
	for u in vertices_prime:
		J[u] = 0
		for j in xrange(1, ceiling(log(n))):
			if not (get_indexes_below_value(distances[u], R) & V[j]) and (S[u] | S_star) <= V[j]:
				J[u] = j
				break
		if J[u] != 0 and tau[(J[u], u)] != 0:
			spins[u] = tau[(J[u], u)]
	# print spins
	for i in xrange(n):
		if spins[i] == 0:
			spins[i] = randspin()

	return spins

def cluster_correlation(spins):
	n = len(spins)
	true_spins = [1.0]*(n/2) + [-1.0]*(n/2)
	correlation = sum([true_spins[i]*spins[i] for i in xrange(n)])/n
	return abs(correlation)

def test():
	n, a, b = 100, 9.0, 1.0
	G = planted_partition_graph(n, a, b)

	# print G, dijsktra(G, 0, n)

	spins = belief_propagation(G, n, a, b)
	print spins

	# print 'Correct clustering: ', (100.0*abs(sum(spins[:(N/2)]))/N+50), '\b%'
	print 'Correlation: ', cluster_correlation(spins)

test()