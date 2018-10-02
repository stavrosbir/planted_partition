import random as r
from math import exp, log, sqrt, ceil
import sys

def print_progress(pr):
	sys.stdout.write("\r%.1f%%" % pr)
	sys.stdout.flush()

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


def randspins(n):
	spins = [1]*(n/2) + [-1]*(n/2)
	r.shuffle(spins)
	return spins

# def randspin():
# 	return 1 if r.random() < 0.5 else -1

# def take_two2(n, spins):
# 	first = int(r.random() * n)
# 	temp = int(r.random() * (n/2))

# 	l = [i for i in xrange(n) if spins[i] != spins[first]]
# 	second = l[temp]
# 	print first, second

# 	# second = 0
# 	# while second < n and temp >= 0:
# 	# 	if spins[second] != spins[first]:
# 	# 		temp -= 1
# 	# 	second += 1

# 	return first, second

# def take_two(n):
# 	first = int(r.random() * n)
# 	second = int(r.random() * (n-1))
# 	return first, second if second < first else (second+1)

# def delta_energy(G, spins, p):
# 	return 2 * spins[p] * sum([spins[j] for j in G[p]])

def eval_clustering(spins):
	n = len(spins)
	return str(100.0*abs(sum(spins[:(n/2)]))/n+50)+'%'

def cluster_correlation(spins):
	n = len(spins)
	true_spins = [1.0]*(n/2) + [-1.0]*(n/2)
	correlation = sum([true_spins[i]*spins[i] for i in xrange(n)])/n
	return abs(correlation)

def ising_on_graph(G):

	n = len(G)
	iterations = 1000*int(ceil(log(n)))
	spins = randspins(n)
	count = [0, n/2, n/2]

	temperature = 5.0

	for i in xrange(iterations):

		for particle in xrange(n):

			if (i*n + particle) % (iterations*n/100) == 0:
				# print temperature
				# print eval_clustering(spins)
				# print_progress(100.0*i/iterations)
				temperature /= 1.0233

			spin = spins[particle]
			energy = 2 * spin * sum([spins[j] for j in G[particle]])
			if count[spin] >= (0.5 - temperature/100.0)*n and r.random() < exp(-energy/temperature):
				count[spin] -= 1
				count[-spin] += 1
				spins[particle] = -spin			

	# print "\r100% "
	return spins

# import matplotlib.pyplot as plt

def scaling(N):
	
	points = 50
	samples = 20
	step = exp(0.022)
	start = 4.0

	b = 1.0

	xs = []
	ys = []

	d = start
	for p in xrange(points):

		sum_y = 0
		
		for s in xrange(samples):
			print_progress(100.0*(p*samples+s)/(points*samples))
			y = cluster_correlation(ising_on_graph(planted_partition_graph(N, b + d, b)))
			sum_y += y

		xs.append(log(d))
		ys.append(log(sum_y/samples))

		d *= step

	print "\r100% "
	print xs, ys

def sampling(N):

	points = 25
	samples = 50
	step = 0.5

	a = 3.0
	b = 1.0

	xs = []
	ys = []

	for p in xrange(points):

		# sum_y = 0
		min_y = 1.0
		
		for s in xrange(samples):
			print_progress(100.0*(p*samples+s)/(points*samples))
			y = cluster_correlation(ising_on_graph(planted_partition_graph(N, a, b)))
			# sum_y += y
			min_y = min(min_y, y)

		xs.append(a)
		# ys.append(sum_y/samples)
		ys.append(min_y)

		a += step

	print "\r100% "
	print xs, ys	

def show_distribution(N, a, b, samples):

	xs = []

	for s in xrange(samples):

		print_progress(100.0*s/samples)
		sample = cluster_correlation(ising_on_graph(planted_partition_graph(N, a, b)))
		xs.append(sample)

	print "\r100% "
	print xs

def variance(a, b):

	points = 101
	samples = 10
	step = 100

	N = 100

	xs = []
	ys = []

	for p in xrange(points):

		zs = []
		
		for s in xrange(samples):
			print_progress(100.0*(p*samples+s)/(points*samples))
			y = cluster_correlation(ising_on_graph(planted_partition_graph(N, a, b)))
			zs.append(y)

		xs.append(N)

		mean = sum(zs)/samples
		sd = sqrt(sum([(z-mean)**2 for z in zs])/(samples-1))
		ys.append(sd)

		N += step

	print "\r100% "
	print xs, ys

def main():

	N = 10000

	print 'Constructing graph...'
	G = planted_partition_graph(N, 5.0, 1.0)

	print 'Simulating ising...'
	spins = ising_on_graph(G)

	print 'Correct clustering: ', (100.0*abs(sum(spins[:(N/2)]))/N+50), '\b%'
	print 'Correlation: ', cluster_correlation(spins)
	# print sum(spins)

main()

# show_distribution(1000, 8.0, 2.0, 1000)
# scaling(1000)
# sampling(2000)
# variance(8.0, 2.0)