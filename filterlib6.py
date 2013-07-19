import networkx as nx
import matplotlib.pyplot as plt
import copy
import json
from networkx.readwrite import json_graph
"""
def load():
	graph=nx.Graph()
	graph=nx.read_gpickle("reducedgraph.gpickle")
	f1=open("path.txt","r")
	seq=f1.readline().strip()
	path={}
	path_lenght={}
	while(seq!=""):
		if seq.startswith(">"):
			seq=seq.split()
			couple=(seq[0][1:],seq[1])
		if path.has_key(couple):	
			#seq=seq.split()
			path[couple].append(seq.split())
			path_lenght[couple]=len(seq.split())
		else:
			path[couple]=[]
		seq=f1.readline().strip()
	return graph,path,path_lenght
load()
"""
def check(graph,path_lenght,removable,protlist,path):
	rem=[]
	ess=[]

	for i in removable:
		count=0
		flag=0
		rem.append(i)
		H=graph.copy()
		H.remove_nodes_from(rem)
		for j in path:
			try:
				lenght=nx.shortest_path_length(H, j[0], j[1])
			except:
				lenght=-1
			
			if lenght==-1 or (lenght+1)!=path_lenght[j]:
				ess.append(i)
				flag=1
				break
			else:
				count=count+1
		if count==len(path):
			rem.append(i)
		elif flag==1:
			rem.remove(i)
	
	graph.add_nodes_from(protlist)
	graph.remove_nodes_from(rem)
	print nx.info(graph)
	#drawgraph(graph,protlist)	
	return graph
	
def drawgraph(graph,fixed):
	pos=nx.spring_layout(graph)
	nx.draw_networkx_nodes(graph,pos,with_labels=True,node_size=100)
	nx.draw_networkx_nodes(graph,pos,with_labels=True,nodelist=fixed,node_size=100,node_color="yellow")
	nx.draw_networkx_edges(graph,pos,with_labels=True,width=0.3)
	nx.draw_networkx_labels(graph,pos,fontsize=10)
	plt.show()
	
def computemaxbw(graph,path,protlist,path_lenght,alone):
	print "------Starting Graph------"
	print nx.info(graph)
	betweennees=nx.betweenness_centrality(graph,normalized=True)
	bwsort = sorted(betweennees, key=betweennees.get)
	removable=set(bwsort).difference(set(protlist))
	graphred=check(graph,path_lenght,list(removable),protlist,path)
	nx.write_gpickle(graphred,"bwfilter.gpickle")
	f1=open("bwproteins.txt","w")
	for i in graphred.nodes():
		if i in alone:
			pass
		else:	
			f1.write(i+"\n")
	#nx.write_gpickle(graphred,"intactgraphremaxbw.gpickle")

def computemaxdegree(graph,path,protlist,path_lenght,alone):
	print "------Starting Graph------"
	print nx.info(graph)
	removable=set(graph.nodes()).difference(set(protlist))
	degree=nx.degree(graph,removable)
	removable = list(sorted(degree, key=degree.get))
	graphred=check(graph,path_lenght,list(removable),protlist,path)
	nx.write_gpickle(graphred,"degreefiltered.gpickle")
	f1=open("degreeproteins.txt","w")
	for i in graphred.nodes():
		if i in alone:
			pass
		else:	
			f1.write(i+"\n")
def computemaxweight(graph,path,protlist,path_lenght,alone):
	elements=[]
	nodes=[]
	ess=[]
	print "------Starting Graph------"
	print nx.info(graph)
	
	for i in path:
		max=0
		for j in path[i]:
			count=0
			for k in range(0,len(j)-1,1):
				count=count+float(graph.edge[j[k]][j[k+1]]["weight"])
			if count>max:
				max=count
				elements=j
		
		ess.extend(elements[1:len(elements)-1])
	ess=list(set(ess))
	H=graph.subgraph(ess+protlist)
	#H.add_nodes_from(protlist)
	graphred=check(H,path_lenght,ess,protlist,path)
	nx.write_gpickle(graphred,"weightmaxfilter.gpickle")
	f1=open("weightproteins.txt","w")
	for i in graphred.nodes():
		if i in alone:
			pass
		else:	
			f1.write(i+"\n")
def computebwmod(graph,path,protlist,path_lenght,alone):
	import itertools 
	fixed=[]
	combination=[]
	betweennees=nx.betweenness_centrality(graph,normalized=True)
	print "------Starting Graph------"
	print nx.info(graph)

	d = json_graph.node_link_data(graph) # node-link format to serialize
# write json
	json.dump(d, open('mcn.json','w'))
	count={}
	for i in graph.nodes():
		if i in protlist:
			continue
		else:
			count[i]=0
	combination=list(itertools.combinations(path,2))
	pa={}
	for i in path:
		pa[i]=[]
		for j in path[i]:
			
			if path.has_key(i):
				pa[i].extend(j[1:len(j)-1])
			else:
				pa[i].extend(j[1:len(j)-1])

	for i in path:
		pa[i]=list(set(sum(path[i],[])))


	for i in pa:
		for j in list(set(pa[i])):
			if j in protlist:
				continue
			else:
				count[j]=count[j]+1
	countsort = sorted(count, key=count.get)
	removable=set(countsort).difference(set(protlist))
	print len(protlist)
	graphred=check(graph,path_lenght,removable,protlist,path)
	for i in graphred.nodes():
	
		if i in protlist:
			graphred.node[i]["group"]=5
		else:
			graphred.node[i]["group"]=10
	f1=open("bwmodproteins.txt","w")
	for i in graphred.nodes():	
			f1.write(i+"\n")
	d = json_graph.node_link_data(graphred) # node-link format to serialize
# write json
	json.dump(d, open('filteredgraph.json','w'))
	nx.write_gpickle(graphred,"bwmodfiltgraph.gpickle")