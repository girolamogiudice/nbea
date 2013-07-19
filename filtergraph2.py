###################################################
#this function load the intact graph and filter it#
#on the basis of degree  or betweenness or intact-#
#miscore.										  #
###################################################
import networkx as nx
import itertools,sys
import filterlib6
#load graph and nodes
def load(graph,nodes):
	f=open(nodes,"r")
	seq=f.readline()
	protlist=[]
	while(seq!=""):
		protlist.append(seq.strip())
		seq=f.readline()
	#graph in gpickle format
	graph=nx.read_gpickle(graph)
	
	print "graph loaded"
	print nx.info(graph)
	print "nodes loaded"
	print "number of nodes loaded: ",len(protlist)
	return graph,protlist

#return all possible path given a couple of nodes
def couple(graph,protlist):
	path={}
	path_lenght={}
	nodes=[]
	graph_reduced=nx.Graph()
	combination=list(itertools.combinations(protlist,2))
	for i in range(len(combination)):
		try:	 
			path[combination[i]]= list(nx.all_shortest_paths(graph,combination[i][0],combination[i][1]))
			path_lenght[combination[i]]=len(path[combination[i]][0])
			nodes.extend(list(sum(path[combination[i]],[])))
		except: 
			pass
	
	nodes=list(set(nodes))
	alone=list(set(protlist).difference(set(nodes)))
	print alone
	for i in protlist:
		graph_reduced.add_node(i)
		graph_reduced.node[i]["group"]=10
	for i in graph_reduced.nodes():
		graph_reduced.node[i]["group"]=20
	
	graph_reduced=graph.subgraph(nodes)
	graph_reduced.add_nodes_from(alone)
	for i in list(set(alone)):
		graph_reduced.node[i]["pathway"]=[]
		graph_reduced.node[i]["gotermcomponent"]=[]
		graph_reduced.node[i]["gotermprocess"]=[]
		graph_reduced.node[i]["gotermfunction"]=[]
	print nx.info(graph_reduced)
	f1=open("subgraphproteins.txt","w")
	for i in graph_reduced.nodes():
		if i in alone:
			pass
		else:
			f1.write(i+"\n")
	return graph_reduced,path,protlist,path_lenght,alone

def main():
	try:
		graph,protlist=load(sys.argv[1],sys.argv[2])
	
	except IOError:
		print "cannon open",sys.argv
    
	
	if len(sys.argv)!=4:
		print "error"
		print "usage python filtergraph.py <graph> <list of nodes> <#filter>"
		print "#1 for degree based filter"
		print "#2 for betweenness based filter"
		print "#3 for intact-miscore filter"
		print "#4 for bwmod filter"
		sys.exit(2)

	graph_reduced,path,protlist,path_lenght,alone=couple(graph,protlist)
	
	if sys.argv[3]=="1":
		nx.write_gpickle(graph_reduced,"subraph.gpickle")
		filterlib6.computemaxdegree(graph_reduced,path,protlist,path_lenght,alone)
	if sys.argv[3]=="2":
		nx.write_gpickle(graph_reduced,"subraph.gpickle")
		filterlib6.computemaxbw(graph_reduced,path,protlist,path_lenght,alone)
	if sys.argv[3]=="3":
		nx.write_gpickle(graph_reduced,"subraph.gpickle")
		filterlib6.computemaxweight(graph_reduced,path,protlist,path_lenght,alone)
	if sys.argv[3]=="4":
		nx.write_gpickle(graph_reduced,"subraph.gpickle")
		filterlib6.computebwmod(graph_reduced,path,protlist,path_lenght,alone)

if __name__ == "__main__":
    main()
