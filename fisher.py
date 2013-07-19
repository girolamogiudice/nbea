from scipy.stats import fisher_exact
import sys
def load(file):
	component={}
	process={}
	function={}
	reaction={}
	protein=[]
	processcount={}
	functioncount={}
	componentcount={}
	reactioncount={}
	gotermname={}
	reactioname={}
	f1=open("9606reactome.txt","r")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split("|")
		reaction[seq[0]]=seq[1].split()
		seq=f1.readline()
	f1.close()
	
	f1=open("protprocess.txt","r")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split()
		process[seq[0]]=seq[1:]
		seq=f1.readline()
	f1.close()	
	f1=open("protfunction.txt","r")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split()
		function[seq[0]]=seq[1:]
		seq=f1.readline()
	f1.close()
	f1=open("protcomponent.txt","r")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split()
		component[seq[0]]=seq[1:]
		seq=f1.readline()
	f1.close()
	f1=open(file,"r")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip()
		protein.append(seq)
		seq=f1.readline()

	protein=list(set(protein))
	
	f1.close()
	f1=open("processcount.txt","r")
	seq=f1.readline()
	totalprocess=int(seq.strip())
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split()
		processcount[seq[0]]=int(seq[1])
		seq=f1.readline()
	f1.close()
	f1=open("functioncount.txt","r")
	seq=f1.readline()
	totalfunction=int(seq.strip())
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split()
		functioncount[seq[0]]=int(seq[1])
		seq=f1.readline()
	f1.close()
	f1=open("componentcount.txt","r")
	seq=f1.readline()
	totalcomponent=int(seq.strip())
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split()
		componentcount[seq[0]]=int(seq[1])
		seq=f1.readline()
	f1.close()
	f1=open("9606totalreactcount.txt","r")
	
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split("|")
		reactioncount[seq[0]]=int(seq[1])
		seq=f1.readline()
	f1.close()
	totalreaction=len(reaction)
	f1=open("gotermname.txt","r")
	seq=f1.readline()
	while (seq!=""):
		gotermname[seq[:10]]=seq[11:].strip()
		seq=f1.readline()
	f1.close()
	f1=open("9606nameidmap.txt","r")
	seq=f1.readline()
	while (seq!=""):
		seq=seq.strip().split("|")
		reactioname[seq[0]]=seq[1]
		seq=f1.readline()
	f1.close()
	return 	component,process,function,protein,processcount,functioncount,componentcount,gotermname,totalprocess,totalcomponent,totalfunction,reaction,totalreaction,reactioncount,reactioname
def fisher():
	try:
		component,process,function,protein,processcount,functioncount,componentcount,gotermname,totalprocess,totalcomponent,totalfunction,reaction,totalreaction,reactioncount,reactioname=load(sys.argv[1])
	
	except IOError:
		print "cannon open",sys.argv
	if len(sys.argv)!=2:
		print "error"
		print "usage python fisher.py <proteinlist>"
	componentset={}
	functionset={}
	processset={}
	reactionset={}
	compsetcount=[]
	funcsetcount=[]
	procsetcount=[]
	reactsetcount=[]

	for i in protein:
		try:
			for j in reaction[i]:
				
				if reactionset.has_key(j):
					reactionset[j]=reactionset[j]+1
				else:
					reactionset[j]=1
		except:
			reactsetcount.append(i)
			pass
		try:
			for j in component[i]:
				
				if componentset.has_key(j):
					componentset[j]=componentset[j]+1
				else:
					componentset[j]=1
		except:
			compsetcount.append(i)
			pass
		try:
			for j in function[i]:
				if functionset.has_key(j):
					functionset[j]=functionset[j]+1
				else:
					functionset[j]=1
		except:
			funcsetcount.append(i)
			pass
		try:
			for j in list(set(process[i])):
			
				if processset.has_key(j):
					processset[j]=processset[j]+1
				else:
					processset[j]=1
		except:
			procsetcount.append(i)
			pass

	numberofproteins=len(protein)-len(list(set(procsetcount)))
	f2=open("gotermprocess.txt","w")
	goprocessfisher={}
	for i in processset:
		a=processset[i]
		b=numberofproteins-a
		c=processcount[i]-a
		d=totalprocess-a-b-c
		table=[[a,b],[c,d]]
		goprocessfisher["GO:"+i[2:]]=fisher_exact(table)[1]
	goprocessfishersort=sorted(goprocessfisher,key=goprocessfisher.get)
	count=0
	for i in goprocessfishersort:
		"""
		count=count+1
		if goprocessfisher[i]*count<(0.001):
			
		
			try:
				name=gotermname[i]
			except:
				name="unknown"
			f2.write(i+"|"+str(goprocessfisher[i]*count)+"|"+name+"\n")
		"""
		if goprocessfisher[i]<(0.05/totalprocess):
			try:
				name=gotermname[i]
			except:
				name="unknown"
	
			f2.write(i+"|"+str(goprocessfisher[i])+"|"+name+"\n")
		
	f2.close()
	
	numberofproteins=len(protein)-len(list(set(funcsetcount)))
	f2=open("gotermfunction.txt","w")
	gofunctionfisher={}
	for i in functionset:
		a=functionset[i]
		b=numberofproteins-a
		c=functioncount[i]-a
		d=totalfunction-a-b-c
		table=[[a,b],[c,d]]
		gofunctionfisher["GO:"+i[2:]]=fisher_exact(table)[1]
	gofunctionfishersort=sorted(gofunctionfisher,key=gofunctionfisher.get)
	for i in gofunctionfishersort:
		if gofunctionfisher[i]<(0.05/totalfunction):
			try:
				name=gotermname[i]
			except:
				name="unknown"
	
			f2.write(i+"|"+str(gofunctionfisher[i])+"|"+name+"\n")
	f2.close()

	numberofproteins=len(protein)-len(list(set(compsetcount)))
	f2=open("gotermcomponent.txt","w")
	gocomponentfisher={}
	for i in componentset:
		a=componentset[i]
		b=numberofproteins-a
		c=componentcount[i]-a
		d=totalcomponent-a-b-c
		
		table=[[a,b],[c,d]]
		gocomponentfisher["GO:"+i[2:]]=fisher_exact(table)[1]
	gocomponentfishersort=sorted(gocomponentfisher,key=gocomponentfisher.get)
	for i in gocomponentfishersort:
		if gocomponentfisher[i]<(0.05/totalcomponent):
			try:
				name=gotermname[i]
			except:
				name="unknown"
	
			f2.write(i+"|"+str(gocomponentfisher[i])+"|"+name+"\n")
	f2.close()
	
	
	numberofproteins=len(protein)-len(list(set(reactsetcount)))
	f2=open("react.txt","w")
	reactfisher={}
	for i in reactionset:
		a=reactionset[i]
		b=numberofproteins-a
		c=reactioncount[i]-a
		d=totalreaction-a-b-c
		table=[[a,b],[c,d]]
		reactfisher[i]=fisher_exact(table)[1]
	reactfishersort=sorted(reactfisher,key=reactfisher.get)
	count=0
	for i in reactfishersort:
		"""
		count=count+1
		if reactfisher[i]*count<(0.05):
			
		
			try:
				name=reactioname[i]
			except:
				name="unknown"
			f2.write(i+"|"+str(reactfisher[i]*count)+"|"+name+"\n")
		"""
		if reactfisher[i]<(0.05/totalreaction):
			try:
				name=reactioname[i]
			except:
				name="unknown"

			f2.write(i+"|"+str(reactfisher[i])+"|"+name+"\n")
	
	f2.close()
	
fisher()

