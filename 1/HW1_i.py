# 3-D KINETIC MONTE CARLO SIMULATION
import numpy
import random
import math
import matplotlib.pyplot as plt


# *************** constants

Pon1=(2.8)*pow(10,-5) 
Poff1=(5.7)*pow(10,-7)
Pdiff1=0.1

# *****************  STEP 1. Create a 3-D lattice

lattice=numpy.zeros((60,60,60))

apaf_cnt=50

#print " initial lattice = ",lattice

# create a molecule array
molecule=[[[0,0,0],'cytoc',None]]   # [x,y,z],flag,complex_no
#print molecule
for j in range(1,100+apaf_cnt):
	if j<100:
		molecule.append([[0,0,0],'cytoc',None])
	else:
		molecule.append([[0,0,0],'apaf',None])




# ****************** STEP 2. 10 Monte Carlo Runs

for num1 in range(1,11):
	print "********************************  MC RUN =",num1
	filename="Complex_MC_Run_"+str(num1)+".txt"
	fo = open(filename, "w+")



	num_apopt=0
	complex_no=0
 	ans=[]

	lattice[:][:][:]=None
	# placing the cytoc and apaf molecules in the lattice
	initialval=random.sample([[x,y,z] for x in xrange(60) for y in xrange(60) for z in xrange(60)], 100+apaf_cnt)
	
	
	for j in range(0,100+apaf_cnt):
		molecule[j][2]=None
		molecule[j][0]=initialval[j]	
		lattice[initialval[j][0]][initialval[j][1]][initialval[j][2]]=j
		#print "molecule=",molecule
	

	#STEP 3. ****************************  Monte Carlo Steps
	for num2 in range(0,100000000): 
		print " MC Step = ",num2	
		ans.append(num_apopt)

		if(num2%100000==0):
			for item in ans:
				fo.write(str(item))
				#print "num_apopt=",num_apopt
				fo.write("\n")
			ans=[]
		
		
		
		for u in range(0,100+apaf_cnt):
			num3=random.randint(0,100+apaf_cnt-1)
				

			if(molecule[num3][2]!=None): 
				
				r=round(random.uniform(0, 1),6)
				if(r<=Poff1):
					tmp=molecule[num3][2]
					num_apopt=num_apopt-1
					
					for k in range(0,100+apaf_cnt):
						if(molecule[k][2]!=None and molecule[k][2]==tmp):
							molecule[k][2]=None

				else:
					rnew=round(random.uniform(0, 1),6)
					if(rnew<Pdiff1):

						pdown=0.166666
						pup=2*0.166666
						pleft=3*0.166666
						pright=4*0.166666
						pback=5*0.166666
						pfront=1


						x=molecule[num3][0][0]
						y=molecule[num3][0][1]
						z=molecule[num3][0][2]

						tmp=molecule[num3][2]
						
					
						for k in range(0,100+apaf_cnt):   #xnew , ynew, znew are coord of the other molecule in that complex
							if(molecule[k][2]!=None and molecule[k][2]==tmp):
								xnew=molecule[k][0][0]   
								ynew=molecule[k][0][1]
								znew=molecule[k][0][2]
								knew=k
								break;

						

						r=round(random.uniform(0, 1),6)
						if(r<pdown and z!=0 and znew!=0):
							if(math.isnan(lattice[x][y][z-1])==True and math.isnan(lattice[xnew][ynew][znew-1])==True):  
								molecule[num3][0][2]=z-1
								molecule[knew][0][2]=znew-1
								lattice[x][y][z]=None
								lattice[xnew][ynew][znew]=None
								lattice[x][y][z-1]=num3
								lattice[xnew][ynew][znew-1]=knew

						elif(r<pup and z!=59 and znew!=59):
							if(math.isnan(lattice[x][y][z+1])==True and math.isnan(lattice[xnew][ynew][znew+1])==True):
								molecule[num3][0][2]=z+1
								molecule[knew][0][2]=znew+1
								lattice[x][y][z]=None
								lattice[xnew][ynew][znew]=None
								lattice[x][y][z+1]=num3
								lattice[xnew][ynew][znew+1]=knew
							
						elif(r<pleft and x!=0 and xnew!=0):
							if(math.isnan(lattice[x-1][y][z])==True and math.isnan(lattice[xnew-1][ynew][znew])==True):
								molecule[num3][0][0]=x-1
								molecule[knew][0][0]=xnew-1
								lattice[x][y][z]=None
								lattice[xnew][ynew][znew]=None
								lattice[x-1][y][z]=num3
								lattice[xnew-1][ynew][znew]=knew

							
						elif(r<pright and x!=59 and xnew!=59):
							if(math.isnan(lattice[x+1][y][z])==True and math.isnan(lattice[xnew+1][ynew][znew])==True):
								molecule[num3][0][0]=x+1
								molecule[knew][0][0]=xnew+1
								lattice[x][y][z]=None
								lattice[xnew][ynew][znew]=None
								lattice[x+1][y][z]=num3
								lattice[xnew+1][ynew][znew]=knew
							
						elif(r<pback and y!=0 and ynew!=0):
							if(math.isnan(lattice[x][y-1][z])==True and math.isnan(lattice[xnew][ynew-1][znew])==True):
								molecule[num3][0][1]=y-1
								molecule[knew][0][1]=ynew-1
								lattice[x][y][z]=None
								lattice[xnew][ynew][znew]=None
								lattice[x][y-1][z]=num3
								lattice[xnew][ynew-1][znew]=knew
							
						elif(r<pback and y!=59 and ynew!=59):
							if(math.isnan(lattice[x][y+1][z])==True and math.isnan(lattice[xnew][ynew+1][znew])==True):
								molecule[num3][0][1]=y+1
								molecule[knew][0][1]=ynew+1
								lattice[x][y][z]=None
								lattice[xnew][ynew][znew]=None
								lattice[x][y+1][z]=num3
								lattice[xnew][ynew+1][znew]=knew
							



			else: 


				
				pdown=0.166666
				pup=2*0.166666
				pleft=3*0.166666
				pright=4*0.166666
				pback=5*0.166666
				pfront=1


				x=molecule[num3][0][0]
				y=molecule[num3][0][1]
				z=molecule[num3][0][2]

				

				r=round(random.uniform(0, 1),6)
				if(r<pdown and z!=0):
					if(math.isnan(lattice[x][y][z-1])==True):  
						molecule[num3][0][2]=z-1
						lattice[x][y][z]=None
						lattice[x][y][z-1]=num3
					elif(molecule[int(lattice[x][y][z-1])][2]==None):
						if((molecule[num3][1]=='cytoc' and molecule[int(lattice[x][y][z-1])][1]=='apaf')or (molecule[num3][1]=='apaf' and molecule[int(lattice[x][y][z-1])][1]=='cytoc')) :
							r1=round(random.uniform(0, 1),6)
							if(r1>Poff1 and r1<(Poff1+Pon1)):
								num_apopt=num_apopt+1
								complex_no=complex_no+1
								molecule[num3][2]= complex_no
								molecule[int(lattice[x][y][z-1])][2]=complex_no


				elif(r<pup and z!=59):
					if(math.isnan(lattice[x][y][z+1])==True):
						molecule[num3][0][2]=z+1
						lattice[x][y][z]=None
						lattice[x][y][z+1]=num3
					elif(molecule[int(lattice[x][y][z+1])][2]==None):	
						if((molecule[num3][1]=='cytoc' and molecule[int(lattice[x][y][z+1])][1]=='apaf')or (molecule[num3][1]=='apaf' and molecule[int(lattice[x][y][z+1])][1]=='cytoc')) :
							r1=round(random.uniform(0, 1),6)
							if(r1>Poff1 and r1<(Poff1+Pon1)):
								num_apopt=num_apopt+1
								complex_no=complex_no+1
								molecule[num3][2]= complex_no
								molecule[int(lattice[x][y][z+1])][2]=complex_no

				elif(r<pleft and x!=0):
					if(math.isnan(lattice[x-1][y][z])==True):
						molecule[num3][0][0]=x-1
						lattice[x][y][z]=None
						lattice[x-1][y][z]=num3
					elif(molecule[int(lattice[x-1][y][z])][2]==None):
						if((molecule[num3][1]=='cytoc' and molecule[int(lattice[x-1][y][z])][1]=='apaf')or (molecule[num3][1]=='apaf' and molecule[int(lattice[x-1][y][z])][1]=='cytoc')) :
							r1=round(random.uniform(0, 1),6)
							if(r1>Poff1 and r1<(Poff1+Pon1)):
								num_apopt=num_apopt+1
								complex_no=complex_no+1
								molecule[num3][2]= complex_no
								molecule[int(lattice[x-1][y][z])][2]=complex_no
				elif(r<pright and x!=59):
					if(math.isnan(lattice[x+1][y][z])==True):
						molecule[num3][0][0]=x+1
						lattice[x][y][z]=None
						lattice[x+1][y][z]=num3
					elif(molecule[int(lattice[x+1][y][z])][2]==None):
						if((molecule[num3][1]=='cytoc' and molecule[int(lattice[x+1][y][z])][1]=='apaf')or (molecule[num3][1]=='apaf' and molecule[int(lattice[x+1][y][z])][1]=='cytoc')) :
							r1=round(random.uniform(0, 1),6)
							if(r1>Poff1 and r1<(Poff1+Pon1)):
								num_apopt=num_apopt+1
								complex_no=complex_no+1
								molecule[num3][2]= complex_no
								molecule[int(lattice[x+1][y][z])][2]=complex_no
				elif(r<pback and y!=0):
					if(math.isnan(lattice[x][y-1][z])==True):
						molecule[num3][0][1]=y-1
						lattice[x][y][z]=None
						lattice[x][y-1][z]=num3
					elif(molecule[int(lattice[x][y-1][z])][2]==None):
						if((molecule[num3][1]=='cytoc' and molecule[int(lattice[x][y-1][z])][1]=='apaf')or (molecule[num3][1]=='apaf' and molecule[int(lattice[x][y-1][z])][1]=='cytoc')) :
							r1=round(random.uniform(0, 1),6)
							if(r1>Poff1 and r1<(Poff1+Pon1)):
								num_apopt=num_apopt+1
								complex_no=complex_no+1
								molecule[num3][2]= complex_no
								molecule[int(lattice[x][y-1][z])][2]=complex_no
				elif(r<pback and y!=59):
					if(math.isnan(lattice[x][y+1][z])==True):
						molecule[num3][0][1]=y+1
						lattice[x][y][z]=None
						lattice[x][y+1][z]=num3
					elif(molecule[int(lattice[x][y+1][z])][2]==None):
						if((molecule[num3][1]=='cytoc' and molecule[int(lattice[x][y+1][z])][1]=='apaf')or (molecule[num3][1]=='apaf' and molecule[int(lattice[x][y+1][z])][1]=='cytoc')) :
							r1=round(random.uniform(0, 1),6)
							if(r1>Poff1 and r1<(Poff1+Pon1)):
								num_apopt=num_apopt+1
								complex_no=complex_no+1
								molecule[num3][2]= complex_no
								molecule[int(lattice[x][y+1][z])][2]=complex_no

	#colors = ['red', 'blue', 'green','cyan','magenta','yellow','black','red', 'blue', 'green']


	#plt.plot(range(0,100000),ans,colors[num1-1])
	
	fo.close()	
#plt.ylim([0,100])
#plt.show()	
	





