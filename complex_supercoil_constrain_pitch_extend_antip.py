#!/usr/bin/env python
from os import system,popen
import string
from sys import argv
import math
from math import sin,cos,atan,tan
#import xyzMath
from xyzMath import *


def enumerate_combinations(list):
    num = len(list)
    ncomb=1
    counter=[]
    for i in range(num):
       ncomb=ncomb*list[i]
       counter.append(0)
    add = 1

    combs=[]
    for j in range(ncomb):
        
      combs.append(string.join(map(lambda x: str(counter[x]), range(num))))  
      for i in range(num):  
        counter[i]=counter[i]+add
        if counter[i]==list[i]:
           counter[i]=0
        else:
           break

    return(combs)

def Generate_Pdbs(Nresl,num_to_output,wl,Rl,orientation,P0,P1,delta_z,output_file_name,chain_name, chain_order):
 
 deg_to_rad= math.pi/180.
# chain parameters
 
##  phase=[]
##  for i in range(num_chain):
##     phase.append(i*ph)
#    phase.append(0)

 chain_set=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

 chain_num=chain_set[0:num_to_output]
 
 R1=2.26
 w1 = 720./7.  # 102.85  for a 2-layer heptad repeat
# w1 = 1080./11.
# w1 = 1800./18.
#rise per residue d  fixed. this constrains pitch (alpha)
 d=1.51
 z1=0.
 z2=0.
 line='ATOM      7  CA  GLY A   2       5.520   2.352   1.361  1.00 20.00     '
 last=line[54:-1]


 for bg in range(1,2,1):  # shift start for loop closure options

   atom_num=1
   res_num=1
   Res_id=[]
   CA_list=[]
   for iter in range(num_to_output):
    #  iter=int(chain_order[counter])   
       CA_list.append([])
       Res_id.append([])
    #  chain=chain_num[iter]
       chain = chain_name[iter]
       orient = orientation[iter]
       Nres=Nresl[iter]
       if orient == 1:
         res_num=0
       else:
         res_num=Nres+1
       w=wl[iter] ##  force twist of 2nd helix to be same as first  to preserve distances    w=wl[iter]
       R=Rl[iter]
       alpha=math.asin(R*w*deg_to_rad/d)
    #  print P0[iter]
       supercoil_phase=P0[iter]+delta_z[iter]*math.tan(alpha)/(R*deg_to_rad)
 
       if orient == 1:
        for t in range(Nres+2):      ## need two extra residues to guide placement of 1st and last residue 

         t_mod=t
         a0=(w*(t_mod-1)+supercoil_phase)*deg_to_rad
         if orient==1: 
          a1=(w1*(t_mod-1)+P1[iter])*deg_to_rad   # set ref point for phase to be along supercoil radius

         x=R*math.cos(a0) + R1*math.cos(a0) * cos(a1) - R1*cos(alpha)*sin(a0)*sin(a1)
         y=R*math.sin(a0) + R1*sin(a0)*cos(a1) + R1*cos(alpha)*cos(a0)*sin(a1)
         if w==0:
          z=d*t_mod+delta_z[iter]
         else:
          z= R*w*t_mod*deg_to_rad/math.tan(alpha)-R1*sin(alpha)*sin(a1)+delta_z[iter]


    #    pdb_lines.append( (res_num,'ATOM   %4d  CA  GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,x,y,z,last)) )

   #     print iter,res_num
         CA_list[iter].append( (res_num,Vec(x,y,z)) )
         Res_id[iter].append(res_num)
         atom_num=atom_num+1

         res_num=res_num+1
       else:
        for t in range(bg,Nres+2+bg):      ## need two extra residues to guide placement of 1st and last residue 

         t_mod=t
         #t_mod=(t+delta_z[iter]/1.5)   ## this makes helix ends flush (assumes only antiparallel helix has z_offset nonzero
         a0=(w*(t_mod-1)+supercoil_phase)*deg_to_rad
         a1=(w1*t_mod-w1*Nresl[iter]-P1[iter])*deg_to_rad

         x=R*math.cos(a0) + R1*math.cos(a0) * cos(a1) - R1*cos(alpha)*sin(a0)*sin(a1)
         y=R*math.sin(a0) + R1*sin(a0)*cos(a1) + R1*cos(alpha)*cos(a0)*sin(a1)
         if w==0:
          z=d*t_mod+delta_z[iter]
         else:
          z= R*w*t_mod*deg_to_rad/math.tan(alpha)-R1*sin(alpha)*sin(a1)+delta_z[iter]


    #    pdb_lines.append( (res_num,'ATOM   %4d  CA  GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,x,y,z,last)) )

   #     print iter,res_num
         CA_list[iter].append( (res_num,Vec(x,y,z)) )
         Res_id[iter].append(res_num)
         atom_num=atom_num+1

         res_num=res_num-1
    # convert CA trace to full backbone model by superimposing on ideal template
    # by matching 3 consecutive CA atoms
    # set up ideal template
   stub_file=map(string.split,open('ideal.pdb','r').readlines())
   atom=[]
   for line in stub_file:
    atom.append( (Vec(float(line[6]),float(line[7]),float(line[8]))))

   ideal_stub=stub(atom[6],atom[1],atom[11])
##  print atom[6]
##  print atom[8]
##  print Vec.distance(atom[6],atom[8]),'dist'
#now make full backbone pdb
#   print bg
   full_pdb=open('%s%s%s%s'%("Ext_",bg,"_",output_file_name),'w')
   atom_num=1

   res_num=0 
   for counter in range(num_to_output):
      iter=int(chain_order[counter])
      chain=chain_name[iter]
#      print iter
      CA_chain_u=CA_list[iter]
      CA_chain = sorted(CA_chain_u, key = lambda res: res[0])
      for res in range(1,Nres+1):
        res_num=res_num+1
    #    print res, res_num
    #    print res_num
    #    print CA_chain[res][1]
        actual_stub=stub(CA_chain[res][1],CA_chain[res-1][1],CA_chain[res+1][1])
        transform=actual_stub * ~ideal_stub

    #    start_d=Vec.distance(atom[5],atom[6])
    #    end_d = Vec.distance(coords,transform*atom[6])
    #    print start_d,end_d,'dist'
    #    print CA_list[res],'ori',coords

    # N
        coords=transform*atom[5]
        full_pdb.write('ATOM %6d  N   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
        atom_num=atom_num+1

    # CA   (use actual CA from trace rather than superimposed one)
        coords=CA_chain[res][1]
        tcoords=transform*atom[6]
    #    print coords,tcoords,'CA'

        full_pdb.write('ATOM %6d  CA  GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
        atom_num=atom_num+1

    #  NH
        coords=transform*atom[7]
        full_pdb.write('ATOM %6d  H   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
        atom_num=atom_num+1

    #  C
        coords=transform*atom[8]
        full_pdb.write('ATOM %6d  C   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
        atom_num=atom_num+1

    # O
        coords=transform*atom[9]
        full_pdb.write('ATOM %6d  O   GLY %s %3d    %8.3f%8.3f%8.3f%s\n'%(atom_num,chain,res_num,coords.x,coords.y,coords.z,last))
        atom_num=atom_num+1

        start_d=Vec.distance(atom[8],atom[6])
        end_d = Vec.distance(transform*atom[8],transform*atom[6])
    #    print start_d,end_d,'dist C  CA '

 return()

def generate_values(data_range):
    
        inc=(data_range[1]-data_range[0])/max(data_range[2]-1.,1.)
        vals=[]          
        for j in range(data_range[2]):
              vals.append(data_range[0]+inc*j)
        return(vals)
    
def input_params(input):

    output_file_prefix=input[0][0]

    Nres=[]
    w0l=[]
    R0l= []
    P0l =[]
    orientation=[]
    P1l = []
    Zl = []
    num_to_output=int(input[1][0])
    chain_name=[]
    chain_order=[]
    for i in range(num_to_output):
        chain_name.append(input[2][i])
        chain_order.append(input[3][i])
    num_vars=7
    for i in range(num_to_output):
        offset=num_vars*(i)+4
        Nres.append(int(input[0+offset][0]))
        w0l.append( (float(input[1+offset][0]),float(input[1+offset][1]),int(input[1+offset][2])))
        R0l.append( (float(input[2+offset][0]),float(input[2+offset][1]),int(input[2+offset][2])))
        P0l.append( (float(input[3+offset][0]),float(input[3+offset][1]),int(input[3+offset][2])))
        orientation.append( int(input[4+offset][0] ))
        P1l.append( (float(input[5+offset][0]),float(input[5+offset][1]),int(input[5+offset][2])))
        Zl.append( (float(input[6+offset][0]),float(input[6+offset][1]),int(input[6+offset][2])))

        
    print ' ###############    SUPERHELIX PARAMS   ############## \n'
    print ' output file prefix:  %s \n'%(output_file_prefix)
    print ' number of chains:  %s \n'%(num_to_output)
    for i in range(num_to_output):
       print 'HELIX %s PARAMS ############################################### \n'%(i)           
       print ' helix length, orientation, and id:  %s %s %s \n'%(Nres[i],orientation[i],chain_name[i])
       print ' starting twist: %s  ending twist: %s  number of samples: %s \n'%(w0l[i][0],w0l[i][1],w0l[i][2]) 
       print ' starting R0:   %s   ending R0:  %s   number of samples:  %s \n'%(R0l[i][0],R0l[i][1],R0l[i][2])
       print ' starting P0:   %s   ending P0:  %s   number of samples:  %s \n'%(P0l[i][0],P0l[i][1],P0l[i][2])
       print ' starting P1:   %s   ending P1:  %s   number of samples:  %s \n'%(P1l[i][0],P1l[i][1],P1l[i][2])                 
       print ' starting Z:   %s   ending Z:  %s   number of samples:  %s \n'%(Zl[i][0],Zl[i][1],Zl[i][2])                 

    print ' ##########  HELIX ORIENTATION CHAIN CHAIN_ORDER ###############  \n '
    for i in range(num_to_output):
        print '%s %s %s %s  \n'%(i, orientation[i],chain_name[i],chain_order[i])
    # sample p1 evenly, w and R around the input values

    ## z_list=[1.5,2.0,2.5,3.0]
    ## p1_s=125.
    ## p2_s=235.


    number_of_combinations=1
    for i in range(num_to_output):
        number_of_combinations=number_of_combinations*R0l[i][2]*P0l[i][2]*P1l[i][2]*Zl[i][2]
        if w0l[i][2] > 0. : number_of_combinations=number_of_combinations*w0l[i][2]

    print  '  NUMBER OF PDBS TO GENERATE:  %s \n'%number_of_combinations

# mofify below
    combinations=[]
    w0=[]
    R0=[]
    P0=[]
    P1=[]
    Z =[]
    for i in range(num_to_output):
        if w0l[i][2] > 0.:
            w0.append(generate_values(w0l[i]))
        else:
            w0.append(["constrained"])
        R0.append(generate_values(R0l[i]))
        P0.append(generate_values(P0l[i]))
        P1.append(generate_values(P1l[i]))
        Z.append(generate_values(Zl[i]))
                   
    return(Nres,num_to_output,orientation,R0,w0,P0,P1,Z,output_file_prefix,chain_name,chain_order)

####################################################                      

input_file=argv[1]
tag = argv[2]
input =map(string.split,open(input_file,'r').readlines())
Nres,num_to_output,orientation,R0,w0,P0,P1,Z,output_file_prefix,chain_name,chain_order = input_params(input)
#print w0,R0,chain_order
items=[]
for i in range(num_to_output):
  items.append(len(w0[i]))
  items.append(len(R0[i]))
  items.append(len(P0[i]))
  items.append(len(P1[i]))
  items.append(len(Z[i]))

combos=enumerate_combinations(items)

for combo in combos:
#    print combo
    id=map(int,string.split(combo))
#    print id

    w0v=[]
    R0v=[]
    P0v=[]
    P1v=[]
    Zv=[]
    
    for i in range(num_to_output):
        if len(w0[i])==1 and w0[i][0]=="constrained":
         w0v.append(1.51*sin(atan(R0[i][id[1+i*5]]*tan(asin(R0[0][id[1]]*w0[0][id[0]]*math.pi/180./1.51))/R0[0][id[1]]))/R0[i][id[1+i*5]]/math.pi*180.)
        else:
         w0v.append(w0[i][id[0+i*5]])
        R0v.append(R0[i][id[1+i*5]])
        P0v.append(P0[i][id[2+i*5]])
        P1v.append(P1[i][id[3+i*5]])
        Zv.append(Z[i][id[4+i*5]])           

    print w0v
 #   print R0v
 #   print P0v
 #   print P1v
 #   print Zv

        
##     ## to get close to C2 symmetry for axis in x-y plane, switch sign of  phase and delta_z of sym mates
##         ## 
## ##     if num_chain==2:
## ## 	helix_phase.append(orientation[1]*helix_phase[0])
## ##         delta_z.append(orientation[1]*delta_z[0])

## ##     if num_chain==4:   ## warning-this is 4 helix bundle centric
## ## 	helix_phase.append(-helix_phase[1])
## ##         helix_phase.append(-helix_phase[0])
## ##         delta_z.append(-delta_z[1])
## ##         delta_z.append(-delta_z[0])

 
    out_file_name='%s'%(tag)
    for i in range(num_to_output):
       
      out_file_name=out_file_name+'_%.2f_%.2f_%.2f_%.2f_%.2f'%(w0v[i],R0v[i],P0v[i],P1v[i],Zv[i])  ##keep w0 fixed!!


    
    out_file_name=out_file_name+'.pdb'
    Generate_Pdbs(Nres,num_to_output,w0v,R0v,orientation,P0v,P1v,Zv,out_file_name,chain_name,chain_order)
    
## for i in range(3):
##     p1=p1_s+(i-1)*5
##     for ii in range(3):
##      p2=p2_s+(ii-1)*5   
##      for j in range(3):
##         R=R0+(j-1)*0.1
##         for k in range(5):
##             w=w0+(k-2)*0.15
##             for l in range(4):
##                 delta_z=z_list[l]
##                 output_file='fine2_%s_%s_%s_%s_%s_%s'%(w,R,p1,p2,delta_z,output_file_name)
## #    print output_file_name,Nres,w,R,p1
##                 Generate_Pdbs(Nres,w,R,num_chain,p1,p2,delta_z,output_file,num_to_output)
## #print atom_list
