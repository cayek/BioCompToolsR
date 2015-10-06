#!/home/cayek/anaconda/bin/python

import sys
import numpy
import os

if (len(sys.argv) != 5) :
    sys.exit(" Use : msTomat.py nbIndiv nbLoci Nm outputFile")
   
   
#ms dir 
ms_dir = "/home/cayek/Projects/UpdateTess/tools/msdir/"

#parameters
nbIndiv = int(sys.argv[1])
nbLoci = int(sys.argv[2])
Nm = float(sys.argv[3])
outputFile = sys.argv[4]

print 'call ms : ' + ms_dir + "ms " \
+ str(nbIndiv) + " " + str(nbLoci) + " -s 1 -I 2 " \
+ str(nbIndiv / 2) + " " + str(nbIndiv / 2) + " -ma x " \
+ str(Nm) + " " + str(Nm) + " x > aux.txt"

os.system(ms_dir + "ms " \
+ str(nbIndiv) + " " + str(nbLoci) + " -s 1 -I 2 " \
+ str(nbIndiv / 2) + " " + str(nbIndiv / 2) + " -ma x " \
+ str(Nm) + " " + str(Nm) + " x > aux.txt")


with open( "aux.txt" , 'r') as file:
    
    #2 line useless
   file.next()
   file.next()

   loci = numpy.empty( (nbIndiv, nbLoci), dtype=int )
    
   for p in range(nbLoci) :
       #skip 3 line
       file.next()
       file.next()
       file.next()
       file.next()
       
       for i in range(nbIndiv) :
           loci[i,p] =int(file.next())
       
   print "export data into " + outputFile
   numpy.savetxt(outputFile, loci, fmt='%i')
   
   
#clean aux file
os.system("rm -f aux.txt seedms")

