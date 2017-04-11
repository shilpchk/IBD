#!/usr/bin/env python
import os
import sys
args=sys.argv
#print type(args)
#print args
herb=open(args[1],'rU')
gene=open(args[2],'rU')

i=1
for Herb in herb:
   gene.seek(0)
   for Gene in gene:
      statement = 'esearch -db pubmed -query "'+Herb+'  [TIAB] AND '+Gene+' [GENE] AND human [ORGN]" | xtract -pattern ENTREZ_DIRECT -element Count >>'+str(i)+'.txt '
      statement = statement.replace("\n"," ")
      os.system(statement)

   i+=1
herb.close()
gene.close()	
