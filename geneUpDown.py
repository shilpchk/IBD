#!/usr/bin/env python
import os
import sys
args=sys.argv
herb=open(args[1],'rU')
geneup=open(args[2],'rU')
genedown=open(args[3],'rU')
i=1
for Herb in herb:
   geneup.seek(0)
   for Gene in geneup:
      statement = 'esearch -db pubmed -query "'+Herb+'  [TIAB] AND Inhibited [TIAB] AND '+Gene+' [GENE] AND human [ORGN]" | xtract -pattern ENTREZ_DIRECT -element Count >>'+str(i)+'.txt '
      statement = statement.replace("\n"," ")
      os.system(statement)
   genedown.seek(0)
   for Gene in genedown:
      statement = 'esearch -db pubmed -query "'+Herb+'  [TIAB] AND Activated [TIAB] AND '+Gene+' [GENE] AND human [ORGN]" | xtract -pattern ENTREZ_DIRECT -element Count >>'+str(i)+'.txt '
      statement = statement.replace("\n"," ")
      os.system(statement)

   i+=1
herb.close()
geneup.close()	
genedown.close()	
