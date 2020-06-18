#!/usr/bin/env python

from optparse import OptionParser
import sys
import os

def main():
  parser = OptionParser()
  parser.add_option("-s",dest="snvFile",default=None)
  parser.add_option("-i",dest="indelFile",default=None)
  (options,args)=parser.parse_args()
  
  snv_file=open(options.snvFile)
  print "Chrom:Pos\tRef\tAlt\tGene\tFunction\tVAF"
  i=1
  for record in snv_file:
    if(i==1):
      i=i+1
      continue
    entry=record.strip().split('\t')
    chrom=entry[0].strip('chr')
    pos=entry[1]
    ref=entry[3]
    alt=entry[4]
    gene=entry[6]
    func=entry[5]
    data=entry[70].split(':')
    geno=entry[69].split(':')
    if(len(geno) != 5):
       print chrom+":"+pos+"\t"+ref+"\t"+alt+"\t"+gene+"\t"+func+"\t"+". ( / )"
    else:
       depth=int(data[2])
       counts=data[1].split(',')
       altD=int(counts[1])
       if(float(depth) != 0):
            freq=round(float(altD)/float(depth),4)
       else:
            freq=0.0
       print chrom+":"+pos+"\t"+ref+"\t"+alt+"\t"+gene+"\t"+func+"\t"+str(freq)+" ("+str(altD)+"/"+str(depth)+")"
  snv_file.close()

  indel_file=open(options.indelFile)
  for record in indel_file:
    if(record=="#"):
      continue
    entry=record.strip().split('\t')
    chrom=entry[0].strip('chr')
    pos=entry[1]
    ref=entry[3]
    alt=entry[4]
    info=entry[7].split('|')
    #pred=info[0].split(';')
    #predVal=pred[0]
    anno=pred[2].split('=')
    gene=anno[1]
    funcInf=info[3].split(',')
    func=funcInf[0]
    ad=entry[9].split(',')
    ref=int(ad[0])
    alt=int(ad[1])
    depth=ref+alt
    freq=round(float(alt)/float(depth),4)

    print chrom+":"+pos+"\t"+ref+"\t"+alt+"\t"+gene+"\t"+func+"\t"+str(freq)+" ("+str(alt)+"/"+str(depth)+")"
    
  indel_file.close()
  
if __name__=="__main__":
  main()
