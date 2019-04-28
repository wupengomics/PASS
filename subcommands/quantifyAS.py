#
## this program quantify exon inclusion level of Alternative Splicing 
#

### import necessary libraries
import re,os,sys,time,datetime,commands;
import scipy,math;
from scipy import stats;

### checking out the number of arguments
if (len(sys.argv)<3): 
  print('Not enough arguments!!');
  print ('It takes 2 argument');
  print ('Usage:\n\tpython quantifyAS.py SAMfile filefolder');
  print ('Example\n\tpython quantifyAS.py hc.sam ./PASS_out');
  sys.exit();

def listToString(x):
  rVal = '';
  for a in x:
    rVal += a+' ';
  return rVal;


#### global variables  ################
#### default configuration values #####
readLength=50;
junctionLength=84;
#######################################

input_1 = sys.argv[1];
outDir = sys.argv[2];

SE = outDir+'/ASevents.SE.txt';
MXE = outDir+'/ASevents.MXE.txt';
A5SS = outDir+'/ASevents.A5SS.txt';
A3SS = outDir+'/ASevents.A3SS.txt';
RI = outDir+'/ASevents.RI.txt';
AFE = outDir+'/ASevents.AFE.txt';
ALE = outDir+'/ASevents.ALE.txt';


#
commands.getstatusoutput('mkdir '+outDir);
#



##### Getting Start Time ######
startTime = time.time();

sample_1 = input_1.split(',');
#
numRep_1 = len(sample_1);
#
ejLength = junctionLength-readLength+1; ## effective junction length
#
#
CT1,CT2=0,1; ### count type 1,2
S1=0; ## sample 1
I,S=0,1; ## inclusion isoform or skipping form
#
#
######### more functions #############
#



def getInitialCounts(): ## getting initial counts for each AS event
  rValue = [[[[],[]],[[],[]]],[[[],[]],[[],[]]],[[[],[]],[[],[]]]]; ## count type 1, 2, and 3
  for i in range(0,numRep_1): ## for sample 1
    rValue[CT1][S1][I].append(0);
    rValue[CT1][S1][S].append(0);
    rValue[CT2][S1][I].append(0);
    rValue[CT2][S1][S].append(0);
  return rValue;
#### end of getInitialCounts()

#
### open AS event files..
#
seFile = open(SE); ## skipped exon event file
mxeFile = open(MXE); ## mxe event file
a5ssFile = open(A5SS); ## a5ss event file
a3ssFile = open(A3SS); ## a3ss event file
afeFile = open(AFE); ## afe event file
aleFile = open(ALE); ## ale event file
riFile = open(RI); ## ri event file
#
### open output files here...
#
JCEC_seFile = open(outDir+'/AS.'+'SE.txt', 'w');
JCEC_mxeFile = open(outDir+'/AS.'+'MXE.txt', 'w');
JCEC_a5ssFile = open(outDir+'/AS.'+'A5SS.txt', 'w');
JCEC_a3ssFile = open(outDir+'/AS.'+'A3SS.txt', 'w');
JCEC_afeFile = open(outDir+'/AS.'+'AFE.txt', 'w');
JCEC_aleFile = open(outDir+'/AS.'+'ALE.txt', 'w');
JCEC_riFile = open(outDir+'/AS.'+'RI.txt', 'w');
#
chunk=1000; ## to speed up the sam file processing
#
se={};mxe={};a5ss={};a3ss={};afe={};ale={};ri={}; ## 7 dictionaries
e_se={};e_mxe={};e_a5ss={};e_a3ss={};e_afe={};e_ale={};e_ri={}; ## exons dictionaries
c_se={};c_mxe={};c_a5ss={};c_a3ss={};c_afe={};c_ale={};c_ri={}; ## count dictionaries
#
#### SE ######
#
c=0;     ## count
numSE=0; ## number of SE
numSEDup=0; ## duplicate SE id 
line=seFile.readline(); ## skipping header
for line in seFile: ## process skipped exon events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3]; 
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4];
  tS = int(ele[5]); tE = int(ele[6]); ## target exon coord
  uS = int(ele[7]); uE = int(ele[8]); ## upstream exon coord
  dS = int(ele[9]); dE = int(ele[10]); ## downstream exon coord


  e_se[id] = [tS,tE,uS,uE,dS,dE];
  c_se[id] = getInitialCounts(); 

  group = range(uS/chunk, uE/chunk+1) +  range(tS/chunk, tE/chunk+1) +  range(dS/chunk, dE/chunk+1); ## groups this event could belong 
  group = list(set(group));  ## remove duplicate groups

  if chr in se: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in se[chr]: ## this group is already there
        if id in se[chr][i]: ## not likely but this group already has the id
          numSEDup+=1;
        else: ## new SE ID
          se[chr][i][id] = [tS,tE,uS,uE,dS,dE]; ## skipping event with coords
          numSE+=1;
      else: ## new group to this chromosome
        se[chr][i]={};
        se[chr][i][id] = [tS,tE,uS,uE,dS,dE]; ## skipping event with coords
        numSE+=1;
  else: ## first time accesing this chromosome
    se[chr]={};
    for i in group: ## for each possible group
      se[chr][i]={};
      se[chr][i][id] = [tS,tE,uS,uE,dS,dE]; ## skipping event with coords
      numSE+=1;
#
#### MXE ####
#
c=0;     ## count
numMXE=0; ## number of MXE
numMXEDup=0; ## duplicate MXE id
line=mxeFile.readline(); ## mxe header
for line in mxeFile: ## process mxe events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3];
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4]; ## '+' or '-'
  tS = int(ele[5]); tE = int(ele[6]); ## target exon coord
  sS = int(ele[7]); sE = int(ele[8]); ## second exon coord
  if strand=='-': ## negative strand, switch target exon and second exon
    tS = int(ele[7]); tE = int(ele[8]); ## target exon coord
    sS = int(ele[5]); sE = int(ele[6]); ## second exon coord    
  uS = int(ele[9]); uE = int(ele[10]); ## upstream exon coord (samller coord)
  dS = int(ele[11]); dE = int(ele[12]); ## downstream exon coord (bigger coord)


  e_mxe[id] = [tS,tE,sS,sE,uS,uE,dS,dE]; ## target, second, up, down (This is different from the input file)
  c_mxe[id] = getInitialCounts();

  group = range(uS/chunk,uE/chunk+1)+range(tS/chunk,tE/chunk+1);
  group = group + range(sS/chunk,sE/chunk+1)+range(dS/chunk,dE/chunk+1);## groups this event could belong
  group = list(set(group));  ## remove duplicate groups

  if chr in mxe: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in mxe[chr]: ## this group is already there
        if id in mxe[chr][i]: ## not likely but this group already has the id
          numMXEDup+=1;
          logging.debug("Duplicate MXE ID: %d" % id);
        else: ## new MXE ID
          mxe[chr][i][id] = [tS,tE,sS,sE,uS,uE,dS,dE]; ## mxe event with coords
          numMXE += 1;
      else: ## new group to this chromosome
        mxe[chr][i]={};
        mxe[chr][i][id] = [tS,tE,sS,sE,uS,uE,dS,dE]; ## mxe event with coords
        numMXE += 1;
  else: ## first time accesing this chromosome
    mxe[chr]={};
    for i in group: ## for each possible group
      mxe[chr][i]={};
      mxe[chr][i][id] = [tS,tE,sS,sE,uS,uE,dS,dE]; ## mxe event with coords
      numMXE+=1;
#
#### A5SS ####
#
c=0;     ## count
numA5SS=0; ## number of A5SS
numA5SSDup=0; ## duplicate A5SS id
line=a5ssFile.readline(); ## skipping header
for line in a5ssFile: ## process a5ss events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3];
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4];
  lS = int(ele[5]); lE = int(ele[6]); ## long exon coord
  sS = int(ele[7]); sE = int(ele[8]); ## short exon coord
  fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord

  e_a5ss[id] = [lS,lE,sS,sE,fS,fE];
  c_a5ss[id] = getInitialCounts();

  group = range(lS/chunk, lE/chunk+1) + range(fS/chunk, fE/chunk+1); ## groups this event could belong
  group = list(set(group));  ## remove duplicate groups

  if chr in a5ss: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in a5ss[chr]: ## this group is already there
        if id in a5ss[chr][i]: ## not likely but this group already has the id
          numA5SSDup+=1;
          logging.debug("Duplicate A5SS ID: %d" % id);
        else: ## new A5SS ID
          a5ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a5ss event with coords
          numA5SS+=1;
      else: ## new group to this chromosome
        a5ss[chr][i]={};
        a5ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a5ss event with coords
        numA5SS+=1;
  else: ## first time accesing this chromosome
    a5ss[chr]={};
    for i in group: ## for each possible group
      a5ss[chr][i]={};
      a5ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a5ss event with coords
      numA5SS+=1;
#
#### A3SS ####
#
c=0;     ## count
numA3SS=0; ## number of A3SS
numA3SSDup=0; ## duplicate A3SS id
line=a3ssFile.readline(); ## skipping header
for line in a3ssFile: ## process a3ss events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3];
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4];
  lS = int(ele[5]); lE = int(ele[6]); ## long exon coord
  sS = int(ele[7]); sE = int(ele[8]); ## short exon coord
  fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord


  e_a3ss[id] = [lS,lE,sS,sE,fS,fE];
  c_a3ss[id] = getInitialCounts();

  group = range(lS/chunk, lE/chunk+1) + range(fS/chunk, fE/chunk+1); ## groups this event could belong
  group = list(set(group));  ## remove duplicate groups

  if chr in a3ss: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in a3ss[chr]: ## this group is already there
        if id in a3ss[chr][i]: ## not likely but this group already has the id
          numA3SSDup+=1;
          logging.debug("Duplicate A3SS ID: %d" % id);
        else: ## new A3SS ID
          a3ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a3ss event with coords
          numA3SS+=1;
      else: ## new group to this chromosome
        a3ss[chr][i]={};
        a3ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a3ss event with coords
        numA3SS+=1;
  else: ## first time accesing this chromosome
    a3ss[chr]={};
    for i in group: ## for each possible group
      a3ss[chr][i]={};
      a3ss[chr][i][id] = [lS,lE,sS,sE,fS,fE]; ## a3ss event with coords
      numA3SS+=1;
#
#### AFE ####
#
c=0;     ## count
numAFE=0; ## number of AFE
numAFEDup=0; ## duplicate AFE id
line=afeFile.readline(); ## skipping header
for line in afeFile: ## process afe events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3];
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4];
  dS = int(ele[5]); dE = int(ele[6]); ## distal exon coord
  pS = int(ele[7]); pE = int(ele[8]); ## proximal exon coord
  fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord


  e_afe[id] = [dS,dE,pS,pE,fS,fE]; ## afe event with coord
  c_afe[id] = getInitialCounts();

  group = range(dS/chunk, dE/chunk+1) + range(pS/chunk, pE/chunk+1) + range(fS/chunk, fE/chunk+1); ## groups this event could belong
  group = list(set(group));  ## remove duplicate groups

  if chr in afe: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in afe[chr]: ## this group is already there
        if id in afe[chr][i]: ## not likely but this group already has the id
          numAFEDup+=1;
          logging.debug("Duplicate AFE ID: %d" % id);
        else: ## new AFE ID
          afe[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## afe event with coord
          numAFE+=1;
      else: ## new group to this chromosome
        afe[chr][i]={};
        afe[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## afe event with coord
        numAFE+=1;
  else: ## first time accesing this chromosome
    afe[chr]={};
    for i in group: ## for each possible group
      afe[chr][i]={};
      afe[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## afe event with coord
      numAFE+=1;
##
##### ALE ####
##
c=0;     ## count
numALE=0; ## number of ALE
numALEDup=0; ## duplicate ALE id
line=aleFile.readline(); ## skipping header
for line in aleFile: ## process ale events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3];
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4];
  dS = int(ele[5]); dE = int(ele[6]); ## distal exon coord
  pS = int(ele[7]); pE = int(ele[8]); ## proximal exon coord
  fS = int(ele[9]); fE = int(ele[10]); ## flanking exon coord


  e_ale[id] = [dS,dE,pS,pE,fS,fE]; ## ale event with coord
  c_ale[id] = getInitialCounts();

  group = range(dS/chunk, dE/chunk+1) + range(pS/chunk, pE/chunk+1) + range(fS/chunk, fE/chunk+1); ## groups this event could belong
  group = list(set(group));  ## remove duplicate groups

  if chr in ale: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in ale[chr]: ## this group is already there
        if id in ale[chr][i]: ## not likely but this group already has the id
          numALEDup+=1;
          logging.debug("Duplicate ALE ID: %d" % id);
        else: ## new ALE ID
          ale[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## ale event with coord
          numALE+=1;
      else: ## new group to this chromosome
        ale[chr][i]={};
        ale[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## ale event with coord
        numALE+=1;
  else: ## first time accesing this chromosome
    ale[chr]={};
    for i in group: ## for each possible group
      ale[chr][i]={};
      ale[chr][i][id] = [dS,dE,pS,pE,fS,fE]; ## ale event with coord
      numALE+=1;
#
#### RI ####
#
c=0;     ## count
numRI=0; ## number of RI
numRIDup=0; ## duplicate RI id
line=riFile.readline(); ## skipping header
for line in riFile: ## process ri events file
  c+=1;
  ele = line.strip().split('\t');
  id = int(ele[0]);
  chr = ele[3];
  if chr[0:3]!='chr': ## X instead of chrX, add 'chr'
    chr = 'chr'+chr;
  strand = ele[4];
  rS = int(ele[5]); rE = int(ele[6]); ## ri exon coord (including up- and down-stream exons)
  uS = int(ele[7]); uE = int(ele[8]); ## upstream exon coord
  dS = int(ele[9]); dE = int(ele[10]); ## downstream exon coord


  e_ri[id] = [rS,rE,uS,uE,dS,dE]; ## ri event with coord
  c_ri[id] = getInitialCounts();

  group = range(rS/chunk, rE/chunk+1); ## groups this event could belong
  group = list(set(group));  ## remove duplicate groups

  if chr in ri: ## already processed this chromosome
    for i in group: ## for each possible group
      if i in ri[chr]: ## this group is already there
        if id in ri[chr][i]: ## not likely but this group already has the id
          numRIDup+=1;
          logging.debug("Duplicate RI ID: %d" % id);
        else: ## new RI ID
          ri[chr][i][id] = [rS,rE,uS,uE,dS,dE]; ## ri event with coord
          numRI+=1;
      else: ## new group to this chromosome
        ri[chr][i]={};
        ri[chr][i][id] = [rS,rE,uS,uE,dS,dE]; ## ri event with coord
        numRI+=1;
  else: ## first time accesing this chromosome
    ri[chr]={};
    for i in group: ## for each possible group
      ri[chr][i]={};
      ri[chr][i][id] = [rS,rE,uS,uE,dS,dE]; ## ri event with coord
      numRI+=1;
#
## sys.exit(0);
#
#
def processSample(sample, sInd): ## call it with processSample(sample_1, S1) something like this
  
  ### process the given sample ###
  for s1 in sample: ## for each sam file 
    rep = sample.index(s1);
    if len(s1.strip())<1: ## incorrect split. a user might accidently put comma at the end of input sam file list
      continue; ### just skip this entry, probably the last one though
    sFile = open(s1.strip()); ## open sam file
    e1 = {}; ## edge count here
    for line in sFile: ## process each line 
      if len(line.strip().split('\t'))<5 or line[0]=='#' or line[0]=='@' : ## blank line or comment
        continue;  ## go to next line
      ele = line.strip().split('\t');
      chr = ele[2];
      if chr[0:3]!='chr':
        chr = 'chr'+chr;
      mc = int(ele[3]); ## 1 base, mapping coordinate
      mString = ele[5]; ## mapping string, 50M or aMbNcM format
      group = mc/chunk; ## group does not change, it's okay to check only one group for a junction read
      if 'D' in mString or 'I' in mString or 'S' in mString or 'H' in mString or 'P' in mString or 'X' in mString or '=' in mString: ## skip
        continue; ## go to next line
      
      ### check to see if the line is either exonic read or junction read
      split_mString = mString.split('M');
      tor = 0; ## type of read, 0 nothing, 1 exonic read, 2 junction read
      if len(split_mString)==2:
        tor = 1; ############ exonic read ######
        rL = int(split_mString[0]); ## read length specified in this mapping string
        mec = mc+rL-1; ## mapping end coord 
 
        ## SE ###
        if chr in se: ## this chromosome has se event(s)
          if group in se[chr]: ## this group has skipped exon event(s)
            for c in se[chr][group]: ## for each skipped exon event in this group
              if (mc>se[chr][group][c][0] and mec<=se[chr][group][c][1]): ## read on the target
                c_se[c][CT2][sInd][I][rep]+=1; 
        ### end of SE ###

        ### MXE ####
        if chr in mxe: ## this chromosome has mxe event(s)
          if group in mxe[chr]: ## this group has mxe event(s)
            for c in mxe[chr][group]: ## for each mxe event in this group
              if (mc>mxe[chr][group][c][0] and mec<=mxe[chr][group][c][1]): ## read on the target exon
                c_mxe[c][CT2][sInd][I][rep]+=1;
              elif (mc>mxe[chr][group][c][2] and mec<=mxe[chr][group][c][3]): ## read on the second exon
                c_mxe[c][CT2][sInd][S][rep]+=1;
        ## end of MXE ###

        ## A5SS ##
        if chr in a5ss: ## this chromosome has a5ss event(s)
          if group in a5ss[chr]: ## this group has a5ss event(s)
            for c in a5ss[chr][group]: ## for each a5ss event in this group

              if a5ss[chr][group][c][4]>a5ss[chr][group][c][1]: ## positive strand
                if (mc>a5ss[chr][group][c][0] and mc<=(a5ss[chr][group][c][3]-(rL-junctionLength/2)+1) and mec<=a5ss[chr][group][c][1] and mec>=(a5ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                  c_a5ss[c][CT1][sInd][I][rep]+=1;
                  c_a5ss[c][CT2][sInd][I][rep]+=1;
                if (mc>a5ss[chr][group][c][3] and mec<=a5ss[chr][group][c][1]): ## exon read supporting target
                  c_a5ss[c][CT2][sInd][I][rep]+=1;
                
              else: ## negative strand
                if (mc>a5ss[chr][group][c][0] and mc<=(a5ss[chr][group][c][2]-(rL-junctionLength/2)+1) and mec<=a5ss[chr][group][c][1] and mec>=(a5ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                  c_a5ss[c][CT1][sInd][I][rep]+=1;
                  c_a5ss[c][CT2][sInd][I][rep]+=1;
                if (mc>a5ss[chr][group][c][0] and mec<=a5ss[chr][group][c][2]): ## exon read supporting target
                  c_a5ss[c][CT2][sInd][I][rep]+=1;

        ## end of A5SS ###
       
        ## A3SS ##
        if chr in a3ss: ## this chromosome has a3ss event(s)
          if group in a3ss[chr]: ## this group has a3ss event(s)
            for c in a3ss[chr][group]: ## for each a3ss event in this group

              if a3ss[chr][group][c][4]>a3ss[chr][group][c][1]: ## negative strand
                if (mc>a3ss[chr][group][c][0] and mc<=(a3ss[chr][group][c][3]-(rL-junctionLength/2)+1) and mec<=a3ss[chr][group][c][1] and mec>=(a3ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                  c_a3ss[c][CT1][sInd][I][rep]+=1;
                  c_a3ss[c][CT2][sInd][I][rep]+=1;
                if (mc>a3ss[chr][group][c][3] and mec<=a3ss[chr][group][c][1]): ## exon read supporting target
                  c_a3ss[c][CT2][sInd][I][rep]+=1;

              else: ## positive strand
                if (mc>a3ss[chr][group][c][0] and mc<=(a3ss[chr][group][c][2]-(rL-junctionLength/2)+1) and mec<=a3ss[chr][group][c][1] and mec>=(a3ss[chr][group][c][3]+(rL-junctionLength/2))): ## multi-exon read supporting target
                  c_a3ss[c][CT1][sInd][I][rep]+=1;
                  c_a3ss[c][CT2][sInd][I][rep]+=1;
                if (mc>a3ss[chr][group][c][0] and mec<=a3ss[chr][group][c][2]): ## exon read supporting target
                  c_a3ss[c][CT2][sInd][I][rep]+=1;

        ## end of A3SS ###

        ## AFE ##
        if chr in afe: ## this chromosome has afe event(s)
          if group in afe[chr]: ## this group has afe event(s)
            for c in afe[chr][group]: ## for each afe event in this group
              if afe[chr][group][c][4]>afe[chr][group][c][1]: ## positive strand
                if (mc>afe[chr][group][c][0] and mec<=afe[chr][group][c][1]): ## genome read supporting inclusion (first exon)
                  c_afe[c][CT2][sInd][I][rep]+=1;
                elif (mc>afe[chr][group][c][2] and mec<=afe[chr][group][c][3]): ## genome read supporting skipping (second exon)
                  c_afe[c][CT2][sInd][S][rep]+=1;
              else: ## negative strand (no difference for genome reads)
                if (mc>afe[chr][group][c][0] and mec<=afe[chr][group][c][1]): ## genome read supporting inclusion (first exon)
                  c_afe[c][CT2][sInd][I][rep]+=1;
                elif (mc>afe[chr][group][c][2] and mec<=afe[chr][group][c][3]): ## genome read supporting skipping (second exon)
                  c_afe[c][CT2][sInd][S][rep]+=1;
        ## end of AFE ###

        ## ALE ##
        if chr in ale: ## this chromosome has ale event(s)
          if group in ale[chr]: ## this group has ale event(s)
            for c in ale[chr][group]: ## for each ale event in this group
              if ale[chr][group][c][4]>ale[chr][group][c][1]: ## negative strand
                if (mc>ale[chr][group][c][2] and mec<=ale[chr][group][c][3]): ## genome read supporting inclusion (first exon)
                  c_ale[c][CT2][sInd][I][rep]+=1;
                elif (mc>ale[chr][group][c][0] and mec<=ale[chr][group][c][1]): ## genome read supporting skipping (second exon)
                  c_ale[c][CT2][sInd][S][rep]+=1;
              else: ## positive strand (no difference for genome reads)
                if (mc>ale[chr][group][c][2] and mec<=ale[chr][group][c][3]): ## genome read supporting inclusion (first exon)
                  c_ale[c][CT2][sInd][I][rep]+=1;
                elif (mc>ale[chr][group][c][0] and mec<=ale[chr][group][c][1]): ## genome read supporting skipping (second exon)
                  c_ale[c][CT2][sInd][S][rep]+=1;
        ## end of ALE ###

        ## RI ##
        if chr in ri: ## this chromosome has ale event(s)
          if group in ri[chr]: ## this group has ri event(s)
            for c in ri[chr][group]: ## for each ri event in this group, strand does not matter for ri events
              if (mc>ri[chr][group][c][0] and mc<=(ri[chr][group][c][3]-(rL-junctionLength/2)+1) and mec<=ri[chr][group][c][4] and mec>=(ri[chr][group][c][3]+(rL-junctionLength/2))) or (mc>ri[chr][group][c][3] and mc<=(ri[chr][group][c][4]-(rL-junctionLength/2)+1) and mec<=ri[chr][group][c][5] and mec>=(ri[chr][group][c][4]+(rL-junctionLength/2))): ## multi-exon read supporting target
                c_ri[c][CT1][sInd][I][rep]+=1;
                c_ri[c][CT2][sInd][I][rep]+=1;
              if (mc>ri[chr][group][c][3] and mec<=ri[chr][group][c][4]): ## exon read supporting target
                c_ri[c][CT2][sInd][I][rep]+=1;
        ## end of RI ###


      elif len(split_mString)==3: ###### junction read ###########
        tor = 2; ## junction read
        jS = mc+int(split_mString[0])-1; ## 1-base
        jE = mc+ int(split_mString[0])+ int(split_mString[1].split('N')[0])  -1; ## 0-base
        key = chr+'_'+str(jS)+'_'+str(jE)+'_0';
        if key in e1: ## exist!
          e1[key] = e1[key]+1;
        else: ## new junction
          e1[key] = 1;

        ## SE ###
        if chr in se: ## this chromosome has se event(s)
          if group in se[chr]: ## this group has skipped exon event(s)
            for c in se[chr][group]: ## for each skipped exon event in this group, examine if the given junction is part of it
              if (jS==se[chr][group][c][3] and jE==se[chr][group][c][0]) or (jS==se[chr][group][c][1] and jE==se[chr][group][c][4]): ## IJC
                c_se[c][CT1][sInd][I][rep]+=1; 
                c_se[c][CT2][sInd][I][rep]+=1; 
              elif jS==se[chr][group][c][3] and jE==se[chr][group][c][4]: ## SJC
                c_se[c][CT1][sInd][S][rep]+=1; 
                c_se[c][CT2][sInd][S][rep]+=1; 
        ### end of SE ###  

        ## MXE ###
        if chr in mxe: ## this chromosome has mxe event(s)
          if group in mxe[chr]: ## this group has mxe event(s)
            for c in mxe[chr][group]: ## for each mxe event in this group, examine if the given junction is part of it
              if (jS==mxe[chr][group][c][5] and jE==mxe[chr][group][c][0]) or (jS==mxe[chr][group][c][1] and jE==mxe[chr][group][c][6]): ## IJC
                c_mxe[c][CT1][sInd][I][rep]+=1;
                c_mxe[c][CT2][sInd][I][rep]+=1;
              elif (jS==mxe[chr][group][c][5] and jE==mxe[chr][group][c][2]) or (jS==mxe[chr][group][c][3] and jE==mxe[chr][group][c][6]): ## SJC
                c_mxe[c][CT1][sInd][S][rep]+=1;
                c_mxe[c][CT2][sInd][S][rep]+=1;
        ### end of MXE ###  

        ## A5SS ###
        if chr in a5ss: ## this chromosome has a5ss event(s)
          if group in a5ss[chr]: ## this group has a5ss event(s)
            for c in a5ss[chr][group]: ## for each a5ss event in this group, examine if the given junction is part of it
              if a5ss[chr][group][c][4]>a5ss[chr][group][c][1]: ## positive strand
                if jS==a5ss[chr][group][c][1] and jE==a5ss[chr][group][c][4]: ## IJC
                  c_a5ss[c][CT1][sInd][I][rep]+=1;
                  c_a5ss[c][CT2][sInd][I][rep]+=1;
                elif jS==a5ss[chr][group][c][3] and jE==a5ss[chr][group][c][4]: ## SJC
                  c_a5ss[c][CT1][sInd][S][rep]+=1;
                  c_a5ss[c][CT2][sInd][S][rep]+=1;
              else: ## negative strand
                if jS==a5ss[chr][group][c][5] and jE==a5ss[chr][group][c][0]: ## IJC
                  c_a5ss[c][CT1][sInd][I][rep]+=1;
                  c_a5ss[c][CT2][sInd][I][rep]+=1;
                elif jS==a5ss[chr][group][c][5] and jE==a5ss[chr][group][c][2]: ## SJC
                  c_a5ss[c][CT1][sInd][S][rep]+=1;
                  c_a5ss[c][CT2][sInd][S][rep]+=1;
        ### end of A5SS ###

        ## A3SS ###
        if chr in a3ss: ## this chromosome has a3ss event(s)
          if group in a3ss[chr]: ## this group has a3ss event(s)
            for c in a3ss[chr][group]: ## for each a3ss event in this group, examine if the given junction is part of it
              if a3ss[chr][group][c][4]>a3ss[chr][group][c][1]: ## negative strand
                if jS==a3ss[chr][group][c][1] and jE==a3ss[chr][group][c][4]: ## IJC
                  c_a3ss[c][CT1][sInd][I][rep]+=1;
                  c_a3ss[c][CT2][sInd][I][rep]+=1;
                elif jS==a3ss[chr][group][c][3] and jE==a3ss[chr][group][c][4]: ## SJC
                  c_a3ss[c][CT1][sInd][S][rep]+=1;
                  c_a3ss[c][CT2][sInd][S][rep]+=1;
              else: ## positive strand
                if jS==a3ss[chr][group][c][5] and jE==a3ss[chr][group][c][0]: ## IJC
                  c_a3ss[c][CT1][sInd][I][rep]+=1;
                  c_a3ss[c][CT2][sInd][I][rep]+=1;
                elif jS==a3ss[chr][group][c][5] and jE==a3ss[chr][group][c][2]: ## SJC
                  c_a3ss[c][CT1][sInd][S][rep]+=1;
                  c_a3ss[c][CT2][sInd][S][rep]+=1;
        ### end of A3SS ###

        ## AFE ###
        if chr in afe: ## this chromosome has afe event(s)
          if group in afe[chr]: ## this group has afe event(s)
            for c in afe[chr][group]: ## for each afe event in this group, examine if the given junction is part of it
              if afe[chr][group][c][4]>afe[chr][group][c][1]: ## positive strand
                if jS==afe[chr][group][c][3] and jE==afe[chr][group][c][4]: ## IJC
                  c_afe[c][CT1][sInd][I][rep]+=1;
                  c_afe[c][CT2][sInd][I][rep]+=1;
                elif jS==afe[chr][group][c][1] and jE==afe[chr][group][c][4]: ## SJC
                  c_afe[c][CT1][sInd][S][rep]+=1;
                  c_afe[c][CT2][sInd][S][rep]+=1;
              else: ## negative strand
                if jS==afe[chr][group][c][5] and jE==afe[chr][group][c][2]: ## IJC
                  c_afe[c][CT1][sInd][I][rep]+=1;
                  c_afe[c][CT2][sInd][I][rep]+=1;
                elif jS==afe[chr][group][c][5] and jE==afe[chr][group][c][0]: ## SJC
                  c_afe[c][CT1][sInd][S][rep]+=1;
                  c_afe[c][CT2][sInd][S][rep]+=1;
        ### end of AFE ###
#
#        ## ALE ###
        if chr in ale: ## this chromosome has ale event(s)
          if group in ale[chr]: ## this group has ale event(s)
            for c in ale[chr][group]: ## for each ale event in this group, examine if the given junction is part of it
              if ale[chr][group][c][4]>ale[chr][group][c][1]: ## negative strand
                if jS==ale[chr][group][c][3] and jE==ale[chr][group][c][4]: ## IJC
                  c_ale[c][CT1][sInd][I][rep]+=1;
                  c_ale[c][CT2][sInd][I][rep]+=1;
                elif jS==ale[chr][group][c][1] and jE==ale[chr][group][c][4]: ## SJC
                  c_ale[c][CT1][sInd][S][rep]+=1;
                  c_ale[c][CT2][sInd][S][rep]+=1;
              else: ## positive strand
                if jS==ale[chr][group][c][5] and jE==ale[chr][group][c][2]: ## IJC
                  c_ale[c][CT1][sInd][I][rep]+=1;
                  c_ale[c][CT2][sInd][I][rep]+=1;
                elif jS==ale[chr][group][c][5] and jE==ale[chr][group][c][0]: ## SJC
                  c_ale[c][CT1][sInd][S][rep]+=1;
                  c_ale[c][CT2][sInd][S][rep]+=1;
        ### end of ALE ###

        ## RI ###
        if chr in ri: ## this chromosome has ri event(s)
          if group in ri[chr]: ## this group has ri event(s)
            for c in ri[chr][group]: ## for each ri event in this group, examine if the given junction is part of it
              if jS==ri[chr][group][c][3] and jE==ri[chr][group][c][4]: ## SJC
                c_ri[c][CT1][sInd][S][rep]+=1;
                c_ri[c][CT2][sInd][S][rep]+=1;
        ### end of RI ###


      else: ## it is not exonic nor junction read. proceed to the next line
        continue;

    sFile.close();  

##### end of processSample ######

processSample(sample_1, S1);

###
def writeInputFile(h1,h2,f2,cnt): ## header 1,2,3, file 1,2,3, count dict)
  ## print header first
  f2.write(h2+'\n');
  
  for k in sorted(cnt.keys()): 
    f2.write(str(k)+'\t'+','.join(map(str,cnt[k][CT2][S1][I]))+'\t'+','.join(map(str,cnt[k][CT2][S1][S]))+'\n');

##### end of writeInputFile function #######

#

#
CT1_header = 'ID\tIJC'+'\tSJC'+'';
CT2_header = 'ID\tIC'+'\tSC'+'';
#
## SE ##
writeInputFile(CT1_header,CT2_header,JCEC_seFile,c_se);
## MXE ##
writeInputFile(CT1_header,CT2_header,JCEC_mxeFile,c_mxe);
## A5SS ##...
writeInputFile(CT1_header,CT2_header,JCEC_a5ssFile,c_a5ss);
## A3SS ##...
writeInputFile(CT1_header,CT2_header,JCEC_a3ssFile,c_a3ss);
### AFE ##
writeInputFile(CT1_header,CT2_header,JCEC_afeFile,c_afe);
### ALE ##...
writeInputFile(CT1_header,CT2_header,JCEC_aleFile,c_ale);
## RI ##...
writeInputFile(CT1_header,CT2_header,JCEC_riFile,c_ri);
#
#### close all files here ##########
seFile.close()
mxeFile.close() 
a5ssFile.close()
a3ssFile.close()
afeFile.close() 
aleFile.close() 
riFile.close()
#
JCEC_seFile.close()
JCEC_mxeFile.close() 
JCEC_a5ssFile.close()
JCEC_a3ssFile.close()
JCEC_afeFile.close() 
JCEC_aleFile.close() 
JCEC_riFile.close()
#
######################################

#############
## calculate total running time
#############
currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
#
sys.exit(0);
#
