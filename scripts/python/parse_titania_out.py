import numpy as np
from pathlib import Path
import re
from itertools import islice

class ParseInfoTNA:
   def __init__ (self, *argv):
      if not len(argv): 
         self.set_list=True
         self.rdcs=True
         self.SECONDA=True
         self.het_csum=True
         self.qfactor=True
         self.chiral_vol=True
         self.polar_angles=True
         self.dynamics=True
         self.damping=True
         self.stop_crit=True
         self.coordinates=True
         return
      else:
         self.set_list=False
         self.rdcs=False
         self.SECONDA=False
         self.het_csum=False
         self.qfactor=False
         self.chiral_vol=False
         self.polar_angles=False
         self.dynamics=False
         self.damping=False
         self.stop_crit=False
         self.coordinates=False

      for arg in argv:  
         if (arg=='all'):
            self.set_list=True
            self.rdcs=True
            self.SECONDA=True
            self.het_csum=True
            self.qfactor=True
            self.chiral_vol=True
            self.polar_angles=True
            self.dynamics=True
            self.damping=True
            self.stop_crit=True
            self.coordinates=True
         elif (arg=='RDCs'):
            self.rdcs=True
         elif (arg=='SECONDA'):
            self.SECONDA=True
            self.het_csum=True
            self.set_list=True
         elif (arg=='chiralVolume'):
            self.chiral_vol=True
         elif (arg=='Qfactor'):
            self.qfactor=True
         elif (arg=='PolarAngles'):
            self.polar_angles=True
         elif (arg=='Damping'):
            self.damping=True
         elif (arg=='StopCrit'):
            self.stop_crit=True
            self.dynamics=True
         elif (arg=='Coordinates'):
            self.coordinates=True
         else:
            print("Argument {0} not known...".format(arg))
            print("available:")
            print('   all')
            print('   SECONDA')
            print('   chiralVolume')
            print('   Qfactor')
            print('   PolarAngles')
            print('   Damping')
            print('   StopCrit')
            print('   Coordinates')
            print('Using arugment: all...')
            self.set_list=True
            self.rdcs=True
            self.SECONDA=True
            self.het_csum=True
            self.qfactor=True
            self.chiral_vol=True
            self.polar_angles=True
            self.dynamics=True
            self.damping=True
            self.stop_crit=True
            self.coordinates=True

def parse_tna_out(tnafile,silent,parser_info):
   tnafiledata = dict()
   tnafiledata['FilePath'] = tnafile
   tnafiledata['Filetrj']  = tnafile + ".trj"
   tnafiledata['Title'] = Path(tnafile).name.split('.')[0]
   if not silent: print("Parsing {}".format(tnafiledata['FilePath']))
   with open(tnafile, "r") as file:
      for line in file:
         if re.search("Number of RDCs:",line):
            tnafiledata['NOR'] = int(line.split()[3])
         elif re.search("Number of atoms:",line):
            tnafiledata['NOA'] = int(line.split()[3])
         elif re.search("Number of RDC sets:",line):
            tnafiledata['NOS'] = int(line.split()[4])
         elif parser_info.set_list and re.search("List of RDC sets",line):
            if not silent: print("Searching {0} RDC sets...".format(tnafiledata['NOS']))
            tnafiledata['Sets'] = read_set_list(file,line,tnafiledata['NOS'])
         elif parser_info.rdcs and re.search("List of RDCs",line):
            if not silent: print("Searching {0} RDCs...".format(tnafiledata['NOR']))
            tnafiledata['RDCs'] = read_rdcs(file,line,tnafiledata['NOR'])
         elif parser_info.SECONDA and re.search("Index\s+|\s+Eigenvalues\s+Kappa\s+CSUM",line):
            if not silent: print("Collecting SECONDA information...")
            tnafiledata['SECONDA']=read_SECONDA(file,line)
         elif parser_info.het_csum and re.search("SECONDA cumulative sum of heterogeneous modes:",line):
            if not silent: print("Collecting SECONDA cumulative sums of heterogeneous modes...")
            tnafiledata['Het_csum'] = red_het_csum(file,line,tnafiledata['NOR'])
         elif parser_info.qfactor and re.search("Q-factors of individual RDC sets during optimization:",line):
            if not silent: print("Collecting Q-factors...")
            tnafiledata['Qfactor']=read_Q_factors(file,line,tnafiledata['NOS'],silent)
         elif parser_info.chiral_vol and re.search("Chiral volumes of the",line):
            if not silent: print("Collecting chiral Volumes...")
            tnafiledata['ChiralVolume']=read_Vc_factors(file,line)
         elif parser_info.polar_angles and re.search("Polar angles of RDC vectors", line):
            if not silent: print("Collecting spherical coordinates...")
            tnafiledata['PolarAngles']=read_polar_angles(file,line,tnafiledata['NOR'])
         elif parser_info.dynamics and re.search("Motional analysis of RDC vectors:", line):
            if not silent: print("Collecting motional information...")
            tnafiledata['Dynamics']=read_motion_information(file,line,tnafiledata['NOR'])
   with open(tnafiledata['Filetrj'], "r") as file:
      for line in file:
         if re.search("Iteration #",line):
            tnafiledata['Iteration']=int(line.strip().split()[3])
         elif parser_info.SECONDA and re.search("Full SECONDA information on initial structure:",line):
            if not silent: print("Collecting SECONDA eigenmodes...")
            read_eigenmodes(file,line,tnafiledata['SECONDA'],tnafiledata['NOR'],silent)
            if not ( parser_info.damping and parser_info.stop_crit ): break
         elif parser_info.damping and re.search("Information on redundant internal coordinate optimizer:", line):
            if not silent: print("Reading redundant information of iteration #{0}...".format(tnafiledata['Iteration']))
            if tnafiledata['Iteration']==1: 
               tnafiledata['Redundants']=dict()
               tnafiledata['Redundants']['Damping']=dict()
               tnafiledata['Redundants']['Damping']['effectiv']=dict()
               tnafiledata['Redundants']['Damping']['effectiv']['Overall']=[]
               tnafiledata['Redundants']['Damping']['effectiv']['Bond']=[]
               tnafiledata['Redundants']['Damping']['effectiv']['Angle']=[]
               tnafiledata['Redundants']['Damping']['effectiv']['Dihedral']=[]
               tnafiledata['Redundants']['Damping']['effectiv']['RDC']=[]
               tnafiledata['Redundants']['Damping']['effectiv']['ChiralVolume']=[]
               tnafiledata['Redundants']['Damping']['static']=dict()
               tnafiledata['Redundants']['Iteration']=[]
               tnafiledata['Redundants']['Delta_x']=[]
               tnafiledata['Redundants']['Delta_q']=[]
               tnafiledata['Redundants']['StopCrit']=[]
            read_damping(file,line,tnafiledata['Redundants'],tnafiledata['Iteration'])
         elif parser_info.stop_crit and re.search("Stop criteria:", line):
            if tnafiledata['Iteration']==1: 
               tnafiledata['StopCrit']=dict()
               tnafiledata['StopCrit']['Values']=dict()
               tnafiledata['StopCrit']['Values']['SaupeMean']=[]
               tnafiledata['StopCrit']['Values']['SaupeSigma']=[]
               tnafiledata['StopCrit']['Values']['PolarMean']=[]
               tnafiledata['StopCrit']['Values']['PolarSigma']=[]
               tnafiledata['StopCrit']['Values']['VectorLength']=[]
               tnafiledata['StopCrit']['Values']['Qfactor']=[]
               tnafiledata['StopCrit']['Threshholds']=dict()
            read_stop_crit(file,line,tnafiledata['StopCrit'],tnafiledata['Iteration'])
         elif parser_info.stop_crit and re.search("overall", line):
            if tnafiledata['Iteration']==1:
               tnafiledata['Soverall']=[]
            tnafiledata['Soverall'].append(float(line.strip().split()[2]))
         elif parser_info.coordinates and re.search("coordinates", line):
            if 'Coordinates' not in tnafiledata:
               tnafiledata['Coordinates'] = []
            tnafiledata['Coordinates'].append(parse_coordinates(file, line, tnafiledata['NOA']))

   return tnafiledata

def read_set_list(file, line, NOS):
   next(file)
   SetLabels=[]
   for Set in islice(file, NOS):
      SetLabels.append(Set.strip().split()[2])
   return SetLabels

def read_rdcs(file,line,NOR):
   next(file)
   RDCs=[]
   for RDC in islice(file, NOR):
      tmp=[RDC.strip().split()[0],RDC.strip().split()[1]]
      tmp.sort()
      RDCs.append("{0} {1}".format( tmp[0], tmp[1]))
   return RDCs

def read_SECONDA(file,line):
   line=file.readline()
   SECONDA_dat=dict()
   SECONDA_dat['SECONDA']=[]
   while re.search("[0-9]+\s+\|",line):
      line=line.strip()
      data=dict()
      data['Index']=int(line.split()[0])
      values=line.split("|")[1]
      lead_index = 0
      if not ( re.search("[0-9]", values.split()[lead_index]) ) : lead_index += 1
      data['EigVal']=float(values.split()[lead_index])
      if ( lead_index == 1 ): 
         data['Kappa']= 0
         SECONDA_dat['Cutoff']=data['EigVal']
      else:
         data['Kappa']= float(values.split()[lead_index+1])
      data['csum']=  float(values.split()[lead_index+2])
      SECONDA_dat['SECONDA'].append(data)
      line=file.readline()
   return SECONDA_dat

def red_het_csum(file,line,NOR):
   next(file)
   het_csum=[]
   for hcsum in islice(file, NOR):
      het_csum.append(float(hcsum.strip().split()[3]))
   return het_csum

def read_Q_factors(file,line,NOS,silent):
   Qfactors=[0]*NOS
   set_index=0
   while set_index < NOS:
      Qfactors[set_index]=[]
      set_index+=1
   while True:
      while not ( re.search("^Block", line) or re.search("^Iter", line) ):
         line=file.readline().strip()
      if re.search("^Block", line):
         set_start=int(line.split()[4])
         set_end  =int(line.split()[6])
         next(file)
      else:
         set_start=1
         set_end=NOS
      if not silent: print ("   * reading set #{0} - #{1}".format(set_start, set_end))
      line=file.readline()
      while line.strip():
         iteration=line.strip().split("|")[0]
         qline=line.strip().split("|")[1]
         set_index=set_start-1
         q_index=0
         while set_index < set_end:
            Qfactors[set_index].append(float(qline.split()[q_index]))
            set_index += 1
            q_index += 1
         line=file.readline()
      if ( set_end == NOS ): break
   return Qfactors

def read_Vc_factors(file,line):
   ChiralVolumes=dict()
   ChiralVolumes['Centers']=[]
   ChiralVolumes['Volumes']=dict()
   while True:
      while not re.search("^Block", line):
         if re.search("Chiral volume analysis:", line):
             read_chiral_volume_analysis(file, ChiralVolumes, line) 
             return ChiralVolumes
         line=file.readline().strip()
      center_start=int(line.split()[4])
      center_end  =int(line.split()[6])
      center_num  =center_end - center_start
      line=file.readline().split("|")[1].strip()
      centers = re.compile("\[\s*?[0-9]*?\]").split(line)
      for i,center in enumerate(centers):
         if( i > center_num ): break
         ChiralVolumes['Centers'].append(center.strip())
         ChiralVolumes['Volumes'][center.strip()]=[]
      line=file.readline()
      while line.strip():
         vcline=line.strip().split("|")[1]
         center_index=0
         for volume in vcline.split():
            ChiralVolumes['Volumes'][ChiralVolumes['Centers'][center_index+center_start-1]].append(float(vcline.split()[center_index]))
            center_index += 1
         line=file.readline()
   return ChiralVolumes

def read_chiral_volume_analysis(file, ChiralVolumes, line):
    ChiralVolumes['Analysis']=dict()
    counter=0
    while True:
        while not re.search("^Block", line):
            line=file.readline().strip()
            if re.search("=======", line):
                return 0
        center_start=int(line.split()[4])
        center_end  =int(line.split()[6])
        center_num  =center_end - center_start
        while not re.search("rmsd", line):
            line=file.readline()
        center_index=0
        vcline=line.split("|")[1]
        for volume in vcline.split():
            ChiralVolumes['Analysis'][ChiralVolumes["Centers"][counter]]=(float(volume))
            center_index += 1
            counter+=1
        line=file.readline()

def read_polar_angles(file,line,NOR):
   PolarAngles=dict()
   while not re.search("Nuc 1 *Nuc 2",line):
      line=file.readline().strip()
   line=file.readline().strip()
   while re.search("([A-Z][a-z]?\S*\s*)([A-Z][a-z]?\S*\s*)",line):
      sline=line.strip().split()
      RDC="{0} {1}".format(sline[0],sline[1] )
      PolarAngles[RDC]=dict()
      PolarAngles[RDC]['initial']=dict()
      PolarAngles[RDC]['final']  =dict()
      PolarAngles[RDC]['initial']['theta']=float(sline[3])
      PolarAngles[RDC]['initial']['phi']  =float(sline[4])
      PolarAngles[RDC]['final']['theta']  =float(sline[6])
      PolarAngles[RDC]['final']['phi']    =float(sline[7])
      line=file.readline().strip()
   return PolarAngles

def read_motion_information(file,line,NOR):
   MotionalInformation=dict()
   while not re.search("^S\(overall\)", line):
      line=file.readline().strip()
   MotionalInformation['Soverall']=float(line.split()[2])
   while not re.search("Nuc 1 *Nuc 2",line):
      line=file.readline().strip()
   line=file.readline().strip()
   while re.search("([A-Z][a-z]?\S*\s*)([A-Z][a-z]?\S*\s*)",line):
      sline=line.strip().split()
      RDC="{0} {1}".format(sline[0],sline[1] )
      MotionalInformation[RDC]=dict()
      MotionalInformation[RDC]['Srdc'] = float(sline[2])
      MotionalInformation[RDC]['Sax']  = float(sline[3])
      MotionalInformation[RDC]['eta']  = float(sline[4])
      MotionalInformation[RDC]['phi']  = float(sline[5])
      MotionalInformation[RDC]['chi2'] = float(sline[6])
      MotionalInformation[RDC]['Dmax'] = float(sline[7])
      line=file.readline().strip()
   return MotionalInformation

def read_eigenmodes(file,line,SECONDA_dict,NOR,silent):
   Eigenmodes=[0]*NOR
   RDC_index=0
   while RDC_index < NOR:
      Eigenmodes[RDC_index]=[]
      RDC_index+=1
   while True:
      while not re.search("^Block", line):
         line=file.readline().strip()
      vec_start=int(line.split()[5])
      vec_end  =int(line.split()[7])
      if not silent: print ("   * reading eigenmode components #{0} - #{1}".format(vec_start, vec_end))
      next(file)
      line=file.readline()
      while line.strip():
         RDC_index=int(line.strip().split("|")[0])-1
         eline=line.strip().split("|")[2]
         vec_index=0
         while vec_index <= (vec_end-vec_start):
            Eigenmodes[RDC_index].append(float(eline.split()[vec_index]))
            vec_index += 1
         line=file.readline()
      if ( vec_end == NOR ): break
   SECONDA_dict['Eigenmodes']=Eigenmodes

def read_damping(file,line,Redundants_dict,iter):
   while not re.search("^Iteration", line):
      line=file.readline().strip()
   while True:
      if re.search("Iteration", line):
         Redundants_dict['Iteration'].append(int(line.split(":")[1].split()[0]))
      elif re.search("^Overall", line): 
         Redundants_dict['Damping']['effectiv']['Overall'].append(float(line.split(":")[1]))
      elif re.search("Bond lengths",line) : 
         Redundants_dict['Damping']['effectiv']['Bond'].append(float(line.split(":")[1].split()[2][:-1]))
         if iter == 1: Redundants_dict['Damping']['static']['Bond']=float(line.split(":")[1].split()[0])
      elif re.search("Bond angles", line):
         Redundants_dict['Damping']['effectiv']['Angle'].append(float(line.split(":")[1].split()[2][:-1]))
         if iter == 1: Redundants_dict['Damping']['static']['Angle']=float(line.split(":")[1].split()[0])
      elif re.search("Torsion angles", line):
         Redundants_dict['Damping']['effectiv']['Dihedral'].append(float(line.split(":")[1].split()[2][:-1]))
         if iter == 1: Redundants_dict['Damping']['static']['Dihedral']=float(line.split(":")[1].split()[0])
      elif re.search("RDC angles", line):
         Redundants_dict['Damping']['effectiv']['RDC'].append(float(line.split(":")[1].split()[2][:-1]))
         if iter == 1: Redundants_dict['Damping']['static']['RDC']=float(line.split(":")[1].split()[0])
      elif re.search("Chiral volumes", line):
         Redundants_dict['Damping']['effectiv']['ChiralVolume'].append(float(line.split(":")[1].split()[2][:-1]))
         if iter == 1: Redundants_dict['Damping']['static']['ChiralVolume']=float(line.split(":")[1].split()[0])
      elif re.search("Distances", line):
         if iter == 1: Redundants_dict['Damping']['static']['Distance']=float(line.split(":")[1])
      elif re.search("Delta x", line):
         Redundants_dict['Delta_x'].append(float(line.split(":")[1]))
      elif re.search("Delta q", line):
         Redundants_dict['Delta_q'].append(float(line.split(":")[1]))
      elif re.search("^Stop crit", line):
         Redundants_dict['StopCrit'].append(line.split(":")[1].split()[0]) 
         break 
      line=file.readline().strip()

def read_stop_crit(file,line,StopCrit_dict,iter):
   while not re.search("parameter", line):
      line=file.readline().strip()
   while True:
      if re.search("rmsd\(S\)", line): 
         if iter==1: StopCrit_dict['Threshholds']['SaupeMean']=float(line.split()[3][:-2])
         StopCrit_dict['Values']['SaupeMean'].append(float(line.split(":")[1].split()[0]))
      elif re.search("rmsd\(sigm\[S\]\)", line):
         if iter==1: StopCrit_dict['Threshholds']['SaupeSigma']=float(line.split()[3][:-2])
         StopCrit_dict['Values']['SaupeSigma'].append(float(line.split(":")[1].split()[0]))
      elif re.search("rmsd\(p\)", line): 
         if iter==1: StopCrit_dict['Threshholds']['PolarMean']=float(line.split()[3][:-2])
         StopCrit_dict['Values']['PolarMean'].append(float(line.split(":")[1].split()[0]))
      elif re.search("rmsd\(sigm\[p\]\)", line):
         if iter==1: StopCrit_dict['Threshholds']['PolarSigma']=float(line.split()[3][:-2])
         StopCrit_dict['Values']['PolarSigma'].append(float(line.split(":")[1].split()[0]))
      elif re.search("rmsd\(R\)", line): 
         if iter==1: StopCrit_dict['Threshholds']['VectorLength']=float(line.split()[3][:-2])
         StopCrit_dict['Values']['VectorLength'].append(float(line.split(":")[1].split()[0]))
      elif re.search("delta\(Q\)", line):
         if iter==1: StopCrit_dict['Threshholds']['Qfactor']=float(line.split()[3][:-2])
         StopCrit_dict['Values']['Qfactor'].append(float(line.split(":")[1].split()[0]))
      elif re.search("max iter", line): break
      line=file.readline().strip()

def parse_coordinates(file, line, noa):
   # The current line is just the identifier
   coordinates=[]
   for atom in range(noa):
       line = file.readline().strip().split()
       coordinates.append([])
       for coordinate in (1, 2, 3):
         coordinates[atom].append( float(line[coordinate]))
   return coordinates

