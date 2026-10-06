print("import")

from ROOT import gROOT,gStyle,TH1F,TChain,TCanvas,TLegend,TH2F,kLightTemperature
import ROOT
import array, math
import numpy as np
import os
from math import *
from DataFormats.FWLite import Events, Handle
import pickle
import sys
import pandas as pd
import matplotlib.pyplot as plt
from array import array

import re
import math
gStyle.SetPalette(kLightTemperature)

print("start")

paf = str(sys.argv[1])     # path where the data is located: root files generated after cmsRun
enename = str(sys.argv[2]) # energy of the particle

path2ring = "<need to specify!>" # path to the Rings.csv file that has the information about rings, radii and areas of all tiles; 


# No longer used TODO: check
#z42 = 4707.71 #Layer 42 Z where there are molded tiles starting from D8 
#y44m= 1537.05 #Radius of the the molded tiles in layers 42+
#area_ring={"R00": 5.3, "R02": 5.78, "R04":6.31, "R06":6.89, "R08":7.51, "R10":8.20, "R12": 8.95, "R14":9.76 , "R16":10.65, "R18":11.62, "R20":12.68, "R22":13.84, "R24":15.10, "R26":16.48, "R28": 17.98 , "R30":19.62, "R32": 21.41, "R34": 23.36, "R36": 25.49, "R38":27.82, "R40":30.36}

layer_z={34: 406.0, 35:412.0, 36: 419.0, 37: 425.0,38:431.0, 39:440.0, 40:448.0, 41:456.0, 42:464.0, 43:472.0, 44:481.0, 45:489.0, 46:497.0, 47:505.0}
board_r={"A": 110.0,"B":135.0, "D":170.0, "E":200.0, "G":240.0}
rings=pd.read_csv(path2rings+'Rings.csv')


# Masks and offsets for the DetID:
PhiOffset = 0
PhiMask = 0x1FF

RadiusOffset = 9
RadiusMask = 0xFF

LayerOffset = 17
LayerMask = 0x1F

TriggerOffset = 22
TriggerMask = 0x1

SiPMOffset = 23
SiPMMask = 0x1
SiPMMask0 = 0xFF7FFFFF

GranularityOffset = 24
GranularityMask = 0x1
GranularityMask0 = 0xFEFFFFFF

ZsideOffset = 25
ZsideMask = 0x1

TypeOffset = 26
TypeMask = 0x3
TypeMask0 = 0xF3FFFFFF

PositionMask = 0xF13FFFFF


# function to switch between ring numbers as ints and naming like "R02" etc.
def ri(x):
	try:
		int(x)
		str_r =  "R"+str(int(x)).zfill(2)
	except:
		str_r =  str(x)
	return str_r

gain = 0.7 # SiPM gain in ADC specified in the simulation; gain = 0.7 ADC is for 9 mm^2, ConvGain 4, OV 2V
thr = 0.25 # Quater of a MIP threshold, same as in simulation

# Functions for Cast and Molded light yield v size for 9mm^2 SiPM and OV 2V
def cast_ly(area):
    sqa = 10*np.sqrt(area)
    ac = 3101.7
    bc = 72.76
    cly4V = (ac/sqa) + bc
    cly2V = cly4V * ((-0.72*4 + 10*2 + 4.41)/(-0.72*16 + 10*4 + 4.41))
    return cly2V

def mold_ly(area):
    sqa = 10*np.sqrt(area)
    am = 4649.7
    bm = -17.1
    mly4V = (am/sqa) + bm
    mly2V = mly4V * ((-0.72*4 + 10*2 + 4.41)/(-0.72*16 + 10*4 + 4.41))
    return mly2V


# Function to identify which type of function to use base on ring and tile type
def which_ly2(r, t):
	#print("R"+str(int(r)).zfill(2))
	a = float(rings[rings.R=="R"+str(int(r)).zfill(2)].area.iloc[0])
	if t==2:
		patt=mold_ly(a)
		name = 'mold'
	else:
		patt=cast_ly(a)
		name = 'cast'
	return patt
# after layer 11 (cmssw nomenclature, in HGCAL terms it's layer 38) have low density tiles; layers 8-11 (or HGCAL's 34-38) are with high density
def which_rad(x, l):
	if l > 11:
		ans="R"+str(int(x)).zfill(2)
	else:
		ans="H"+str(int(x))
	return ans

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

# converting phi numbers (taken from det id) to degrees
def mydeg2(x, g):
	if g == 0:
		ans1 = x*1.25
	else:
		ans1 = x*0.8333333333
	if ans1<181:
		ans=ans1
	else:
		ans=ans1 #-360 # if wanted -180 to 180; instead 0 - 360
	return ans

# Function to print the mask name and value if detid provided; 
# Useful for debugging
def ident_s(deti):
	masks = [PhiMask, RadiusMask, LayerMask, TriggerMask, SiPMMask, GranularityMask, ZsideMask, TypeMask]#, PositionMask]
	shifts = [PhiOffset, RadiusOffset, LayerOffset, TriggerOffset, SiPMOffset, GranularityOffset, ZsideOffset, TypeOffset]
	masks_n = ['PhiMask', 'RadiusMask', 'LayerMask', 'TriggerMask', 'SiPMMask', 'GranularityMask', 'ZsideMask', 'TypeMask']
	hdet = deti #hex(deti)
	print('DET ID: {}'.format(deti))
	print('HEX DET ID: {}'.format(hex(hdet)))
	for i, m in enumerate(masks):
		hsh = hdet >> shifts[i]
		print('{}: {}'.format(masks_n[i], (hsh&m)))


# converting layer number to lambda, from the starting layer
# lambda = nuclear interaction length
# For high density layers lambda per layer = 0.314
# For low density layers lambda pwer layer = 0.428
# Calculation can be found on slide 20 of https://indico.cern.ch/event/1629299/contributions/7254179/attachments/3335925/5977362/HGCALdev_CMSSW.pdf
def lam(init_lay0, lay_i):
	init_lay = init_lay0 - 26
	if (init_lay < 12) & (lay_i < 12):
		lamb = 0.314 + 0.314*(lay_i-init_lay)
	elif (init_lay < 12) & (lay_i > 11):
		lamb0 = 0.314 + 0.314*(11-init_lay)
		lamb = lamb0 + 0.428 + 0.428*(lay_i-12)
	elif (init_lay > 11) & (lay_i > 11):
		lamb = 0.428 + 0.428*(lay_i-init_lay)
	else:
		lamb = 99
	return lamb



# regex pattern for the filename:
pc_patt = re.compile(r'gensimdigireco_pion-plus_(\d+)_board_(\w+)_z_(\d+)_seed_(\d+).\d+_sigma_(\d+).root')
#gensimdigireco_
#pion-plus_(\d+)_ 	# particle name and energy
#board_(\w+) 		# module letter, e.g. D
#_z_(\d+) 		# z position as the layer number, e.g. 34
#_seed_(\d+).\d+_	# seed number, usually of format 1.0, so inly the integer is parsed
#sigma_(\d+)		# ly sigma in %, e.g. 5

nevperseed = 100 # number of events per seed specified for the simulation
# later will be used with the extracted seed number to compute event number when combining files from different seeds

files=[] # a list of all files in the data directory to be filled
b=paf
for d in os.listdir(b):
	#print(d)
	if d.endswith(".root"):
		if pc_patt.match(os.path.basename(d)).group(1)==str(enename): # choose only files that have the requested energy of the particle
			files+=[b+d]
		else:
			continue
	else:
		continue
#print(files)

#lweight = 428/314 # ratios of lambda
lweight = 82.42/60.04   #from https://cms-hgcal-validation.web.cern.ch/cms-hgcal-validation/HGCAL_Geometry_v19_D120_dEdx_EMWeights/index.html#dedx
#to be multiplied with #MIPs, only for layers 38 and above

# new root file into which all the analysis will be written: some descriptive words about what is done + the info about particle, board, z, etc copied from the original root file
f1 = ROOT.TFile('sum_layerweights-em-event_sel_l10MT3_analysis_{}'.format(os.path.basename(files[0]).split("gensimdigireco_")[1]), 'RECREATE')
tree = ROOT.TTree("digis", "digis") 		# digis tree 
tree2 = ROOT.TTree("sums", "sums")		# tree with sums of MIPs for all events
tree3 = ROOT.TTree("layersums", "layersums")	# tree for sums per layer for all events
#treep = ROOT.TTree("pcalo", "pcalo")		# PCaloHits were also written before, but takes a lot longer and no longer needed

'''
eventidp = array('d', [0])
treep.Branch('eventid', eventidp, 'eventid/D')
pcalo_did = array('d', [0])
treep.Branch('detid', pcalo_did, 'detid/D')
layerp = array('d', [0])
treep.Branch('layer', layerp, 'layer/D')
phip = array('d', [0])
treep.Branch('phi', phip, 'phi/D')
radiusp = array('d', [0])
treep.Branch('radius', radiusp, 'radius/D')
granularityp = array('d', [0])
treep.Branch('granularity', granularityp, 'granularity/D')
typeidp = array('d', [0])
treep.Branch('typeid', typeidp, 'typeid/D')
initzp = array('d', [0])
treep.Branch('initz', initzp, 'initz/D')
initrp = array('d', [0])
treep.Branch('initr', initrp, 'initr/D')
initep = array('d', [0])
treep.Branch('inite', initep, 'inite/D')

#tiletype = array('u', ['0'])
#tree.Branch('tiletype', tiletype, 'tiletype/U')
phi_degp = array('d', [0])
treep.Branch('phi_deg', phi_degp, 'phi_deg/D')
xp = array('d', [0])
treep.Branch('x', xp, 'x/D')
yp = array('d', [0])
treep.Branch('y', yp, 'y/D')

'''


eventid = array('d', [0])
tree.Branch('eventid', eventid, 'eventid/D')
adc = array('d', [0])
tree.Branch('adc', adc, 'adc/D')
toaflag = array('d', [0])
tree.Branch('toaflag', toaflag, 'toaflag/D')	# setToAValid flag
detid = array('d', [0])
tree.Branch('detid', detid, 'detid/D')
layer = array('d', [0])
tree.Branch('layer', layer, 'layer/D')
phi = array('d', [0])
tree.Branch('phi', phi, 'phi/D')
radius = array('d', [0])
tree.Branch('radius', radius, 'radius/D')
granularity = array('d', [0])
tree.Branch('granularity', granularity, 'granularity/D')
typeid = array('d', [0])
tree.Branch('typeid', typeid, 'typeid/D')
lambdas = array('d', [0])
tree.Branch('lambdas', lambdas, 'lambdas/D')

#tiletype = array('u', ['0'])
#tree.Branch('tiletype', tiletype, 'tiletype/U')
phi_deg = array('d', [0])
tree.Branch('phi_deg', phi_deg, 'phi_deg/D')
x = array('d', [0])
tree.Branch('x', x, 'x/D')
y = array('d', [0])
tree.Branch('y', y, 'y/D')

ADC_thresh = array('d', [0])
tree.Branch('ADC_thresh', ADC_thresh, 'ADC_thresh/D')
ExpectedADC = array('d', [0])
tree.Branch('ExpectedADC', ExpectedADC, 'ExpectedADC/D')
LY = array('d', [0])
tree.Branch('LY', LY, 'LY/D')
initz = array('d', [0])
tree.Branch('initz', initz, 'initz/D')
initr = array('d', [0])
tree.Branch('initr', initr, 'initr/D')
inite = array('d', [0])
tree.Branch('inite', inite, 'inite/D')
satmips = array('d', [0])
tree.Branch('satmips', satmips, 'satmips/D')
mips = array('d', [0])
tree.Branch('mips', mips, 'mips/D')

unweightmips = array('d', [0])
tree.Branch('unweightmips', unweightmips, 'unweightmips/D')

lysigma = array('d', [0])
tree.Branch('lysigma', lysigma, 'lysigma/D')

eventid2 = array('d', [0])
tree2.Branch('eventid', eventid2, 'eventid/D')

summips = array('d', [0])
tree2.Branch('summips', summips, 'summips/D')

numdigis = array('d', [0])
tree2.Branch('numdigis', numdigis, 'numdigis/D')

sumsatmips = array('d', [0])
tree2.Branch('sumsatmips', sumsatmips, 'sumsatmips/D')

sumunwmips = array('d', [0])
tree2.Branch('sumunwmips', sumunwmips, 'sumunwmips/D')

suminite = array('d', [0])
tree2.Branch('suminite', suminite, 'suminite/D')
suminitmips = array('d', [0])
tree2.Branch('suminitmips', suminitmips, 'suminitmips/D')

suminitunwmips = array('d', [0])
tree2.Branch('suminitunwmips', suminitunwmips, 'suminitunwmips/D')


summipsL8 = array('d', [0])
tree2.Branch('summipsL8', summipsL8, 'summipsL8/D')

summipsL9 = array('d', [0])
tree2.Branch('summipsL9', summipsL9, 'summipsL9/D')

summipsL10 = array('d', [0])
tree2.Branch('summipsL10', summipsL10, 'summipsL10/D')

summipsL11 = array('d', [0])
tree2.Branch('summipsL11', summipsL11, 'summipsL11/D')

summipsL12 = array('d', [0])
tree2.Branch('summipsL12', summipsL12, 'summipsL12/D')

summipsL13 = array('d', [0])
tree2.Branch('summipsL13', summipsL13, 'summipsL13/D')

summipsL14 = array('d', [0])
tree2.Branch('summipsL14', summipsL14, 'summipsL14/D')

summipsL15 = array('d', [0])
tree2.Branch('summipsL15', summipsL15, 'summipsL15/D')

summipsL16 = array('d', [0])
tree2.Branch('summipsL16', summipsL16, 'summipsL16/D')

summipsL17 = array('d', [0])
tree2.Branch('summipsL17', summipsL17, 'summipsL17/D')

summipsL18 = array('d', [0])
tree2.Branch('summipsL18', summipsL18, 'summipsL18/D')

summipsL19 = array('d', [0])
tree2.Branch('summipsL19', summipsL19, 'summipsL19/D')

summipsL20 = array('d', [0])
tree2.Branch('summipsL20', summipsL20, 'summipsL20/D')

summipsL21 = array('d', [0])
tree2.Branch('summipsL21', summipsL21, 'summipsL21/D')



numdigisL8 = array('d', [0])
tree2.Branch('numdigisL8', numdigisL8, 'numdigisL8/D')

numdigisL9 = array('d', [0])
tree2.Branch('numdigisL9', numdigisL9, 'numdigisL9/D')

numdigisL10 = array('d', [0])
tree2.Branch('numdigisL10', numdigisL10, 'numdigisL10/D')

numdigisL11 = array('d', [0])
tree2.Branch('numdigisL11', numdigisL11, 'numdigisL11/D')

numdigisL12 = array('d', [0])
tree2.Branch('numdigisL12', numdigisL12, 'numdigisL12/D')

numdigisL13 = array('d', [0])
tree2.Branch('numdigisL13', numdigisL13, 'numdigisL13/D')

numdigisL14 = array('d', [0])
tree2.Branch('numdigisL14', numdigisL14, 'numdigisL14/D')

numdigisL15 = array('d', [0])
tree2.Branch('numdigisL15', numdigisL15, 'numdigisL15/D')

numdigisL16 = array('d', [0])
tree2.Branch('numdigisL16', numdigisL16, 'numdigisL16/D')

numdigisL17 = array('d', [0])
tree2.Branch('numdigisL17', numdigisL17, 'numdigisL17/D')

numdigisL18 = array('d', [0])
tree2.Branch('numdigisL18', numdigisL18, 'numdigisL18/D')

numdigisL19 = array('d', [0])
tree2.Branch('numdigisL19', numdigisL19, 'numdigisL19/D')

numdigisL20 = array('d', [0])
tree2.Branch('numdigisL20', numdigisL20, 'numdigisL20/D')

numdigisL21 = array('d', [0])
tree2.Branch('numdigisL21', numdigisL21, 'numdigisL21/D')

layers = array('d', [0])
tree3.Branch('layers', layers, 'layers/D')
eventid3 = array('d', [0])
tree3.Branch('eventid', eventid3, 'eventid/D')
sumMlayers = array('d', [0])
tree3.Branch('sumMlayers', sumMlayers, 'sumMlayers/D')
numHitslayers = array('d', [0])
tree3.Branch('numHitslayers', numHitslayers, 'numHitslayers/D')
lambdas3 = array('d', [0])
tree3.Branch('lambdas', lambdas3, 'lambdas/D')

layr_sums_map = {8: summipsL8, 9: summipsL9, 10: summipsL10, 11: summipsL11, 12: summipsL12, 13: summipsL13, 14: summipsL14, 15: summipsL15, 16: summipsL16, 17: summipsL17, 18: summipsL18, 19: summipsL19, 20: summipsL20, 21: summipsL21}

layr_nums_map = {8: numdigisL8, 9: numdigisL9, 10: numdigisL10, 11: numdigisL11, 12: numdigisL12, 13: numdigisL13, 14: numdigisL14, 15: numdigisL15, 16: numdigisL16, 17: numdigisL17, 18: numdigisL18, 19: numdigisL19, 20: numdigisL20, 21: numdigisL21}

'''
inimips = array('d', [0])
treep.Branch('inimips', inimips, 'inimips/D')
'''


# SiPM de-saturation function; unitless, to be applied to MIPs
# !! It is possible that the number of fired pixels would be higher than number of SiPM pixels, which will produce error in the logarithm
# TODO: need to adjust to avoid this issue
# Has not been a problem so far (10/2026)
def desat(sig, g):#, ly):
	m=39998
	n = sig/g
	m1 = m/n
	lo = math.log(1-(n/m))
	ans = -m1*math.log(1-(n/m))
	return ans



for f in files:
	#print(f)
	fname=os.path.basename(f)
	lay = int(pc_patt.match(fname).group(3)) 	# parsing the start layer from the filename
	boa = pc_patt.match(fname).group(2)		# parsing the module from the filename
	ene = float(pc_patt.match(fname).group(1))	# parsing particle energy from the filename
	sig = float(pc_patt.match(fname).group(5))	# parsing LY sigma from the filename
	see = float(pc_patt.match(fname).group(4)) - 1	# parsing the seed from the filename, since seed 1 is first, to use in function calculating event number subtracting 1
	events=TChain("Events")
	events.Add(f,-1)
	# Event selection criteria
	# l1mt = layer 1 more than: sum of MIPs in the first layer more than some threshold
	# l2mt, l3mt - conditions for subsequent 2 layers; depending on which layer number is the one that has insident particles (the z layer specified)
	# used to select events with shower containment
	# To avoid event selection, set all to 0
	# !!! needs to be specified
	if str(boa) == 'D':
		l1mt = 2
		l2mt = 2
		l3mt = 3
		lsummt = 20 # lsummt = layer sum more than; sum of mips in all three layers
		l1i = lay-26 # which layer is the first layer in cmssw nomenclature (hgcal layer - 26)
		l2i = l1i+1
		l3i = l1i+2
	else:
		l1mt = 2 # 0
		l2mt = 2 # 0
		l3mt = 3 # 0
		lsummt = 20 # 0
		l1i = lay-26 # from layer 38: 12
		l2i = l1i+1 # 13
		l3i = l1i+2 # 14

	for i, event in enumerate(events):
		#print("event {}".format(i))
		li = event.GetListOfBranches()
		le = event.GetListOfLeaves()
		summips1=0
		sumsatmips1=0
		sumunwmips1=0
		suminitmips1=0
		suminitunwmips1=0
		sumpcalo1=0
		sumdigis1=0
		layers1=0
		summipsL8i=summipsL9i=summipsL10i=summipsL11i=summipsL12i=summipsL13i=0
		summipsL14i=summipsL15i=summipsL16i=summipsL17i=summipsL18i=summipsL19i=summipsL20i=summipsL21i=0
		numdigisL8i = numdigisL9i = numdigisL10i = numdigisL11i = numdigisL12i = numdigisL13i = numdigisL14i = 0
		numdigisL15i = numdigisL16i = numdigisL17i = numdigisL18i = numdigisL19i = numdigisL20i = numdigisL21i = 0
		layr_sums_mapi = {8: summipsL8i, 9: summipsL9i, 10: summipsL10i, 11: summipsL11i, 12: summipsL12i, 13: summipsL13i, 14: summipsL14i, 15: summipsL15i, 16: summipsL16i, 17: summipsL17i, 18: summipsL18i, 19: summipsL19i, 20: summipsL20i, 21: summipsL21i}
		layr_nums_mapi = {8: numdigisL8i, 9: numdigisL9i, 10: numdigisL10i, 11: numdigisL11i, 12: numdigisL12i, 13: numdigisL13i, 14: numdigisL14i, 15: numdigisL15i, 16: numdigisL16i, 17: numdigisL17i, 18: numdigisL18i, 19: numdigisL19i, 20: numdigisL20i, 21: numdigisL21i}
		pid=[]
		for l in list(li):
			# PCaloHits are not written to the new root file, but still accessed to get the initial mips for checking the digitiser
			if "PCaloHits_g4SimHits_HGCHitsHEback_GENSIMDIGIRECO" in  l.GetName():
				for p in getattr(event, l.GetName()).product():
					pcalo_did1 = p.id()
					layr = (pcalo_did1>>LayerOffset)&LayerMask
					#pcalo_did[0] = pcalo_did1
					#layerp[0] = (pcalo_did1>>LayerOffset)&LayerMask
					#phip[0] = (pcalo_did1>>PhiOffset)&PhiMask
					#radiusp[0] = (pcalo_did1>>RadiusOffset)&RadiusMask
					#granularityp[0] = (pcalo_did1>>GranularityOffset)&GranularityMask
					#typeidp[0] = (pcalo_did1>>TypeOffset)&TypeMask
					#phi_degp[0] = mydeg2((pcalo_did1>>PhiOffset)&PhiMask, (pcalo_did1>>GranularityOffset)&GranularityMask)
					#xp[0] = pol2cart(rings[rings.R==which_rad(radius[0], layer[0])].center.iloc[0]/10, math.radians(phi_deg[0]))[0]
					#yp[0] = pol2cart(rings[rings.R==which_rad(radius[0], layer[0])].center.iloc[0]/10, math.radians(phi_deg[0]))[1]
					#initzp[0] = layer_z[lay]
					#initrp[0] = board_r[boa]
					#initep[0] = ene
					enr = p.energy() * 1000000
					initmipsi0 = enr/455.0
					if(layr > 11):
						initmipsi = lweight*initmipsi0
					else:
						initmipsi = initmipsi0


					#inimips[0] = initmipsi
					suminitmips1+=initmipsi
					suminitunwmips1+=initmipsi0
					#treep.Fill()
			if "DetIdHGCSampleHGCDataFramesSorted_mix_HGCDigisHEback_GENSIMDIGIRECO" in  l.GetName():
				for p in getattr(event, l.GetName()).product():
					eventid[0] = i+(nevperseed*see)
					did = (p.id()).rawId() # detector id
					toaflag[0] = p.data().data().toa() # Not used, TODO: remove
					detid[0] = did
					# !!! Very important to have this in order to properly interpret the data value from the digis!
					toaval = p.data()[2].getToAValid() # ToAValid was set to True if adc value was in the TOT range
					adcs1 = p.data()[2].data()
					if (toaval==True):
						adcs=3*adcs1 # emulating the TOT to ADC conversion; a very basic implementation
					else:
						adcs=adcs1
					adc[0] = adcs
					layr = (did>>LayerOffset)&LayerMask
					layer[0] = layr
					lambd = lam(lay, layr)
					lambdas[0] = lambd
					phi[0] = (did>>PhiOffset)&PhiMask
					radius[0] = (did>>RadiusOffset)&RadiusMask
					granularity[0] = (did>>GranularityOffset)&GranularityMask
					typeid[0] = (did>>TypeOffset)&TypeMask
					phi_deg[0] = mydeg2((did>>PhiOffset)&PhiMask, (did>>GranularityOffset)&GranularityMask)
					x[0] = pol2cart(rings[rings.R==which_rad(radius[0], layer[0])].center.iloc[0]/10, math.radians(phi_deg[0]))[0]
					y[0] = pol2cart(rings[rings.R==which_rad(radius[0], layer[0])].center.iloc[0]/10, math.radians(phi_deg[0]))[1]
					ADC_thresh[0] = thr*gain*which_ly2(radius[0], typeid[0])
					expadc=gain*which_ly2(radius[0], typeid[0])
					ExpectedADC[0] = expadc
					ly = which_ly2(radius[0], typeid[0])
					LY[0] = ly
					initz[0] = layer_z[lay]
					initr[0] = board_r[boa]
					inite[0] = ene
					lysigma[0] = sig
					sumdigis1+=1
					satmipss = adcs/expadc
					# applying weights
					unweightmips1=satmipss*desat(adcs, gain) # filling unweighted branches to compare later
					if(layr > 11):
						mipss = lweight*satmipss*desat(adcs, gain)
					else:
						mipss = satmipss*desat(adcs, gain)

					satmips[0] = satmipss
					sumsatmips1+=satmipss
					mips[0] = mipss
					summips1+=mipss
					unweightmips[0]=unweightmips1
					sumunwmips1+=unweightmips1
					layr_sums_mapi[layr]+=mipss
					layr_nums_mapi[layr]+=1
					tree.Fill() 	# digis tree; saving all digis regardless of event selection
		if ((layr_sums_mapi[l1i]>l1mt) or (layr_sums_mapi[l2i]>l2mt) or (layr_sums_mapi[l3i]>l3mt)):	# event selection for sums and layer sums
			if (layr_sums_mapi[l1i]+layr_sums_mapi[l2i]+layr_sums_mapi[l3i]>lsummt) and (ene > 10.0):
				sumsatmips[0]=sumsatmips1
				sumunwmips[0]=sumunwmips1
				suminitmips[0]=suminitmips1
				suminitunwmips[0]=suminitunwmips1
				summips[0]=summips1
				numdigis[0]=sumdigis1
				suminite[0]=ene
				eventid2[0]=i+(nevperseed*see)
				for k, la in enumerate(np.arange(lay-26, 22)):
					layr_sums_map[la][0]=layr_sums_mapi[la]
					layr_nums_map[la][0]=layr_nums_mapi[la]
					sumMlayers[0]=layr_sums_mapi[la]
					numHitslayers[0]=layr_nums_mapi[la]
					layers[0]=la
					lambdas3[0]=lam(lay, la)
					eventid3[0]=i+(nevperseed*see)
					tree3.Fill()	# filling layersums tree
				tree2.Fill() # filling sums tree
			elif (ene <=10.0):	# no event selection done for energies below 10GeV (assume reasonable shower containment, since energy is small)
				sumsatmips[0]=sumsatmips1
				sumunwmips[0]=sumunwmips1
				suminitmips[0]=suminitmips1
				suminitunwmips[0]=suminitunwmips1
				summips[0]=summips1
				numdigis[0]=sumdigis1
				suminite[0]=ene
				eventid2[0]=i+(nevperseed*see)
				for k, la in enumerate(np.arange(lay-26, 22)):
					layr_sums_map[la][0]=layr_sums_mapi[la]
					layr_nums_map[la][0]=layr_nums_mapi[la]
					sumMlayers[0]=layr_sums_mapi[la]
					numHitslayers[0]=layr_nums_mapi[la]
					eventid3[0]=i+(nevperseed*see)
					layers[0]=la
					lambdas3[0]=lam(lay, la)
					tree3.Fill()
				tree2.Fill()
			else:
				continue
		else:
			continue
tree.Write()
#treep.Write()
tree2.Write()
tree3.Write()
f1.Close()
sys.exit(0)
