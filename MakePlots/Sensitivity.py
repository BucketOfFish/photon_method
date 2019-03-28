import ROOT
import sys,os
import array
from math import *
from ROOT import *
from array import *

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat("0.2f")

#source = "/afs/cern.ch/work/r/rshang/public/photon_method_output/2016v02/ew_susy/"
#source = "/afs/cern.ch/work/r/rshang/public/photon_method_output/2016v02/v02_multijets/"
source = "/afs/cern.ch/work/r/rshang/public/photon_method_output/2016v02/fast/"
year = '2016winter'
sr = ''
#sr = 'SRmedium'
#sr = 'SRhigh'
#sr = 'SRlow'
#sr = 'VRmedium'
#sr = 'VRhigh'
sr = 'VRlow'
#sr = 'oslo'
#sr = 'SRcompress'
#sr = 'SRJigsawMed'

isTest = False

def HighMassSignalRegionCut(tree):  # jet pT > 30 GeV
	if tree.bjet_n>0: return False
	if tree.lep_n!=2: return False
	if tree.jet_n<2: return False
	if tree.MET<250: return False
	if tree.lep_pT[0]<25: return False
	if tree.lep_pT[1]<25: return False
	if abs(tree.lep_eta[0])>2.5: return False
	if abs(tree.lep_eta[1])>2.5: return False
	if tree.jet_pT[0]<30: return False
	if tree.jet_pT[1]<30: return False
	if tree.Z_pt<80.: return False
	if tree.W01_pt<100.: return False
	if tree.mll<81: return False
	if tree.mll>101: return False
	if tree.mj0j1<70: return False
	if tree.mj0j1>100: return False
	if tree.MT2W<100: return False
	if tree.DPhi_METW01<0.5: return False
	if tree.DPhi_METW01>3.0: return False
	if tree.DR_J0J1>1.5: return False
	if tree.DR_2Lep>1.8: return False
	return True
def HighMassValidationRegionCut(tree):  # jet pT > 30 GeV
	if tree.bjet_n>0: return False
	if tree.lep_n!=2: return False
	if tree.jet_n<2: return False
	if tree.MET<250: return False
	if tree.lep_pT[0]<25: return False
	if tree.lep_pT[1]<25: return False
	if abs(tree.lep_eta[0])>2.5: return False
	if abs(tree.lep_eta[1])>2.5: return False
	if tree.jet_pT[0]<30: return False
	if tree.jet_pT[1]<30: return False
	if tree.Z_pt<80.: return False
	if tree.W01_pt<100.: return False
	if tree.mll<81: return False
	if tree.mll>101: return False
	if tree.mj0j1>60 and tree.mj0j1<100: return False
	if tree.MT2W<100: return False
	if tree.DPhi_METW01<0.5: return False
	if tree.DPhi_METW01>3.0: return False
	if tree.DR_J0J1>1.5: return False
	if tree.DR_2Lep>1.8: return False
	return True
def MediumMassSignalRegionCut(tree):  # jet pT > 30 GeV
	if tree.bjet_n>0: return False
	if tree.lep_n!=2: return False
	if tree.jet_n<2: return False
	if tree.MET<150: return False
	if tree.lep_pT[0]<25: return False
	if tree.lep_pT[1]<25: return False
	if abs(tree.lep_eta[0])>2.5: return False
	if abs(tree.lep_eta[1])>2.5: return False
	if tree.jet_pT[0]<30: return False
	if tree.jet_pT[1]<30: return False
	if tree.Z_pt<80.: return False
	if tree.W01_pt<100.: return False
	if tree.mll<81: return False
	if tree.mll>101: return False
	if tree.mj0j1<70: return False
	if tree.mj0j1>100: return False
	if tree.MT2W<100: return False
	if tree.DPhi_METW01<0.5: return False
	if tree.DPhi_METW01>3.0: return False
	if tree.DR_J0J1>1.5: return False
	if tree.DR_2Lep>1.8: return False
	return True
def MediumMassValidationRegionCut(tree):  # jet pT > 30 GeV
	if tree.bjet_n>0: return False
	if tree.lep_n!=2: return False
	if tree.jet_n<2: return False
	if tree.MET<150: return False
	if tree.lep_pT[0]<25: return False
	if tree.lep_pT[1]<25: return False
	if abs(tree.lep_eta[0])>2.5: return False
	if abs(tree.lep_eta[1])>2.5: return False
	if tree.jet_pT[0]<30: return False
	if tree.jet_pT[1]<30: return False
	if tree.Z_pt<80.: return False
	if tree.W01_pt<100.: return False
	if tree.mll<81: return False
	if tree.mll>101: return False
	if tree.mj0j1>60 and tree.mj0j1<100: return False
	if tree.MT2W<100: return False
	if tree.DPhi_METW01<0.5: return False
	if tree.DPhi_METW01>3.0: return False
	if tree.DR_J0J1>1.5: return False
	if tree.DR_2Lep>1.8: return False
	return True
def LowMassSignalRegionCut(tree): #jet pT>30 GeV
	if tree.bjet_n>0: return False
	if tree.lep_n!=2: return False
	if tree.jet_n<2: return False
	if tree.lep_pT[0]<25: return False
	if tree.lep_pT[1]<25: return False
	if abs(tree.lep_eta[0])>2.5: return False
	if abs(tree.lep_eta[1])>2.5: return False
	if tree.jet_pT[0]<30: return False
	if tree.jet_pT[1]<30: return False
	if tree.mll<81: return False
	if tree.mll>101: return False
	if tree.mWmin<70: return False
	if tree.mWmin>90: return False
	if tree.MET<100.: return False
	DPhi_METJetMin = min(tree.DPhi_METJetLeading,tree.DPhi_METJetSecond)
	DPhi_METLepMin = min(tree.DPhi_METLepLeading,tree.DPhi_METLepSecond)
	if tree.jet_n==2:
		if tree.Z_pt<60.: return False
		if tree.DPhi_METPhoton>0.8: return False
		if tree.DPhi_METWmin<1.5: return False
		if tree.Wmin_pt>0: MET_2_Wpt = tree.MET/tree.Wmin_pt
		if MET_2_Wpt>0.8: return False
		if tree.Z_pt>0: MET_2_Zpt = tree.MET/tree.Z_pt
		if MET_2_Zpt<0.6: return False
		if MET_2_Zpt>1.6: return False
	if tree.jet_n>=3:
		if tree.jet_pT[2]<30: return False
		if tree.jet_n>=6: return False
		if abs(tree.Z_eta)>1.6: return False
		if tree.Z_pt<40.: return False
		if tree.mll<86: return False
		if tree.mll>96: return False
		if tree.DPhi_METNonWminJet<2.4: return False
		if tree.DPhi_METJetLeading<2.6: return False
		if tree.DPhi_METWmin>2.2: return False
		if tree.NonWminJet_pT>0: MET_2_NonWJet = tree.MET/tree.NonWminJet_pT
		if tree.jet_pT[0]>0: MET_2_J0 = tree.MET/tree.jet_pT[0]
		if MET_2_NonWJet<0.4: return False
		if MET_2_NonWJet>0.8: return False
		if tree.DR_Wmin2Jet>2.2: return False
	return True
def LowMassValidationRegionCut(tree): #jet pT>30 GeV
	if tree.bjet_n>0: return False
	if tree.lep_n!=2: return False
	if tree.jet_n<2: return False
	if tree.lep_pT[0]<25: return False
	if tree.lep_pT[1]<25: return False
	if abs(tree.lep_eta[0])>2.5: return False
	if abs(tree.lep_eta[1])>2.5: return False
	if tree.jet_pT[0]<30: return False
	if tree.jet_pT[1]<30: return False
	if tree.mll<81: return False
	if tree.mll>101: return False
	if tree.mWmin>60 and tree.mWmin<100: return False
	if tree.MET<100.: return False
	DPhi_METJetMin = min(tree.DPhi_METJetLeading,tree.DPhi_METJetSecond)
	DPhi_METLepMin = min(tree.DPhi_METLepLeading,tree.DPhi_METLepSecond)
	if tree.jet_n==2:
		if tree.Z_pt<60.: return False
		if tree.DPhi_METPhoton>0.8: return False
		if tree.DPhi_METWmin<1.5: return False
		if tree.Wmin_pt>0: MET_2_Wpt = tree.MET/tree.Wmin_pt
		if MET_2_Wpt>0.8: return False
		if tree.Z_pt>0: MET_2_Zpt = tree.MET/tree.Z_pt
		if MET_2_Zpt<0.6: return False
		if MET_2_Zpt>1.6: return False
	if tree.jet_n>=3:
		if tree.jet_pT[2]<30: return False
		if tree.jet_n>=6: return False
		if abs(tree.Z_eta)>1.6: return False
		if tree.Z_pt<40.: return False
		if tree.mll<86: return False
		if tree.mll>96: return False
		if tree.DPhi_METNonWminJet<2.4: return False
		if tree.DPhi_METJetLeading<2.6: return False
		if tree.DPhi_METWmin>2.2: return False
		if tree.NonWminJet_pT>0: MET_2_NonWJet = tree.MET/tree.NonWminJet_pT
		if tree.jet_pT[0]>0: MET_2_J0 = tree.MET/tree.jet_pT[0]
		if MET_2_NonWJet<0.4: return False
		if MET_2_NonWJet>0.8: return False
		if tree.DR_Wmin2Jet>2.2: return False
	return True
#def JigsawMedSignalRegionCut(tree): #jet pT>30 GeV
#	if tree.bjet_n>0: return False
#	if tree.lep_n!=2: return False
#	if tree.jet_n<2: return False
#	if tree.lep_pT[0]<25: return False
#	if tree.lep_pT[1]<25: return False
#	if abs(tree.lep_eta[0])>2.5: return False
#	if abs(tree.lep_eta[1])>2.5: return False
#	if tree.jet_pT[0]<30: return False
#	if tree.jet_pT[1]<30: return False
#	if tree.mll<80.: return False
#	if tree.mll>100.: return False
#	if tree.mjj<60.: return False
#	if tree.mjj>100.: return False
#	if tree.H5PP<600.: return False
#	if tree.dphiVP>2.8: return False
#				if dphiVP<0.3: continue
#				if RPT_HT5PP<0.05: continue
#				if minH2P_to_minH3P<0.6: continue
def OsloSignalRegionCut(tree):  # jet pT > 30 GeV
	if tree.MET_rel<100: return False
	if tree.bjet_n>0: return False
	if tree.jet_n<2: return False
	if tree.lep_n!=2: return False
	if tree.mll<81: return False
	if tree.mll>101: return False
	if tree.mj0j1<60: return False
	if tree.mj0j1>100: return False
	if tree.MT2W<100: return False
	if tree.DR_2Lep<0.3: return False
	if tree.DR_2Lep>1.5: return False
	return True
def CompressedSignalRegionCut(tree): #jet pT>20 GeV
	if tree.lep_n!=2: return False
	if tree.lep_pT[0]<20: return False
	if tree.lep_pT[1]<20: return False
	if abs(tree.lep_eta[0])>2.5: return False
	if abs(tree.lep_eta[1])>2.5: return False
	if tree.jet_n<3: return False
	#if tree.jet_n>4: return False
	if tree.bjet_n>0: return False
	if tree.jet_pT[0]<20: return False
	if tree.jet_pT[1]<20: return False
	if tree.mll>50: return False
	if tree.mWmin<=0: return False
	if tree.mWmin>50: return False
	if tree.MET/tree.NonWminJet_pT>1.0: return False
	if tree.DPhi_METNonWminJet<2.9: return False
	if tree.jet_n==3:
		#if tree.jet_pT[1]<30: return False
		if tree.MET<150: return False
	if tree.jet_n==4:
		#if tree.MT2Top<300: return False
		if tree.MET<250: return False
	return True
def SignalRegion(tree):
	if sr=='SRmedium':
		if not MediumMassSignalRegionCut(tree): return False
	if sr=='SRhigh':
		if not HighMassSignalRegionCut(tree): return False
	if sr=='SRlow':
		if not LowMassSignalRegionCut(tree): return False
	if sr=='VRmedium':
		if not MediumMassValidationRegionCut(tree): return False
	if sr=='VRhigh':
		if not HighMassValidationRegionCut(tree): return False
	if sr=='VRlow':
		if not LowMassValidationRegionCut(tree): return False
	if sr=='oslo':
		if not OsloSignalRegionCut(tree): return False
	if sr=='SRcompress':
		if not CompressedSignalRegionCut(tree): return False
	if sr=='SRJigsawMed':
		if not JigsawMedSignalRegionCut(tree): return False
	return True

tag = ''
if year=='2015':
	lumi = 3200.
	table_dataset = '2015 data 3.2/fb'
	tag += '_2015'
	run_cut = ROOT.TCut('RunNumber<=301973')
if year=='2016':
	lumi = 11500.
	table_dataset = '2016 data 11.5/fb'
	tag += '_2016'
	run_cut = ROOT.TCut('RunNumber<=301973')
if year=='2016winter':
	lumi = 33000. + 3500
	table_dataset = '2016 data 36.5/fb'
	tag += '_2016winter'
	run_cut = ROOT.TCut('RunNumber<=301973')
if year=='all':
	lumi = 3200.+11500.
	table_dataset = '2015+2016 data 14.7/fb'
	tag += '_all'
	run_cut = ROOT.TCut('')

tt_n = 0
tt_raw = 0
tt_ee_file = ROOT.TFile.Open(source+'tt/ttee.root')
tt_ee_tree = tt_ee_file.Get("BaselineTree")
tt_ee_file.cd('')
tt_ee_n = 0
for entry in range(0,tt_ee_tree.GetEntries()):
	if isTest: continue
	tt_ee_tree.GetEntry(entry)
	if not SignalRegion(tt_ee_tree): continue
	tt_ee_n += 0.5*lumi*tt_ee_tree.totalWeight
	tt_raw += 1
tt_mm_file = ROOT.TFile.Open(source+'tt/ttmm.root')
tt_mm_tree = tt_mm_file.Get("BaselineTree")
tt_mm_file.cd('')
tt_mm_n = 0
for entry in range(0,tt_mm_tree.GetEntries()):
	if isTest: continue
	tt_mm_tree.GetEntry(entry)
	if not SignalRegion(tt_mm_tree): continue
	tt_mm_n += 0.5*lumi*tt_mm_tree.totalWeight
	tt_raw += 1
tt_n = tt_ee_n + tt_mm_n
print 'tt = %0.2f, raw = %0.2f'%(tt_n,tt_raw)

vv_n = 0
vv_ee_file = ROOT.TFile.Open(source+'vv/vvee.root')
vv_ee_tree = vv_ee_file.Get("BaselineTree")
vv_ee_file.cd('')
vv_ee_n = 0
for entry in range(0,vv_ee_tree.GetEntries()):
	if isTest: continue
	vv_ee_tree.GetEntry(entry)
	if not SignalRegion(vv_ee_tree): continue
	vv_ee_n += 0.5*lumi*vv_ee_tree.totalWeight
vv_mm_file = ROOT.TFile.Open(source+'vv/vvmm.root')
vv_mm_tree = vv_mm_file.Get("BaselineTree")
vv_mm_file.cd('')
vv_mm_n = 0
for entry in range(0,vv_mm_tree.GetEntries()):
	if isTest: continue
	vv_mm_tree.GetEntry(entry)
	if not SignalRegion(vv_mm_tree): continue
	vv_mm_n += 0.5*lumi*vv_mm_tree.totalWeight
vv_n = vv_ee_n + vv_mm_n
print 'vv = %0.2f'%(vv_n)

z_n = 0
#z_ee_file = ROOT.TFile.Open(source+'zjets/zee.root')
z_ee_file = ROOT.TFile.Open(source+'gdata_2016/data_2016_ee_mcmetl.root')
z_ee_tree = z_ee_file.Get("BaselineTree")
z_ee_file.cd('')
z_ee_n = 0
for entry in range(0,z_ee_tree.GetEntries()):
	if isTest: continue
	z_ee_tree.GetEntry(entry)
	if not SignalRegion(z_ee_tree): continue
	#z_ee_n += lumi*z_ee_tree.totalWeight
	if sr=='SRmedium' or sr=='SRhigh':
		z_ee_n += z_ee_tree.totalWeight*z_ee_tree.ptrw_2j30
	else:
		z_ee_n += z_ee_tree.totalWeight*z_ee_tree.ptrw_2j30
#z_mm_file = ROOT.TFile.Open(source+'zjets/zmm.root')
z_mm_file = ROOT.TFile.Open(source+'gdata_2016/data_2016_mm_mcmetl.root')
z_mm_tree = z_mm_file.Get("BaselineTree")
z_mm_file.cd('')
z_mm_n = 0
for entry in range(0,z_mm_tree.GetEntries()):
	if isTest: continue
	z_mm_tree.GetEntry(entry)
	if not SignalRegion(z_mm_tree): continue
	#z_mm_n += lumi*z_mm_tree.totalWeight
	if sr=='SRmedium' or sr=='SRhigh':
		z_mm_n += z_mm_tree.totalWeight*z_mm_tree.ptrw_2j30
	else:
		z_mm_n += z_mm_tree.totalWeight*z_mm_tree.ptrw_2j30
z_n = z_ee_n + z_mm_n
print 'z = %0.2f'%(z_n)

if tt_n<0: tt_n = 0
if vv_n<0: vv_n = 0
if z_n<0: z_n = 0
bkg_n = tt_n+vv_n+z_n

if isTest: bkg_n = 3.5

hist_weight = ROOT.TH1D('hist_weight',"",100,-1,1)
hist_weight.Sumw2()
hist_gen = ROOT.TH2D('hist_gen',"",14,100,800,14,0,700)
hist_req = ROOT.TH2D('hist_req',"",14,100,800,14,0,700)
hist_raw = ROOT.TH2D('hist_raw',"",14,100,800,14,0,700)
hist_sig = ROOT.TH2D('hist_sig',"",14,100,800,14,0,700)
hist_bkg = ROOT.TH2D('hist_bkg',"",14,100,800,14,0,700)
hist_fom = ROOT.TH2D('hist_fom',"",14,100,800,14,0,700)
hist_zn = ROOT.TH2D('hist_zn',"",14,100,800,14,0,700)

#dsid = ['392330']
#m_n2 = ['200']
#m_n1 = ['100']

#dsid = ['392359', '392327', '392319', '392305', '392336', '392364', '392321', '392353', '392311', '392326', '392334', '392307', '392351', '392318', '392312', '392315', '392348', '392335', '392324', '392341', '392303', '392352', '392350', '392302', '392328', '392310', '392363', '392325', '392361', '392365', '392358', '392355', '392306', '392320', '392317', '392300', '392337', '392309', '392343', '392304', '392345', '392340', '392362', '392354', '392347', '392313', '392333', '392360', '392314', '392308', '392330', '392322', '392342', '392344', '392349', '392339', '392338', '392331', '392329', '392357', '392301', '392346', '392316', '392332', '392356', '392323']
#m_n2 = ['650', '300', '450', '350', '500', '700', '400', '600', '350', '100', '250', '150', '600', '500', '150', '250', '600', '500', '400', '550', '450', '600', '600', '500', '450', '350', '700', '500', '700', '700', '650', '600', '350', '350', '400', '350', '500', '400', '550', '300', '550', '550', '700', '600', '550', '250', '250', '650', '450', '300', '200', '200', '550', '550', '600', '550', '500', '350', '100', '650', '300', '550', '450', '500', '600', '500']
#m_n1 = ['150', '150', '400', '100', '250', '100', '200', '150', '200', '0', '50', '100', '250', '400', '50', '150', '400', '350', '300', '300', '150', '200', '300', '100', '350', '50', '200', '200', '400', '0', '250', '50', '300', '0', '0', '250', '150', '100', '200', '100', '100', '350', '300', '100', '0', '100', '200', '50', '250', '200', '100', '150', '250', '150', '350', '400', '50', '150', '50', '350', '250', '50', '50', '300', '0', '0']

dsid = ['392359', '392327', '392305', '392336', '392364', '392321', '392353', '392311', '392326', '392334', '392307', '392351', '392318', '392312', '392315', '392348', '392335', '392324', '392341', '392303', '392350', '392302', '392328', '392310', '392363', '392325', '392361', '392365', '392358', '392355', '392306', '392320', '392317', '392300', '392337', '392309', '392343', '392304', '392345', '392340', '392362', '392354', '392347', '392313', '392333', '392360', '392314', '392308', '392330', '392322', '392342', '392344', '392349', '392338', '392331', '392329', '392357', '392301', '392346', '392316', '392332', '392356', '392323']
m_n2 = ['650', '300', '350', '500', '700', '400', '600', '350', '100', '250', '150', '600', '500', '150', '250', '600', '500', '400', '550', '450', '600', '500', '450', '350', '700', '500', '700', '700', '650', '600', '350', '350', '400', '350', '500', '400', '550', '300', '550', '550', '700', '600', '550', '250', '250', '650', '450', '300', '200', '200', '550', '550', '600', '500', '350', '100', '650', '300', '550', '450', '500', '600', '500']
m_n1 = ['150', '150', '100', '250', '100', '200', '150', '200', '0', '50', '100', '250', '400', '50', '150', '400', '350', '300', '300', '150', '300', '100', '350', '50', '200', '200', '400', '0', '250', '50', '300', '0', '0', '250', '150', '100', '200', '100', '100', '350', '300', '100', '0', '100', '200', '50', '250', '200', '100', '150', '250', '150', '350', '50', '150', '50', '350', '250', '50', '50', '300', '0', '0']

for s in range(0,len(dsid)):
	sig_n = 0
	susy_ee_file = ROOT.TFile.Open(source+'susy/%s_ee.root'%(dsid[s]))
	susy_ee_tree = susy_ee_file.Get("BaselineTree")
	susy_ee_gen = susy_ee_file.Get("hist_EventCount")
	susy_ee_file.cd('')
	sig_ee_n = 0
	sig_ee_n_raw = 0
	for entry in range(0,susy_ee_tree.GetEntries()):
		susy_ee_tree.GetEntry(entry)
		if not SignalRegion(susy_ee_tree): continue
		sig_ee_n += lumi*susy_ee_tree.totalWeight
		sig_ee_n_raw += 1
		if '392330' in dsid[s]: hist_weight.Fill(lumi*susy_ee_tree.totalWeight)
		if '392330' in dsid[s]: hist_weight.Fill(-0.5,lumi*susy_ee_tree.totalWeight)
	susy_mm_file = ROOT.TFile.Open(source+'susy/%s_mm.root'%(dsid[s]))
	susy_mm_tree = susy_mm_file.Get("BaselineTree")
	susy_mm_gen = susy_mm_file.Get("hist_EventCount")
	susy_mm_file.cd('')
	sig_mm_n = 0
	sig_mm_n_raw = 0
	for entry in range(0,susy_mm_tree.GetEntries()):
		susy_mm_tree.GetEntry(entry)
		if not SignalRegion(susy_mm_tree): continue
		sig_mm_n += lumi*susy_mm_tree.totalWeight
		sig_mm_n_raw += 1
		if '392330' in dsid[s]: hist_weight.Fill(lumi*susy_mm_tree.totalWeight)
		if '392330' in dsid[s]: hist_weight.Fill(-0.5,lumi*susy_mm_tree.totalWeight)
	sig_n = sig_ee_n + sig_mm_n
	sig_n_raw = sig_ee_n_raw + sig_mm_n_raw
	FOM = 0
	zn = 0
	if not bkg_n==0: 
		FOM = sig_n/pow(bkg_n+pow(0.1*tt_n,2)+pow(0.3*vv_n,2)+pow(0.5*z_n,2),0.5)
		#zn = ROOT.RooStats.NumberCountingUtils.BinomialExpZ(sig_n,bkg_n,0.3)
		zn = ROOT.RooStats.NumberCountingUtils.BinomialExpZ(sig_n,bkg_n,2.7/3.5)
	binx = hist_fom.GetXaxis().FindBin(float(m_n2[s]))
	biny = hist_fom.GetYaxis().FindBin(float(m_n1[s]))
	hist_fom.SetBinContent(binx,biny,FOM)
	hist_zn.SetBinContent(binx,biny,zn)
	hist_sig.SetBinContent(binx,biny,sig_n)
	hist_bkg.SetBinContent(binx,biny,bkg_n)
	hist_raw.SetBinContent(binx,biny,sig_n_raw)
	hist_gen.SetBinContent(binx,biny,susy_ee_gen.GetBinContent(1))
	reqest_number = 0
	if not sig_n_raw==0 and not FOM<0.1:
		reqest_number = susy_ee_gen.GetBinContent(1)*susy_ee_gen.GetBinContent(3)/susy_ee_gen.GetBinContent(2)*100./sig_n_raw
	hist_req.SetBinContent(binx,biny,reqest_number)
	print 'DSID: %s, N2 %s\t, N1 %s\t, sig = %0.2f\t, raw = %0.0f\t, FOM = %0.2f, Zn = %0.2f'%(dsid[s],m_n2[s],m_n1[s],sig_n,sig_n_raw,FOM,zn)

workingDir = '/afs/cern.ch/user/r/rshang/workarea/SmallExamples/NtupleTreeReader/Susy2LTree_svn/Susy2LTree'
hist_weight.SetStats(0)
hist_fom.SetStats(0)
hist_zn.SetStats(0)
hist_sig.SetStats(0)
hist_bkg.SetStats(0)
hist_raw.SetStats(0)
hist_gen.SetStats(0)
hist_req.SetStats(0)
canvas = ROOT.TCanvas("canvas","canvas", 200, 10, 600, 600)
hist_fom.SetMaximum(3)
hist_fom.GetYaxis().SetTitle('N1 mass [GeV]')
hist_fom.GetXaxis().SetTitle('N2/C1 mass [GeV]')
hist_fom.Draw("COL4Z")
hist_fom.Draw("TEXT45 same")
canvas.SaveAs('%s/output/'%(workingDir)+'SUSY_'+sr+'_FOM.pdf')
hist_zn.SetMaximum(1.6)
hist_zn.GetYaxis().SetTitle('N1 mass [GeV]')
hist_zn.GetXaxis().SetTitle('N2/C1 mass [GeV]')
hist_zn.Draw("COL4Z")
hist_zn.Draw("TEXT45 same")
canvas.SaveAs('%s/output/'%(workingDir)+'SUSY_'+sr+'_ZN.pdf')
hist_sig.SetMaximum(10.0)
hist_sig.GetYaxis().SetTitle('N1 mass [GeV]')
hist_sig.GetXaxis().SetTitle('N2/C1 mass [GeV]')
hist_sig.Draw("COL4Z")
hist_sig.Draw("TEXT45 same")
canvas.SaveAs('%s/output/'%(workingDir)+'SUSY_'+sr+'_SIG.pdf')
hist_bkg.Draw("COL4Z")
hist_bkg.Draw("TEXT45 same")
canvas.SaveAs('%s/output/'%(workingDir)+'SUSY_'+sr+'_BKG.pdf')
hist_raw.Draw("COL4Z")
hist_raw.Draw("TEXT45 same")
canvas.SaveAs('%s/output/'%(workingDir)+'SUSY_'+sr+'_RAW.pdf')
hist_gen.Draw("COL4Z")
hist_gen.Draw("TEXT45 same")
canvas.SaveAs('%s/output/'%(workingDir)+'SUSY_'+sr+'_GEN.pdf')
hist_req.Draw("COL4Z")
hist_req.Draw("TEXT45 same")
canvas.SaveAs('%s/output/'%(workingDir)+'SUSY_'+sr+'_REQ.pdf')
outputFile = ROOT.TFile('%s/output/'%(workingDir)+'SUSY_'+sr+'_FOM.root',"recreate")
hist_weight.Write()
hist_fom.Write()
hist_zn.Write()
hist_sig.Write()
hist_bkg.Write()
hist_raw.Write()
hist_gen.Write()
hist_req.Write()
