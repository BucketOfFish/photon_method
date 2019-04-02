#
# make run scripts for v20 ntuples
# the input ntuples should be named as DSID.root
#

import ROOT
from ROOT import TFile, TIter, TKey, TTree, TH1, TText

def getBigDataFilesFromOthers(eos_path, filelist, files):
	if filelist=="": return
	inputFile = open(filelist, 'r')
	line = 1
	n = 0
	while line:
		line = inputFile.readline()
		if not line: break
		line = line.strip('\n')
		dsid = line.split('.')[1]
		files += [dsid]
		n += 1
def getMassPoints(eos_path, filelist, mass_n2, mass_n1):
	inputFile = open(filelist, 'r')
	line = 1
	n = 0
	while line:
		line = inputFile.readline()
		if not line: break
		line = line.strip('\n')
		m2 = line.split('_')[5]
		m1 = line.split('_')[6]
		mass_n2 += [m2]
		mass_n1 += [m1]
		n += 1

script_dir = '../Scripts'
exec_dir = '../Root'


#sm_path = "/eos/atlas/user/b/benhoob/forRuo_APR27/highpt_bkg/"
#data_path = "/eos/atlas/user/b/benhoob/forRuo_APR27/highpt_data16/"

sm_path = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/MC/"
data_path = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Data/"

#sm_path = "/eos/atlas/user/b/benhoob/forRuo_APR27/highpt_bkg/"
#data_path = "/eos/atlas/user/r/rshang/bphy_test/data/"

susy_mass_n2 = []
susy_mass_n1 = []


gmc_path = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Ruo/photon_method/gmc/"
#gmc_path = "/eos/atlas/user/b/benhoob/forRuo2/gmc/"
gmc_list = "list_photon_mc15c_jetm4.txt"
gmc_files = []

Vg_path = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Ruo/photon_method/gV/"
#Vg_path = "/eos/atlas/user/b/benhoob/forRuo2/gX/"
Vg_list = "list_Vgamma_mc15c_jetm4.txt"
Vg_files = []

gdata_path_2016 = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Ruo/photon_method/gdata16/"
#gdata_path_2016 = "/eos/atlas/user/b/benhoob/forRuo2/gdata16/"
gdata_list_2016 = "list_photon_data16_p2840.txt"
gdata_files_2016 = []

gdata_path_2015 = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Ruo/photon_method/gdata15/"
#gdata_path_2015 = "/eos/atlas/user/b/benhoob/forRuo2/gdata15/"
gdata_list_2015 = "list_photon_data15_p2667.txt"
gdata_files_2015 = []

bphys_path = "/eos/atlas/user/r/rshang/bphy_test/bdata/"
#bphys_path = data_path
bphys_list = ""
bphys_files = ['bdata']

data_ee_path_2016 = data_path
data_ee_list_2016 = ""
data_ee_files_2016 = ['data']
#data_ee_files_2016 = ['periodA','periodB','periodC','periodD','periodE','periodF','periodG','periodI','periodK','periodL']
data_mm_path_2016 = data_path
data_mm_list_2016 = ""
data_mm_files_2016 = ['data']
#data_mm_files_2016 = ['periodA','periodB','periodC','periodD','periodE','periodF','periodG','periodI','periodK','periodL']
data_em_path_2016 = data_path
data_em_list_2016 = ""
data_em_files_2016 = ['data']
#data_em_files_2016 = ['periodA','periodB','periodC','periodD','periodE','periodF','periodG','periodI','periodK','periodL']

#data_ee_path_2016 = "/afs/cern.ch/work/r/rshang/private/Frameworks/SusyZMETjetsWorkArea_test/SusyZMETjetsWorkArea/Run/SusyZMETjetsOutput_susy2/data-outputTree/"
#data_ee_list_2016 = ""
#data_ee_files_2016 = ['data16_13TeV.00308047.physics_Main.merge.AOD.r9264_p3083']
#data_mm_path_2016 = "/afs/cern.ch/work/r/rshang/private/Frameworks/SusyZMETjetsWorkArea_test/SusyZMETjetsWorkArea/Run/SusyZMETjetsOutput_susy2/data-outputTree/"
#data_mm_list_2016 = ""
#data_mm_files_2016 = ['data16_13TeV.00308047.physics_Main.merge.AOD.r9264_p3083']

vvee_path = sm_path
vvee_list = "list_vv_susy2.txt"
vvee_files = []

vvmm_path = sm_path
vvmm_list = "list_vv_susy2.txt"
vvmm_files = []

vvem_path = sm_path
vvem_list = "list_vv_susy2.txt"
vvem_files = []

zee_path = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Ruo/photon_method/zmc/"
#zee_path = "/eos/atlas/user/b/benhoob/forRuo2/zmc/"
zee_list = "list_ZJets_ee_25ns.txt"
zee_files = []

zmm_path = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/Ruo/photon_method/zmc/"
#zmm_path = "/eos/atlas/user/b/benhoob/forRuo2/zmc/"
zmm_list = "list_ZJets_mm_25ns.txt"
zmm_files = []

wee_path = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/MC/"
wee_list = "list_WJets_Sherpa221_ee.txt"
wee_files = []

wmm_path = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04/MC/"
wmm_list = "list_WJets_Sherpa221_mm.txt"
wmm_files = []

z221ee_path = sm_path
z221ee_list = "list_ZJets_Sherpa221_ee.txt"
z221ee_files = []

z221mm_path = sm_path
z221mm_list = "list_ZJets_Sherpa221_mm.txt"
z221mm_files = []

zttee_path = sm_path
zttee_list = "list_ZJets_Sherpa221_tt.txt"
zttee_files = []

zttmm_path = sm_path
zttmm_list = "list_ZJets_Sherpa221_tt.txt"
zttmm_files = []

ttee_path = sm_path
ttee_list = "list_top_susy2.txt"
ttee_files = []
ttmm_path = sm_path
ttmm_list = "list_top_susy2.txt"
ttmm_files = []
ttem_path = sm_path
ttem_list = "list_top_susy2.txt"
ttem_files = []

susy_ee_path = "root://eosatlas//eos/atlas/user/r/rshang/ew_susy_highpt/"
susy_ee_list = ""
susy_ee_files = ['392354','392304','392330','392318','392308'] #392354=(600,100); 392304=(300,100); 392330=(200,100); 392318=(500,400); 392308=(300,200)
susy_mm_path = "root://eosatlas//eos/atlas/user/r/rshang/ew_susy_highpt/"
susy_mm_list = ""
susy_mm_files = ['392354','392304','392330','392318','392308'] #392354=(600,100); 392304=(300,100); 392330=(200,100); 392318=(500,400); 392308=(300,200)

#susy_ee_path = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04_lowpt_new_signal/"
#susy_ee_list = ""
#susy_ee_files = ['374226','374231','374238','374239','374241','372449'] #374226=(1000,970), 374231=(1200,1160), 374238=(400,370), 374239=(400,376), 374241=(400,384), m(gluino,N1)
#susy_mm_path = "/eos/atlas/atlascerngroupdisk/phys-susy/strong2L/v02-04_lowpt_new_signal/"
#susy_mm_list = ""
#susy_mm_files = ['374226','374231','374238','374239','374241','372449'] #374226=(1000,970), 374231=(1200,1160), 374238=(400,370), 374239=(400,376), 374241=(400,384), m(gluino,N1)

getBigDataFilesFromOthers(zee_path,zee_list,zee_files)
getBigDataFilesFromOthers(zmm_path,zmm_list,zmm_files)
getBigDataFilesFromOthers(wee_path,wee_list,wee_files)
getBigDataFilesFromOthers(wmm_path,wmm_list,wmm_files)
getBigDataFilesFromOthers(z221ee_path,z221ee_list,z221ee_files)
getBigDataFilesFromOthers(z221mm_path,z221mm_list,z221mm_files)
getBigDataFilesFromOthers(zttee_path,zttee_list,zttee_files)
getBigDataFilesFromOthers(zttmm_path,zttmm_list,zttmm_files)
getBigDataFilesFromOthers(vvee_path,vvee_list,vvee_files)
getBigDataFilesFromOthers(vvmm_path,vvmm_list,vvmm_files)
getBigDataFilesFromOthers(vvem_path,vvem_list,vvem_files)
getBigDataFilesFromOthers(gmc_path,gmc_list,gmc_files)
getBigDataFilesFromOthers(Vg_path,Vg_list,Vg_files)
getBigDataFilesFromOthers(gdata_path_2016,gdata_list_2016,gdata_files_2016)
getBigDataFilesFromOthers(gdata_path_2015,gdata_list_2015,gdata_files_2015)
getBigDataFilesFromOthers(data_ee_path_2016,data_ee_list_2016,data_ee_files_2016)
getBigDataFilesFromOthers(data_mm_path_2016,data_mm_list_2016,data_mm_files_2016)
getBigDataFilesFromOthers(data_em_path_2016,data_em_list_2016,data_em_files_2016)
getBigDataFilesFromOthers(ttee_path,ttee_list,ttee_files)
getBigDataFilesFromOthers(ttmm_path,ttmm_list,ttmm_files)
getBigDataFilesFromOthers(ttem_path,ttem_list,ttem_files)


bsub_data_File = open(script_dir+'/bsub_data.sh', 'w')
local_data_File = open(script_dir+'/local_data.sh', 'w')
bsub_tt_File = open(script_dir+'/bsub_tt.sh', 'w')
local_tt_File = open(script_dir+'/local_tt.sh', 'w')
bsub_vv_File = open(script_dir+'/bsub_vv.sh', 'w')
local_vv_File = open(script_dir+'/local_vv.sh', 'w')
bsub_z_File = open(script_dir+'/bsub_z.sh', 'w')
local_z_File = open(script_dir+'/local_z.sh', 'w')
local_w_File = open(script_dir+'/local_w.sh', 'w')
bsub_z221_File = open(script_dir+'/bsub_z221.sh', 'w')
local_z221_File = open(script_dir+'/local_z221.sh', 'w')
bsub_gdata_File = open(script_dir+'/bsub_gdata.sh', 'w')
local_gdata_File = open(script_dir+'/local_gdata.sh', 'w')
bsub_gmc_File = open(script_dir+'/bsub_gmc.sh', 'w')
local_gmc_File = open(script_dir+'/local_gmc.sh', 'w')
bsub_Vg_File = open(script_dir+'/bsub_Vg.sh', 'w')
local_Vg_File = open(script_dir+'/local_Vg.sh', 'w')
local_bphys_File = open(script_dir+'/local_bphys.sh', 'w')
bsub_susy_File = open(script_dir+'/bsub_susy.sh', 'w')
local_susy_File = open(script_dir+'/local_susy.sh', 'w')

bsub_data_File.write("bsub -q 8nh -J other < local_data.sh\n")
bsub_z_File.write("bsub -q 8nh -J zmc < local_z.sh\n")
bsub_z221_File.write("bsub -q 8nh -J z221 < local_z221.sh\n")
bsub_gmc_File.write("bsub -q 8nh -J gmc < local_gmc.sh\n")
bsub_Vg_File.write("bsub -q 8nh -J Vg < local_Vg.sh\n")
bsub_gdata_File.write("bsub -q 8nh -J gdata < local_gdata.sh\n")


for i in range(0,len(data_em_files_2016)):
	outputFile = open(script_dir+'/run_%s_data_em_2016.sh'%(data_em_files_2016[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_data_File.write("sh run_%s_data_em_2016.sh\n"%(data_em_files_2016[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+data_em_files_2016[i]+'"','"data"','"'+data_em_path_2016+'"','"em"','true'))

for i in range(0,len(data_mm_files_2016)):
	outputFile = open(script_dir+'/run_%s_data_mm_2016.sh'%(data_mm_files_2016[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_data_File.write("sh run_%s_data_mm_2016.sh\n"%(data_mm_files_2016[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+data_mm_files_2016[i]+'"','"data"','"'+data_mm_path_2016+'"','"mm"','true'))

for i in range(0,len(data_ee_files_2016)):
	outputFile = open(script_dir+'/run_%s_data_ee_2016.sh'%(data_ee_files_2016[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_data_File.write("sh run_%s_data_ee_2016.sh\n"%(data_ee_files_2016[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+data_ee_files_2016[i]+'"','"data"','"'+data_ee_path_2016+'"','"ee"','true'))


for i in range(0,len(zmm_files)):
	outputFile = open(script_dir+'/run_%s_zmm.sh'%(zmm_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_z_File.write("sh run_%s_zmm.sh\n"%(zmm_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+zmm_files[i]+'"','"zjets"','"'+zmm_path+'"','"mm"','false'))

for i in range(0,len(zee_files)):
	outputFile = open(script_dir+'/run_%s_zee.sh'%(zee_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_z_File.write("sh run_%s_zee.sh\n"%(zee_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+zee_files[i]+'"','"zjets"','"'+zee_path+'"','"ee"','false'))

for i in range(0,len(wmm_files)):
	outputFile = open(script_dir+'/run_%s_wmm.sh'%(wmm_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_w_File.write("sh run_%s_wmm.sh\n"%(wmm_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+wmm_files[i]+'"','"wjets"','"'+wmm_path+'"','"mm"','false'))

for i in range(0,len(wee_files)):
	outputFile = open(script_dir+'/run_%s_wee.sh'%(wee_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_w_File.write("sh run_%s_wee.sh\n"%(wee_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+wee_files[i]+'"','"wjets"','"'+wee_path+'"','"ee"','false'))

for i in range(0,len(zttee_files)):
	outputFile = open(script_dir+'/run_%s_zttee.sh'%(zttee_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_z221_File.write("sh run_%s_zttee.sh\n"%(zttee_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+zttee_files[i]+'"','"zjets_221"','"'+zttee_path+'"','"ee"','false'))

for i in range(0,len(zttmm_files)):
	outputFile = open(script_dir+'/run_%s_zttmm.sh'%(zttmm_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_z221_File.write("sh run_%s_zttmm.sh\n"%(zttmm_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+zttmm_files[i]+'"','"zjets_221"','"'+zttmm_path+'"','"mm"','false'))

for i in range(0,len(z221mm_files)):
	outputFile = open(script_dir+'/run_%s_z221mm.sh'%(z221mm_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_z221_File.write("sh run_%s_z221mm.sh\n"%(z221mm_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+z221mm_files[i]+'"','"zjets_221"','"'+z221mm_path+'"','"mm"','false'))

for i in range(0,len(z221ee_files)):
	outputFile = open(script_dir+'/run_%s_z221ee.sh'%(z221ee_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_z221_File.write("sh run_%s_z221ee.sh\n"%(z221ee_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+z221ee_files[i]+'"','"zjets_221"','"'+z221ee_path+'"','"ee"','false'))


for i in range(0,len(ttem_files)):
	outputFile = open(script_dir+'/run_%s_ttem.sh'%(ttem_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_tt_File.write("sh run_%s_ttem.sh\n"%(ttem_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+ttem_files[i]+'"','"tt"','"'+ttem_path+'"','"em"','false'))

for i in range(0,len(ttmm_files)):
	outputFile = open(script_dir+'/run_%s_ttmm.sh'%(ttmm_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_tt_File.write("sh run_%s_ttmm.sh\n"%(ttmm_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+ttmm_files[i]+'"','"tt"','"'+ttmm_path+'"','"mm"','false'))

for i in range(0,len(ttee_files)):
	outputFile = open(script_dir+'/run_%s_ttee.sh'%(ttee_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_tt_File.write("sh run_%s_ttee.sh\n"%(ttee_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+ttee_files[i]+'"','"tt"','"'+ttee_path+'"','"ee"','false'))

for i in range(0,len(vvem_files)):
	outputFile = open(script_dir+'/run_%s_vvem.sh'%(vvem_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_vv_File.write("sh run_%s_vvem.sh\n"%(vvem_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+vvem_files[i]+'"','"vv"','"'+vvem_path+'"','"em"','false'))

for i in range(0,len(vvmm_files)):
	outputFile = open(script_dir+'/run_%s_vvmm.sh'%(vvmm_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_vv_File.write("sh run_%s_vvmm.sh\n"%(vvmm_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+vvmm_files[i]+'"','"vv"','"'+vvmm_path+'"','"mm"','false'))

for i in range(0,len(vvee_files)):
	outputFile = open(script_dir+'/run_%s_vvee.sh'%(vvee_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_vv_File.write("sh run_%s_vvee.sh\n"%(vvee_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+vvee_files[i]+'"','"vv"','"'+vvee_path+'"','"ee"','false'))

for i in range(0,len(susy_mm_files)):
	outputFile = open(script_dir+'/run_%s_susy_mm.sh'%(susy_mm_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_susy_File.write("sh run_%s_susy_mm.sh\n"%(susy_mm_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+susy_mm_files[i]+'"','"susy"','"'+susy_mm_path+'"','"mm"','false'))

for i in range(0,len(susy_ee_files)):
	outputFile = open(script_dir+'/run_%s_susy_ee.sh'%(susy_ee_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_susy_File.write("sh run_%s_susy_ee.sh\n"%(susy_ee_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBaseLineEvents.C+(%s,%s,%s,%s,%s)' \n" %('"'+susy_ee_files[i]+'"','"susy"','"'+susy_ee_path+'"','"ee"','false'))

for i in range(0,len(gdata_files_2016)):
	outputFile = open(script_dir+'/run_%s_gdata16.sh'%(gdata_files_2016[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_gdata_File.write("sh run_%s_gdata16.sh\n"%(gdata_files_2016[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetPhotonEvents.C+(%s,%s,%s,%s)' \n" %('"'+gdata_files_2016[i]+'"','"gdata"','"'+gdata_path_2016+'"','1'))

#for i in range(0,len(gdata_files_2015)):
#	outputFile = open(script_dir+'/run_%s_gdata15.sh'%(gdata_files_2015[i]), 'w')
#	outputFile.write("cd %s \n"%(exec_dir))
#	local_gdata_File.write("sh run_%s_gdata15.sh\n"%(gdata_files_2015[i]))
#	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
#	outputFile.write("lsetup ROOT\n")
#	outputFile.write("root -l -b -q 'GetPhotonEvents.C+(%s,%s,%s,%s)' \n" %('"'+gdata_files_2015[i]+'"','"gdata"','"'+gdata_path_2015+'"','1'))

for i in range(0,len(gmc_files)):
	outputFile = open(script_dir+'/run_%s_gmc.sh'%(gmc_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_gmc_File.write("sh run_%s_gmc.sh\n"%(gmc_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetPhotonEvents.C+(%s,%s,%s,%s)' \n" %('"'+gmc_files[i]+'"','"gmc"','"'+gmc_path+'"','0'))

for i in range(0,len(Vg_files)):
	outputFile = open(script_dir+'/run_%s_Vg.sh'%(Vg_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_Vg_File.write("sh run_%s_Vg.sh\n"%(Vg_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetPhotonEvents.C+(%s,%s,%s,%s)' \n" %('"'+Vg_files[i]+'"','"Vg"','"'+Vg_path+'"','2'))

for i in range(0,len(bphys_files)):
	outputFile = open(script_dir+'/run_%s_bphys.sh'%(bphys_files[i]), 'w')
	outputFile.write("cd %s \n"%(exec_dir))
	local_bphys_File.write("sh run_%s_bphys.sh\n"%(bphys_files[i]))
	outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
	outputFile.write("lsetup ROOT\n")
	outputFile.write("root -l -b -q 'GetBPhysicsEvents.C+(%s,%s,%s,%s)' \n" %('"'+bphys_files[i]+'"','"bphys"','"'+bphys_path+'"','1'))


outputFile = open(script_dir+'/run_gmc_sm.sh', 'w')
outputFile.write("cd %s \n"%(exec_dir))
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
outputFile.write("lsetup ROOT\n")
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/utilities/oldAliasSetup.sh root\n")
outputFile.write("root -l -b -q 'GetPhotonSmearing.C+(%s,%s,%s)' \n" %('"ee"','0','"2016"'))
outputFile.write("root -l -b -q 'GetPhotonSmearing.C+(%s,%s,%s)' \n" %('"mm"','0','"2016"'))

outputFile = open(script_dir+'/run_gdata_sm.sh', 'w')
outputFile.write("cd %s \n"%(exec_dir))
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
outputFile.write("lsetup ROOT\n")
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/utilities/oldAliasSetup.sh root\n")
outputFile.write("root -l -b -q 'GetPhotonSmearing.C+(%s,%s,%s)' \n" %('"ee"','1','"2016"'))
outputFile.write("root -l -b -q 'GetPhotonSmearing.C+(%s,%s,%s)' \n" %('"mm"','1','"2016"'))

outputFile = open(script_dir+'/run_bphys_sm.sh', 'w')
outputFile.write("cd %s \n"%(exec_dir))
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
outputFile.write("lsetup ROOT\n")
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/utilities/oldAliasSetup.sh root\n")
outputFile.write("root -l -b -q 'GetJpsiSmearing.C+(%s,%s,%s)' \n" %('"ee"','0','"2016"'))
outputFile.write("root -l -b -q 'GetJpsiSmearing.C+(%s,%s,%s)' \n" %('"mm"','0','"2016"'))
outputFile.write("root -l -b -q 'GetJpsiSmearing.C+(%s,%s,%s)' \n" %('"ee"','1','"2016"'))
outputFile.write("root -l -b -q 'GetJpsiSmearing.C+(%s,%s,%s)' \n" %('"mm"','1','"2016"'))
outputFile.write("root -l -b -q 'GetJpsiSmearing.C+(%s,%s,%s)' \n" %('"ee"','2','"2016"'))
outputFile.write("root -l -b -q 'GetJpsiSmearing.C+(%s,%s,%s)' \n" %('"mm"','2','"2016"'))

outputFile = open(script_dir+'/run_gmc_rw.sh', 'w')
outputFile.write("cd %s \n"%(exec_dir))
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
outputFile.write("lsetup ROOT\n")
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/utilities/oldAliasSetup.sh root\n")
outputFile.write("root -l -b -q 'GetPhotonReweighting.C+(%s,%s,%s)' \n" %('"ee"','0','"2016"'))
outputFile.write("root -l -b -q 'GetPhotonReweighting.C+(%s,%s,%s)' \n" %('"mm"','0','"2016"'))

outputFile = open(script_dir+'/run_gdata_rw.sh', 'w')
outputFile.write("cd %s \n"%(exec_dir))
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
outputFile.write("lsetup ROOT\n")
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/utilities/oldAliasSetup.sh root\n")
outputFile.write("root -l -b -q 'GetPhotonReweighting.C+(%s,%s,%s)' \n" %('"ee"','1','"2016"'))
outputFile.write("root -l -b -q 'GetPhotonReweighting.C+(%s,%s,%s)' \n" %('"mm"','1','"2016"'))

outputFile = open(script_dir+'/run_bphys_rw.sh', 'w')
outputFile.write("cd %s \n"%(exec_dir))
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh\n")
outputFile.write("lsetup ROOT\n")
outputFile.write("source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/utilities/oldAliasSetup.sh root\n")
outputFile.write("root -l -b -q 'GetJpsiReweighting.C+(%s,%s,%s)' \n" %('"ee"','0','"2016"'))
outputFile.write("root -l -b -q 'GetJpsiReweighting.C+(%s,%s,%s)' \n" %('"mm"','0','"2016"'))
outputFile.write("root -l -b -q 'GetJpsiReweighting.C+(%s,%s,%s)' \n" %('"ee"','1','"2016"'))
outputFile.write("root -l -b -q 'GetJpsiReweighting.C+(%s,%s,%s)' \n" %('"mm"','1','"2016"'))
outputFile.write("root -l -b -q 'GetJpsiReweighting.C+(%s,%s,%s)' \n" %('"ee"','2','"2016"'))
outputFile.write("root -l -b -q 'GetJpsiReweighting.C+(%s,%s,%s)' \n" %('"mm"','2','"2016"'))

