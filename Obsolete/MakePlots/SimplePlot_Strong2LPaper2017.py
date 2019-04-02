import sys,ROOT
from inspect import currentframe, getframeinfo
import sys, os
sys.path.insert(0, 'python/')

ROOT.gROOT.Reset()
ROOT.gROOT.SetBatch(True)

from setAtlasStyle import *
from histUtils import *

hist_tag = []
histoNames = ''
xLabel = ''
tag = ''
flavor = ''
label = ''

def MakeThePlots():

  doMC = 0
  #doMC = 1
  #doMC = 2
  #doMC = 3 # ee v.s. mm
  doRatio = True
  #doRatio = False
  #isLog = True
  isLog = False
  if 'Hist_MET' in histoNames: isLog = True
  PaperStyle = False

  histoFiles = []
  legendTitles = []
  color_code = []
  if doMC==0:
    histoFiles += ["../OutputHistogram/Data"]
    histoFiles += ["../OutputHistogram/VV"]
    histoFiles += ["../OutputHistogram/Top"]
    #histoFiles += ["../OutputHistogram/DataFS"]
    #histoFiles += ["../OutputHistogram/GData_Syst"]
    histoFiles += ["../OutputHistogram/GData"]
    legendTitles += ['Data']
    legendTitles += ['VV']
    legendTitles += ['Top']
    legendTitles += ['Z+jets(DD)']
    color_code += [1]
    color_code += [46]
    color_code += [30]
    color_code += [38]
  elif doMC==1:
    histoFiles += ["../OutputHistogram/ZMC"]
    #histoFiles += ["../OutputHistogram/GMC_Syst"]
    histoFiles += ["../OutputHistogram/GMC"]
    legendTitles += ['Z/#gamma*+jets MC']
    legendTitles += ['Z/#gamma*+jets (from #gamma+jets MC)']
    color_code += [1]
    color_code += [kAzure+1]
  elif doMC==2:
    histoFiles += ["../OutputHistogram/Data"]
    histoFiles += ["../OutputHistogram/VV"]
    histoFiles += ["../OutputHistogram/Top"]
    histoFiles += ["../OutputHistogram/ZMC"]
    legendTitles += ['Data']
    legendTitles += ['VV']
    legendTitles += ['Top']
    legendTitles += ['Z+jets(MC)']
    color_code += [1]
    color_code += [46]
    color_code += [30]
    color_code += [38]
  elif doMC==3:
    histoFiles += ["../OutputHistogram/ZMC"]
    legendTitles += ['Z+jets (ee)']
    legendTitles += ['Z+jets (mm)']
    color_code += [1]
    color_code += [38]
  signalFiles = []
  #signalFiles += ['../OutputHistogram/Susy_392330']
  signalTitles = []
  #signalTitles += ['m(N2,N1)=(200,100)']
  if doMC==0:
    #saveName = '../OutputPlots/%s_%s_Data'%(tag,histoNames)
    saveName = '/afs/cern.ch/user/r/rshang/workarea/SmallExamples/NtupleTreeReader/Susy2LTree_svn/Susy2LTree/Giulia_plots/output/%s_%s_%s_Data'%(tag,flavor,histoNames)
  elif doMC==1:
    #saveName = '../OutputPlots/%s_%s_MC'%(tag,histoNames)
    saveName = '/afs/cern.ch/user/r/rshang/workarea/SmallExamples/NtupleTreeReader/Susy2LTree_svn/Susy2LTree/Giulia_plots/output/%s_%s_%s_MC'%(tag,flavor,histoNames)
  elif doMC==2:
    #saveName = '../OutputPlots/%s_%s_221'%(tag,histoNames)
    saveName = '/afs/cern.ch/user/r/rshang/workarea/SmallExamples/NtupleTreeReader/Susy2LTree_svn/Susy2LTree/Giulia_plots/output/%s_%s_%s_MC'%(tag,flavor,histoNames)
  elif doMC==3:
    #saveName = '../OutputPlots/%s_%s_eemm'%(tag,histoNames)
    saveName = '/afs/cern.ch/user/r/rshang/workarea/SmallExamples/NtupleTreeReader/Susy2LTree_svn/Susy2LTree/Giulia_plots/output/%s_%s_eemm'%(tag,histoNames)
  histoList = []
  signalList = []
  for h in range(0,len(histoFiles)):
    if 'ee' in flavor or 'mm' in flavor:
        File=ROOT.TFile(histoFiles[h]+'_'+flavor+'_'+tag+'.root')
        print 'Getting %s'%(histoNames)
        histo=File.Get(histoNames)
        histo.SetDirectory(0)
        File.Close()
        histoList += [histo]
    elif 'em' in flavor:
        File=ROOT.TFile(histoFiles[h]+'_'+flavor+'_'+tag+'.root')
        print 'Getting %s'%(histoNames)
        histo=File.Get(histoNames)
        histo.SetDirectory(0)
        histo.Scale(scale_em[h])
        File.Close()
        histoList += [histo]
    else:
        File=ROOT.TFile(histoFiles[h]+'_ee_'+tag+'.root')
        print 'Getting %s'%(histoNames)
        histo=File.Get(histoNames)
        histo.SetDirectory(0)
        File.Close()
        File=ROOT.TFile(histoFiles[h]+'_mm_'+tag+'.root')
        print 'Getting %s'%(histoNames)
        histo.Add(File.Get(histoNames))
        histo.SetDirectory(0)
        File.Close()
        histoList += [histo]
    if doMC==3:
        histoList = []
        File=ROOT.TFile(histoFiles[h]+'_ee_'+tag+'.root')
        print 'Getting %s'%(histoNames)
        histo=File.Get(histoNames)
        histo.SetDirectory(0)
        File.Close()
        histoList += [histo]
        File=ROOT.TFile(histoFiles[h]+'_mm_'+tag+'.root')
        print 'Getting %s'%(histoNames)
        histo=File.Get(histoNames)
        histo.SetDirectory(0)
        File.Close()
        histoList += [histo]
  for s in range(0,len(signalFiles)):
    File=ROOT.TFile(signalFiles[s]+'_'+flavor+'_'+tag+'.root')
    signal=File.Get(histoNames)
    signal.SetDirectory(0)
    File.Close()
    signalList += [signal]
  
  
  # set negative bins to zero
  for h in range(0,len(histoList)):
      for b in range(1,histoList[h].GetNbinsX()+1):
          histoList[h].SetBinContent( b, max(histoList[h].GetBinContent(b),0.) )
  
  dataH=histoList[0]
  
  #dataGraph = ROOT.RooHist(histoList[0],ROOT.RooAbsData.Poisson)
  ##for binx in range(1,dataH.GetNbinsX()+1):
  ##  dataGraph.SetPointEXhigh(binx-1,0)
  ##  dataGraph.SetPointEXlow(binx-1,0)
  #bin_this = 0
  #for binx in range(1,dataH.GetNbinsX()+1):
  #  bin_this += 1
  #  if dataH.GetBinContent(binx)==0: 
  #    dataGraph.RemovePoint(bin_this-1)
  #    bin_this += -1


  signal_color_code = [kBlue+2,kGreen-4]
  signal_style_code = [5,2]
  for s in range(0,len(signalList)):
  	set2L_histStyle_signal( signalList[s], signal_color_code[s], signal_style_code[s] )
  
  # set styles
  for h in range(0,len(histoList)):
  	set2L_histStyle( histoList[h], color_code[h] )
  
  
  # make stack histogram
  stack = ROOT.THStack("stack", "stack")
  if len(histoList)>1:
  	for h in range(1,len(histoList)): stack.Add( histoList[h] )
  else:
  	stack.Add( histoList[0] )
  
  # add photon systematics
  #for h in range(0,len(histoFiles)):
  #    if 'GData' in histoFiles[h]: ResetErrorBand(histoList[h],range_syst,rel_syst)
  
  # make total SM histogram
  if len(histoList)>1:
  	totH = histoList[1].Clone('tot')
  	setICHEP_histStyle_totH(totH)
  	for h in range(2,len(histoList)): totH.Add( histoList[h] )
  	hsave = totH.Clone("hsave")
  else:
  	totH = histoList[0].Clone('tot')
  	hsave = totH.Clone("hsave")
  
  # set style of total SM MC histogram
  setICHEP_histStyle_totSM(totH) 
  hsave.SetLineColor(1)
  hsave.SetLineWidth(3)
  
  # set up the canvas
  if not doRatio:
      c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
      pad1 = ROOT.TPad("pad1","pad1",0,0.0,1,1)
      pad1.SetBottomMargin(0.13)
  else:
      c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 600)
      #c_both = ROOT.TCanvas("c_both","c both", 200, 10, 600, 400)
      pad1 = ROOT.TPad("pad1","pad1",0,0.3,1,1)
      pad1.SetBottomMargin(0.0)
      pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.3)
      pad2.SetBottomMargin(0.39)
      pad2.SetTopMargin(0.0)
      pad2.SetBorderMode(0)
      pad2.Draw()
       
  maxContent=max(totH.GetMaximum()+totH.GetBinError(totH.GetMaximumBin()),dataH.GetMaximum()+dataH.GetBinError(dataH.GetMaximumBin()))
  minContent = min(maxContent*0.0001,0.1)
  
  yMax=2.0 #times maximum bin content
  if isLog: 
  	if minContent>0: yMax = pow(maxContent,0.2)/minContent
  	else: isLog = False
  if isLog:
      pad1.SetLogy()  
  pad1.SetTopMargin(0.03)
  pad1.SetBorderMode(0)
  pad1.Draw()
  pad1.cd()
  
  dataH.Draw("pe")
  #dataGraph.SetLineColor(2)
  #dataGraph.Draw("p")
  
  dataH.SetMaximum( yMax*maxContent)
  #dataH.SetMinimum( 0.8*minContent)
  dataH.SetMinimum( 0.015 )
  #if doRatio: stack.SetMinimum( 0.001 )
  
  dataH.GetXaxis().SetTitle( xLabel )
  dataH.GetYaxis().SetTitle( "Events / 20 GeV")
  
  dataH.GetYaxis().SetTitleOffset(1.3)
  dataH.GetXaxis().SetTitleOffset(1.1)
  dataH.GetXaxis().SetTitleSize(0.064)
  dataH.GetYaxis().SetTitleSize(0.05)
  dataH.GetYaxis().SetLabelSize(0.058)
  dataH.GetXaxis().SetTickLength(0.030)
  dataH.SetMarkerStyle(ROOT.kFullCircle)
  dataH.GetXaxis().SetNdivisions(505)
  
  # make the uncertainty band
  #makeErrorBand(mcList,relUncList,statIncluded,totH)
  
  # draw total background and data
  if len(histoList)>1: 
  	stack.Draw("hist same")
  	totH.Draw("e2 same")
  	hsave.Draw("hist same")
  	dataH.Draw("pe same")
        #dataGraph.Draw("p same")
  for s in range(0,len(signalList)):
      signalList[s].Draw("hist same")
  
  lines={}
  #for i in [150,250]:
  #    lines[i] = ROOT.TLine(i,minContent,i,maxContent)
  #    lines[i].SetLineStyle(2)
  #    lines[i].SetLineColor(2)
  #    lines[i].SetLineWidth(2)
  #    lines[i].Draw()
  
  # make the legend
  legend = ROOT.TLegend(0.55,0.64,0.94,0.89)
  legend.SetBorderSize(0)
  legend.SetTextFont(42)
  if doRatio:
      legend.SetTextSize(0.035)
  else: 
      legend.SetTextSize(0.035)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.SetLineColor(0)
  
  #legend.AddEntry(dataH,legendTitles[0]+", %0.1f#pm%0.1f"%(getYields(dataH,0),getErrors(dataH,0)),"pl")
  #legend.AddEntry(dataH,legendTitles[0]+", %0.1f"%(getYields(dataH,0)),"pl")
  legend.AddEntry(dataH,legendTitles[0],"pl")
  #for h in range(1,len(histoList)): legend.AddEntry(histoList[h],legendTitles[h]+", %0.1f#pm%0.1f"%(getYields(histoList[h],0),getErrors(histoList[h],0)),"f")
  #for h in range(1,len(histoList)): legend.AddEntry(histoList[h],legendTitles[h]+", %0.1f"%(getYields(histoList[h],0)),"f")
  for h in range(1,len(histoList)): legend.AddEntry(histoList[len(histoList)-h],legendTitles[len(histoList)-h],"f")
  for s in range(0,len(signalList)): legend.AddEntry(signalList[s],signalTitles[s],"l")
  legend.Draw("SAME")
  
  # add labels to plot
  lumilab = ROOT.TLatex(0.2,0.84,'#sqrt{s} = 13 TeV, 36.1 fb^{-1}' )
  lumilab.SetNDC()
  lumilab.SetTextSize(0.038)
  lumilab.Draw()
  rs.myText( 0.2,0.78, 1,label, size=0.05)  
  label_2 = ''
  if 'LowMET' in histoNames:
    label_2 += 'E_{T}^{miss} #in [100,200]'
  if not 'mll' in histoNames and not 'MET' in histoNames:
    label_2 += ' m_{ll} #in [81,101]'
  rs.myText( 0.2,0.68, 1,label_2, size=0.05)  
  rs.ATLAS_LABEL(0.2,0.90, tsize=0.06)
  rs.myText( 0.37,0.90, 1,"Simulation Internal", size=0.055)
  
  c_both.Update()
  pad1.cd()
  
  #make ratio plot
  if doRatio:
      ratioH=makeRatio(dataH,totH)
      bandH=makeRelErrorBand(totH)
      mymin=0.25
      mymax=2.2
      pad2.cd()
  
      bandH.SetFillColor(kBlack)
      bandH.SetFillStyle(3354)
      bandH.SetMarkerSize(0)
      
      frameH = ROOT.TH2D("frame",'',100,totH.GetBinLowEdge(1),totH.GetBinLowEdge(totH.GetNbinsX())+totH.GetBinWidth(totH.GetNbinsX()),2,mymin,mymax)
      frameH.GetXaxis().SetTitle( dataH.GetXaxis().GetTitle() )
      frameH.GetYaxis().SetTitle( 'Z MC/#gamma+jets est.' )
      frameH.GetYaxis().SetTitleOffset(0.5)
      frameH.GetXaxis().SetTitleOffset(1.1)
      frameH.GetYaxis().SetTitleSize(0.13)
      frameH.GetYaxis().SetLabelSize(0.13)
      frameH.GetXaxis().SetTitleSize(0.13)
      frameH.GetXaxis().SetLabelSize(0.13)
      frameH.GetXaxis().SetTickLength(0.030)
      frameH.GetXaxis().SetTickLength(0.030)
      frameH.GetXaxis().SetNdivisions(505)
  
      frameH.GetXaxis().SetTitle( xLabel )
      frameH.GetYaxis().SetNdivisions(int(mymax)+4)
      frameH.Draw()
      ratioH.Draw('pesame')
      bandH.Draw('e2same')
  
      lines={}
      #for i in [0.6,0.8,1.0,1.2,1.4]:
      for i in [1.0]:
          lines[i] = ROOT.TLine(totH.GetBinLowEdge(1),i,totH.GetBinLowEdge(totH.GetNbinsX())+totH.GetBinWidth(totH.GetNbinsX()),i)
          if i==1: lines[i].SetLineStyle(1)
          else: lines[i].SetLineStyle(2)
          lines[i].Draw()
   
      marker = ROOT.TArrow()
      marker.SetLineWidth(1)
      marker.SetLineColor(2)
      marker.SetFillColor(2)
      marker.SetAngle(45)
      plrange = mymax - mymin
      arrowOffset = 0.08
      arrowLength = 0.24
      arrowHeadSize = 0.06
      for i in range(0,ratioH.GetNbinsX()):
          if ratioH.GetBinContent(i+1)>mymax:
              marker.DrawArrow(ratioH.GetBinCenter(i+1),mymax - (arrowOffset+arrowLength)*plrange,ratioH.GetBinCenter(i+1),mymax - arrowOffset*plrange,arrowHeadSize*0.35,"|>")
  
      c_both.Update()
  
  pad1.cd()
  hAxis=stack.Clone()
  hAxis.Draw("AXISSAME")
  c_both.Update()
  
  # save
  #c_both.SaveAs(saveName+".eps")
  c_both.SaveAs(saveName+".pdf")
  
  # clean up
  del(c_both)
  if doRatio:
      del(frameH)
  #del(frameH2)
  #result = [corr,getYields(dataH,0)]
  result = 0
  return result




# import histograms
#tag = 'Strong2L_Strong2LSRlow_PTRW_McSmear'
#label = 'SR-low ee+#mu#mu'
#tag = 'Strong2L_Strong2LSRmedium_PTRW_McSmear'
#label = 'SR-medium ee+#mu#mu'
#tag = 'Strong2L_Strong2LSRhigh_PTRW_McSmear'
#label = 'SR-high ee+#mu#mu'
#tag = 'Strong2L_Strong2LVRlow_PTRW_McSmear'
#label = 'VR-low-#Delta#phi ee+#mu#mu'
#tag = 'Strong2L_Strong2LVRmedium_PTRW_McSmear'
#label = 'VR-med-#Delta#phi ee+#mu#mu'
tag = 'Strong2L_Strong2LVRhigh_PTRW_McSmear'
label = 'VR-high-#Delta#phi ee+#mu#mu'
#flavor = 'ee'
#flavor = 'mm'
flavor = 'sf'

histoNames = 'Hist_MET'
xLabel = 'E_{T}^{miss} [GeV]'
result = MakeThePlots()
histoNames = 'Hist_MET_20GeV'
xLabel = 'E_{T}^{miss} [GeV]'
result = MakeThePlots()
histoNames = 'Hist_mll'
xLabel = 'dilepton mass [GeV]'
result = MakeThePlots()
histoNames = 'Hist_DPhiMETJetMin'
xLabel = 'min. #Delta#phi(E_{T}^{miss},jets)'
result = MakeThePlots()
histoNames = 'Hist_nlep'
xLabel = 'N. leptons'
result = MakeThePlots()
histoNames = 'Hist_LowMET_mll'
xLabel = 'dilepton mass [GeV]'
result = MakeThePlots()
histoNames = 'Hist_LowMET_DPhiMETJetMin'
xLabel = 'min. #Delta#phi(E_{T}^{miss},jets)'
result = MakeThePlots()
histoNames = 'Hist_LowMET_nlep'
xLabel = 'N. leptons'
result = MakeThePlots()
