from math import *
from ROOT import *

gROOT.SetBatch(True)
gStyle.SetOptStat(0)

blind = True
setLog = True



texRegionNameLine = "textbf"


class Region:
    def __init__(self, name):
        self.name=name

        self.obs = 0
        self.exp = 0
        self.unc = 0
        self.sig = 0
        self.bkgs = dict()


def getSig(reg):
    if reg.sig !=0: return reg.sig
    if "sr" in reg.name.lower() and blind: return 0

    #return (reg.obs - reg.exp) / reg.unc

    nbObs = reg.obs
    nbExp = reg.exp
    nbExpEr = reg.unc
    
    factor1 = nbObs and nbObs*log( (nbObs*(nbExp+nbExpEr**2))/(nbExp**2+nbObs*nbExpEr**2) )
    factor2 = (nbExp**2/nbExpEr**2)*log( 1 + (nbExpEr**2*(nbObs-nbExp))/(nbExp*(nbExp+nbExpEr**2)) )
    sig = sqrt(2*max(0., factor1 - factor2))
    if nbObs < nbExp:
        sig = -sig

    return sig


    
    
        
regions = ["VRZ", "VRZ_MET0_50", "VRZ_MET50_100", "VRZ_MET100_150", "VRZ_MET150_200", "VRC", "VRLow", "VRMed", "VRHigh", "VRLowZ", "VRMedZ", "VRHighZ", "SRC", "SRLow", "SRMed", "SRHigh", "SRLowZ", "SRMedZ", "SRHighZ"]
bkgOrder = ['Zjets', 'other']
bkgColors = {'other':TColor.GetColor("#FFFB00"), 'Zjets':TColor.GetColor("#41ab5d")}

yielddict = {
    "nobs": [3072267, 2760892, 282322, 21134, 4940, 301, 325, 362, 85, 62, 54, 19],
    "TOTAL_FITTED_bkg_events": [3064670.0, 2758693.8, 279187.5, 21419.6, 4886.6, 400.0, 327.0, 405.2, 99.6, 73.2, 66.3, 21.6, 39.3, 128.8, 94.9, 36.9, 33.8, 18.9, 8.1],
    "TOTAL_FITTED_bkg_events_err": [3919.1, 3791.2, 1016.0, 169.8, 37.0, 37.3, 5.9, 5.3, 1.9, 1.8, 1.8, 0.5, 1.1, 1.7, 1.7, 1.1, 0.6, 0.5, 0.3],
    "Fitted_events_Zjets": [2918632.4, 2676506.5, 237585.5, 6828.4, 329.9, 82.5, 40.7, 28.1, 15.4, 12.5, 11.8, 8.0, 0.7, 4.6, 2.3, 1.9, 2.4, 1.4, 0.7],
    "Fitted_events_other": [146037.6, 82187.3, 41602.0, 14591.2, 4556.7, 317.5, 286.3, 377.1, 84.2, 60.7, 54.5, 13.6, 38.6, 124.2, 92.6, 35.0, 31.4, 17.4, 7.4],
}

resultsDict = dict()

for i,reg in enumerate(regions):
    newReg = Region(reg)

    for key,val in yielddict.items():
        if key == "nobs":
            if blind and "sr" in reg.lower():
                newReg.obs = 0
            else:
                newReg.obs = val[i]
        if key == "TOTAL_FITTED_bkg_events_err":
            newReg.unc = val[i]
            
        if key == "TOTAL_FITTED_bkg_events":
            newReg.exp = val[i]
            
        if "Fitted_events" in key:
            newReg.bkgs[key.replace("Fitted_events_","")] = val[i]
            
    resultsDict[reg] = newReg



nbins = len(regions)

Bkgs = resultsDict[regions[0]].bkgs.keys()
nBkgs = len(resultsDict[regions[0]].bkgs)

print(Bkgs)


# make all the histograms
hist_dict = dict()
for b in Bkgs + ["exp","obs","sig"]:
    hist_dict[b]=TH1F(b,b,nbins,0,nbins)
    hist_dict[b].SetTitle("")
    hist_dict[b].SetLineWidth(1)
    hist_dict[b].SetLineColor(1)
    
    if b not in ["exp","obs","sig"]:
        hist_dict[b].SetFillColor(bkgColors[b])

    for i,r in enumerate(regions):
        ibin = i+1
        reg = resultsDict[r]
        if b == "exp":
            hist_dict[b].SetBinContent(ibin, reg.exp)
            hist_dict[b].SetBinError(ibin, reg.unc)
            #hist_dict[b].GetXaxis().SetBinLabel(ibin, r)            
        elif b=="obs":
            hist_dict[b].SetBinContent(ibin, reg.obs)
            hist_dict[b].SetBinError(ibin, 0.)
            #hist_dict[b].GetXaxis().SetBinLabel(ibin, r)            
        elif b=="sig":
            hist_dict[b].SetBinContent(ibin, getSig(reg) )
            hist_dict[b].SetBinError(ibin, 0.)
            hist_dict[b].GetXaxis().SetBinLabel(ibin, r)            
        else:
            hist_dict[b].SetBinContent(ibin, reg.bkgs[b])
            hist_dict[b].SetBinError(ibin, 0.)
            #hist_dict[b].GetXaxis().SetBinLabel(ibin, r)
            

# stack of backgrounds
stack = THStack("stack","stack")
leg = TLegend(0.51,0.53,0.88,0.88)
leg.SetNColumns(2)
leg.SetBorderSize(0)

for b in bkgOrder:
    stack.Add(hist_dict[b])
    leg.AddEntry(hist_dict[b], b, "f")


# tgraph for errrs
tge = TGraphErrors(hist_dict["exp"])
tge.SetFillStyle(3254)
tge.SetFillColor(kGray+2)


    
hist_dict["obs"].SetMarkerSize(1.0)
hist_dict["obs"].SetMarkerStyle(20)
hist_dict["obs"].SetLineColor(1)
hist_dict["obs"].SetBinErrorOption(TH1.kPoisson)
hist_dict["obs"].GetYaxis().SetTitle("Events")
hist_dict["obs"].GetYaxis().SetTitleSize(0.07)
hist_dict["obs"].GetYaxis().SetLabelSize(0.07)
hist_dict["obs"].GetYaxis().SetTitleOffset(0.7)

hist_dict["obs"].SetMinimum( 1 if setLog else 0 )
hist_dict["obs"].SetMaximum( 30000000 if setLog else hist_dict["obs"].GetMaximum() * 1.5   ) 
# do some y axis scaling todo

leg.AddEntry(hist_dict["obs"], "Data", "ep")

hist_dict["exp"].SetFillColor(2)
hist_dict["exp"].SetLineColor(1)
hist_dict["exp"].SetFillColor(kGray+2)
hist_dict["exp"].SetFillStyle(3254)
hist_dict["exp"].SetMarkerSize(0)
leg.AddEntry(hist_dict["exp"], "Standard Model (SM)", "lf")

c1 = TCanvas("c1","c1",800,600)

upperPad = TPad("upperPad","upperPad",0.001,0.4,0.995,0.995)
lowerPad = TPad("lowerPad","lowerPad",0.001,0.001,0.995,0.405)

upperPad.SetFillColor(0)
upperPad.SetBorderMode(0)
upperPad.SetBorderSize(2)
upperPad.SetTicks()
upperPad.SetTopMargin   ( 0.1 )
upperPad.SetRightMargin ( 0.1 )
upperPad.SetBottomMargin( 0.0025 )
upperPad.SetLeftMargin( 0.10 )
upperPad.SetFrameBorderMode(0)
upperPad.SetFrameBorderMode(0)
if( setLog ):  upperPad.SetLogy(1)
upperPad.Draw()

lowerPad.SetGridy()
lowerPad.SetFillColor(0)
lowerPad.SetBorderMode(0)
lowerPad.SetBorderSize(2)
lowerPad.SetTickx(1)
lowerPad.SetTicky(1)
lowerPad.SetTopMargin   ( 0.003 )
lowerPad.SetRightMargin ( 0.1 )
lowerPad.SetBottomMargin( 0.6 )
lowerPad.SetLeftMargin( 0.10 )
lowerPad.Draw()

c1.SetFrameFillColor(kWhite)




# draw top stuff
upperPad.cd()

hist_dict["obs"].Draw()
stack.Draw("same")
tge.Draw("2")
hist_dict["obs"].Draw("esame")


atlasLabel = TLatex()
atlasLabel.SetNDC()
atlasLabel.SetTextFont( 72 )
atlasLabel.SetTextColor( 1 )
atlasLabel.SetTextSize( 0.09 )

tex = TLatex()
tex.SetNDC()
tex.SetTextSize( 0.07 )

atlasLabel.DrawLatex(0.15,0.8, "ATLAS")
atlasLabel.SetTextFont(42)
atlasLabel.DrawLatex(0.3,0.8, "Internal")
tex.DrawLatex(0.15,0.7,"#sqrt{s}=13 TeV, 139.0 fb^{-1}")


leg.Draw()

# draw lower stuff
lowerPad.cd()

hist_dict["sig"].SetFillColor(kOrange)
hist_dict["sig"].SetLineColor(kOrange)

hist_dict["sig"].SetBarWidth(0.8)
hist_dict["sig"].SetBarOffset(0.1)
hist_dict["sig"].GetYaxis().SetNdivisions(5)
hist_dict["sig"].GetXaxis().LabelsOption("v")
hist_dict["sig"].GetXaxis().SetLabelSize(0.1)
hist_dict["sig"].GetXaxis().SetLabelOffset(0.02)
hist_dict["sig"].GetYaxis().SetLabelSize(0.1)

hist_dict["sig"].GetYaxis().SetTitleSize(0.11)
hist_dict["sig"].GetYaxis().SetTitleOffset(0.4)

hist_dict["sig"].GetYaxis().SetTitle("Significance   ")
#hist_dict["sig"].GetYaxis().SetTitle("(Obs - Exp)/Unc")
hist_dict["sig"].GetXaxis().SetTitle("Region")
hist_dict["sig"].GetXaxis().SetTitleSize(0.12)
hist_dict["sig"].GetXaxis().SetTitleOffset(2)
hist_dict["sig"].SetMinimum(-2)
hist_dict["sig"].SetMaximum(2)
hist_dict["sig"].Draw("hb")





c1.Print("plots/summaryPlot.pdf")
