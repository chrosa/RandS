#!/usr/bin/python

from ROOT import TH1F, TFile, TCanvas, TGraphAsymmErrors

filename = "out_noJVT.txt"
fobj_run = open(filename)
line = fobj_run.readline()

MHT_reco = TH1F("MHT_reco","MHT_reco",50,0,500)
MHT_reco.Sumw2();
MHT_reco_triggered = TH1F("MHT_reco","MHT_reco",50,0,500)
MHT_reco_triggered.Sumw2();

while line:

	MHT = float(line.split(" ")[4])
	triggered = int(line.split(" ")[5])
	hasBadJet = int(line.split(" ")[6])
	VBFremoved = int(line.split(" ")[7])
	
	if ( hasBadJet == 0 ):
		print MHT, triggered, hasBadJet, VBFremoved
		MHT_reco.Fill(MHT)
		if ( triggered == 1 ):
			MHT_reco_triggered.Fill(MHT)
	
	line = fobj_run.readline()
	
eff = TGraphAsymmErrors();
eff.Divide( MHT_reco_triggered, MHT_reco );

#eff = TH1F(MHT_reco_triggered)
#eff.Divide(MHT_reco)

c = TCanvas("c", "c", 800, 800 )
c.cd()

eff.GetHistogram().SetTitle("")
eff.GetHistogram().SetXTitle("MHT (no JVT)")
eff.GetHistogram().SetYTitle("Efficiency")

eff.GetHistogram().SetStats(0)

eff.SetMarkerStyle(20)
eff.SetMarkerColor(1)
eff.SetMarkerSize(1)

eff.SetLineColor(1)

eff.SetMaximum(1.1)
eff.Draw("AP")
c.Update()
c.SaveAs("TurnOn_MHT_noJVT.pdf")

