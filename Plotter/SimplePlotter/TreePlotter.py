from PlotterTools import *
from ROOT import *
gROOT.SetBatch(1)
gStyle.SetOptStat(0)

for v in Vars:
   hists = []
   leg = TLegend(0.6, 0.7, 0.89, 0.89)
   leg.SetBorderSize(0)
   Max = -0.
   for fi,f in enumerate(files):
      ch = TChain()
      ch.Add(f[0])
      hname = v[1]+'_'+str(fi)
      h = TH1F(hname, v[2], v[3], v[4], v[5])

      ch.Draw(v[0]+'>>'+hname,TCut(Cut))
      total = h.Integral()
      h.Scale(float(f[4])/float(total))
      h.SetLineColor(f[2])
      h.SetLineWidth(2)
      print h.GetMaximum()
      h.SetMaximum(h.GetMaximum()+0.8*h.GetMaximum())
      h.GetYaxis().SetTitleOffset(1.5)

      hists.append([h,f[1]])

   c0 = TCanvas('a', 'a', 800, 600)
   h.Draw()
   for fi,hh in enumerate(hists):
      leg.AddEntry(hh[0],hh[1], 'lf')
      leg.SetTextSize(0.02)
      if fi == 0:
         hh[0].Draw('hist')
      if fi > 0:
         hh[0].Draw('hist same')

   leg.SetFillStyle(0)
   leg.SetTextSize(0.02)
   leg.Draw('same')

   c0.Update()
   c0.SaveAs(outputLoc+v[1]+'.pdf')
   c0.SaveAs(outputLoc+v[1]+'.png')
   c0.SetLogy()
   c0.SaveAs(outputLoc+v[1]+'.pdf')
   c0.SaveAs(outputLoc+v[1]+'_log.png')
