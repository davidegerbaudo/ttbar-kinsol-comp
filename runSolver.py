#!/bin/env python

from array import array
import collections
import sys
from math import sin, cos, pi
import numpy 
LinAlgError = numpy.linalg.linalg.LinAlgError
import time
import ROOT as r
r.gROOT.SetBatch(1)
r.gROOT.LoadMacro('vecUtils.h'+'+')
lv = r.Math.LorentzVector(r.Math.PtEtaPhiE4D('float'))
tlv = r.TLorentzVector
DeltaPhi = r.Math.VectorUtil.DeltaPhi
r.gROOT.LoadMacro('dileptonSolver/DileptonAnalyticalSolver.cc+')
solver = r.llsolver.DileptonAnalyticalSolver()

sys.path.append('./analytic-nu')
import nuSolutions as ns
doubleNeutrinoSolutions = ns.doubleNeutrinoSolutions

verbose = False
inputFileName = sys.argv[1]
inputFile = r.TFile.Open(inputFileName)
treeName = 'ToyNt'
inputTree = inputFile.Get(treeName)

nEntries = inputTree.GetEntries()
if verbose : print "%d entries"%nEntries

#nEntries = 10 # just for a quick test

nEventsPerNsol_Betchart = collections.defaultdict(int)
nEventsPerNsol_Sonnensc = collections.defaultdict(int)
nCallsB, timeCallsB = 0, 0
nCallsS, timeCallsS = 0, 0

maxNsol = 8
hnamePre, htitlePre = 'h_Dphi_', 'N solutions vs. #Delta#phi'
titleX, titleY = '#Delta#phi(MET_{reco}, MET_{fit})', 'N of kinematic solutions'
h_Dphi_B = r.TH2D(hnamePre+'B', htitlePre+' Betchart;'+titleX+';'+titleY,
                  50, 0.0, pi, maxNsol+1, -0.5, maxNsol+0.5)
h_Dphi_S = r.TH2D(hnamePre+'S', htitlePre+' Sonnenschein;'+titleX+';'+titleY,
                  50, 0.0, pi, maxNsol+1, -0.5, maxNsol+0.5)

def phi_0_twopi(phi) :
    while phi <  0.0    : phi = phi+2.0*pi
    while phi >= 2.0*pi : phi = phi-2.0*pi
    return phi
def lv2str(l) :
    return "Pt : %f Eta : %f Phi : %f E : %f"%(l.Pt(), l.Eta(), l.Phi(), l.E())
for iEntry in xrange(9100, 9104): #xrange(nEntries) :
    inputTree.GetEntry(iEntry)
    it = inputTree
    lepsPt, lepsEta, lepsPhi, lepsE = it.l_pt, it.l_eta, it.l_phi, it.l_e
    jetsPt, jetsEta, jetsPhi, jetsE = it.j_pt, it.j_eta, it.j_phi, it.j_e
    if len(jetsPt) < 2 : continue # cannot solve without two jets
    metP, metPhi = it.met, it.met_phi
    l0q, l1q = it.l_q[0], it.l_q[1]
    l0 = lv(lepsPt[0], lepsEta[0], lepsPhi[0], lepsE[0])
    l1 = lv(lepsPt[1], lepsEta[1], lepsPhi[1], lepsE[1])
    j0 = lv(jetsPt[0], jetsEta[0], jetsPhi[0], jetsE[0])
    j1 = lv(jetsPt[1], jetsEta[1], jetsPhi[1], jetsE[1])
    metx, mety = metP*cos(metPhi), metP*sin(metPhi)
    metRec =  tlv()
    metRec.SetXYZM(metx, mety, 0.0, 0.0)
    metRec = lv(metRec.Pt(), metRec.Eta(), metRec.Phi(), metRec.E())
    try :
        t1 = time.time()
        dns = doubleNeutrinoSolutions((j0, j1), (l0, l1), (metx, mety))
        solutions = dns.nunu_s
        t2 = time.time()
        nCallsB += 1
        timeCallsB += (t2-t1)
        nSolB = len(solutions)
        if verbose or nSolB==6 :
            print "entry[%d] found %d solutions"%(iEntry, len(solutions))
            print '-'*4+" input "+'-'*4
            for p in ['l0', 'l1', 'j0', 'j1', 'metRec'] :
                print "%6s"%p,': ',lv2str(eval(p))
            for i,ss in enumerate(dns.solutionSets) :
                print "solutionSets[%d].Z %f"%(i,ss.Z)

        for i,sol in enumerate(solutions) :
            n0, n1 = sol[0], sol[1]
            nu, nu_ = tlv(), tlv()
            nu.SetXYZM(n0[0], n0[1], n0[2], 0.0)
            nu_.SetXYZM(n1[0], n1[1], n1[2], 0.0)
            metFit = nu + nu_
            metFit = lv(metFit.Pt(), metFit.Eta(), metFit.Phi(), metFit.E())
            h_Dphi_B.Fill(phi_0_twopi(DeltaPhi(metFit, metRec)), nSolB)
            if verbose or nSolB==6:
                print '-'*4+" solution %d "%i+'-'*4
                print "n0: ",n0
                print "n1: ",n1
                print "met fit: Pt %.1f, phi %.2f"%(metFit.Pt(), metFit.Phi())
                print "met rec  Pt %.1f, phi %.2f"%(metRec.Pt(), metRec.Phi())
        nEventsPerNsol_Betchart[nSolB] += 1
    except LinAlgError :
        nEventsPerNsol_Betchart[0] += 1
        #print "skipping singular matrix"
    
    lp = l0 if l0q>0 else l1
    ln = l0 if l1q<0 else l1
    
    lp = array('d',[lp.E(), lp.px(), lp.py(), lp.pz()])
    lm = array('d',[ln.E(), ln.px(), ln.py(), ln.pz()])
    b  = array('d',[j0.E(), j0.px(), j0.py(), j0.pz()])
    bb = array('d',[j1.E(), j1.px(), j1.py(), j1.pz()])
    ETmiss = array('d', [metx, mety])
    nu     = array('d',4*[0.])
    nub    = array('d',4*[0.])

    pnux = r.std.vector('double')()
    pnuy = r.std.vector('double')()
    pnuz = r.std.vector('double')()
    pnubx = r.std.vector('double')()
    pnuby = r.std.vector('double')()
    pnubz = r.std.vector('double')()
    mt = mtb = 175.
    mWp = mWm = 80.41
    mnu = mnub = 0.0
    cd_diff = r.std.vector('double')()
    cubic_single_root_cmplx = r.Long()

    t1 = time.time()    
    solver.solve(ETmiss, b, bb, lp, lm, mWp, mWm, mt, mtb, mnu, mnub,
                 pnux, pnuy, pnuz, pnubx, pnuby, pnubz,
                 cd_diff, cubic_single_root_cmplx)
    t2 = time.time()
    nCallsS += 1
    timeCallsS += (t2-t1)
    nEventsPerNsol_Sonnensc[int(pnux.size())] += 1
    nSolS = pnux.size()
    if verbose : print "solver: %d solutions"%nSolS
    for nx, ny, nz, Nx, Ny, Nz in zip(pnux, pnuy, pnuz, pnubx, pnuby, pnubz) :
        nu, nu_ = tlv(), tlv()
        nu.SetXYZM( nx, ny, nz, 0.0)
        nu_.SetXYZM(Nx, Ny, Nz, 0.0)
        metFit = nu + nu_
        metFit = lv(metFit.Pt(), metFit.Eta(), metFit.Phi(), metFit.E())
        h_Dphi_S.Fill(phi_0_twopi(DeltaPhi(metFit, metRec)), nSolS)
    

print "number of solutions Betchart     : ", nEventsPerNsol_Betchart
print "time : ",(timeCallsB/nCallsB)
print "number of solutions Sonnenschein : ", nEventsPerNsol_Sonnensc
print "time : ",(timeCallsS/nCallsS)

c = r.TCanvas('c_solver_comparison','solver comparison')
opt = 'colz'
c.Divide(2)
c.cd(1)
h_Dphi_B.Draw(opt)
c.cd(2)
h_Dphi_S.Draw(opt)
for ext in ['png'] : c.SaveAs(c.GetName()+'.'+ext)
