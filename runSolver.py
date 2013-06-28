#!/bin/env python

from array import array
import collections
import sys
from math import sin, cos
import numpy 
LinAlgError = numpy.linalg.linalg.LinAlgError
import time
import ROOT as r
r.gROOT.LoadMacro('vecUtils.h'+'+')
lv = r.Math.LorentzVector(r.Math.PtEtaPhiE4D('float'))
tlv = r.TLorentzVector
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

    #nEntries = 10

nEventsPerNsol_Betchart = collections.defaultdict(int)
nEventsPerNsol_Sonnensc = collections.defaultdict(int)
nCallsB, timeCallsB = 0, 0
nCallsS, timeCallsS = 0, 0

for iEntry in xrange(nEntries) :
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
    t1 = time.time()
    try :
        dns = doubleNeutrinoSolutions((j0, j1), (l0, l1), (metx, mety))
        #print solutions.n_
        solutions = dns.nunu_s
        if verbose : print "found %d solutions"%len(solutions)
        for sol in solutions :
            n0, n1 = sol[0], sol[1]
            nu, nu_ = tlv(), tlv()
            nu.SetXYZM(n0[0], n0[1], n0[2], 0.0)
            nu_.SetXYZM(n1[0], n1[1], n1[2], 0.0)
            totMet = nu + nu_
            if verbose :
                print "met: magnitude %.1f, phi %.2f"%(totMet.Pt(), totMet.Phi())
                print "          reco %.1f,     %.2f"%(metP, metPhi)
        nEventsPerNsol_Betchart[len(sol)] += 1
    except LinAlgError :
        nEventsPerNsol_Betchart[0] += 1
        #print "skipping singular matrix"
    t2 = time.time()
    nCallsB += 1
    timeCallsB += (t2-t1)
    
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
    nEventsPerNsol_Sonnensc[int(pnux.size())] += 1
    if verbose : print "solver: %d solutions"%pnux.size()
    t2 = time.time()
    nCallsS += 1
    timeCallsS += (t2-t1)
    
print "number of solutions Betchart     : ", nEventsPerNsol_Betchart
print "time : ",(timeCallsB/nCallsB)
print "number of solutions Sonnenschein : ", nEventsPerNsol_Sonnensc
print "time : ",(timeCallsS/nCallsS)
