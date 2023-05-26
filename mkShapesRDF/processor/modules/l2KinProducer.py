import ROOT
from mkShapesRDF.processor.framework.Module import Module


class l2KinProducer(Module):
    def __init__(self):
        super().__init__("l2KinProducer")

    def runModule(self, df, values):
        #df = df.Define(
        #    "Lepton_4DV",
        #    "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>"
        #    "(Lepton_pt, Lepton_eta, Lepton_phi, "
        #    "ROOT::RVecF(Lepton_pt.size(), 0))",
        #)

        #df = df.Define(
        #    "CleanJet_4DV",
        #    "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>"
        #    "(CleanJet_pt, CleanJet_eta, CleanJet_phi, "
        #    "Take(Jet_mass, CleanJet_jetIdx))",
        #)

        df = df.Define(
            "MET_4DV", "ROOT::Math::PtEtaPhiMVector" "(PuppiMET_pt, 0, PuppiMET_phi, 0)"
        )

        df = df.Define(
            "TkMET_4DV", "ROOT::Math::PtEtaPhiMVector" "(TkMET_pt, 0, TkMET_phi, 0)"
        )

        df = df.Define("_lepOk", "Lepton_pt[Lepton_pt > 0].size() >= 2 ")
        df = df.Define("_metOk", "MET_4DV.E() > 0")
        df = df.Define("_tkMetOk", "TkMET_4DV.E() > 0")
        df = df.Define("_jetOk", "CleanJet_pt[CleanJet_pt > 0].size() >= 1")

        # FIXME complete l2kin module!
        #df = df.Define(
        #    "mjj", 
        #    "ROOT::VecOps::InvariantMasses"
        #    "(ROOT::RVecF(CleanJet_pt[0]), ROOT::RVecF(CleanJet_eta[0]), ROOT::RVecF(CleanJet_phi[0]), ROOT::RVecF(CleanJet_mass[0]), "
        #    " ROOT::RVecF(CleanJet_pt[1]), ROOT::RVecF(CleanJet_eta[1]), ROOT::RVecF(CleanJet_phi[1]), ROOT::RVecF(CleanJet_mass[1]))"
        #)

        ROOT.gInterpreter.Declare("""
        float compute_mjj(float jetpt1, float jeteta1, float jetphi1, float jetmass1, float jetpt2, float jeteta2, float jetphi2, float jetmass2)
        {
            #include <TLorentzVector.h>

            TLorentzVector J1,J2;
            if (jetpt1 > 0 && jetpt2 > 0){
              J1.SetPtEtaPhiM(jetpt1, jeteta1, jetphi1, jetmass1);
              J2.SetPtEtaPhiM(jetpt2, jeteta2, jetphi2, jetmass2);
              return (J1+J2).M();
            }
            else return -9999.9;
        }
        """)

        df = df.Define(
            "mjj", 
            "compute_mjj(CleanJet_pt[0], CleanJet_eta[0], CleanJet_phi[0], CleanJet_mass[0], CleanJet_pt[1], CleanJet_eta[1], CleanJet_phi[1], CleanJet_mass[1])" 
        )


        #df.DropColumns("Lepton_4DV")
        #df.DropColumns("CleanJet_4DV")
        df.DropColumns("MET_4DV")
        df.DropColumns("TkMET_4DV")

        df.DropColumns("_lepOk")
        df.DropColumns("_metOk")
        df.DropColumns("_tkMetOk")
        df.DropColumns("_jetOk")

        return df
