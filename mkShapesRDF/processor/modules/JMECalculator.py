import ROOT
from mkShapesRDF.processor.framework.Module import Module


class JMECalculator(Module):
    def __init__(self, JEC_era, JER_era, p_object, do_JER=True):
        super().__init__("JMECalculator")
        self.JEC_era = JEC_era
        self.JER_era = JER_era
        self.p_object = p_object
        self.do_JER = do_JER

    def runModule(self, df, values):
        from CMSJMECalculators.jetdatabasecache import JetDatabaseCache

        jecDBCache = JetDatabaseCache("JECDatabase", repository="cms-jet/JECDatabase")
        jrDBCache = JetDatabaseCache("JRDatabase", repository="cms-jet/JRDatabase")
        # usage example, returns the local path
        JEC_era = self.JEC_era
        JER_era = self.JER_era
        p_object = self.p_object

        txtJECs = []
        txtJECs.append(jecDBCache.getPayload(JEC_era, "L1FastJet", p_object))
        txtJECs.append(jecDBCache.getPayload(JEC_era, "L2Relative", p_object))

        txtUnc = jecDBCache.getPayload(
            JEC_era, "UncertaintySources", p_object, "RegroupedV2_"
        )
        # print(txtUnc)
        txtPtRes = jrDBCache.getPayload(JER_era, "PtResolution", p_object)
        txtSF = jrDBCache.getPayload(JER_era, "SF", p_object)
        print("Path for SF:", txtSF)
        # print(pl)

        from CMSJMECalculators import loadJMESystematicsCalculators

        loadJMESystematicsCalculators()
        ROOT.gInterpreter.ProcessLine("JetVariationsCalculator myJetVarCalc{}")
        calc = getattr(ROOT, "myJetVarCalc")

        ROOT.gInterpreter.ProcessLine("Type1METVariationsCalculator myMETVarCalc{}")
        metcalc = getattr(ROOT, "myMETVarCalc")


        # redo JEC, push_back corrector parameters for different levels
        jecParams = getattr(ROOT, "std::vector<JetCorrectorParameters>")()
        for txtJEC in txtJECs:
            jecParams.push_back(ROOT.JetCorrectorParameters(txtJEC))
        calc.setJEC(jecParams)
        metcalc.setJEC(jecParams)

        jecL1Params = getattr(ROOT, "std::vector<JetCorrectorParameters>")()
        txtL1JECs = []
        txtL1JECs.append(jecDBCache.getPayload(JEC_era, "L1FastJet", p_object))
        for txtJEC in txtL1JECs:
            jecL1Params.push_back(ROOT.JetCorrectorParameters(txtJEC))
        metcalc.setL1JEC(jecL1Params)

        metcalc.setUnclusteredEnergyTreshold(20.)
        metcalc.setIsT1SmearedMET(True)

        # calculate JES uncertainties (repeat for all sources)

        # uncert_sources = ['Total']
        with open(txtUnc) as f:
            lines = f.read().split("\n")
            sources = [x for x in lines if x.startswith("[") and x.endswith("]")]
            sources = [x[1:-1] for x in sources]
            # sources = list(filter(lambda source: source in , sources))
        # print(sources)

        for s in sources:
            jcp_unc = ROOT.JetCorrectorParameters(txtUnc, s)
            calc.addJESUncertainty(s, jcp_unc)
            metcalc.addJESUncertainty(s, jcp_unc)

        if self.do_JER:
            # Smear jets, with JER uncertainty
            calc.setSmearing(
                txtPtRes,
                txtSF,
                True,
                True,
                0.2,
                3.0,  # decorrelate for different regions
            )  # use hybrid recipe, matching parameters

                        # Smear jets, with JER uncertainty
            metcalc.setSmearing(
                txtPtRes,
                txtSF,
                True,
                True,
                0.2,
                3.0,  # decorrelate for different regions
            )  # use hybrid recipe, matching parameters

        cols = []
        # reco jet coll
        cols.append("CleanJet_pt")
        cols.append("CleanJet_eta")
        cols.append("CleanJet_phi")
        cols.append("Take(Jet_mass, CleanJet_jetIdx)")
        cols.append("Take(Jet_rawFactor, CleanJet_jetIdx)")
        cols.append("Take(Jet_area, CleanJet_jetIdx)")
        cols.append("Take(Jet_jetId, CleanJet_jetIdx)")
        # rho
        cols.append("fixedGridRhoFastjetAll")

        cols.append("Take(Jet_partonFlavour, CleanJet_jetIdx)")

        # seed
        cols.append(
            "(run<<20) + (luminosityBlock<<10) + event + 1 + int(Jet_eta[0]/.01)"
        )

        # gen jet coll
        cols.append("GenJet_pt")
        cols.append("GenJet_eta")
        cols.append("GenJet_phi")
        cols.append("GenJet_mass")

        df = df.Define("jetVars", f'myJetVarCalc.produce({", ".join(cols)})')


        metcols = []
        metcols.append("CleanJet_pt")
        metcols.append("CleanJet_eta")
        metcols.append("CleanJet_phi")
        metcols.append("Take(Jet_mass, CleanJet_jetIdx)")
        metcols.append("Take(Jet_rawFactor, CleanJet_jetIdx)")
        metcols.append("Take(Jet_area, CleanJet_jetIdx)")
        metcols.append("Take(Jet_muonSubtrFactor, CleanJet_jetIdx)")
        metcols.append("Take(Jet_neEmEF, CleanJet_jetIdx)")
        metcols.append("Take(Jet_chEmEF, CleanJet_jetIdx)")
        metcols.append("Take(Jet_jetId, CleanJet_jetIdx)")
        metcols.append("fixedGridRhoFastjetAll")

        metcols.append("Take(Jet_partonFlavour, CleanJet_jetIdx)")

        # seed
        metcols.append(
            "(run<<20) + (luminosityBlock<<10) + event + 1 + int(Jet_eta[0]/.01)"
        )

        # gen jet coll
        metcols.append("GenJet_pt")
        metcols.append("GenJet_eta")
        metcols.append("GenJet_phi")
        metcols.append("GenJet_mass")

        metcols.append("RawMET_phi")
        metcols.append("RawMET_pt")
        metcols.append("MET_MetUnclustEnUpDeltaX")
        metcols.append("MET_MetUnclustEnUpDeltaY")
        metcols.append("Take(CorrT1METJet_rawPt, CleanJet_jetIdx)")
        metcols.append("Take(CorrT1METJet_eta, CleanJet_jetIdx)")
        metcols.append("Take(CorrT1METJet_phi, CleanJet_jetIdx)")
        metcols.append("Take(CorrT1METJet_area, CleanJet_jetIdx)")
        metcols.append("Take(CorrT1METJet_muonSubtrFactor, CleanJet_jetIdx)")
        metcols.append("ROOT::VecOps::RVec<float>()")
        metcols.append("ROOT::VecOps::RVec<float>()")

        df = df.Define("metVars", f'myMETVarCalc.produce({", ".join(metcols)})')



        # resortCols = getColumnsNamesRegex(df, 'Jet_*')

        df = df.Define("CleanJet_pt_raw", "CleanJet_pt")
        df = df.Define("CleanJet_pt", "jetVars.pt(0)")

        df = df.Define(
            "CleanJet_sorting",
            "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(CleanJet_pt))",
        )

        df = df.Redefine("CleanJet_pt", "Take(CleanJet_pt, CleanJet_sorting)")

        # df = df.Redefine('CleanJet_pt', 'Jet_pt_nominal')

        resortCols = ["CleanJet_" + prop for prop in ["pt", "eta", "phi", "jetIdx"]]
        print(resortCols)
        for col in resortCols:
            df.Redefine(col, f"Take({col}, CleanJet_sorting)")

        df = df.Define("Jet_mass_raw", "Jet_mass")

        ROOT.gInterpreter.Declare(
            """
            using namespace ROOT;
            RVecF propagateVector(RVecI jetIdx, RVecF jetVar, RVecF jetVar_raw) {
                RVecF out(jetVar_raw);
                for (int i = 0; i < jetIdx.size(); i++) {
                    out[jetIdx[i]] = jetVar[i];
                }
                return out;
            }
            """
        )
        df = df.Define(
            "Jet_mass",
            "propagateVector(CleanJet_jetIdx, Take(jetVars.mass(0), CleanJet_sorting), Jet_mass_raw)",
        )

        sources = list(map(lambda k: "JES_" + k, sources))
        _sources = []
        if self.do_JER:
            _sources = [f"JER_{i}" for i in range(6)]
        _sources += sources
        sources = _sources.copy()


        for variable in ["CleanJet_pt", "Jet_mass"]:
            for i, source in enumerate(sources):
                up = f"Take(jetVars.{variable.split('_')[-1]}({2*i+1}), CleanJet_sorting)"
                do = f"Take(jetVars.{variable.split('_')[-1]}({2*i+1+1}), CleanJet_sorting)"
                df = df.Vary(variable, "ROOT::RVec<ROOT::RVecF>{" + up + ", " + do + "}", ["up", "down"], source) 

        for variable in ["PuppiMET_pt"]:
            for i, source in enumerate(sources):
                up = f"metVars.{variable.split('_')[-1]}({2*i+1})"
                do = f"metVars.{variable.split('_')[-1]}({2*i+1+1})"
                df = df.Vary(variable, "ROOT::RVec<ROOT::RVecF>{" + up + ", " + do + "}", ["up", "down"], source)

        df.DropColumns("jetVars")
        df.DropColumns("metVars")
        df.DropColumns("CleanJet_sorting")

        return df
