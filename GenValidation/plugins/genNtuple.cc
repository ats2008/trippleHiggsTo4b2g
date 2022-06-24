#include <memory>
/***************************************************************/
/********************HHbbggEvent selection**********************/
/***************************************************************/
/************************TIFR August20202***********************/



// class declaration
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"


#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "TLorentzVector.h"
#include "FWCore/Framework/interface/Run.h"

#include "TTree.h"
#include "TH1.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include <fastjet/GhostedAreaSpec.hh>
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/Selector.hh"

#define N_ITEMS_MAX 100

using namespace fastjet;
using namespace std;
//using namespace fastjet::contrib;
const Int_t kMaxVertices = 300;
const Int_t kMaxGenJet = 300;
/*
//corrections for xy-corrections//
void Corr_met(float met, float met_phi, int npv, float* corr_met, float* corr_met_phi, float* corr_met_x, float* corr_met_y)
{

  double METxcorr = -(0.296713*npv -0.141506);
  double METycorr = -(0.115685*npv +0.0128193);
   
  double CorrectedMET_x = met*cos(met_phi)+METxcorr;
  double CorrectedMET_y = met*sin(met_phi)+METycorr;

  double CorrectedMET = sqrt(CorrectedMET_x*CorrectedMET_x+CorrectedMET_y*CorrectedMET_y);
  double CorrectedMETPhi;
  if(CorrectedMET_x==0 && CorrectedMET_y>0) CorrectedMETPhi = TMath::Pi();
  else if(CorrectedMET_x==0 && CorrectedMET_y<0 )CorrectedMETPhi = -TMath::Pi();
  else if(CorrectedMET_x >0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x);
  else if(CorrectedMET_x <0&& CorrectedMET_y>0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) + TMath::Pi();
  else if(CorrectedMET_x <0&& CorrectedMET_y<0) CorrectedMETPhi = TMath::ATan(CorrectedMET_y/CorrectedMET_x) - TMath::Pi();
  else CorrectedMETPhi =0;

  *corr_met= CorrectedMET;
  *corr_met_phi = CorrectedMETPhi;
  *corr_met_x = CorrectedMET_x;
  *corr_met_y = CorrectedMET_y;
}
*/
//redefinition of delta_phi variable between -pi < delta_phi < pi
double PhiInRange(const double& phi) {
  double phiout = phi;

  if( phiout > 2*M_PI || phiout < -2*M_PI) {
    phiout = fmod( phiout, 2*M_PI);
  }
  if (phiout <= -M_PI) phiout += 2*M_PI;
  else if (phiout >  M_PI) phiout -= 2*M_PI;

  return phiout;
}
bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
//particle is already the ancestor
        if(ancestor == particle ) return true;

//otherwise loop on mothers, if any and return true if the ancestor is found
        for(size_t i=0;i< particle->numberOfMothers();i++)
        {   
                if(particle->mother(i)->pdgId() != ancestor->pdgId()) continue;
                if(isAncestor(ancestor,particle->mother(i))) return true;
        }
//if we did not return yet, then particle and ancestor are not relatives
        return false;
}






class genNtuple : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit genNtuple(const edm::ParameterSet&);
      ~genNtuple();
      
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      int getNPU(edm::Handle <std::vector <PileupSummaryInfo> >  puInfo);
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      // ----------member data ---------------------------
       //edm::EDGetTokenT <reco::GenParticleCollection> genparticlesToken;
       edm::EDGetTokenT <edm::View<pat::PackedGenParticle>> genparticlesToken;
       edm::EDGetTokenT <pat::ElectronCollection> electronsToken;
       edm::EDGetTokenT <edm::View<pat::MET>> metToken_ ;
       edm::EDGetTokenT <GenEventInfoProduct> genInfoToken_;
       edm::EDGetTokenT <double > rhoToken;
       edm::EDGetTokenT <edm::View<reco::GenJet>> genJetToken_;
       edm::EDGetTokenT <edm::View<pat::Jet>> recJetToken_;
       edm::EDGetTokenT <edm::View<pat::Photon>> photonsToken_ ;
       edm::EDGetTokenT <std::vector<reco::Vertex>>verticesToken_;
       edm::EDGetTokenT <GenRunInfoProduct>   genRunInfoToken_;
       edm::EDGetTokenT <edm::View<reco::GenParticle>> prunedGenToken_;



       TTree* m_tree;
       std::map<TString,Int_t>storageMapInt;
       std::map<TString,Float_t*>storageMapFloatArray;

       void addGenObjectBranches();
       void fillGenPrticle(TString tag,const reco::Candidate *part);
       void fillGenPrticle(TString tag,const reco::GenParticle &part);
       
       void addJetBranches();
       void fillJetBranches(const edm::Event & iEvent);
       
       void addGenJetBranches();
       void fillGenJetBranches(const edm::Event & iEvent);
       
       void addPhotonBranches();
       void fillPhotonBranches(const edm::Event & iEvent);


       static const int njetmx = 300;
       int ngenjetantikt;
       float genjetantiktpt[njetmx], genjetantikty[njetmx], genjetantiktphi[njetmx], genjetantiktmass[njetmx]; 
       
       int n_PU, Tnpv;
       float rho;
       double weight,xsec;
       //*******************variables for GsFelectron**********************//
       float leadJet_pt, leadJet_eta, leadJet_phi, leadJet_e, leadJet_m, leadJet_btg, leadJet_dj, leadJet_hflv, subleadJet_hflv,subleadJet_pt, subleadJet_eta, subleadJet_phi, subleadJet_e, subleadJet_m, subleadJet_btg, subleadJet_dj,ratio_mjj, mjj, gen_mjj, dphi1, dphi2, corr_dphi1, corr_dphi2, met, met_phi, corrected_met, corrected_met_phi, corrmet_px, corrmet_py, ht_tot, vbf_mjj, vbf_deta ;int vtx_size;
       float vtx_x[kMaxVertices], vtx_y[kMaxVertices], vtx_z[kMaxVertices];
       float vtx_pt2[kMaxVertices];
       int nj, nj_m; float ptgj[kMaxGenJet], ptgj_m[kMaxGenJet];
       float pt, eta, phi, en;
       float pt_V, eta_V, phi_V, en_V;
       Float_t q1_pt, q1_eta, q1_phi, q1_E, q2_pt, q2_eta, q2_phi, q2_E; Int_t q2_id, q1_id;
       Float_t b1_pt, b1_eta, b1_phi, b1_E, b2_pt, b2_eta, b2_phi, b2_E; Int_t b2_id, b1_id;
       Float_t Hm;
//******************************************************************//
};

genNtuple::genNtuple(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   prunedGenToken_      = consumes <edm::View<reco::GenParticle>> (iConfig.getParameter<edm::InputTag>("pruned"));
   photonsToken_       = consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons"));
   genJetToken_=consumes<edm::View<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJet"));
   recJetToken_=consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("recJet"));

   edm::Service<TFileService> fs;
    m_tree = fs->make<TTree>("tree", "");

    addGenObjectBranches();
    addGenJetBranches();
    addPhotonBranches();
    addJetBranches();
}


genNtuple::~genNtuple()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

                 
/*
int genNtuple::getNPU(edm::Handle <std::vector<PileupSummaryInfo>>puInfo)
{
    int nPU = -1;
    for(int iPileUp = (int) puInfo->size() - 1; iPileUp >= 0; iPileUp--)
    {
        PileupSummaryInfo  pileUpInfo = puInfo->at(iPileUp);
        
        int BX = pileUpInfo.getBunchCrossing();
        
        nPU = pileUpInfo.getPU_NumInteractions();
        if(BX == 0)
        {
            break;
        }
    }
    return nPU;
}
*/

void genNtuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
    
    edm::Handle<edm::View<reco::GenParticle>> pruned;
    iEvent.getByToken(prunedGenToken_, pruned);

    TLorentzVector b1;
    TLorentzVector b2;
    int nHiggs=0;
    size_t higgsIdx[10];
    Float_t higgsPt[10];
    for(size_t i=0; i< pruned->size();i++)
    {
        const reco::GenParticle * part_pru = (&(*pruned)[i]);
          if(not (&(*pruned)[i])->isHardProcess()) continue;
    //      if((&(*pruned)[i])->isHardProcess())
    //       std::cout<<i<<" pdgid : "<< part_pru->pdgId() << "    status : " << part_pru->status() 
    //                << "   isHard " << part_pru->isHardProcess()  
    //                << "   mother " << part_pru->mother()->pdgId()   << std::endl;
          if( part_pru->pdgId()==25 )
          {
                higgsPt[nHiggs]=part_pru->pt();
                higgsIdx[nHiggs]=i;
             //   std::cout<<" got higgs : "<<part_pru->pdgId()<<" , "
             //             <<higgsIdx[nHiggs]<<" , "<<higgsPt[nHiggs]<<" , "<<part_pru->numberOfDaughters()<<"\n";
                nHiggs++;
          }

    }

    if(nHiggs!=3)
    {
        std::cout<<"No Higgs found in event !\n";
        return;
    }

    size_t h1Idx(0),h2Idx(1),h3Idx(2);
    
    if(higgsPt[h1Idx] < higgsPt[h2Idx]) { auto x=h2Idx ; h2Idx=h1Idx ; h1Idx=x ;}
    if(higgsPt[h2Idx] < higgsPt[h3Idx]) { auto x=h3Idx ; h3Idx=h2Idx ; h2Idx=x ;}
    if(higgsPt[h1Idx] < higgsPt[h2Idx]) { auto x=h2Idx ; h2Idx=h1Idx ; h1Idx=x ;}
    
    h1Idx=higgsIdx[h1Idx];
    h2Idx=higgsIdx[h2Idx];
    h3Idx=higgsIdx[h3Idx];
    
    fillGenPrticle("H1",pruned->at(h1Idx));
    const reco::GenParticle *  hig= &(pruned->at(h1Idx));
    const reco::GenParticle * dau[2];
    Int_t dIdx(0);
    for(size_t i=0; i< pruned->size();i++)
    {
        const reco::GenParticle * part_pru = (&(*pruned)[i]);
        if(isAncestor(hig,part_pru) and  hig->pdgId()!=part_pru->pdgId())
        {
                dau[dIdx]=part_pru;
                dIdx++;
        }
    }
    if(dIdx ==2) {
        fillGenPrticle("H1_dau1",dau[0]);
        fillGenPrticle("H1_dau2",dau[1]);
    }
    else std::cout<<" null here \n";
    
    fillGenPrticle("H2",pruned->at(h2Idx));
    hig= &(pruned->at(h2Idx));
    dau[0]=nullptr; dau[1]=nullptr;
    dIdx=0;
    for(size_t i=0; i< pruned->size();i++)
    {
        const reco::GenParticle * part_pru = (&(*pruned)[i]);
        if(isAncestor(hig,part_pru) and  hig->pdgId()!=part_pru->pdgId())
        {
                dau[dIdx]=part_pru;
                dIdx++;
        }
    }
    if(dIdx ==2) {
        fillGenPrticle("H2_dau1",dau[0]);
        fillGenPrticle("H2_dau2",dau[1]);
    }
    else std::cout<<" null here \n";
    
    fillGenPrticle("H3",pruned->at(h3Idx));
    hig= &(pruned->at(h3Idx));
    dau[0]=nullptr; dau[1]=nullptr;
    dIdx=0;
    for(size_t i=0; i< pruned->size();i++)
    {
        const reco::GenParticle * part_pru = (&(*pruned)[i]);
        if(isAncestor(hig,part_pru) and  hig->pdgId()!=part_pru->pdgId())
        {
                dau[dIdx]=part_pru;
                dIdx++;
        }
    }
    if(dIdx ==2) {
        fillGenPrticle("H3_dau1",dau[0]);
        fillGenPrticle("H3_dau2",dau[1]);
    }
    else std::cout<<" null here \n";
    fillJetBranches(iEvent);
    fillGenJetBranches(iEvent);
    fillPhotonBranches(iEvent);


   
    m_tree->Fill();
}


void genNtuple::addGenObjectBranches()
{
    for(TString tag : {"H1","H2","H3"} )
    {
      
      std::cout<<" setting tag : "<<tag<<" \n";
      storageMapFloatArray["gen_"+tag+"_pt"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_pt", storageMapFloatArray["gen_"+tag+"_pt"]);
      storageMapFloatArray["gen_"+tag+"_pdgId"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_pdgId", storageMapFloatArray["gen_"+tag+"_pdgId"]);
      storageMapFloatArray["gen_"+tag+"_y"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_y", storageMapFloatArray["gen_"+tag+"_y"]);
      storageMapFloatArray["gen_"+tag+"_eta"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_eta", storageMapFloatArray["gen_"+tag+"_eta"]);
      storageMapFloatArray["gen_"+tag+"_phi"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_phi", storageMapFloatArray["gen_"+tag+"_phi"]);
      storageMapFloatArray["gen_"+tag+"_e"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_e", storageMapFloatArray["gen_"+tag+"_e"]);
      storageMapFloatArray["gen_"+tag+"_numberOfDaughters"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_numberOfDaughters", storageMapFloatArray["gen_"+tag+"_numberOfDaughters"]);
    }

    for(TString tag : {"H1_dau1","H2_dau1","H3_dau1","H1_dau2","H2_dau2","H3_dau2"} )
    {
      std::cout<<" setting tag : "<<tag<<" \n";
      storageMapFloatArray["gen_"+tag+"_pt"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_pt", storageMapFloatArray["gen_"+tag+"_pt"]);
      storageMapFloatArray["gen_"+tag+"_pdgId"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_pdgId", storageMapFloatArray["gen_"+tag+"_pdgId"]);
      storageMapFloatArray["gen_"+tag+"_y"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_y", storageMapFloatArray["gen_"+tag+"_y"]);
      storageMapFloatArray["gen_"+tag+"_eta"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_eta", storageMapFloatArray["gen_"+tag+"_eta"]);
      storageMapFloatArray["gen_"+tag+"_phi"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_phi", storageMapFloatArray["gen_"+tag+"_phi"]);
      storageMapFloatArray["gen_"+tag+"_e"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_e", storageMapFloatArray["gen_"+tag+"_e"]);
      storageMapFloatArray["gen_"+tag+"_numberOfDaughters"] = new Float_t;
      m_tree->Branch("gen_"+tag+"_numberOfDaughters", storageMapFloatArray["gen_"+tag+"_numberOfDaughters"]);
    }
}
void genNtuple::fillGenPrticle(TString tag, const reco::Candidate* particle)
{
      storageMapFloatArray["gen_"+tag+"_pt"][0]   = particle->pt() ;
      storageMapFloatArray["gen_"+tag+"_pdgId"][0]   = particle->pdgId() ;
      storageMapFloatArray["gen_"+tag+"_y"][0]   = particle->y() ;
      storageMapFloatArray["gen_"+tag+"_eta"][0]   = particle->eta() ;
      storageMapFloatArray["gen_"+tag+"_phi"][0]   = particle->phi() ;
      storageMapFloatArray["gen_"+tag+"_e"][0]   = particle->energy() ;
      storageMapFloatArray["gen_"+tag+"_numberOfDaughters"][0]   = particle->numberOfDaughters() ;
}

void genNtuple::fillGenPrticle(TString tag, const reco::GenParticle &particle)
{
      storageMapFloatArray["gen_"+tag+"_pt"][0]   = particle.pt() ;
      storageMapFloatArray["gen_"+tag+"_pdgId"][0]   = particle.pdgId() ;
      storageMapFloatArray["gen_"+tag+"_y"][0]   = particle.y() ;
      storageMapFloatArray["gen_"+tag+"_eta"][0]   = particle.eta() ;
      storageMapFloatArray["gen_"+tag+"_phi"][0]   = particle.phi() ;
      storageMapFloatArray["gen_"+tag+"_e"][0]   = particle.energy() ;
      storageMapFloatArray["gen_"+tag+"_numberOfDaughters"][0]   = particle.numberOfDaughters() ;
}



void genNtuple::addPhotonBranches()
{
      storageMapInt["nPhotons"]=0;
      m_tree->Branch("nPhotons",   &storageMapInt["nPhotons"]);
      storageMapFloatArray["photons_pt"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("photons_pt",   storageMapFloatArray["photons_pt"],"photons_pt[nPhotons]/F");
      storageMapFloatArray["photons_eta"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("photons_eta",   storageMapFloatArray["photons_eta"],"photons_eta[nPhotons]/F");
      storageMapFloatArray["photons_phi"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("photons_phi",   storageMapFloatArray["photons_phi"],"photons_phi[nPhotons]/F");
      storageMapFloatArray["photons_e"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("photons_e",   storageMapFloatArray["photons_e"],"photons_e[nPhotons]/F");
}

void genNtuple::fillGenJetBranches(const edm::Event &iEvent)
{
      edm::Handle<edm::View<reco::GenJet>> genjets;
      iEvent.getByToken(genJetToken_,genjets);
      storageMapInt["nGenJets"]=0;
      Int_t idx=0;

      for(auto const& gjet : *genjets) {
         storageMapFloatArray["genJets_pt"][idx]  = gjet.pt(); 
         storageMapFloatArray["genJets_eta"][idx] = gjet.eta();
         storageMapFloatArray["genJets_phi"][idx] = gjet.phi();
         storageMapFloatArray["genJets_phi"][idx] = gjet.phi();
         storageMapFloatArray["genJets_mass"][idx]   = gjet.mass();  
         storageMapFloatArray["genJets_e"][idx]   = gjet.energy();  
         idx++;
      }
      storageMapInt["nGenJets"]=idx;
}


void genNtuple::fillPhotonBranches(const edm::Event & iEvent)
{
   edm::Handle<edm::View<pat::Photon>>  photonsHandle ; 
   iEvent.getByToken(photonsToken_,photonsHandle);
   storageMapInt["nPhotons"]=0;
   Int_t idx=0;
    for (auto const& pho : *photonsHandle){
        storageMapFloatArray["photons_pt"][idx]       =             pho.et();
        storageMapFloatArray["photons_eta"][idx]      =             pho.eta();
        storageMapFloatArray["photons_phi"][idx]      =             pho.phi();
        storageMapFloatArray["photons_e"][idx]        =             pho.energy();
        idx++;
    }
    storageMapInt["nPhotons"]=idx;
}

void genNtuple::addGenJetBranches()
{
      storageMapInt["nGenJets"]=0;
      m_tree->Branch("nGenJets",   &storageMapInt["nGenJets"]);
      storageMapFloatArray["genJets_pt"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("genJets_pt",   storageMapFloatArray["genJets_pt"],"genJets_pt[nGenJets]/F");
      storageMapFloatArray["genJets_eta"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("genJets_eta",   storageMapFloatArray["genJets_eta"],"genJets_eta[nGenJets]/F");
      storageMapFloatArray["genJets_phi"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("genJets_phi",   storageMapFloatArray["genJets_phi"],"genJets_phi[nGenJets]/F");
      storageMapFloatArray["genJets_mass"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("genJets_mass",   storageMapFloatArray["genJets_mass"],"genJets_mass[nGenJets]/F");
      storageMapFloatArray["genJets_e"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("genJets_e",   storageMapFloatArray["genJets_e"],"genJets_e[nGenJets]/F");
}


void genNtuple::addJetBranches()
{
      storageMapInt["nJets"]=0;
      m_tree->Branch("nJets",   &storageMapInt["nJets"]);
      storageMapFloatArray["jets_pt"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("jets_pt",   storageMapFloatArray["jets_pt"],"jets_pt[nJets]/F");
      storageMapFloatArray["jets_eta"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("jets_eta",   storageMapFloatArray["jets_eta"],"jets_eta[nJets]/F");
      storageMapFloatArray["jets_phi"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("jets_phi",   storageMapFloatArray["jets_phi"],"jets_phi[nJets]/F");
      storageMapFloatArray["jets_mass"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("jets_mass",   storageMapFloatArray["jets_mass"],"jets_mass[nJets]/F");
      storageMapFloatArray["jets_csvScore"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("jets_csvScore",   storageMapFloatArray["jets_csvScore"],"jets_csvScore[nJets]/F");
      storageMapFloatArray["jets_deepCSVScore"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("jets_deepCSVScore",   storageMapFloatArray["jets_deepCSVScore"],"jets_deepCSVScore[nJets]/F");
      storageMapFloatArray["jets_deepJetScore"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("jets_deepJetScore",   storageMapFloatArray["jets_deepJetScore"],"jets_deepJetScore[nJets]/F");
      storageMapFloatArray["jets_flavour"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("jets_flavour",   storageMapFloatArray["jets_flavour"],"jets_flavour[nJets]/F");
      storageMapFloatArray["jets_pFlavour"] = new Float_t[N_ITEMS_MAX];
      m_tree->Branch("jets_pFlavour",   storageMapFloatArray["jets_pFlavour"],"jets_pFlavour[nJets]/F");
}


void genNtuple::fillJetBranches( const edm::Event& iEvent)
{
     
      edm::Handle<edm::View<pat::Jet>> recjets;
      iEvent.getByToken(recJetToken_, recjets);
      
      storageMapInt["nJets"]=0;
      Int_t idx=0;
      for(auto const & recJet : *recjets) {
            storageMapFloatArray["jets_pt"][idx]        = recJet.pt();
            storageMapFloatArray["jets_eta"][idx]       = recJet.eta();
            storageMapFloatArray["jets_phi"][idx]       = recJet.phi();
            storageMapFloatArray["jets_mass"][idx]         = recJet.mass();
            storageMapFloatArray["jets_csvScore"][idx]  = recJet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
            storageMapFloatArray["jets_deepCSVScore"][idx]   = recJet.bDiscriminator("pfDeepCSVJetTags:probb") + recJet.bDiscriminator("pfDeepCSVJetTags:probbb");
            storageMapFloatArray["jets_deepJetScore"][idx]   = recJet.bDiscriminator("pfDeepFlavourJetTags:probb")+
                                                               recJet.bDiscriminator("pfDeepFlavourJetTags:probbb")+
                                                               recJet.bDiscriminator("pfDeepFlavourJetTags:problepb");
            storageMapFloatArray["jets_flavour"][idx]        = recJet.hadronFlavour() ;
            storageMapFloatArray["jets_pFlavour"][idx]       = recJet.partonFlavour();
        idx++;
    }
    storageMapInt["nJets"]=idx;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
genNtuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(genNtuple);
