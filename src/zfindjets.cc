#include <iomanip>
#include "zfindjets.h"
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoContainerv4.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexMapv1.h>
#include <globalvertex/GlobalVertex.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertexMapv1.h>
#include <globalvertex/MbdVertex.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <mbd/MbdPmtHit.h>
#include <jetbackground/TowerBackgroundv1.h>
#include <cmath>
#include <mbd/MbdOut.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/Jet.h>
#include <calobase/RawTowerv1.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos
#include <phool/recoConsts.h>
#include <iostream>
#include <centrality/CentralityInfo.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <TH2D.h>
#include <TColor.h>
#include <TMarker.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TH1D.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/MinimumBiasInfov1.h>
#include <calotrigger/MinimumBiasClassifier.h>
#include <ffarawobjects/Gl1Packetv3.h>
#include <TLorentzVector.h>
#include <ffaobjects/EventHeader.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <calobase/TowerInfoDefs.h>
#include <ffarawobjects/CaloPacketv1.h>
#include <globalvertex/Vertex.h>
using namespace std;

//____________________________________________________________________________..
zfindjets::zfindjets(const std::string &filename, const std::string &name, const int debug, int isdat):
  SubsysReco(name)
{
  _name = name;
  _debug = debug;
  _filename = filename;
  _nprocessed = 0;
  jet_tree = new TTree("jet_tree","a persevering date tree");
  _f = new TFile(_filename.c_str(), "RECREATE");
  _f->cd();
  jet_tree->SetDirectory(_f);
  _isdat = isdat;
}

//____________________________________________________________________________..
zfindjets::~zfindjets()
{

}

//____________________________________________________________________________..
int zfindjets::Init(PHCompositeNode *topNode)
{
  
  if(_debug > 1) cout << "Begin init: " << endl;
  jet_tree->Branch("zvtx",_zvtx,"zvtx[3]/F");

  jet_tree->Branch("dzjet_n",&_dzjet_n,"dzjet_n/I");
  jet_tree->Branch("dzjet_pt",_dzjet_pt,"dzjet_pt[dzjet_n]/F");
  jet_tree->Branch("dzjet_eta",_dzjet_eta,"dzjet_eta[dzjet_n]/F");
  jet_tree->Branch("dzjet_phi",_dzjet_phi,"dzjet_phi[dzjet_n]/F");
  jet_tree->Branch("dzjet_e",_dzjet_e,"dzjet_e[dzjet_n]/F");
  
  jet_tree->Branch("zzjet_n",&_zzjet_n,"zzjet_n/I");
  jet_tree->Branch("zzjet_pt",_zzjet_pt,"zzjet_pt[zzjet_n]/F");
  jet_tree->Branch("zzjet_eta",_zzjet_eta,"zzjet_eta[zzjet_n]/F");
  jet_tree->Branch("zzjet_phi",_zzjet_phi,"zzjet_phi[zzjet_n]/F");
  jet_tree->Branch("zzjet_e",_zzjet_e,"zzjet_e[zzjet_n]/F");

  jet_tree->Branch("czjet_n",&_czjet_n,"czjet_n/I");
  jet_tree->Branch("czjet_pt",_czjet_pt,"czjet_pt[czjet_n]/F");
  jet_tree->Branch("czjet_eta",_czjet_eta,"czjet_eta[czjet_n]/F");
  jet_tree->Branch("czjet_phi",_czjet_phi,"czjet_phi[czjet_n]/F");
  jet_tree->Branch("czjet_e",_czjet_e,"czjet_e[czjet_n]/F");
  if(!_isdat)
    {
      jet_tree->Branch("tjet_n",&_tjet_n,"tjet_n/I");
      jet_tree->Branch("tjet_e",_tjet_e,"tjet_e[tjet_n]/F");
      jet_tree->Branch("tjet_pt",_tjet_pt,"tjet_pt[tjet_n]/F");
      jet_tree->Branch("tjet_eta",_tjet_eta,"tjet_eta[tjet_n]/F");
      jet_tree->Branch("tjet_phi",_tjet_phi,"tjet_phi[tjet_n]/F");
    }
  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int zfindjets::InitRun(PHCompositeNode *topNode)
{
  if(_debug > 1) cout << "Initializing!" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

void print_debug(float jet_eta, float jet_phi, float tower_eta, float tower_phi, float dphi, float deta)
{
  cout << "printing debug info for dphi deta:" << endl;
  cout << "jet eta/phi: " << jet_eta << " " << jet_phi << endl;
  cout << "tower eta/phi: " << tower_eta << " " << tower_phi << endl;
  cout << "deta dphi:" << deta << " " << dphi << endl;
}

int zfindjets::process_event(PHCompositeNode *topNode)
{

  

  if(_debug > 1) cout << endl << endl << endl << "zfindjets: Beginning event processing" << endl;
  if(_nprocessed % 1000 == 0) cout << "processing event " << _nprocessed << endl;

  ++_nprocessed;
  GlobalVertexMap* gvtxmap = findNode::getClass<GlobalVertexMapv1>(topNode, "GlobalVertexMap");
  _zvtx[0] = 0;
  _zvtx[1] = 0;

  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (truthinfo)
    {
      PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
      _zvtx[2] = gvertex->get_z();
    }
  GlobalVertex::VTXTYPE vtxtype[2] = {GlobalVertex::VTXTYPE::CALO, GlobalVertex::VTXTYPE::MBD};

  if (gvtxmap)
    {
      if (gvtxmap->empty())
	{
	  if (_debug > 0)
	    {
	      std::cout << "gvtxmap empty - aborting event." << std::endl;
	    }
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
      for(int i=0; i<2; ++i)
	{
	  
	  for(auto gvtxit : *gvtxmap)
	    {
	      GlobalVertex* gvtx = gvtxit.second;
	      if (gvtx)
		{
		  auto startIter = gvtx->find_vertexes(vtxtype[i]);
		  auto endIter = gvtx->end_vertexes();
		  for (auto iter = startIter; iter != endIter; ++iter)
		    {
		      const auto &[type, vertexVec] = *iter;
		      if (type != vtxtype[i])
			{
			  continue;
			}
		      for (const auto *vertex : vertexVec)
			{
			  if (!vertex)
			    {
			      continue;
			    }
			  _zvtx[i] = vertex->get_z();
			  if(_zvtx[i] != 0) break;
			}
		    }
		}
	    }
	}
    }


  JetContainer *czjets = findNode::getClass<JetContainerv1>(topNode, "czjets");
  JetContainer *dzjets = findNode::getClass<JetContainerv1>(topNode, "dzjets");
  JetContainer *zzjets = findNode::getClass<JetContainerv1>(topNode, "zzjets");
  JetContainer* truthjets = NULL;
  truthjets = findNode::getClass<JetContainerv1>(topNode, "AntiKt_Truth_r04");
  _tjet_n = 0;
  _czjet_n = 0;
  _dzjet_n = 0;
  _zzjet_n = 0;
  if(truthjets)
    {
      for(int i=0; i<truthjets->size(); ++i)
        {
          Jet* jet = truthjets->get_jet(i);
	  _tjet_pt[_tjet_n] = jet->get_pt();
	  if(_tjet_pt[_tjet_n] < 1) continue;
          _tjet_e[_tjet_n] = jet->get_e();
          _tjet_eta[_tjet_n] = jet->get_eta();
          _tjet_phi[_tjet_n] = jet->get_phi();
          _tjet_n++;
        }
    }
  if(czjets)
    {
      int tocheck = czjets->size();
      if(_debug > 2) cout << "Found " << tocheck << " czjets to check..." << endl;
      for(int i=0; i<tocheck; ++i)
        {
          Jet *czjet = czjets->get_jet(i);
          if(czjet)
            {
	      if(czjet->get_e() < 10) continue;
	      _czjet_e[_czjet_n] = czjet->get_e();
	      _czjet_pt[_czjet_n] = czjet->get_pt();
	      _czjet_eta[_czjet_n] = czjet->get_eta();
	      _czjet_phi[_czjet_n] = czjet->get_phi();
	      ++_czjet_n; 
	    }
	}
    }

    if(dzjets)
    {
      int tocheck = dzjets->size();
      if(_debug > 2) cout << "Found " << tocheck << " dzjets to check..." << endl;
      for(int i=0; i<tocheck; ++i)
        {
          Jet *dzjet = dzjets->get_jet(i);
          if(dzjet)
            {
	      if(dzjet->get_e() < 10) continue;
	      _dzjet_e[_dzjet_n] = dzjet->get_e();
	      _dzjet_pt[_dzjet_n] = dzjet->get_pt();
	      _dzjet_eta[_dzjet_n] = dzjet->get_eta();
	      _dzjet_phi[_dzjet_n] = dzjet->get_phi();
	      ++_dzjet_n;
	    }
	}
    }

    if(zzjets)
    {
      int tocheck = zzjets->size();
      if(_debug > 2) cout << "Found " << tocheck << " zzjets to check..." << endl;
      for(int i=0; i<tocheck; ++i)
        {
          Jet *zzjet = zzjets->get_jet(i);
          if(zzjet)
            {
	      if(zzjet->get_e() < 10) continue;
	      _zzjet_e[_zzjet_n] = zzjet->get_e();
	      _zzjet_pt[_zzjet_n] = zzjet->get_pt();
	      _zzjet_eta[_zzjet_n] = zzjet->get_eta();
	      _zzjet_phi[_zzjet_n] = zzjet->get_phi();
	      ++_zzjet_n;
	    }
	}
    }

    jet_tree->Fill();


  if(_debug > 3) cout << "end event" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
    
}
//____________________________________________________________________________..
int zfindjets::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "zfindjets::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int zfindjets::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "zfindjets::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int zfindjets::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "zfindjets::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }
  cout << "ending run" << endl;
  if(jet_tree->GetEntries() > 0)
    {
      jet_tree->Write();
      cout << "tree written" << endl;
      _f->Write();
      cout << "file written" << endl;
      _f->Close();
      cout << "file closed" << endl;
    }
  cout << "file saved if necessary" << endl;
  //delete jet_tree;
  //delete _f;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int zfindjets::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "zfindjets::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void zfindjets::Print(const std::string &what) const
{
  std::cout << "zfindjets::Print(const std::string &what) const Printing info for " << what << std::endl;
}
