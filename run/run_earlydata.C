#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <caloreco/CaloTowerStatus.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/BeamBackgroundFilterAndQA.h>
#include <fstream>
#include <phool/recoConsts.h>
#include <TSystem.h>
#include <caloreco/CaloTowerCalib.h>
#include <frog/FROG.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <globalvertex/GlobalVertexReco.h>
#include <Calo_Calib.C>
#include <GlobalVertex.h>
#include <TruthJetInput.h>//#include <G4Setup_sPHENIX.C>
#include <MbdDigitization.h>
#include <MbdReco.h>
#include <zfinder/zfinder.h>
#include <zfindjets/zfindjets.h>

using namespace std;

R__LOAD_LIBRARY(libfun4all.so)
//R__LOAD_LIBRARY(libr24earlytreemaker.so)
R__LOAD_LIBRARY(libg4centrality.so)
//R__LOAD_LIBRARY(libFROG.so)
R__LOAD_LIBRARY(libg4vertex.so)
R__LOAD_LIBRARY(libglobalvertex.so);
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libcalo_io.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4mbd.so)
R__LOAD_LIBRARY(libmbd_io.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libtrigzvtxchecker.so)
R__LOAD_LIBRARY(libg4dst.so)
//gSystem->Load("libg4detectors.so");
//gSystem->Load("libg4detectors.so");
R__LOAD_LIBRARY(libzfinder.so)
R__LOAD_LIBRARY(libzfindjets.so)
int run_earlydata(string tag = "", int nproc = 0, int debug = 0, int nevt = 0, int usez = 0, int setz = 1)
{
  
  int verbosity = debug;
  string zfindfilename = "./"+to_string(nproc)+"_zfind/events_"+tag+"_"+to_string(nproc)+"_"+to_string(nproc)+"_"+to_string(nevt)+"_zfindfile.root";
  cout << "test1" << endl;
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc =  recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG","MDC2");
  rc->set_uint64Flag("TIMESTAMP",28);
  
  se->Verbosity(verbosity);
  // just if we set some flags somewhere in this macro


  

  Fun4AllDstInputManager *in_1 = new Fun4AllDstInputManager("DSTin1");
  Fun4AllDstInputManager *in_2 = new Fun4AllDstInputManager("DSTin2");
  Fun4AllDstInputManager *in_3 = new Fun4AllDstInputManager("DSTin3");
  Fun4AllDstInputManager *in_4 = new Fun4AllDstInputManager("DSTin4");
  in_1->Verbosity(debug);
  in_2->Verbosity(debug);
  in_3->Verbosity(debug);
  in_4->Verbosity(debug);
  cout << "get filenames" << endl;
  string line1, line2, line3, line4;
  line1 = "./dsts/"+to_string(nproc)+"/calo_cluster_"+to_string(nproc)+".root";
  line2 = "./dsts/"+to_string(nproc)+"/g4hits_"+to_string(nproc)+".root";
  line3 = "./dsts/"+to_string(nproc)+"/mbd_epd_"+to_string(nproc)+".root";
  line4 = "./dsts/"+to_string(nproc)+"/truth_jet_"+to_string(nproc)+".root";
  in_1->AddFile(line1);
  in_2->AddFile(line2);
  in_3->AddFile(line3);
  in_4->AddFile(line4);
  se->registerInputManager( in_1 );
  se->registerInputManager( in_2 );
  se->registerInputManager( in_3 );
  se->registerInputManager( in_4 );
  
  std::cout << "status setters" << std::endl;
  
  CDBInterface::instance()->Verbosity(0);

  Process_Calo_Calib();

  RetowerCEMC *rcemc = new RetowerCEMC();
  rcemc->set_towerinfo(true);
  rcemc->Verbosity(verbosity);
  se->registerSubsystem(rcemc);
  cout << "set up retower emcal" << endl;


  JetReco *zztowerjetreco = new JetReco("zz");
  TowerJetInput* zzemtji = new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER,"TOWERINFO_CALIB");
  TowerJetInput* zzohtji = new TowerJetInput(Jet::HCALIN_TOWERINFO,"TOWERINFO_CALIB");
  TowerJetInput* zzihtji = new TowerJetInput(Jet::HCALOUT_TOWERINFO,"TOWERINFO_CALIB");
  zzemtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
  zzohtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
  zzihtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);      
  zztowerjetreco->add_input(zzemtji);
  zztowerjetreco->add_input(zzohtji);
  zztowerjetreco->add_input(zzihtji);
  zztowerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "zzjets");
  zztowerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.6), "zzjets06");
  zztowerjetreco->set_algo_node("ANTIKT");
  zztowerjetreco->set_input_node("TOWER");
  zztowerjetreco->Verbosity(verbosity);
  se->registerSubsystem(zztowerjetreco);

  auto mbddigi = new MbdDigitization();                                                                              
  auto mbdreco = new MbdReco();
  GlobalVertexReco* gblvertex = new GlobalVertexReco();
  mbddigi->Verbosity(verbosity);
  se->registerSubsystem(mbddigi);
  mbdreco->Verbosity(verbosity);
  se->registerSubsystem(mbdreco);
  
  gblvertex->Verbosity(verbosity);
  se->registerSubsystem(gblvertex);
  
  se->Print("NODETREE");


  JetReco *dztowerjetreco = new JetReco("dz");
  TowerJetInput* dzemtji = new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER,"TOWERINFO_CALIB");
  TowerJetInput* dzohtji = new TowerJetInput(Jet::HCALIN_TOWERINFO,"TOWERINFO_CALIB");
  TowerJetInput* dzihtji = new TowerJetInput(Jet::HCALOUT_TOWERINFO,"TOWERINFO_CALIB");
  dzemtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
  dzohtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
  dzihtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);      
  dztowerjetreco->add_input(dzemtji);
  dztowerjetreco->add_input(dzohtji);
  dztowerjetreco->add_input(dzihtji);
  dztowerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "dzjets");
  dztowerjetreco->set_algo_node("ANTIKT");
  dztowerjetreco->set_input_node("TOWER");
  dztowerjetreco->Verbosity(verbosity);
  se->registerSubsystem(dztowerjetreco);


  /*
  JetReco *truthjetreco = new JetReco();  
  TruthJetInput *tji = new TruthJetInput(Jet::PARTICLE);
  tji->add_embedding_flag(0);  // changes depending on signal vs. embedded
  truthjetreco->add_input(tji);
  truthjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "AntiKt_Truth_r04");
  truthjetreco->set_algo_node("ANTIKT");
  truthjetreco->set_input_node("TRUTH");
  se->registerSubsystem(truthjetreco);
  */
  
  zfinder* zf = new zfinder("zf",verbosity,(usez==0?false:true),(setz==1?true:false));
  se->registerSubsystem(zf);

  JetReco *cztowerjetreco = new JetReco("cz");
  TowerJetInput* czemtji = new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER,"TOWERINFO_CALIB");
  TowerJetInput* czohtji = new TowerJetInput(Jet::HCALIN_TOWERINFO,"TOWERINFO_CALIB");
  TowerJetInput* czihtji = new TowerJetInput(Jet::HCALOUT_TOWERINFO,"TOWERINFO_CALIB");
  czemtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::CALO);
  czohtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::CALO);
  czihtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::CALO);      
  cztowerjetreco->add_input(czemtji);
  cztowerjetreco->add_input(czohtji);
  cztowerjetreco->add_input(czihtji);
  cztowerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "czjets");
  cztowerjetreco->set_algo_node("ANTIKT");
  cztowerjetreco->set_input_node("TOWER");
  cztowerjetreco->Verbosity(verbosity);
  se->registerSubsystem(cztowerjetreco);
  
  cout << "set up jetreco" << endl;
      //}

  zfindjets* zj = new zfindjets(zfindfilename,"zj",verbosity);
  se->registerSubsystem(zj);
  
  cout << "test4" << endl;
  se->Print("NODETREE");
  cout << "run " << nevt << endl;
  se->run(nevt);
  cout << "ran " << nevt << endl;
  cout << "Ran all events" << endl;
  se->Print("NODETREE");
  cout << "ending" << endl;
  se->End();
  cout << "printing timer" << endl;
  se->PrintTimer();
  cout << "Ended server" << endl;
  //delete se;
  cout << "Deleted server" << endl;
  //cout << "wrote " << filename << endl;
  //gSystem->Exit(0);
  //cout << "Exited gSystem" << endl;
  return 0;

}
