#ifndef ZFINDJETS_H
#define ZFINDJETS_H
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <fun4all/SubsysReco.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include <calobase/RawTowerGeomContainer.h>
#include <phparameter/PHParameters.h>
#include <globalvertex/GlobalVertex.h>
class PHCompositeNode;
class CentralityInfo;
class zfindjets : public SubsysReco
{
 public:

  zfindjets(const std::string &filename = "/sphenix/user/jocl/projects/run2024_earlydata/run/output/temphists/debug.root", const std::string &name = "zfindjets", const int debug = 0, const int isdat = 0);

  virtual ~zfindjets();

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

 
 private:
  TTree* jet_tree;
  int _debug;
  std::string _filename;
  TFile* _f;
  std::string _name;
  float _zvtx[3];

  int _dzjet_n;
  float _dzjet_e[100];
  float _dzjet_pt[100];
  float _dzjet_eta[100];
  float _dzjet_phi[100];

  int _czjet_n;
  float _czjet_e[100];
  float _czjet_pt[100];
  float _czjet_eta[100];
  float _czjet_phi[100];

  int _zzjet_n;
  float _zzjet_e[100];
  float _zzjet_pt[100];
  float _zzjet_eta[100];
  float _zzjet_phi[100];

  int _tjet_n;
  float _tjet_e[100];
  float _tjet_pt[100];
  float _tjet_eta[100];
  float _tjet_phi[100];
  
  int _nprocessed;
  int _isdat;
};

#endif // ZFINDJETS
