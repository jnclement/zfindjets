double lorentz(double* x, double* par)
{
  double num = par[0];
  double den = par[2]+(x[0]-par[1])*(x[0]-par[1])/par[2];
  return num/den;
}

double gaussian(double* x, double* par)
{
  double inexp = -(x[0]-par[1])*(x[0]-par[1])/(2*par[2]*par[2]);
  return par[0]*exp(inexp);
}

void fit_calib(int calib = 0, int isdat = 0, float thresh = 10, const char* filename = "hadded_20251014_eta.root",
	       const char* treename = "jet_tree")
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
    // Open the ROOT file
    TFile *f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        Error("analyze_zvtx", "Cannot open file %s", filename);
        return;
    }

    // Get the TTree
    TTree *t = (TTree*)f->Get(treename);
    if (!t) {
        Error("analyze_zvtx", "Tree %s not found in file", treename);
        return;
    }

    t->SetBranchStatus("*",0);
    t->SetBranchStatus("zvtx",1);
    t->SetBranchStatus("zzjet_n",1);
    t->SetBranchStatus("zzjet_e",1);
    
    // Set up the zvtx array
    float zvtx[3] = {0};
    int zzjet_n;
    float zzjet_e[100];
    t->SetBranchAddress("zvtx", zvtx);
    t->SetBranchAddress("zzjet_n",&zzjet_n);
    t->SetBranchAddress("zzjet_e",zzjet_e);
    float scalefactor = 1; //1.406 for unmodified
    float wasscaledby = 1.406;
    if(!calib)
      {
	wasscaledby = 1.0;
      }
    // Define histogram and profile
    TH2D *h2 = new TH2D("h2", ";z_{calo} [cm];z_{truth} [cm]",
                        61,-152.15*wasscaledby,152.5*wasscaledby,121,-302.5,302.5);

    Long64_t nentries = t->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);
	//cout << zvtx[0]/wasscaledby << endl;
	if(isdat && zzjet_e[zzjet_n-1] < thresh) continue;
        h2->Fill(zvtx[0]*scalefactor, zvtx[2]);
    }


    // Make the TProfile along X (average zvtx[2] vs zvtx[0])
    TProfile *prof = h2->ProfileX("prof");
    prof->SetMarkerStyle(20);

    TCanvas *c = new TCanvas("c", "", 600, 1200); 
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.15);
    c->SetTopMargin(0.15);
    c->SetBottomMargin(0.15);
    
    for(int i=1; i<h2->GetNbinsX()+1; ++i)
      {
	//if(prof->GetBinContent(i) != 0) cout << prof->GetBinLowEdge(i) << ":" << prof->GetBinContent(i) << endl;
	TH1D* proj = h2->ProjectionY("h2py",i,i,"e");
	proj->SetMarkerStyle(20);
	proj->Draw("PE");
	c->SaveAs(("proj_"+to_string(i)+".pdf").c_str());
	delete proj;
      }
    cout << endl;
    
    // Fit with a straight line
    TF1 *fitfunc = new TF1("fitfunc", "[0]*x", -150, 150);
    fitfunc->SetParameters(1);
    prof->Fit(fitfunc, "R");

    // Draw everything

    h2->Draw("COLZ");
    prof->SetLineColor(kRed);
    prof->Draw("SAME");

    c->SaveAs("zvtx_profile_fit"+(isdat?"data_thresh"+to_string(thresh):"pythia50")+".pdf");
    gPad->SetLogz();
    c->SaveAs("zvtx_profile_log"+(isdat?"data_thresh"+to_string(thresh):"pythia50")+".pdf");
    float slope = fitfunc->GetParameter(0);
    TCanvas* d = new TCanvas("","",1000,1000);
    d->SetLeftMargin(0.15);
    d->SetRightMargin(0.15);
    d->SetTopMargin(0.15);
    d->SetBottomMargin(0.15);
    d->cd();
    std::string s = (std::ostringstream() << std::fixed << std::setprecision(3) << slope).str();
    TH1D* h1_diff = new TH1D("h1_diff",(";z_{truth} - "+s+"*z_{calo} [cm]").c_str(),400,-200,200);

    for(int i=0; i<nentries; ++i)
      {
	t->GetEntry(i);
	h1_diff->Fill(zvtx[2] - slope*scalefactor*zvtx[0]);
      }
    TF1* f1 = new TF1("f1",gaussian,-50,50,3);
    f1->SetParameters(3e5,0,25);
    f1->SetParLimits(2,0,0);
    h1_diff->Fit(f1,"R");
    h1_diff->SetMarkerStyle(20);
    h1_diff->SetLineColor(kRed);
    h1_diff->Draw("PE");
    d->SaveAs("fitted_diff"+(isdat?"data_thresh"+to_string(thresh):"pythia50")+".pdf");
    gPad->SetLogy();
    d->SaveAs("fitted_diff_log"+(isdat?"data_thresh"+to_string(thresh):"pythia50")+".pdf");
    f->Close();
}
