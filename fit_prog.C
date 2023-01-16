#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooWorkspace.h"
#include "RooBinning.h"
//roostat
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/PointSetInterval.h"
#include "RooMinuit.h"
#include "RooProfileLL.h"


#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <fstream>
using namespace std;

#include "fitf.h"
#include "TF1.h"
#include "TLegend.h"
#include "functions.C"
#include "RooTFnBinding.h"
#include "TGraphErrors.h"
#include "TMath.h"

using namespace RooFit;
using namespace RooStats;

void fit_program(TString isoSym, TString dauSym, TString gdauSym, Int_t AtNo, Int_t NeuNo, char *decayFile, char* bkgFile, double st = 5, double et = 10000, Int_t neuEff_par = 0){
  Int_t massA = AtNo + NeuNo;
  TString isoName = isoSym + massA;
    cout<<"I am going to fit the isotope which has"<<endl;
    cout<<"Z = "<<AtNo<<endl;
    cout<<"N = "<<NeuNo<<endl;
    cout<<"Neu Eff parameters = "<<neuEff_par<<endl;
    //0 for normal neutron efficeincy, 1 for lower band and 2 for upper band of neutron detection efficiency
    
    //to pull the information of parameters
    Int_t isoA=9999, isoB=9999, isoC=9999, isoD=9999, isoE=9999, isoF=9999, isoG=9999, isoH=9999, isoI=9999, isoJ=9999, isoK=9999, isoL=9999; //crossponding isotope line numbers, isotope A, B, C etc

    Double_t a, b, c, d, e, f, g, h, l, n, o;
    Int_t m = 157;
    Double_t Z[m], N[m], T[m], TL[m], TU[m], Pn[m], PnL[m], PnU[m], Pnn[m], PnnL[m], PnnU[m], neff1[m], neff2[m], neff1L[m], neff2L[m], neff1H[m];
    Double_t n2H, neff2H[m];
    Double_t n1, n2, n1L, n2L, n1H; //for the neutron efficiency from paramter file
    ifstream infile;
    infile.open("parameters.txt"); //open nuclear data file or parameters
    //to check wheather the file is properly open or not
    if( !infile ){
      cout<<"Unable to open the file!!! check the file name and try again"<<endl;
      exit(1); //terminate with error
    }

    Int_t i = 0;
    while (!infile.eof()) {
      infile>>a>>b>>c>>d>>e>>f>>g>>h>>l>>n>>o>>n1>>n2>>n1L>>n2L>>n1H>>n2H;
      Z[i] = a;
      N[i] = b;
      T[i] = c;
      TL[i] = d;
      TU[i] = e;
      Pn[i] = f;
      PnL[i] = g;
      PnU[i] = h;
      Pnn[i] = l;
      PnnL[i] = n;
      PnnU[i] = o;
      neff1[i] = n1;
      neff2[i] = n2;
      neff1L[i] = n1L;
      neff2L[i] = n2L;
      neff1H[i] = n1H;
      neff2H[i] = n2H;
      i++;
      cout<<n1<<"\t"<<n2<<"\t"<<n1L<<"\t"<<n2L<<"\t"<<n1H<<"\t"<<n2H<<endl;
    }

    //Let's print one random row
    // cout<<Z[4]<<"\t"<<N[4]<<"\t"<<T[4]<<"\t"<<TL[4]<<"\t"<<TU[4]<<"\t"<<Pn[4]<<"\t"<<PnL[4]<<"\t"<<PnU[4]<<"\t"<<Pnn[4]<<"\t"<<PnnL[4]<<"\t"<<PnnU[4]<<endl;

    
    //funcA or isotope A
    //parent nuclei
    for(Int_t j=0; j<m; j++){
        if( Z[j] == AtNo && N[j] == NeuNo) isoA = j;
    }
    Double_t actA = log(2)/T[isoA];
    Double_t actA_L = (log(2)/T[isoA])+(log(2)/TL[isoA]);
    Double_t actA_U = (log(2)/T[isoA])+(log(2)/TU[isoA]);
    Double_t pn1A = Pn[isoA]/100.;
    Double_t pn2A = Pnn[isoA]/100.;
    Double_t n1A = neff1[isoA]/100.;
    Double_t n2A = neff2[isoA]/100.;
    Double_t n1AL = neff1L[isoA]/100.;
    Double_t n2AL = neff2L[isoA]/100.;
    Double_t n1AH = neff1H[isoA]/100.;
    Double_t n2AH = neff2H[isoA]/100.;
    cout<<"n2AH = "<<n2AH<<endl;
    
    //funcB or isotope B
    for(Int_t j=0; j<m; j++){
        if( Z[j] == Int_t(AtNo+1) && N[j] == Int_t(NeuNo-1)) isoB = j;
    }
    Double_t actB = log(2)/T[isoB];
    Double_t actB_L = (log(2)/T[isoB])+(log(2)/TL[isoB]);
    Double_t actB_U = (log(2)/T[isoB])+(log(2)/TU[isoB]);
    Double_t pn1B = Pn[isoB]/100.;
    Double_t pn1B_L = pn1B + (PnL[isoB]/100.);
    Double_t pn1B_U = pn1B + (PnU[isoB]/100.);
    Double_t pn2B = Pnn[isoB]/100.;
    Double_t pn2B_L = pn2B + (PnnL[isoB]/100.);
    Double_t pn2B_U = pn2B + (PnnU[isoB]/100.);

    //funcC or isotope C
    for(Int_t j=0; j<m; j++){
        if( Z[j] == Int_t(AtNo+2) && N[j] == Int_t(NeuNo-2)) isoC = j;
    }
    Double_t actC = log(2)/T[isoC];
    Double_t actC_L = (log(2)/T[isoC])+(log(2)/TL[isoC]);
    Double_t actC_U = (log(2)/T[isoC])+(log(2)/TU[isoC]);
    Double_t pn1C = Pn[isoC]/100.;
    Double_t pn1C_L = pn1C + (PnL[isoC]/100.);
    Double_t pn1C_U = pn1C + (PnU[isoC]/100.);
    Double_t pn2C = Pnn[isoC]/100.;
    Double_t pn2C_L = pn2C + (PnnL[isoC]/100.);
    Double_t pn2C_U = pn2C + (PnnU[isoC]/100.);

    //funcD or isotope D
    for(Int_t j=0; j<m; j++){
        if( Z[j] == Int_t(AtNo+3) && N[j] == Int_t(NeuNo-3)) isoD = j;
    }
    Double_t actD = log(2)/T[isoD];
    Double_t actD_L = (log(2)/T[isoD])+(log(2)/TL[isoD]);
    Double_t actD_U = (log(2)/T[isoD])+(log(2)/TU[isoD]);
    Double_t pn1D = Pn[isoD]/100.;
    Double_t pn1D_L = pn1D + (PnL[isoD]/100.);
    Double_t pn1D_U = pn1D + (PnU[isoD]/100.);
    Double_t pn2D = Pnn[isoD]/100.;
    Double_t pn2D_L = pn2D + (PnnL[isoD]/100.);
    Double_t pn2D_U = pn2D + (PnnU[isoD]/100.);

    //funcE or isotope E
    for(Int_t j=0; j<m; j++){
        if( Z[j] == Int_t(AtNo+1) && N[j] == Int_t(NeuNo-2)) isoE = j;
    }
    Double_t actE = log(2)/T[isoE];
    Double_t actE_L = (log(2)/T[isoE])+(log(2)/TL[isoE]);
    Double_t actE_U = (log(2)/T[isoE])+(log(2)/TU[isoE]);
    Double_t pn1E = Pn[isoE]/100.;
    Double_t pn1E_L = pn1E + (PnL[isoE]/100.);
    Double_t pn1E_U = pn1E + (PnU[isoE]/100.);
    Double_t pn2E = Pnn[isoE]/100.;
    Double_t pn2E_L = pn2E + (PnnL[isoE]/100.);
    Double_t pn2E_U = pn2E + (PnnU[isoE]/100.);

    //funcF or isotope F
    for(Int_t j=0; j<m; j++){
        if( Z[j] == Int_t(AtNo+2) && N[j] == Int_t(NeuNo-3)) isoF = j;
    }
    Double_t actF = log(2)/T[isoF];
    Double_t actF_L = (log(2)/T[isoF])+(log(2)/TL[isoF]);
    Double_t actF_U = (log(2)/T[isoF])+(log(2)/TU[isoF]);
    Double_t pn1F = Pn[isoF]/100.;
    Double_t pn1F_L = pn1F + (PnL[isoF]/100.);
    Double_t pn1F_U = pn1F + (PnU[isoF]/100.);
    Double_t pn2F = Pnn[isoF]/100.;
    Double_t pn2F_L = pn2F + (PnnL[isoF]/100.);
    Double_t pn2F_U = pn2F + (PnnU[isoF]/100.);

    //funcG or isotope G
    for(Int_t j=0; j<m; j++){
        if( Z[j] == Int_t(AtNo+3) && N[j] == Int_t(NeuNo-4)) isoG = j;
    }
    Double_t actG = log(2)/T[isoG];
    Double_t actG_L = (log(2)/T[isoG])+(log(2)/TL[isoG]);
    Double_t actG_U = (log(2)/T[isoG])+(log(2)/TU[isoG]);
    Double_t pn1G = Pn[isoG]/100.;
    Double_t pn1G_L = pn1G + (PnL[isoG]/100.);
    Double_t pn1G_U = pn1G + (PnU[isoG]/100.);
    Double_t pn2G = Pnn[isoG]/100.;
    Double_t pn2G_L = pn2G + (PnnL[isoG]/100.);
    Double_t pn2G_U = pn2G + (PnnU[isoG]/100.);

    //funcH or isotope H
    for(Int_t j=0; j<m; j++){
        if( Z[j] == Int_t(AtNo+1) && N[j] == Int_t(NeuNo-3)) isoH = j;
    }
    Double_t actH = log(2)/T[isoH];
    Double_t actH_L = (log(2)/T[isoH])+(log(2)/TL[isoH]);
    Double_t actH_U = (log(2)/T[isoH])+(log(2)/TU[isoH]);
    Double_t pn1H = Pn[isoH]/100.;
    Double_t pn1H_L = pn1H + (PnL[isoH]/100.);
    Double_t pn1H_U = pn1H + (PnU[isoH]/100.);
    Double_t pn2H = Pnn[isoH]/100.;
    Double_t pn2H_L = pn2H + (PnnL[isoH]/100.);
    Double_t pn2H_U = pn2H + (PnnU[isoH]/100.);

    //funcI or isotope I
    for(Int_t j=0; j<m; j++){
        if( Z[j] == Int_t(AtNo+2) && N[j] == Int_t(NeuNo-4)) isoI = j;
    }
    Double_t actI = log(2)/T[isoI];
    Double_t actI_L = (log(2)/T[isoI])+(log(2)/TL[isoI]);
    Double_t actI_U = (log(2)/T[isoI])+(log(2)/TU[isoI]);
    Double_t pn1I = Pn[isoI]/100.;
    Double_t pn1I_L = pn1I + (PnL[isoI]/100.);
    Double_t pn1I_U = pn1I + (PnU[isoI]/100.);
    Double_t pn2I = Pnn[isoI]/100.;
    Double_t pn2I_L = pn2I + (PnnL[isoI]/100.);
    Double_t pn2I_U = pn2I + (PnnU[isoI]/100.);

    //funcJ or isotope J
    for(Int_t j=0; j<m; j++){
        if( Z[j] == Int_t(AtNo+3) && N[j] == Int_t(NeuNo-5)) isoJ = j;
    }
    Double_t actJ = log(2)/T[isoJ];
    Double_t actJ_L = (log(2)/T[isoJ])+(log(2)/TL[isoJ]);
    Double_t actJ_U = (log(2)/T[isoJ])+(log(2)/TU[isoJ]);
    Double_t pn1J = Pn[isoJ]/100.;
    Double_t pn1J_L = pn1J + (PnL[isoJ]/100.);
    Double_t pn1J_U = pn1J + (PnU[isoJ]/100.);
    Double_t pn2J = Pnn[isoJ]/100.;
    Double_t pn2J_L = pn2J + (PnnL[isoJ]/100.);
    Double_t pn2J_U = pn2J + (PnnU[isoJ]/100.);

    //funcK or isotope K
    for(Int_t j=0; j<m; j++){
        if( Z[j] == Int_t(AtNo+2) && N[j] == Int_t(NeuNo-5)) isoK = j;
    }
    Double_t actK = log(2)/T[isoK];
    Double_t actK_L = (log(2)/T[isoK])+(log(2)/TL[isoK]);
    Double_t actK_U = (log(2)/T[isoK])+(log(2)/TU[isoK]);
    Double_t pn1K = Pn[isoK]/100.;
    Double_t pn1K_L = pn1K + (PnL[isoK]/100.);
    Double_t pn1K_U = pn1K + (PnU[isoK]/100.);
    Double_t pn2K = Pnn[isoK]/100.;
    Double_t pn2K_L = pn2K + (PnnL[isoK]/100.);
    Double_t pn2K_U = pn2K + (PnnU[isoK]/100.);

    //funcL or isotope L
    for(Int_t j=0; j<m; j++){
        if( Z[j] == Int_t(AtNo+3) && N[j] == Int_t(NeuNo-6)) isoL = j;
    }
    Double_t actL = log(2)/T[isoL];
    Double_t actL_L = (log(2)/T[isoL])+(log(2)/TL[isoL]);
    Double_t actL_U = (log(2)/T[isoL])+(log(2)/TU[isoL]);
    
    

    //Initializing the backgrounds, neutron efficiency etc ..
    Double_t rN0, rN0L, rN0H;
    Double_t rbkg1, rbkg1L, rbkg1H;
    Double_t rbkg2, rbkg2L, rbkg2H;
    Double_t rbkg3, rbkg3L, rbkg3H;
    Double_t rneuEff, rneudEff;

    //open text file which has information of histogram backgrounds and neutron efficeincy
    cout<<"Hey I am going to opent background txt file"<<endl;

    std::ifstream ifs(bkgFile);
    if(ifs.is_open()){
      cout<<"I am opening txt file, which has information of histogram backgrounds"<<endl;
      ifs>>rN0>>rN0L>>rN0H;
      ifs>>rbkg1>>rbkg1L>>rbkg1H;
      ifs>>rbkg2>>rbkg2L>>rbkg2H;
      ifs>>rbkg3>>rbkg3L>>rbkg3H;
      //      ifs>>rneuEff>>rneudEff; no  need of these since neutron efficiencies are defined from parameters.txt
    }
    else cout<<"Can't open background file"<<endl;
    
    if(neuEff_par == 0) {
      rneuEff = n1A;
      rneudEff = n2A;
    }
    else if(neuEff_par == 1) {
      rneuEff = n1AL;
      rneudEff = n2AL;
    }
    else if(neuEff_par == 2) {
      rneuEff = n1AH;
      rneudEff = n2AH;
    }
    else {
      cout<<"DEFINE NEUTRON EFFICEINCY FIRST"<<endl;
    }
      
    
    //root file
    TFile *in = new TFile(decayFile);
    TTree *tree = (TTree*)in->Get("bdtree"); //NEED TO CHANGE LATER
    
    RooRealVar decay_time("decay_time","decay_time",st,et);
    RooCategory neu_mult("neu_mult","neu_mult");
    neu_mult.defineType("0neutron",0);
    neu_mult.defineType("1neutron",1);
    neu_mult.defineType("2neutron",2);
    
    RooDataSet *decaydata = new RooDataSet("decaydata","decaydata",RooArgSet(decay_time,neu_mult),Import(*tree));
    
    //calculate r1 and r2 using histogram or make histogram from tree
    TH1F *decaycurve = (TH1F*)in->Get("decayall");
    TH1F *n1decaycurve = (TH1F*)in->Get("decay1n");
    TH1F *n2decaycurve = (TH1F*)in->Get("decay2n");
    TH1F *n1decay_bkg = (TH1F*)in->Get("decay1n_bkg");
    TH1F *n2decay_bkg = (TH1F*)in->Get("decay2n_bkg");

    Double_t id = decaycurve->Integral();
    Double_t id1n = n1decaycurve->Integral();
    Double_t id2n = n2decaycurve->Integral();
    Double_t id1nb = n1decay_bkg->Integral();
    Double_t id2nb = n2decay_bkg->Integral();
    Double_t r1 = id1nb/id;
    Double_t r2 = id2nb/id;
    cout<<"THE VALUE OF r1 = "<<r1<<endl;
    cout<<"THE VALUE OF r2 = "<<r2<<endl;
    
    

    //TO INITILISING FITTING PARAMETERS
    RooRealVar N0("N0","N0",rN0,rN0L,rN0H);
    //RooRealVar N0("N0","N0",rN0);

    RooRealVar lamA("lamA","lamA",0.0693,0.00069314718,0.69314718);
    //RooRealVar lamA("lamA","lamA",actA);
    
    RooRealVar p1A("p1A","p1A",pn1A,0.0,1.0);
    //RooRealVar p1A("p1A","p1A",pn1A);
//    p1A.setVal(0.5);
    //RooRealVar p2A("p2A","p2A",0.0);    
     RooRealVar p2A("p2A","p2A",pn2A,0.0,1.0);
    //if(parent_p2n=0){p2A.setVal(0); p2A.setRange(0,0);}
    //RooRealVar p2A("p2A","p2A",pn2A);
    //    p2A.setVal(0.5);
    
    RooRealVar lamB("lamB","lamB",actB);
    RooRealVar p1B("p1B","p1B",pn1B);
    RooRealVar p2B("p2B","p2B",pn2B);
    
    RooRealVar lamC("lamC","lamC",actC);
    RooRealVar p1C("p1C","p1C",pn1C);
    RooRealVar p2C("p2C","p2C",pn2C);
    
    RooRealVar lamD("lamD","lamD",actD);
    RooRealVar p1D("p1D","p1D",pn1D);
    RooRealVar p2D("p2D","p2D",pn2D);

    RooRealVar lamE("lamE","lamE",actE);
    RooRealVar p1E("p1E","p1E",pn1E);
    RooRealVar p2E("p2E","p2E",pn2E);
    
    RooRealVar lamF("lamF","lamF",actF);
    RooRealVar p1F("p1F","p1F",pn1F);
    RooRealVar p2F("p2F","p2F",pn2F);
    
    RooRealVar lamG("lamG","lamG",actG);
    RooRealVar p1G("p1G","p1G",pn1G);
    RooRealVar p2G("p2G","p2G",pn2G);
    
    RooRealVar lamH("lamH","lamH",actH);
    RooRealVar p1H("p1H","p1H",pn1H);
    RooRealVar p2H("p2H","p2H",pn2H);

    RooRealVar lamI("lamI","lamI",actI);
    RooRealVar p1I("p1I","p1I",pn1I);
    RooRealVar p2I("p2I","p2I",pn2I);
    
    RooRealVar lamJ("lamJ","lamJ",actJ);
    RooRealVar p1J("p1J","p1J",pn1J);
    RooRealVar p2J("p2J","p2J",pn2J);
    
    RooRealVar lamK("lamK","lamK",actK);
    RooRealVar p1K("p1K","p1K",pn1K);
    RooRealVar p2K("p2K","p2K",pn2K);
    
    RooRealVar lamL("lamL","lamL",actL);

    RooRealVar neuEff("neuEff","neuEff",rneuEff);

    RooRealVar bkg1("bkg1","bkg1",rbkg1,rbkg1L,rbkg1H);
    //RooRealVar bkg1("bkg1","bkg1",rbkg1);
    
    RooRealVar bkg2("bkg2","bkg2",rbkg2,rbkg2L,rbkg2H);
    //RooRealVar bkg2("bkg2","bkg2",rbkg2);
    
    RooRealVar bkg3("bkg3","bkg3",rbkg3,rbkg3L,rbkg3H);

    //RooRealVar bkg3("bkg3","bkg3",rbkg3);
    
    RooRealVar rc1("rc1","rc1",r1);
    RooRealVar rc2("rc2","rc2",r2);

    RooRealVar neudEff("neudEff","neudEff",rneudEff);
    //added above line for neutron detection efficeincy for each of two delayed neutrons    
    cout<<"------> FITTING STARTS NOW -------->"<<endl;
    
    fitf fitddata("fitddata","fitddata",decay_time,neu_mult,N0,lamA,lamB,lamC,lamD,lamE,lamF,lamG,lamH,lamI,lamJ,lamK,lamL,p1A,p1B,p1C,p1D,p1E,p1F,p1G,p1H,p1I,p1J,p1K,p2A,p2B,p2C,p2D,p2E,p2F,p2G,p2H,p2I,p2J,p2K,neuEff,bkg1,bkg2,bkg3,rc1,rc2,neudEff);

    //TEST
   
    
    RooFitResult *res = fitddata.fitTo(*decaydata,NumCPU(60),Save());
    
    cout<<"DONE WITH EVERYTHING, FOR FITTING PART"<<endl;
    
    //to add the output histograms in a root files
    TFile *fOut;
    TString foutname = "fittingCanvas/"+isoName+"_results.root";
    fOut = new TFile(foutname,"recreate");
    
    //for binning part
    RooBinning tbins(5,10000);
    tbins.addUniform(1999,5,10000);

    //plots and residuals
    RooPlot *xframet = decay_time.frame(Title("Total decay"));
    decaydata->plotOn(xframet, Binning(tbins),DataError(RooAbsData::SumW2),RooFit::Name("totaldecay"));
    fitddata.plotOn(xframet, LineColor(kBlue), RooFit::Name("totaldecay_mod")); 
    RooHist *residualt = xframet->residHist();
    RooPlot *xframet_r = decay_time.frame(Title("Total decay residual"),RooFit::Name("r_totaldecay"));
    xframet_r->addPlotable(residualt,"P* X0");


    //to pull chi2 value from the plot
    //Double_t chi2_td = xframet->chiSquare(42); //this line is commented as this is not proper way to measure the chi2
    
    RooPlot *xframe0 = decay_time.frame(Title("decay with 0n delayed"));
    decaydata->plotOn(xframe0,Cut("neu_mult==neu_mult::0neutron"),Binning(tbins),DataError(RooAbsData::SumW2),RooFit::Name("n0decay"));
    fitddata.plotOn(xframe0,Slice(neu_mult,"0neutron"),LineColor(kBlue), RooFit::Name("n0decay_mod"));
    RooHist *residual0 = xframe0->residHist();
    RooPlot *xframe0_r = decay_time.frame(Title("decay with 0n delayed residual"),RooFit::Name("r_n0decay"));
    xframe0_r->addPlotable(residual0,"P* X0");

    RooPlot *xframe1 = decay_time.frame(Title("decay with 1n delayed"));
    decaydata->plotOn(xframe1,Cut("neu_mult==neu_mult::1neutron"),Binning(tbins),DataError(RooAbsData::SumW2),RooFit::Name("n1decay"));
    fitddata.plotOn(xframe1,Slice(neu_mult,"1neutron"),LineColor(kBlue), RooFit::Name("n1decay_mod"));
    RooHist *residual1 = xframe1->residHist();
    RooPlot *xframe1_r = decay_time.frame(Title("decay with 1n delayed residual"),RooFit::Name("r_n1decay"));
    xframe1_r->addPlotable(residual1,"P* X0");

    RooPlot *xframe2 = decay_time.frame(Title("decay with 2n delayed"));
    decaydata->plotOn(xframe2,Cut("neu_mult==neu_mult::2neutron"),Binning(tbins),DataError(RooAbsData::SumW2),RooFit::Name("n2decay"));
    fitddata.plotOn(xframe2,Slice(neu_mult,"2neutron"),LineColor(kBlue), RooFit::Name("n2decay_mod"));
    RooHist *residual2 = xframe2->residHist();
    RooPlot *xframe2_r = decay_time.frame(Title("decay with 2n delayed residual"),RooFit::Name("r_n2decay"));
    xframe2_r->addPlotable(residual2,"P* X0");
    
    TCanvas *c1 = new TCanvas("c1","c1",1000,800);
    c1->Divide(1,2);
    c1->cd(1);
    xframet->Draw();
    c1->cd(2);
    xframet_r->Draw();
    
    TCanvas *c0n = new TCanvas("c0n","c0n",1000,800);
    c0n->Divide(1,2);
    c0n->cd(1);
    xframe0->Draw();
    c0n->cd(2);
    xframe0_r->Draw();
    
    TCanvas *c1n = new TCanvas("c1n","c1n",1000,800);
    c1n->Divide(1,2);
    c1n->cd(1);
    xframe1->Draw();
    c1n->cd(2);
    xframe1_r->Draw();
    
    TCanvas *c2n = new TCanvas("c2n","c2n",1000,800);
    c2n->Divide(1,2);
    c2n->cd(1);
    xframe2->Draw();
    c2n->cd(2);
    xframe2_r->Draw();
    
    c1->Write();
    c0n->Write();
    c1n->Write();
    c2n->Write();
    //fOut->Close();


    //to calculate chi2 for each of the canvas //todo: need find the paper that describes this method of calculating chi2/ndf

    RooHist *hist_total = (RooHist*)xframet->getHist("totaldecay");
    RooHist *hist_n0 = (RooHist*)xframe0->getHist("n0decay");
    RooHist *hist_n1 = (RooHist*)xframe1->getHist("n1decay");
    RooHist *hist_n2 = (RooHist*)xframe2->getHist("n2decay");

    RooCurve *curve_total = (RooCurve*)xframet->getCurve("totaldecay_mod");
    RooCurve *curve_n0 = (RooCurve*)xframe0->getCurve("n0decay_mod");
    RooCurve *curve_n1 = (RooCurve*)xframe1->getCurve("n1decay_mod");
    RooCurve *curve_n2 = (RooCurve*)xframe2->getCurve("n2decay_mod");
  
    //for the first canvas //for decay spectra with out any neutron gate or conditions
    Double_t chisqua=0;
    Double_t xres[10000];
    Double_t yres[10000];
    Double_t yreserr[10000];
    Double_t chiVal = 0, logVal = 0;
    cout<<"total hist bins"<<hist_total->GetN();
    for (Int_t i=0;i<hist_total->GetN();i++){
        Double_t xi=hist_total->GetX()[i];
        Double_t yi=hist_total->GetY()[i];
        Double_t yeval=curve_total->Eval(xi);
        Double_t reldev=yeval-yi;
	if(yi !=0 && yeval != 0) logVal = TMath::Log(yi/yeval); else logVal = 0;
       	chiVal=yeval-yi+yi*logVal;
        chisqua += chiVal;
    }    
    chisqua=2*chisqua;
    cout<<"chi2 = "<<chisqua<<endl;
    cout<<"ndf="<<res->floatParsFinal().getSize()<<endl;
    cout<<"chisquare/ndf="<<chisqua/(1999-res->floatParsFinal().getSize())<<endl;
    Double_t chi2_over_ndf = chisqua/(1999-res->floatParsFinal().getSize());
    //1999 is total bins from 5 to 10000 ms
    //lets check if GetNbinsX provides the total number of bins with the given range
    Int_t test1 = hist_total->GetNbinsX();
    Int_t test2 = hist_total->GetN();
    cout<<"test1 NbinsX = "<<test1<<endl;
    cout<<"test2 N = "<<test2<<endl;    


    //for decay spectra with 1n emission
    Double_t chisqua_n1=0;
    Double_t chiVal_n1 = 0, logVal_n1 = 0;
    for (Int_t i=0;i<hist_n1->GetN();i++){
        Double_t xi_n1=hist_n1->GetX()[i];
        Double_t yi_n1=hist_n1->GetY()[i];
        Double_t yeval_n1=curve_n1->Eval(xi);
        Double_t reldev_n1=yeval_n1-yi_n1;
	if(yi_n1 !=0 && yeval_n1 != 0) logVal_n1 = TMath::Log(yi_n1/yeval_n1); else logVal_n1 = 0;
       	chiVal_n1=yeval_n1-yi_n1+yi_n1*logVal_n1;
        chisqua_n1 += chiVal_n1;
    }
   chisqua_n1=2*chisqua_n1;
   cout<<"chi2 = "<<chisqua_n1<<endl;
   cout<<"ndf="<<res->floatParsFinal().getSize()<<endl;
   cout<<"chisquare/ndf for 1n emi ="<<chisqua_n1/(1999-res->floatParsFinal().getSize())<<endl;
   Double_t chi2_over_ndf_n1 = chisqua_n1/(1999-res->floatParsFinal().getSize());


   //for decay spectra with 2n emission
   Double_t chisqua_n2=0;
   Double_t chiVal_n2 = 0, logVal_n2 = 0;
   for (Int_t i=0;i<hist_n2->GetN();i++){
     Double_t xi_n2=hist_n2->GetX()[i];
     Double_t yi_n2=hist_n2->GetY()[i];
     Double_t yeval_n2=curve_n2->Eval(xi);
     Double_t reldev_n2=yeval_n2-yi_n2;
     if(yi_n2 !=0 && yeval_n2 != 0) logVal_n2 = TMath::Log(yi_n2/yeval_n2); else logVal_n2 = 0;
     chiVal_n2=yeval_n2-yi_n2+yi_n2*logVal_n2;
     chisqua_n2 += chiVal_n2;
   }
   chisqua_n2=2*chisqua_n2;
   cout<<"chi2 = "<<chisqua_n2<<endl;
   cout<<"ndf="<<res->floatParsFinal().getSize()<<endl;
   cout<<"chisquare/ndf for 2n emi ="<<chisqua_n2/(1999-res->floatParsFinal().getSize())<<endl;
   Double_t chi2_over_ndf_n2 = chisqua_n2/(1999-res->floatParsFinal().getSize());
    
    
    cout<<"///////////////////////////////"<<endl;
    cout<<"\n\n\n\n\n\\n"<<endl;
  //to save the fitting spectra in PNG
  /////////////////////////////////////
  ////////////////////////////////////
  ///////////////////////////////////
  /* ------------------------------------------------------*/
  //added new part to draw isotope by isotope for total decay curve//
    TF1* ftall = new TF1("ftall",FT,5,1000,42);
    TF1* fparent = new TF1("fparent",funcA_d,5,1000,42);
    TF1* fdaughter = new TF1("fdaughter",funcB_d,5,1000,42);
    TF1* fgdaughter = new TF1("fgdaughter",funcC_d,5,1000,42);
    TF1* fdaughter1n = new TF1("fdaughter1n",funcE_d,5,1000,42);
    TF1* fdaughter2n = new TF1("fdaughter2n",funcH_d,5,1000,42);

    //for decay 1n curve
    TF1* ft1n = new TF1("ft1n",F1NT,5,1000,42);
    TF1 *fisoA1 = new TF1("fisoA1",funcA_d1,5,1000,42);
    TF1 *fisoB1 = new TF1("fisoB1",funcB_d1,5,1000,42);
    TF1 *fisoC1 = new TF1("fisoC1",funcC_d1,5,1000,42);
    TF1 *fisoE1 = new TF1("fisoE1",funcE_d1,5,1000,42);
    TF1 *fisoH1 = new TF1("fisoH1",funcH_d1,5,1000,42);

    //for decay 2n curve
    TF1* ft2n = new TF1("ft2n",F2NT,5,1000,42);
    TF1 *fisoA2 = new TF1("fisoA2",funcA_d2,5,1000,42);
    TF1 *fisoB2 = new TF1("fisoB2",funcB_d2,5,1000,42);
    TF1 *fisoC2 = new TF1("fisoC2",funcC_d2,5,1000,42);
    TF1 *fisoE2 = new TF1("fisoE2",funcE_d2,5,1000,42);
    TF1 *fisoH2 = new TF1("fisoH2",funcH_d2,5,1000,42);
    
    /*
    //////work for normalization//////////
    //retrive normailized const
    //RooArgSet nset(decay_time);
    //cout<<"pdf_Norm[decay_time] = "<<fitddata.getVal(&nset) <<endl;
    //RooAbsReal *igx = fitddata.createIntegral(decay_time,NormSet(decay_time));
    //cout<<"pdf_Int[decay_time] = "<<igx->getVal()<<endl;
    //Double_t igx = fitddata.getVal();
    // Double_t pdfint = igx->getValV();
    
    //return 'raw' unnormalized value of pdf
    Double_t gx = fitddata.getVal();
    //return value of pdf normalized in range 5 to 10s
    RooArgSet nset(decay_time);
    Double_t gx2 = fitddata.getVal(&nset);
    //return used norm factor
    RooAbsReal *igx = fitddata.createIntegral(decay_time);
    Double_t gx3 = igx->getVal();

    //return used norm factor 2nd method
    RooAbsReal *igx0 = fitddata.createIntegral(decay_time,NormSet(decay_time));
    Double_t gx4 = igx0->getVal();
   
    cout<<"raw unnormalized pdf (gx) = "<<gx<<endl;
    cout<<"normalized pdf in range (gx2) = "<<gx2<<endl;
    cout<<"used norm factor in roofit (gx3)= gx/gx1 = "<<gx3<<endl;
    cout<<"used norm factor in roofit (gx4)= gx/gx1 = "<<gx4<<endl;
    
    decaycurve->GetXaxis()->SetRangeUser(5,10000);
    Double_t int1 = decaycurve->Integral();

    // Double_t norm = int1/igx;
    cout<<"hist int = "<<int1<<endl;
    // cout<<"pdf int = "<<pdfint<<endl;
    Double_t norm = int1/gx3;
    cout<<"norm factor = "<<norm<<endl;
    */ //work for normalization
    /////////////////////////
    

    //fixed N0 to 1 and background components to zero to draw decay components
    Double_t array[42];
    array[0] = 1.0;
    //    array[0] = N0.getVal()*norm;
    array[1] = lamA.getVal();
    array[4] = lamB.getVal();
    array[7] = lamC.getVal();
    array[10] = lamD.getVal();
    array[13] = lamE.getVal();
    array[16] = lamF.getVal();
    array[19] = lamG.getVal();
    array[22] = lamH.getVal();
    array[25] = lamI.getVal();
    array[28] = lamJ.getVal();
    array[31] = lamK.getVal();
    array[34] = lamL.getVal();
    array[2] = p1A.getVal();
    array[5] = p1B.getVal();
    array[8] = p1C.getVal();
    array[11] = p1D.getVal();
    array[14] = p1E.getVal();
    array[17] = p1F.getVal();
    array[20] = p1G.getVal();
    array[23] = p1H.getVal();
    array[26] = p1I.getVal();
    array[29] = p1J.getVal();
    array[32] = p1K.getVal();
    array[3] = p2A.getVal();
    array[6] = p2B.getVal();
    array[9] = p2C.getVal();
    array[12] = p2D.getVal();
    array[15] = p2E.getVal();
    array[18] = p2F.getVal();
    array[21] = p2G.getVal();
    array[24] = p2H.getVal();
    array[27] = p2I.getVal();
    array[30] = p2J.getVal();
    array[33] = p2K.getVal();	
    array[35] = neuEff.getVal();
    //   array[36] = bkg1.getVal()*norm;
    array[36] = 0;
    //array[37] = bkg2.getVal()*norm;
    array[37] = 0;
    //array[38] = bkg3.getVal()*norm;
    array[38] = 0;
    array[39] = rc1.getVal();
    array[40] = rc2.getVal();
    array[41] = neudEff.getVal();
    
    
    for(Int_t i=0; i<42; i++){
      cout<<i<<"\t"<<array[i]<<endl;
      fparent->SetParameter(i,array[i]);
      fdaughter->SetParameter(i,array[i]);
      fdaughter1n->SetParameter(i,array[i]);
      fdaughter2n->SetParameter(i,array[i]);
      fgdaughter->SetParameter(i,array[i]);
      ftall->SetParameter(i,array[i]);

      ft1n->SetParameter(i,array[i]);
      fisoA1->SetParameter(i,array[i]);
      fisoB1->SetParameter(i,array[i]);
      fisoC1->SetParameter(i,array[i]);
      fisoE1->SetParameter(i,array[i]);
      fisoH1->SetParameter(i,array[i]);

      ft2n->SetParameter(i,array[i]);
      fisoA2->SetParameter(i,array[i]);
      fisoB2->SetParameter(i,array[i]);
      fisoC2->SetParameter(i,array[i]);
      fisoE2->SetParameter(i,array[i]);
      fisoH2->SetParameter(i,array[i]);

  //	  fbkg->SetParameter(i,array[i]); 
	}	  
    

    
    
    RooBinning tbins2(5,1000);
    tbins2.addUniform(199,5,1000);
    
    //for plots and pulls
    RooPlot *xframetp = decay_time.frame(Title("Total decay"));
    decaydata->plotOn(xframetp, Binning(tbins2),DataError(RooAbsData::SumW2));
    fitddata.plotOn(xframetp,LineColor(kBlue));
    RooHist *residualtp = xframet->residHist();
    RooPlot *xframet_rp = decay_time.frame(Title("Residual"));
    xframet_rp->addPlotable(residualtp,"P* X0");

    
    RooPlot *xframe0p = decay_time.frame(Title("Decay with 0 neutron emission"));
    decaydata->plotOn(xframe0p,Cut("neu_mult==neu_mult::0neutron"),Binning(tbins2),DataError(RooAbsData::SumW2));
    fitddata.plotOn(xframe0p,Slice(neu_mult,"0neutron"),LineColor(kBlue));
    RooHist *residual0p = xframe0p->residHist();
    RooPlot *xframe0_rp = decay_time.frame(Title("Residual"));
    xframe0_rp->addPlotable(residual0p,"P* X0");

    RooPlot *xframe1p = decay_time.frame(Title("Decay with 1 neutron emission"));
    decaydata->plotOn(xframe1p,Cut("neu_mult==neu_mult::1neutron"),Binning(tbins2),DataError(RooAbsData::SumW2));
    fitddata.plotOn(xframe1p,Slice(neu_mult,"1neutron"),LineColor(kBlue));
    RooHist *residual1p = xframe1p->residHist();
    RooPlot *xframe1_rp = decay_time.frame(Title("Residual"));
    xframe1_rp->addPlotable(residual1p,"P* X0");

    RooPlot *xframe2p = decay_time.frame(Title("Decay with 2 neutron emission"));
    decaydata->plotOn(xframe2p,Cut("neu_mult==neu_mult::2neutron"),Binning(tbins2),DataError(RooAbsData::SumW2));
    fitddata.plotOn(xframe2p,Slice(neu_mult,"2neutron"),LineColor(kBlue));
    RooHist *residual2p = xframe2->residHist();
    RooPlot *xframe2_rp = decay_time.frame(Title("Residual"));
    xframe2_rp->addPlotable(residual2p,"P* X0");

    xframetp->GetXaxis()->SetLabelSize(0.05);
    xframetp->GetYaxis()->SetLabelSize(0.05);
    xframetp->GetXaxis()->SetTitle("Decay time (ms)");
    xframetp->GetYaxis()->SetTitle("Counts per 5 ms");
    xframetp->GetXaxis()->SetTitleSize(0.05);
    xframetp->GetYaxis()->SetTitleSize(0.05);
    xframetp->GetXaxis()->SetTitleOffset(0.98);
    xframetp->GetYaxis()->SetTitleOffset(0.98);
    xframetp->GetXaxis()->CenterTitle();
    xframetp->GetYaxis()->CenterTitle();
    xframetp->GetXaxis()->SetMaxDigits(2);
    xframetp->GetYaxis()->SetMaxDigits(2);
    xframetp->GetXaxis()->SetRangeUser(5,1000);
    
    xframet_rp->GetXaxis()->SetLabelSize(0.05);
    xframet_rp->GetYaxis()->SetLabelSize(0.05);
    xframet_rp->GetXaxis()->SetTitle("Decay time (ms)");
    xframet_rp->GetYaxis()->SetTitle("Fit - data");
    xframet_rp->GetXaxis()->SetTitleSize(0.05);
    xframet_rp->GetYaxis()->SetTitleSize(0.05);
    xframet_rp->GetXaxis()->SetTitleOffset(0.98);
    xframet_rp->GetYaxis()->SetTitleOffset(0.98);
    xframet_rp->GetXaxis()->CenterTitle();
    xframet_rp->GetYaxis()->CenterTitle();
    xframet_rp->GetXaxis()->SetMaxDigits(2);
    xframet_rp->GetYaxis()->SetMaxDigits(2);
    xframet_rp->GetXaxis()->SetRangeUser(5,1000);
    TCanvas *cdg = new TCanvas("cdg","cdg",1500,600);
    cdg->Divide(3,1);
    cdg->cd(1);
    xframetp->Draw();
    cdg->cd(2);
    xframet_rp->Draw();
    cdg->cd(3);

    
    ftall->SetLineColor(kBlue);
    ftall->SetTitle("Decay Components");
    ftall->GetXaxis()->SetLabelSize(0.05);
    ftall->GetYaxis()->SetLabelSize(0.05);
    ftall->GetXaxis()->SetTitle("Decay time (ms)");
    ftall->GetYaxis()->SetTitle("Normalized counts");
    ftall->GetXaxis()->SetTitleSize(0.05);
    ftall->GetYaxis()->SetTitleSize(0.05);
    ftall->GetXaxis()->SetTitleOffset(0.98);
    ftall->GetYaxis()->SetTitleOffset(0.98);
    ftall->GetXaxis()->CenterTitle();
    ftall->GetYaxis()->CenterTitle();
    ftall->GetXaxis()->SetMaxDigits(2);
    ftall->GetYaxis()->SetMaxDigits(2);

    ftall->Draw();
    fparent->SetLineColor(kRed);
    fparent->SetTitle("parent");
    fparent->SetLineWidth(2);
    fparent->Draw("same");
    fdaughter->SetLineColor(kGreen);
    fdaughter->SetLineWidth(2);
    fdaughter->Draw("same");
    fgdaughter->SetLineColor(kMagenta);
    fgdaughter->SetLineWidth(2);
    fgdaughter->Draw("same");
    fdaughter1n->SetLineColor(kBlack);
    fdaughter1n->SetLineWidth(2);
    fdaughter1n->Draw("same");
    fdaughter2n->SetLineColor(kYellow);
    fdaughter2n->SetLineWidth(2);
    if(p2A.getVal()>0.001)    fdaughter2n->Draw("same");
    //legend
    TLegend *leg = new TLegend(0.50,0.50,0.85,0.85);
    Int_t massNum = AtNo + NeuNo;
    Int_t massNum1 = massNum-1;
    Int_t massNum2 = massNum-2;
    //TString parentNam = isoName;
    //TString daughtNam = dauName + massNum;
    TString parentNam = Form("^{%d}"+isoSym,massNum);
    TString daughtNam = Form("^{%d}"+dauSym,massNum);
    TString daughtNam1n = Form("^{%d}"+dauSym,massNum1);
    //    TString daughtNam1n = dauSym + massNum1;
    TString daughtNam2n = Form("^{%d}"+dauSym,massNum2);
    //TString daughtNam2n = dauSym + massNum2;
    //TString gdaughtNam = gdauSym + massNum;
    TString gdaughtNam = Form("^{%d}"+gdauSym,massNum);
    
    leg->AddEntry(ftall,"Total Fit","l");
    leg->AddEntry(fparent,parentNam,"l");
    leg->AddEntry(fdaughter,daughtNam,"l");
    leg->AddEntry(fgdaughter,gdaughtNam,"l");
    leg->AddEntry(fdaughter1n,daughtNam1n,"l");
    if(p2A.getVal()>0.001) leg->AddEntry(fdaughter2n,daughtNam2n,"l");
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->Draw();

    cdg->SaveAs("savedCanvas/"+isoName+"decayall.png");
    cdg->Write();

    xframe1p->GetXaxis()->SetLabelSize(0.05);
    xframe1p->GetYaxis()->SetLabelSize(0.05);
    xframe1p->GetXaxis()->SetTitle("Decay time (ms)");
    xframe1p->GetYaxis()->SetTitle("Counts per 5 ms");
    xframe1p->GetXaxis()->SetTitleSize(0.05);
    xframe1p->GetYaxis()->SetTitleSize(0.05);
    xframe1p->GetXaxis()->SetTitleOffset(0.98);
    xframe1p->GetYaxis()->SetTitleOffset(0.98);
    xframe1p->GetXaxis()->CenterTitle();
    xframe1p->GetYaxis()->CenterTitle();
    xframe1p->GetXaxis()->SetMaxDigits(2);
    xframe1p->GetYaxis()->SetMaxDigits(2);
    xframe1p->GetXaxis()->SetRangeUser(5,1000);
    
    xframe1_rp->GetXaxis()->SetLabelSize(0.05);
    xframe1_rp->GetYaxis()->SetLabelSize(0.05);
    xframe1_rp->GetXaxis()->SetTitle("Decay time (ms)");
    xframe1_rp->GetYaxis()->SetTitle("Fit - data");
    xframe1_rp->GetXaxis()->SetTitleSize(0.05);
    xframe1_rp->GetYaxis()->SetTitleSize(0.05);
    xframe1_rp->GetXaxis()->SetTitleOffset(0.98);
    xframe1_rp->GetYaxis()->SetTitleOffset(0.98);
    xframe1_rp->GetXaxis()->CenterTitle();
    xframe1_rp->GetYaxis()->CenterTitle();
    xframe1_rp->GetXaxis()->SetMaxDigits(2);
    xframe1_rp->GetYaxis()->SetMaxDigits(2);
    xframe1_rp->GetXaxis()->SetRangeUser(5,1000);
    TCanvas *cdh = new TCanvas("cdh","cdh",1500,600);
    cdh->Divide(3,1);
    cdh->cd(1);
    xframe1p->Draw();
    cdh->cd(2);
    xframe1_rp->Draw();
    cdh->cd(3);
    ft1n->SetLineColor(kBlue);
    ft1n->SetTitle("Decay Components 1n gated");
    ft1n->GetXaxis()->SetLabelSize(0.05);
    ft1n->GetYaxis()->SetLabelSize(0.05);
    ft1n->GetXaxis()->SetTitle("Decay time (ms)");
    ft1n->GetYaxis()->SetTitle("Normalized counts");
    ft1n->GetXaxis()->SetTitleSize(0.05);
    ft1n->GetYaxis()->SetTitleSize(0.05);
    ft1n->GetXaxis()->SetTitleOffset(0.98);
    ft1n->GetYaxis()->SetTitleOffset(0.98);
    ft1n->GetXaxis()->CenterTitle();
    ft1n->GetYaxis()->CenterTitle();
    ft1n->GetXaxis()->SetMaxDigits(2);
    ft1n->GetYaxis()->SetMaxDigits(2);
    ft1n->Draw();

    fisoA1->SetLineColor(kRed);
    fisoA1->SetLineWidth(2);
    fisoA1->Draw("same");

    fisoB1->SetLineColor(kGreen);
    fisoB1->SetLineWidth(2);
    fisoB1->Draw("same");

    fisoC1->SetLineColor(kMagenta);
    fisoC1->SetLineWidth(2);
    fisoC1->Draw("same");

    fisoE1->SetLineColor(kBlack);
    fisoE1->SetLineWidth(2);
    fisoE1->Draw("same");

    fisoH1->SetLineColor(kYellow);
    fisoH1->SetLineWidth(2);
    if(p2A.getVal()>0.001) fisoH1->Draw("same");

    //legend
    TLegend *leg2 = new TLegend(0.50,0.50,0.85,0.85);
    /*   Int_t massNum = AtNo + NeuNo;
    Int_t massNum1 = massNum-1;
    Int_t massNum2 = massNum-2;
    TString parentNam = isoName;
    TString daughtNam = dauName + massNum;
    TString daughtNam1n = dauName + massNum1;
    TString daughtNam2n = dauName + massNum2;
    TString gdaughtNam = gdauName + massNum; */

    leg2->AddEntry(ft1n,"Total Fit","l");
    leg2->AddEntry(fisoA1,parentNam,"l");
    leg2->AddEntry(fisoB1,daughtNam,"l");
    leg2->AddEntry(fisoC1,gdaughtNam,"l");
    leg2->AddEntry(fisoE1,daughtNam1n,"l");
    if(p2A.getVal()>0.001) leg2->AddEntry(fisoH1,daughtNam2n,"l");
    leg2->SetBorderSize(0);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.05);
    leg2->Draw();
    
    cdh->SaveAs("savedCanvas/"+isoName+"decay_1n.png");
    cdh->Write();

    xframe2p->GetXaxis()->SetLabelSize(0.05);
    xframe2p->GetYaxis()->SetLabelSize(0.05);
    xframe2p->GetXaxis()->SetTitle("Decay time (ms)");
    xframe2p->GetYaxis()->SetTitle("Counts per 5 ms");
    xframe2p->GetXaxis()->SetTitleSize(0.05);
    xframe2p->GetYaxis()->SetTitleSize(0.05);
    xframe2p->GetXaxis()->SetTitleOffset(0.98);
    xframe2p->GetYaxis()->SetTitleOffset(0.98);
    xframe2p->GetXaxis()->CenterTitle();
    xframe2p->GetYaxis()->CenterTitle();
    xframe2p->GetXaxis()->SetMaxDigits(2);
    xframe2p->GetYaxis()->SetMaxDigits(2);
    xframe2p->GetXaxis()->SetRangeUser(5,1000);
    
    xframe2_rp->GetXaxis()->SetLabelSize(0.05);
    xframe2_rp->GetYaxis()->SetLabelSize(0.05);
    xframe2_rp->GetXaxis()->SetTitle("Decay time (ms)");
    xframe2_rp->GetYaxis()->SetTitle("Fit - data");
    xframe2_rp->GetXaxis()->SetTitleSize(0.05);
    xframe2_rp->GetYaxis()->SetTitleSize(0.05);
    xframe2_rp->GetXaxis()->SetTitleOffset(0.98);
    xframe2_rp->GetYaxis()->SetTitleOffset(0.98);
    xframe2_rp->GetXaxis()->CenterTitle();
    xframe2_rp->GetYaxis()->CenterTitle();
    xframe2_rp->GetXaxis()->SetMaxDigits(2);
    xframe2_rp->GetYaxis()->SetMaxDigits(2);
    xframe2_rp->GetXaxis()->SetRangeUser(5,1000);
    TCanvas *cdf = new TCanvas("cdf","cdf",1500,600);
    cdf->Divide(3,1);
    cdf->cd(1);
    xframe2p->Draw();
    cdf->cd(2);
    xframe2_rp->Draw();
    cdf->cd(3);

    ft2n->SetLineColor(kBlue);
    ft2n->SetTitle("Decay Components 2n gated");
    ft2n->GetXaxis()->SetLabelSize(0.05);
    ft2n->GetYaxis()->SetLabelSize(0.05);
    ft2n->GetXaxis()->SetTitle("Decay time (ms)");
    ft2n->GetYaxis()->SetTitle("Normalized counts");
    ft2n->GetXaxis()->SetTitleSize(0.05);
    ft2n->GetYaxis()->SetTitleSize(0.05);
    ft2n->GetXaxis()->SetTitleOffset(0.98);
    ft2n->GetYaxis()->SetTitleOffset(0.98);
    ft2n->GetXaxis()->CenterTitle();
    ft2n->GetYaxis()->CenterTitle();
    ft2n->GetXaxis()->SetMaxDigits(2);
    ft2n->GetYaxis()->SetMaxDigits(2);
    ft2n->Draw();

    fisoA2->SetLineColor(kRed);
    fisoA2->SetLineWidth(2);
    fisoA2->Draw("same");

    fisoB2->SetLineColor(kGreen);
    fisoB2->SetLineWidth(2);
    fisoB2->Draw("same");

    fisoC2->SetLineColor(kMagenta);
    fisoC2->SetLineWidth(2);
    fisoC2->Draw("same");

    fisoE2->SetLineColor(kBlack);
    fisoE2->SetLineWidth(2);
    fisoE2->Draw("same");

    fisoH2->SetLineColor(kYellow);
    fisoH2->SetLineWidth(2);
    if(p2A.getVal()>0.001) fisoH2->Draw("same");

    //legend
    TLegend *leg3 = new TLegend(0.50,0.50,0.85,0.85);
    /*   Int_t massNum = AtNo + NeuNo;
    Int_t massNum1 = massNum-1;
    Int_t massNum2 = massNum-2;
    TString parentNam = isoName;
    TString daughtNam = dauName + massNum;
    TString daughtNam1n = dauName + massNum1;
    TString daughtNam2n = dauName + massNum2;
    TString gdaughtNam = gdauName + massNum; */

    leg3->AddEntry(ft2n,"Total Fit","l");
    leg3->AddEntry(fisoA2,parentNam,"l");
    leg3->AddEntry(fisoB2,daughtNam,"l");
    leg3->AddEntry(fisoC2,gdaughtNam,"l");
    leg3->AddEntry(fisoE2,daughtNam1n,"l");
    if(p2A.getVal()>0.001) leg3->AddEntry(fisoH2,daughtNam2n,"l");
    leg3->SetBorderSize(0);
    leg3->SetTextFont(42);
    leg3->SetTextSize(0.05);
    leg3->Draw();

    cdf->SaveAs("savedCanvas/"+isoName+"decay_2n.png");
    cdf->Write();
    
    cout<<"\n\nIsotope A, T_1/2 = "<<T[isoA]<<"\tP1n = "<<Pn[isoA]<<"\tP2n = "<<Pnn[isoA]<<endl;
    cout<<"Isotope B, T_1/2 = "<<T[isoB]<<"\tP1n = "<<Pn[isoB]<<"\tP2n = "<<Pnn[isoB]<<endl;
    cout<<"Isotope B, T_1/2 = "<<T[isoC]<<"\tP1n = "<<Pn[isoC]<<"\tP2n = "<<Pnn[isoC]<<endl;
    cout<<"Isotope D, T_1/2 = "<<T[isoD]<<"\tP1n = "<<Pn[isoD]<<"\tP2n = "<<Pnn[isoD]<<endl;
    cout<<"Isotope E, T_1/2 = "<<T[isoE]<<"\tP1n = "<<Pn[isoE]<<"\tP2n = "<<Pnn[isoE]<<endl;
    cout<<"Isotope F, T_1/2 = "<<T[isoF]<<"\tP1n = "<<Pn[isoF]<<"\tP2n = "<<Pnn[isoF]<<endl;
    cout<<"Isotope G, T_1/2 = "<<T[isoG]<<"\tP1n = "<<Pn[isoG]<<"\tP2n = "<<Pnn[isoG]<<endl;
    cout<<"Isotope H, T_1/2 = "<<T[isoH]<<"\tP1n = "<<Pn[isoH]<<"\tP2n = "<<Pnn[isoH]<<endl;
    cout<<"Isotope I, T_1/2 = "<<T[isoI]<<"\tP1n = "<<Pn[isoI]<<"\tP2n = "<<Pnn[isoI]<<endl;
    cout<<"Isotope J, T_1/2 = "<<T[isoJ]<<"\tP1n = "<<Pn[isoJ]<<"\tP2n = "<<Pnn[isoJ]<<endl;
    cout<<"Isotope K, T_1/2 = "<<T[isoK]<<"\tP1n = "<<Pn[isoK]<<"\tP2n = "<<Pnn[isoK]<<endl;
    cout<<"Isotope L, T_1/2 = "<<T[isoL]<<endl;
    cout<<"r1 = "<<r1<<endl;
    cout<<"r2 = "<<r2<<endl;

    
    cout<<"\n\n******* RESULTS *******\n"<<endl;
    cout<<"N0 = "<<N0.getVal()<<endl;
    cout<<"Half-life ="<<log(2)/lamA.getVal()<<endl;
    cout<<"Half-life error ="<<log(2)/pow(lamA.getVal(),2)*lamA.getErrorHi()<<endl;
    cout<<"Prob of 1 neu emi ="<<p1A.getVal()<<endl;
    cout<<"Prob of 1 neu emi error ="<<p1A.getErrorHi()<<endl;
    cout<<"Prob of 2 neu emi ="<<p2A.getVal()<<endl;
    cout<<"Prob of 2 neu emi error ="<<p2A.getErrorHi()<<endl;
    cout<<"bkg1 = "<<bkg1.getVal()<<endl;
    cout<<"bkg2 = "<<bkg2.getVal()<<endl;
    cout<<"bkg3 = "<<bkg3.getVal()<<endl;
    cout<<"\n"<<endl;
    cout<<"chi2 / ndf "<<chi2_over_ndf<<endl;

    
    //close output root canvas file
    hist_total->Write();
    curve_total->Write();
        
    hist_n0->Write();
    curve_n0->Write();

    hist_n1->Write();
    curve_n1->Write();
    
    hist_n2->Write();
    curve_n2->Write();
    
    fOut->Close();
    //write results in output text file
    ofstream myfile;
    myfile.open("fittingResults/"+isoName+"_results.txt");
    myfile<<"result copy to spreadsheet"<<endl;
    myfile<<"HL\t"<<"HL ErrL\t"<<"HL ErrH\t"<<"P1n\t"<<"P1n ErrL\t"<<"P1n ErrH\t"<<"P2n\t"<<"P2n ErrL\t"<<"P2n ErrH"<<"chi2/ndf"<<endl;
    Double_t hlr = log(2)/lamA.getVal();
    Double_t hlLr = log(2)/pow(lamA.getVal(),2)*lamA.getErrorLo();
    Double_t hlHr = log(2)/pow(lamA.getVal(),2)*lamA.getErrorHi();
    Double_t p1nr = p1A.getVal() * 100.0;
    Double_t p1nrL = p1A.getErrorLo() *100.0;
    Double_t p1nrH = p1A.getErrorHi() *100.0;
    
    Double_t p2nr = p2A.getVal() * 100.0;
    Double_t p2nrL = p2A.getErrorLo() *100.0;
    Double_t p2nrH = p2A.getErrorHi() *100.0;
    
    myfile<<hlr<<"\t"<<hlLr<<"\t"<<hlHr<<"\t"<<p1nr<<"\t"<<p1nrL<<"\t"<<p1nrH<<"\t"<<p2nr<<"\t"<<p2nrL<<"\t"<<p2nrH<<"\t"<<chi2_over_ndf<<endl;
    myfile<<"\nchi2 over ndf from total decay curve = "<<chi2_over_ndf<<endl;    
    myfile<<"\nNeutron Efficiency par for this fit  = "<<neuEff_par<<endl;
    myfile<<"0 means effective neutron eff"<<endl;
    myfile<<"1 means lower band uncertainty of neutron eff"<<endl;
    myfile<<"2 means higher band uncertainty of neutron eff"<<endl;
	
    //myfile<<"\n\nchi2 = "<<chi2_td<<endl; //this is not the correct way to get chi2, use values from manual method //commented on Jan15, 2023
    myfile<<"Neu Eff = "<<neuEff.getVal()<<endl;
    myfile<<"Neu Eff of for two delayed neutrons = "<<neudEff.getVal()<<endl;    
    myfile<<"N0 = "<<N0.getVal()<<endl;
    myfile<<"Half-life = "<<log(2)/lamA.getVal()<<" ms"<<endl;
    myfile<<"Half-life error Hi ="<<log(2)/pow(lamA.getVal(),2)*lamA.getErrorHi()<<endl;
	myfile<<"Half-life error Lo ="<<log(2)/pow(lamA.getVal(),2)*lamA.getErrorLo()<<endl;
    myfile<<"Prob of 1 neu emi ="<<p1A.getVal()<<endl;
    myfile<<"Prob of 1 neu emi error Hi ="<<p1A.getErrorHi()<<endl;
    myfile<<"Prob of 1 neu emi error Lo ="<<p1A.getErrorLo()<<endl;
    myfile<<"Prob of 2 neu emi ="<<p2A.getVal()<<endl;
    myfile<<"Prob of 2 neu emi error Hi ="<<p2A.getErrorHi()<<endl;
    myfile<<"Prob of 2 neu emi error Lo ="<<p2A.getErrorLo()<<endl;	
    myfile<<"bkg1 = "<<bkg1.getVal()<<endl;
    myfile<<"bkg2 = "<<bkg2.getVal()<<endl;
    myfile<<"bkg3 = "<<bkg3.getVal()<<endl;

    myfile<<"test1 NbinsX= "<<test1<<endl;
    myfile<<"test2 N= "<<test2<<endl;

    myfile<<"\n\n\n Used fitting parameters ---->"<<endl;
    myfile<<"\n\nIsotope A, T_1/2 = "<<T[isoA]<<"\tP1n = "<<Pn[isoA]<<"\tP2n = "<<Pnn[isoA]<<endl;
    myfile<<"Isotope B, T_1/2 = "<<T[isoB]<<"\tP1n = "<<Pn[isoB]<<"\tP2n = "<<Pnn[isoB]<<endl;
    myfile<<"Isotope C, T_1/2 = "<<T[isoC]<<"\tP1n = "<<Pn[isoC]<<"\tP2n = "<<Pnn[isoC]<<endl;
    myfile<<"Isotope D, T_1/2 = "<<T[isoD]<<"\tP1n = "<<Pn[isoD]<<"\tP2n = "<<Pnn[isoD]<<endl;
    myfile<<"Isotope E, T_1/2 = "<<T[isoE]<<"\tP1n = "<<Pn[isoE]<<"\tP2n = "<<Pnn[isoE]<<endl;
    myfile<<"Isotope F, T_1/2 = "<<T[isoF]<<"\tP1n = "<<Pn[isoF]<<"\tP2n = "<<Pnn[isoF]<<endl;
    myfile<<"Isotope G, T_1/2 = "<<T[isoG]<<"\tP1n = "<<Pn[isoG]<<"\tP2n = "<<Pnn[isoG]<<endl;
    myfile<<"Isotope H, T_1/2 = "<<T[isoH]<<"\tP1n = "<<Pn[isoH]<<"\tP2n = "<<Pnn[isoH]<<endl;
    myfile<<"Isotope I, T_1/2 = "<<T[isoI]<<"\tP1n = "<<Pn[isoI]<<"\tP2n = "<<Pnn[isoI]<<endl;
    myfile<<"Isotope J, T_1/2 = "<<T[isoJ]<<"\tP1n = "<<Pn[isoJ]<<"\tP2n = "<<Pnn[isoJ]<<endl;
    myfile<<"Isotope K, T_1/2 = "<<T[isoK]<<"\tP1n = "<<Pn[isoK]<<"\tP2n = "<<Pnn[isoK]<<endl;
    myfile<<"Isotope L, T_1/2 = "<<T[isoL]<<endl;
    myfile<<"r1 = "<<r1<<endl;
    myfile<<"r2 = "<<r2<<endl;
	

    
    //closing the output result text file
    myfile.close();
    
    //file about bkg and N0
    ofstream bkgfile;
    bkgfile.open("outputbkg/"+isoName+"_bkgoutput.txt");
    Double_t n90 = N0.getVal();
    Double_t bk1 = bkg1.getVal();
    Double_t bk2 = bkg2.getVal();
    Double_t bk3 = bkg3.getVal();
    bkgfile<<n90<<"\t"<<n90*0.90<<"\t"<<n90*1.10<<endl;
    bkgfile<<bk1<<"\t"<<bk1*0.90<<"\t"<<bk1*1.10<<endl;
    bkgfile<<bk2<<"\t"<<bk2*0.90<<"\t"<<bk2*1.10<<endl;
    bkgfile<<bk3<<"\t"<<bk3*0.90<<"\t"<<bk3*1.10<<endl;
    bkgfile.close();
    //*/new file having all parameters value for the drawing of decay spectras
    ofstream myfile2;
    myfile2.open("fittingResults_forDraw/"+isoName+"_result_forDrawing.txt");
    myfile2<<"0\t"<<N0.getVal()<<"\t"<<N0.getErrorLo()<<"\t"<<N0.getErrorHi()<<endl;
    //    myfile2<<"0\t"<<N0.getValV()<<"\t"<<N0.getErrorLo()<<"\t"<<N0.getErrorHi()<<endl;
    myfile2<<"1\t"<<lamA.getVal()<<"\t"<<lamA.getErrorLo()<<"\t"<<lamA.getErrorHi()<<endl;
    myfile2<<"2\t"<<p1A.getVal()<<"\t"<<p1A.getErrorLo()<<"\t"<<p1A.getErrorHi()<<endl;
    myfile2<<"3\t"<<p2A.getVal()<<"\t"<<p2A.getErrorLo()<<"\t"<<p2A.getErrorHi()<<endl;
    myfile2<<"4\t"<<lamB.getVal()<<"\t"<<lamB.getErrorLo()<<"\t"<<lamB.getErrorHi()<<endl;
    myfile2<<"5\t"<<p1B.getVal()<<"\t"<<p1B.getErrorLo()<<"\t"<<p1B.getErrorHi()<<endl;
    myfile2<<"6\t"<<p2B.getVal()<<"\t"<<p2B.getErrorLo()<<"\t"<<p2B.getErrorHi()<<endl;
    myfile2<<"7\t"<<lamC.getVal()<<"\t"<<lamC.getErrorLo()<<"\t"<<lamC.getErrorHi()<<endl;
    myfile2<<"8\t"<<p1C.getVal()<<"\t"<<p1C.getErrorLo()<<"\t"<<p1C.getErrorHi()<<endl;
    myfile2<<"9\t"<<p2C.getVal()<<"\t"<<p2C.getErrorLo()<<"\t"<<p2C.getErrorHi()<<endl;
    myfile2<<"10\t"<<lamD.getVal()<<"\t"<<lamD.getErrorLo()<<"\t"<<lamD.getErrorHi()<<endl;
    myfile2<<"11\t"<<p1D.getVal()<<"\t"<<p1D.getErrorLo()<<"\t"<<p1D.getErrorHi()<<endl;
    myfile2<<"12\t"<<p2D.getVal()<<"\t"<<p2D.getErrorLo()<<"\t"<<p2D.getErrorHi()<<endl;
    myfile2<<"13\t"<<lamE.getVal()<<"\t"<<lamE.getErrorLo()<<"\t"<<lamE.getErrorHi()<<endl;
    myfile2<<"14\t"<<p1E.getVal()<<"\t"<<p1E.getErrorLo()<<"\t"<<p1E.getErrorHi()<<endl;
    myfile2<<"15\t"<<p2E.getVal()<<"\t"<<p2E.getErrorLo()<<"\t"<<p2E.getErrorHi()<<endl;
    myfile2<<"16\t"<<lamF.getVal()<<"\t"<<lamF.getErrorLo()<<"\t"<<lamF.getErrorHi()<<endl;
    myfile2<<"17\t"<<p1F.getVal()<<"\t"<<p1F.getErrorLo()<<"\t"<<p1F.getErrorHi()<<endl;
    myfile2<<"18\t"<<p2F.getVal()<<"\t"<<p2F.getErrorLo()<<"\t"<<p2F.getErrorHi()<<endl;
    myfile2<<"19\t"<<lamG.getVal()<<"\t"<<lamG.getErrorLo()<<"\t"<<lamG.getErrorHi()<<endl;
    myfile2<<"20\t"<<p1G.getVal()<<"\t"<<p1G.getErrorLo()<<"\t"<<p1G.getErrorHi()<<endl;
    myfile2<<"21\t"<<p2G.getVal()<<"\t"<<p2G.getErrorLo()<<"\t"<<p2G.getErrorHi()<<endl;
    myfile2<<"22\t"<<lamH.getVal()<<"\t"<<lamH.getErrorLo()<<"\t"<<lamH.getErrorHi()<<endl;
    myfile2<<"23\t"<<p1H.getVal()<<"\t"<<p1H.getErrorLo()<<"\t"<<p1H.getErrorHi()<<endl;
    myfile2<<"24\t"<<p2H.getVal()<<"\t"<<p2H.getErrorLo()<<"\t"<<p2H.getErrorHi()<<endl;
    myfile2<<"25\t"<<lamI.getVal()<<"\t"<<lamI.getErrorLo()<<"\t"<<lamI.getErrorHi()<<endl;
    myfile2<<"26\t"<<p1I.getVal()<<"\t"<<p1I.getErrorLo()<<"\t"<<p1I.getErrorHi()<<endl;
    myfile2<<"27\t"<<p2I.getVal()<<"\t"<<p2I.getErrorLo()<<"\t"<<p2I.getErrorHi()<<endl;
    myfile2<<"28\t"<<lamJ.getVal()<<"\t"<<lamJ.getErrorLo()<<"\t"<<lamJ.getErrorHi()<<endl;
    myfile2<<"29\t"<<p1J.getVal()<<"\t"<<p1J.getErrorLo()<<"\t"<<p1J.getErrorHi()<<endl;
    myfile2<<"30\t"<<p2J.getVal()<<"\t"<<p2J.getErrorLo()<<"\t"<<p2J.getErrorHi()<<endl;
    myfile2<<"31\t"<<lamK.getVal()<<"\t"<<lamK.getErrorLo()<<"\t"<<lamK.getErrorHi()<<endl;
    myfile2<<"32\t"<<p1K.getVal()<<"\t"<<p1K.getErrorLo()<<"\t"<<p1K.getErrorHi()<<endl;
    myfile2<<"33\t"<<p2K.getVal()<<"\t"<<p2K.getErrorLo()<<"\t"<<p2K.getErrorHi()<<endl;
    myfile2<<"34\t"<<lamL.getVal()<<"\t"<<lamL.getErrorLo()<<"\t"<<lamL.getErrorHi()<<endl;
    myfile2<<"35\t"<<neuEff.getVal()<<"\t"<<neuEff.getErrorLo()<<"\t"<<neuEff.getErrorHi()<<endl;
    myfile2<<"36\t"<<bkg1.getVal()<<"\t"<<bkg1.getErrorLo()<<"\t"<<bkg1.getErrorHi()<<endl;
    // myfile2<<"36\t"<<bkg1.getValV()<<"\t"<<bkg1.getErrorLo()<<"\t"<<bkg1.getErrorHi()<<endl;
    myfile2<<"37\t"<<bkg2.getVal()<<"\t"<<bkg2.getErrorLo()<<"\t"<<bkg2.getErrorHi()<<endl;
    //myfile2<<"37\t"<<bkg2.getValV()<<"\t"<<bkg2.getErrorLo()<<"\t"<<bkg2.getErrorHi()<<endl;
    myfile2<<"38\t"<<bkg3.getVal()<<"\t"<<bkg3.getErrorLo()<<"\t"<<bkg3.getErrorHi()<<endl;
    myfile2<<"39\t"<<rc1.getVal()<<"\t"<<rc1.getErrorLo()<<"\t"<<rc1.getErrorHi()<<endl;
    myfile2<<"40\t"<<rc2.getVal()<<"\t"<<rc2.getErrorLo()<<"\t"<<rc2.getErrorHi()<<endl;
    myfile2<<"35\t"<<neudEff.getVal()<<"\t"<<neudEff.getErrorLo()<<"\t"<<neudEff.getErrorHi()<<endl;
    myfile2.close();
    
}

int main(){
  
  /*/////done part1
  fit_program("Y","Zr","Nb",39,65,"decayData_final/Y104_decaydata.root","bkg.txt",5,10000);
  fit_program("Y","Zr","Nb",39,66,"decayData_final/Y105_decaydata.root","bkg.txt",5,10000);
  fit_program("Y","Zr","Nb",39,67,"decayData_final/Y106_decaydata.root","bkg.txt",5,10000);
  
  fit_program("Y_neuEffLo","Zr","Nb",39,65,"decayData_final/Y104_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Y_neuEffLo","Zr","Nb",39,66,"decayData_final/Y105_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Y_neuEffLo","Zr","Nb",39,67,"decayData_final/Y106_decaydata.root","bkg.txt",5,10000,1);
  
  fit_program("Y_neuEffHi","Zr","Nb",39,65,"decayData_final/Y104_decaydata.root","bkg.txt",5,10000,2);
  fit_program("Y_neuEffHi","Zr","Nb",39,66,"decayData_final/Y105_decaydata.root","bkg.txt",5,10000,2);
  fit_program("Y_neuEffHi","Zr","Nb",39,67,"decayData_final/Y106_decaydata.root","bkg.txt",5,10000,2);
  
  fit_program("Y","Zr","Nb",39,68,"decayData_final/Y107_decaydata.root","bkg.txt",5,10000);
  fit_program("Y","Zr","Nb",39,69,"decayData_final/Y108_decaydata.root","bkg.txt",5,10000); 
  fit_program("Y","Zr","Nb",39,70,"decayData_final/Y109_decaydata.root","bkg.txt",5,10000); //use bkg3 0

  fit_program("Y_neuEffLo","Zr","Nb",39,68,"decayData_final/Y107_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Y_neuEffLo","Zr","Nb",39,69,"decayData_final/Y108_decaydata.root","bkg.txt",5,10000,1); 
  fit_program("Y_neuEffLo","Zr","Nb",39,70,"decayData_final/Y109_decaydata.root","bkg.txt",5,10000,1); //use p2n or bkg3 0

  fit_program("Y_neuEffHi","Zr","Nb",39,68,"decayData_final/Y107_decaydata.root","bkg.txt",5,10000,2);
  fit_program("Y_neuEffHi","Zr","Nb",39,69,"decayData_final/Y108_decaydata.root","bkg.txt",5,10000,2);
  fit_program("Y_neuEffHi","Zr","Nb",39,70,"decayData_final/Y109_decaydata.root","bkg.txt",5,10000,2);//use p2n bkg3 0
  /////// //////
  fit_program("Sr","Y","Zr",38,63,"decayData_final/Sr101_decaydata.root","bkg.txt",5,10000);
  fit_program("Sr_neuEffLo","Y","Zr",38,63,"decayData_final/Sr101_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Sr_neuEffHi","Y","Zr",38,63,"decayData_final/Sr101_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Sr","Y","Zr",38,64,"decayData_final/Sr102_decaydata.root","bkg.txt",5,10000);
  fit_program("Sr_neuEffLo","Y","Zr",38,64,"decayData_final/Sr102_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Sr_neuEffHi","Y","Zr",38,64,"decayData_final/Sr102_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Sr","Y","Zr",38,65,"decayData_final/Sr103_decaydata.root","bkg.txt",5,10000);
  fit_program("Sr_neuEffLo","Y","Zr",38,65,"decayData_final/Sr103_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Sr_neuEffHi","Y","Zr",38,65,"decayData_final/Sr103_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Rb","Sr","Y",37,62,"decayData_final/Rb99_decaydata.root","bkg.txt",5,10000);
  fit_program("Rb_neuEffLo","Sr","Y",37,62,"decayData_final/Rb99_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Rb_neuEffHi","Sr","Y",37,62,"decayData_final/Rb99_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Rb","Sr","Y",37,63,"decayData_final/Rb100_decaydata.root","bkg.txt",5,10000);
  fit_program("Rb_neuEffLo","Sr","Y",37,63,"decayData_final/Rb100_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Rb_neuEffHi","Sr","Y",37,63,"decayData_final/Rb100_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Rb","Sr","Y",37,64,"decayData_final/Rb101_decaydata.root","bkg.txt",5,10000);
  fit_program("Rb_neuEffLo","Sr","Y",37,64,"decayData_final/Rb101_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Rb_neuEffHi","Sr","Y",37,64,"decayData_final/Rb101_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Br","Kr","Rb",35,59,"decayData_final/Br94_decaydata.root","bkg.txt",5,10000);
  fit_program("Br_neuEffLo","Kr","Rb",35,59,"decayData_final/Br94_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Br_neuEffHi","Kr","Rb",35,59,"decayData_final/Br94_decaydata.root","bkg.txt",5,10000,2);
  
  fit_program("Se","Br","Kr",34,58,"decayData_final/Se92_decaydata.root","bkg.txt",5,10000);
  fit_program("Se_neuEffLo","Br","Kr",34,58,"decayData_final/Se92_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Se_neuEffHi","Br","Kr",34,58,"decayData_final/Se92_decaydata.root","bkg.txt",5,10000,2);
  ////// //part 1 
  
  /// //for Sr isotopes, part2
  fit_program("Sr","Y","Zr",38,66,"decayData_final/Sr104_decaydata.root","bkg.txt",5,10000);
  fit_program("Sr_neuEffLo","Y","Zr",38,66,"decayData_final/Sr104_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Sr_neuEffHi","Y","Zr",38,66,"decayData_final/Sr104_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Sr","Y","Zr",38,67,"decayData_final/Sr105_decaydata.root","bkg.txt",5,10000);
  fit_program("Sr_neuEffLo","Y","Zr",38,67,"decayData_final/Sr105_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Sr_neuEffHi","Y","Zr",38,67,"decayData_final/Sr105_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Sr","Y","Zr",38,68,"decayData_final/Sr106_decaydata.root","bkg.txt",5,10000);
  fit_program("Sr_neuEffLo","Y","Zr",38,68,"decayData_final/Sr106_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Sr_neuEffHi","Y","Zr",38,68,"decayData_final/Sr106_decaydata.root","bkg.txt",5,10000,2);
  / ///done for part2

  /////for Rb isotopes, part3
  fit_program("Rb","Sr","Y",37,65,"decayData_final/Rb102_decaydata.root","bkg.txt",5,10000);
  fit_program("Rb_neuEffLo","Sr","Y",37,65,"decayData_final/Rb102_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Rb_neuEffHi","Sr","Y",37,65,"decayData_final/Rb102_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Rb","Sr","Y",37,66,"decayData_final/Rb103_decaydata.root","bkg.txt",5,10000);
  fit_program("Rb_neuEffLo","Sr","Y",37,66,"decayData_final/Rb103_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Rb_neuEffHi","Sr","Y",37,66,"decayData_final/Rb103_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Rb","Sr","Y",37,67,"decayData_final/Rb104_decaydata.root","bkg.txt",5,10000);
  fit_program("Rb_neuEffLo","Sr","Y",37,67,"decayData_final/Rb104_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Rb_neuEffHi","Sr","Y",37,67,"decayData_final/Rb104_decaydata.root","bkg.txt",5,10000,2);
   ///part 3
  
  ///////for Kr istopes    
  fit_program("Kr","Rb","Sr",36,63,"decayData_final/Kr99_decaydata.root","bkg.txt",5,10000);
  fit_program("Kr_neuEffLo","Rb","Sr",36,63,"decayData_final/Kr99_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Kr_neuEffHi","Rb","Sr",36,63,"decayData_final/Kr99_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Kr","Rb","Sr",36,64,"decayData_final/Kr100_decaydata.root","bkg.txt",5,10000);
  fit_program("Kr_neuEffLo","Rb","Sr",36,64,"decayData_final/Kr100_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Kr_neuEffHi","Rb","Sr",36,64,"decayData_final/Kr100_decaydata.root","bkg.txt",5,10000,2);
  */
  fit_program("tttKr_neuEffHi","Rb","Sr",36,65,"decayData_final/Kr101_decaydata.root","bkg.txt",5,10000,2);
  fit_program("tttKr","Rb","Sr",36,65,"decayData_final/Kr101_decaydata.root","bkg.txt",5,10000);
  fit_program("tttKr_neuEffLo","Rb","Sr",36,65,"decayData_final/Kr101_decaydata.root","bkg.txt",5,10000,1);

  /*

  fit_program("Kr","Rb","Sr",36,66,"decayData_final/Kr102_decaydata.root","bkg.txt",5,10000);
  fit_program("Kr_neuEffLo","Rb","Sr",36,66,"decayData_final/Kr102_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Kr_neuEffHi","Rb","Sr",36,66,"decayData_final/Kr102_decaydata.root","bkg.txt",5,10000,2);
   //////kr -isotopes done
  
  ////Br isotopes
    //fit_program("Br","Kr","Rb",35,62,"decayData_final/Br97_decaydata.root","bkg.txt",5,10000);
    //fit_program("Br_neuEffLo","Kr","Rb",35,62,"decayData_final/Br97_decaydata.root","bkg.txt",5,10000,1);
    //fit_program("Br_neuEffHi","Kr","Rb",35,62,"decayData_final/Br97_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Br","Kr","Rb",35,63,"decayData_final/Br98_decaydata.root","bkg.txt",5,10000);
  fit_program("Br_neuEffLo","Kr","Rb",35,63,"decayData_final/Br98_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Br_neuEffHi","Kr","Rb",35,63,"decayData_final/Br98_decaydata.root","bkg.txt",5,10000,2);
  
  fit_program("Br","Kr","Rb",35,64,"decayData_final/Br99_decaydata.root","bkg.txt",5,10000);
  fit_program("Br_neuEffLo","Kr","Rb",35,64,"decayData_final/Br99_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Br_neuEffHi","Kr","Rb",35,64,"decayData_final/Br99_decaydata.root","bkg.txt",5,10000,2);
   
  *////

  /*///////last part
  fit_program("Se","Br","Kr",34,60,"decayData_final/Se94_decaydata.root","bkg.txt",5,10000);
  fit_program("Se_neuEffLo","Br","Kr",34,60,"decayData_final/Se94_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Se_neuEffHi","Br","Kr",34,60,"decayData_final/Se94_decaydata.root","bkg.txt",5,10000,2);

  fit_program("Se","Br","Kr",34,61,"decayData_final/Se95_decaydata.root","bkg.txt",5,10000);
  fit_program("Se_neuEffLo","Br","Kr",34,61,"decayData_final/Se95_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Se_neuEffHi","Br","Kr",34,61,"decayData_final/Se95_decaydata.root","bkg.txt",5,10000,2);
  
  fit_program("Se","Br","Kr",34,62,"decayData_final/Se96_decaydata.root","bkg.txt",5,10000);
  fit_program("Se_neuEffLo","Br","Kr",34,62,"decayData_final/Se96_decaydata.root","bkg.txt",5,10000,1);
  fit_program("Se_neuEffHi","Br","Kr",34,62,"decayData_final/Se96_decaydata.root","bkg.txt",5,10000,2);
  
  *//////last part


  
  return 0;
  
}
  


