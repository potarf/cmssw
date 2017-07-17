//include TMath;
const int fibermin=0,fibermax=15;
TFile *fileinput;
const int numberofruns=6;
TF1 *g1,*g2;
TH1F *hbreakdown;
TH2F *hbreakdown2d,*hbreakdown2dBV;
Double_t peak[2][6][fibermax-fibermin+1];
TGraphErrors *g[6][fibermax+1][3];
TGraphErrors *echo[6];
TMultiGraph *mg=new TMultiGraph();
TMultiGraph *mg2=new TMultiGraph();
//int runnumber[numberofruns]={4988,4989,4990,4991,4992,4993,4994,4995,4996};
//Double_t BV[numberofruns]={5.,5.05,5.1,5.15,5.2,5.25,5.3,5.35,5.4};
char name[100],label[100];
void runcompare_multishunt(int fiberN){
  // printf("start");
  for(int ifiberch=0;ifiberch<6;ifiberch++){
      echo[ifiberch]=new TGraphErrors;
      sprintf(name,"fiberch_%d",ifiberch);
      echo[ifiberch]->SetName(name);
      echo[ifiberch]->SetTitle(name);
  }
  echo[0]->GetYaxis()->SetRangeUser(0,0.5);
  for (int ifiberch=0;ifiberch<6;ifiberch++){
    sprintf (name,"ana_averagePS_run1000021920_EMAP-kalinin_HTR1_phi.root");
    fileinput= TFile::Open(name);
    if (fileinput==NULL) continue;
 

    //  sprintf(label,"pulse eta%d phi%d led=0",ifiberch,ifiber);
    sprintf(label,"pulse_eta%d_phi%d_led=0",ifiberch,fiberN);
    printf("%s\n",label);
    TGraphErrors *gg1=(TGraphErrors*)fileinput->Get(label);
    /*
    sprintf(label,"pulse_eta%d_phi%d_led=0",ifiberch,fiberN);
    printf("%s\n",label);
    TGraphErrors *gg2=(TGraphErrors*)fileinput->Get(label);
    if ((gg1==NULL)||(gg2==NULL))continue;
    if(gg1->GetN()==0)continue;
    if(gg2->GetN()==0)continue;
    gg2->SetMarkerColor(ifiberch+1);
    gg2->SetLineColor(ifiberch+1);
    gg2->SetMarkerStyle(20+ifiberch);
    mg2->Add(gg2); */
    Double_t x,y1,y2;
    for(int ipoint=1;ipoint<14;ipoint++){
      gg1->GetPoint(ipoint,x,y1);
      //     gg2->GetPoint(ipoint,x,y2);
      
      printf("%f %f %f\n",x,y1,y2);
      printf("%f %f\n",x,y1);
      //      if ((y1==0)||(y2==0))continue;
      //     echo[ifiberch]->SetPoint(echo[ifiberch]->GetN(),x,y1/y2);
    }
    //echo[ifiberch]->Draw("ap");

    fileinput->Close();
  }
  // echo[0]->GetXaxis()->SetTitle("Voltage");
  for (int ifiberch=0;ifiberch<6;ifiberch++){
    echo[ifiberch]->SetMarkerColor(ifiberch+1);
    echo[ifiberch]->SetLineColor(ifiberch+1);
    echo[ifiberch]->SetMarkerStyle(20+ifiberch);
    mg->Add(echo[ifiberch]);
  }
  sprintf (name,"signal 1TS%d",fiberN);
  TCanvas *c1=new TCanvas("c1","c1",0,0,800,800);
  mg2->SetTitle(name);
  /*  mg->GetXaxis()->SetTitle("Voltage");
  mg2->GetXaxis()->SetTitle("Voltage");
  mg->GetYaxis()->SetTitle("Ratio");
  mg2->GetYaxis()->SetTitle("Signal");
  c1->Divide(2);
  c1->cd(1);
  mg2->Draw("apl");
  c1->cd(2);
  mg->Draw("apl");
  */
}
