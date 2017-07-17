//#include "Riostream.h"
#include <string>
#include <stdlib.h>
#include <iostream> 
#include <cstdio>
#include <cstdlib>
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include "TTree.h"
#include "TF1.h"
#include "TDirectory.h"
//#include <gRoot>
const int Numberofchannels=1000;
int Numberofchannels_emap;

	double const adc2fC[256] = {
			1.58, 4.73, 7.88, 11.0, 14.2, 17.3, 20.5, 23.6, 
			26.8, 29.9, 33.1, 36.2, 39.4, 42.5, 45.7, 48.8,
	  		53.6, 60.1, 66.6, 73.0, 79.5, 86.0, 92.5, 98.9,
	    		105, 112, 118, 125, 131, 138, 144, 151,
		 	157, 164, 170, 177, 186, 199, 212, 225,
		   	238, 251, 264, 277, 289, 302, 315, 328,
		     	341, 354, 367, 380, 393, 406, 418, 431,
			444, 464, 490, 516, 542, 568, 594, 620,
	  		569, 594, 619, 645, 670, 695, 720, 745,
	    		771, 796, 821, 846, 871, 897, 922, 947,
		  	960, 1010, 1060, 1120, 1170, 1220, 1270, 1320,
		    	1370, 1430, 1480, 1530, 1580, 1630, 1690, 1740,
	  		1790, 1840, 1890, 1940,  2020, 2120, 2230, 2330,
	    		2430, 2540, 2640, 2740, 2850, 2950, 3050, 3150,
		  	3260, 3360, 3460, 3570, 3670, 3770, 3880, 3980,
		    	4080, 4240, 4450, 4650, 4860, 5070, 5280, 5490,

	  		5080, 5280, 5480, 5680, 5880, 6080, 6280, 6480,
	    		6680, 6890, 7090, 7290, 7490, 7690, 7890, 8090,
		  	8400, 8810, 9220, 9630, 10000, 10400, 10900, 11300,
		    	11700, 12100, 12500, 12900, 13300, 13700, 14100, 14500,
	  		15000, 15400, 15800, 16200, 16800, 17600, 18400, 19300,
	    		20100, 20900, 21700, 22500, 23400, 24200, 25000, 25800,
		  	26600, 27500, 28300, 29100, 29900, 30700, 31600, 32400,
		    	33200, 34400, 36100, 37700, 39400, 41000, 42700, 44300,
	  		41100, 42700, 44300, 45900, 47600, 49200, 50800, 52500,
	    		54100, 55700, 57400, 59000, 60600, 62200, 63900, 65500,
	  		68000, 71300, 74700, 78000, 81400, 84700, 88000, 91400,
	    		94700, 98100, 101000, 105000, 108000, 111000, 115000, 118000,
	  		121000, 125000, 128000, 131000, 137000, 145000, 152000, 160000,
	    		168000, 176000, 183000, 191000, 199000, 206000, 214000, 222000,
	  		230000, 237000, 245000, 253000, 261000, 268000, 276000, 284000,
	    		291000, 302000, 316000, 329000, 343000, 356000, 370000, 384000
		};


TCanvas *c[40];
int fiber,fiberch;
char label[50],name[70];
TFile *fileinput,*fileout;
TH1F *hsinglePS[40][6][32][10][2];//rbx fiberch fiber ts led
TH1F *TSmax[40][6][32][2];//rbx fiberch fiber led
TH1F *hpulse[40][6][32][2];

float pulse[Numberofchannels][10];
int TSn,fibermin=32,fibermax=0,numberofrbx=3;
//int phi;
//int channel[10][6];
TH1F *avgPS[40][6][32][2];
//string listn ;
int ch2eta[Numberofchannels];
int ch2phi[Numberofchannels];
int ch2depth[Numberofchannels];
TH2F *spectr[40][2],*meansignalmap[40][2],*rmssignalmap[40][2];

TH1F *gain1d[40][5];//= new TH1F("gain1d","gain1d",100,0,100);
TH2F *gain2d[40][5];//
int number_TS_to_process=2;
int maxts(float TSvalue[10],int number_TS_to_process);
void DrawPlots(int rbxnumber,int whattodraw);
void DrawPS(int rbxnumber);
void DrawPedSpectra(int rbxnumber);
void DrawLedSpectra(int rbxnumber);
void SaveRBX2file(int rbxnumber,char *filename);
void FitPedSpectra(int rbxnumber);
void tree_analyze(int runnumber,int uhtr,int rbxnumber=0){
  // getting map from file
 ///
 ///
double  low_edge[256];
low_edge[0]=0;
 for(int i=1;i<256;i++){
   low_edge[i]=(adc2fC[i]-adc2fC[i-1])/2+adc2fC[i-1];//printf("%d %2f\t %2f\n",i,low_edge[i],adc2fC[i]);
   if (low_edge[i]<low_edge[i-1]){low_edge[i]=low_edge[i-1];printf("!!!%2f\n",low_edge[i]);}
 }//printf("\n");return;
 FILE * in;
 char label1[50];
 sprintf(name,"tb_chanmap_EMAP-kalinin_HTR%d_phi.py",uhtr);
 in = fopen (name, "rt");
 rewind(in);
 fgets(label1,50,in);
 fgets(label1,50,in);
 // printf("%s\n",label1);
 int etaN,phiN,channelN,led_trigger,depth;
 for(int channel=0;channel<Numberofchannels;channel++){
   
   fgets(label1,50,in);
   if (label1[0]=='\n') {printf("channels %d\n",channel);Numberofchannels_emap=channel;break;}// else {printf("exist%d\n",channel);}
   //  printf("%s\n",label1);
     //  printf("%s\n",listn);
   sscanf(label1,"chanmap[%d,%d,%d] = %d",&phiN,&etaN,&depth,&channelN);
     // channel[etaN][phiN]=channelN;
   //   printf("channel[%d][%d]=%d\n",etaN,phiN,channelN);
     ch2eta[channelN-1]=etaN;
     ch2phi[channelN-1]=phiN;
     
   
 }
 fclose(in);
 //getting fibermin fibermax
for(int channel=0;channel<Numberofchannels_emap;channel++){
  if (ch2phi[channel]>fibermax) fibermax = ch2phi[channel];
  if (ch2phi[channel]<fibermin) fibermin = ch2phi[channel];
 }
 printf("fiber %d - %d\n",fibermin,fibermax);
 ///create histos
 for(int idepth=0;idepth<numberofrbx;idepth++){
   for(int ietan=0;ietan<10;ietan++){
     for(int phin=fibermin;phin<fibermax+1;phin++){
       for(int led=0;led<2;led++){
	 for(int TSnn=0;TSnn<10;TSnn++){
	   
	   sprintf(label,"charge_rbx%d_eta%d_phi%d_TS%d_led=%d",idepth,ietan,phin,TSnn,led);
	   if(gROOT->FindObject(label))  hsinglePS[idepth][ietan][phin][TSnn][led]=(TH1F*)gROOT->FindObject(label);
	   //	 else	 hsinglePS[ietan][phin][TSnn][led]= new TH1F(label,label,1000000,0,1000000);
	   else	 hsinglePS[idepth][ietan][phin][TSnn][led]= new TH1F(label,label,255,low_edge);
	   hsinglePS[idepth][ietan][phin][TSnn][led]->Reset();
	   
	 }
	 
	 
	 sprintf(label,"avgPS_rbx%d_eta%d_phi%d_led=%d",idepth,ietan,phin,led);
	 if(gROOT->FindObject(label))  avgPS[idepth][ietan][phin][led] = (TH1F*)gROOT->FindObject(label);
	 else   avgPS[idepth][ietan][phin][led]= new TH1F(label,label,10,0,10);
	 avgPS[idepth][ietan][phin][led]->Reset();
	 
	 sprintf(label,"TSmax_rbx%d_eta%d_phi%d_led=%d",idepth,ietan,phin,led);
	 if(gROOT->FindObject(label))  TSmax[idepth][ietan][phin][led] = (TH1F*)gROOT->FindObject(label);
	 else   TSmax[idepth][ietan][phin][led]= new TH1F(label,label,10,0,10);
	 TSmax[idepth][ietan][phin][led]->Reset();
	 
	 sprintf(label,"pulse_rbx%d_eta%d_phi%d_led=%d",idepth,ietan,phin,led);
	 if(gROOT->FindObject(label)) hpulse[idepth][ietan][phin][led] = (TH1F*)gROOT->FindObject(label);
	 else   hpulse[idepth][ietan][phin][led]=new TH1F(label,label,1000000/4,0,1000000);
	 hpulse[idepth][ietan][phin][led]->Reset();
	 
       }
       
     }
   }
 }
 printf("next step1\n");
 for (int idepth=0;idepth<numberofrbx;idepth++){
   for (int led=0;led<2;led++){
     printf("idepth=%d led=%d\n",idepth,led);
     sprintf(label,"spectraVS_fiber__rbx%d_led=%d",idepth,led);
     if(gROOT->FindObject(label))  spectr[idepth][led] =  (TH2F*)gROOT->FindObject(label);
     else  spectr[idepth][led] = new TH2F(label,label,10*(fibermax-fibermin+1),fibermin*10,(fibermax+1)*10,1000000/4,0,1000000);
     spectr[idepth][led]->Reset();
     printf("idepth=%d led=%d finished1\n",idepth,led);
     sprintf(label,"mean_rbx%d_led=%d",idepth,led);
     if(gROOT->FindObject(label)) meansignalmap[idepth][led] =  (TH2F*)gROOT->FindObject(label);
     else   meansignalmap[idepth][led] = new TH2F(label,label,fibermax-fibermin+1,fibermin,fibermax+1,6,0,6);
     meansignalmap[idepth][led]->Reset();
     printf("idepth=%d led=%d finished2\n",idepth,led);
     sprintf(label,"rms_rbx%d_led=%d",idepth,led);
     if(gROOT->FindObject(label)) rmssignalmap[idepth][led] = (TH2F*)gROOT->FindObject(label);
     else  rmssignalmap[idepth][led] = new TH2F(label,label,fibermax-fibermin+1,fibermin,fibermax+1,6,0,6);
     rmssignalmap[idepth][led]->Reset();
     printf("idepth=%d led=%d finished3\n",idepth,led);
     
   }
 }
 printf("next step 2\n");
 for (int idepth=0;idepth<numberofrbx;idepth++){
   for (int i=0;i<5;i++){
     sprintf(label,"rbx%d_gain1d[%d]",idepth,i);
     if(gROOT->FindObject(label)) gain1d[idepth][i] = (TH1F*)gROOT->FindObject(label);
     else  gain1d[idepth][i] = new TH1F(label,label,100,0,100);
     gain1d[idepth][i]->Reset();
     
     sprintf(label,"rbx%d_gain2d[%d]",idepth,i);
     if(gROOT->FindObject(label)) gain2d[idepth][i] = (TH2F*)gROOT->FindObject(label);
     else  gain2d[idepth][i] = new TH2F(label,label,fibermax-fibermin+1,fibermin,fibermax+1,6,0,6);
     gain2d[idepth][i]->Reset();
   }
 }
 sprintf(name,"ana_h2_tb_run%.6d_EMAP-kalinin_HTR%d_phi.root",runnumber,uhtr);
 fileinput = TFile::Open(name);
 
 

 //read TTree
 int ieta[Numberofchannels],phi[Numberofchannels],idep[Numberofchannels];
 TTree *Events = (TTree*)fileinput->Get("QIE11Data/Events");
 Events->SetBranchAddress("pulse",&pulse);
 // Events->SetBranchAddress("TSn",&TSn);
 Events->SetBranchAddress("ieta",&ieta);
 Events->SetBranchAddress("iphi",&phi);
 Events->SetBranchAddress("depth",&idep);
 Events->SetBranchAddress("led_trigger",&led_trigger);
 int sumTS,sumall;
 
 //sprintf(cut_criteria,"ieta==%d && phi==%d",2,0);
 Int_t nentries = (Int_t)Events->GetEntries();
 printf("total number of events in tree = %d\n",nentries);
 // starting mapping
Events->GetEntry(0);
 for (channelN=0;channelN<Numberofchannels_emap;channelN++){
   ch2eta[channelN]=phi[channelN];
   ch2phi[channelN]=ieta[channelN];
   //ch2depth[channelN]=idep[channelN];
   ch2depth[channelN]=0;
   printf("channelN= %d ieta=%d   phi=%d\n",channelN,ieta[channelN],phi[channelN]); 
}
 //end mapping



 for (Int_t i=0; i<nentries; i++) {
   if(i%(100+(900*(i>999)))==0)  printf("processing %d event\n",i);
   Events->GetEntry(i);
 
    
   for (channelN=0;channelN<Numberofchannels_emap;channelN++){
     sumTS=0;
     sumall=0;
     //here we will skip unused RBX
     if (idep[channelN]!=rbxnumber) {continue;}


     int maxtsn=maxts(pulse[channelN],number_TS_to_process);
     TSmax[ch2depth[channelN]][ch2eta[channelN]][ch2phi[channelN]][led_trigger]->Fill(maxtsn);

     for(int TS=0;TS<10;TS++){
       
       // if (led_trigger==1){
       // if (led_trigger==1 ) hsinglePS[ch2eta[channelN]][ch2phi[channelN]][TS][led_trigger]->Fill(pulse[channelN][TS]);
       //       if ((led_trigger==0)&&(pulse[channelN][TS]<50)) 
       hsinglePS[ch2depth[channelN]][ch2eta[channelN]][ch2phi[channelN]][TS][led_trigger]->Fill(pulse[channelN][TS]);
     
       if((TS>=maxtsn)&&(TS<maxtsn+number_TS_to_process)) sumTS+=pulse[channelN][TS];
	 sumall+=pulse[channelN][TS];
	 // }
     }
     if (led_trigger==0) {
       spectr[ch2depth[channelN]][led_trigger]->Fill(ch2phi[channelN]*10+ch2eta[channelN],sumTS);
       hpulse[ch2depth[channelN]][ch2eta[channelN]][ch2phi[channelN]][led_trigger]->Fill(sumTS);
     }
      if (led_trigger==1) {
	spectr[ch2depth[channelN]][led_trigger]->Fill(ch2phi[channelN]*10+ch2eta[channelN],sumTS);
       hpulse[ch2depth[channelN]][ch2eta[channelN]][ch2phi[channelN]][led_trigger]->Fill(sumTS);
       }

   }
   }

 
 
 //Filling Pulse shape histo
     for (int idepth=0;idepth<numberofrbx;idepth++){
       for(int ietan=0;ietan<6;ietan++){
	 for(int phin=fibermin;phin<fibermax+1;phin++){
	   for (int led=0;led<2;led++){
	     for(int TSnn=0;TSnn<10;TSnn++){
	       avgPS[idepth][ietan][phin][led]->Fill(TSnn,hsinglePS[idepth][ietan][phin][TSnn][led]->GetMean());
	       avgPS[idepth][ietan][phin][led]->SetBinError(TSnn+1,hsinglePS[idepth][ietan][phin][TSnn][led]->GetRMS());
	     }
	     meansignalmap[idepth][led]->SetBinContent(phin-fibermin+1,ietan+1,hpulse[idepth][ietan][phin][led]->GetMean());
	     rmssignalmap[idepth][led]->SetBinContent(phin-fibermin+1,ietan+1,hpulse[idepth][ietan][phin][led]->GetRMS());
	   }
	 }
       }
     }
     

  fileinput->Close();


  //outrootfile




  //end outrootfile



}


void DrawPS(int rbxnumber){DrawPlots(rbxnumber,0);}
void DrawPedSpectra(int rbxnumber){DrawPlots(rbxnumber,1);}
void DrawLedSpectra(int rbxnumber){DrawPlots(rbxnumber,2);}
void SaveRBX2file(int rbxnumber, char *filename){
  // sprintf(name,"ana_averagePS_run%.6d_EMAP-kalinin_HTR%d_phi.root",runnumber,uhtr);
  //fileout = TFile::Open(name);
  fileout = new TFile(filename,"update");
  fileout->cd();
  sprintf(name,"rbx%d",rbxnumber);
  TDirectory *rbxdir=fileout->mkdir(name);
  rbxdir->cd();
  TDirectory *pulsedir = rbxdir->mkdir("pulse");
  pulsedir->cd();
  for(fiber=fibermin;fiber<fibermax+1;fiber++){
    for(fiberch=0;fiberch<6;fiberch++){
      
      hpulse[rbxnumber][fiberch][fiber][0]->Write();
      hpulse[rbxnumber][fiberch][fiber][1]->Write();
    }
  }

TDirectory *PSdir = rbxdir->mkdir("PS");
 PSdir->cd();
  for(fiber=fibermin;fiber<fibermax+1;fiber++){
    for(fiberch=0;fiberch<6;fiberch++){
      
      avgPS[rbxnumber][fiberch][fiber][0]->Write();
      avgPS[rbxnumber][fiberch][fiber][1]->Write();
    }
  }

  fileout->Close();

}


void FitPedSpectra(int rbxnumber=0){
  for(fiber=fibermin;fiber<fibermax+1;fiber++){
    for(fiberch=0;fiberch<6;fiberch++){
      Double_t par[9];
      int whattodraw=1;
      //hpulse[rbxnumber][fiberch][fiber][0]->Rebin(4);
      int binmax=   hpulse[rbxnumber][fiberch][fiber][0]->GetMaximumBin();
      int maximumvalue=  hpulse[rbxnumber][fiberch][fiber][0]->GetBinCenter(binmax);
      hpulse[rbxnumber][fiberch][fiber][0]->GetXaxis()->SetRangeUser(0,maximumvalue-4*(maximumvalue>4));
      
      while( (hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Integral())>((hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetEntries())*0.05))
	{
	  if ((hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Integral())>=(hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetEntries())) break;
	  
	  binmax=   hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetMaximumBin();
	  maximumvalue=  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetBinCenter(binmax);
	  if (maximumvalue<5) break;
	  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(0,maximumvalue-4*(maximumvalue>4));
	  
	  
	  // printf("%d %d %d %f \n",fiberch,fiber,maximumvalue,hpulse[fiberch][fiber][whattodraw-1]->Integral());
	}
      hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(maximumvalue-20,maximumvalue+20);
      TF1 *g1    = new TF1("g1","gaus",maximumvalue-20,maximumvalue+20);
      hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Fit(g1,"R");
      g1->GetParameters(&par[0]);
      hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(maximumvalue+20,maximumvalue+100);
      binmax=   hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetMaximumBin();
      maximumvalue=  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetBinCenter(binmax);
      TF1  *g2    = new TF1("g2","gaus",maximumvalue-15,maximumvalue+30);
      hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Fit(g2,"R+");
      g2->GetParameters(&par[3]);
      hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(maximumvalue+20,maximumvalue+100);
      binmax=   hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetMaximumBin();
      maximumvalue=  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetBinCenter(binmax);
      TF1  *g3    = new TF1("g3","gaus",maximumvalue-15,maximumvalue+30);
      hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Fit(g3,"R+");
      g3->GetParameters(&par[6]);
    }
  }
}


void DrawPlots(int rbxnumber=0,int whattodraw=0){
  //Drawing PS histos
  for (fiber=fibermin;fiber<fibermax+1+2;fiber++){
    sprintf(label,"fiber %d",fiber);
    if (fiber==fibermax+1) sprintf(label,"spectrVSchannel");
    if (fiber==fibermax+2) sprintf(label,"mean rms map");
    if (fiber==fibermax+3) sprintf(label,"reflectionmap");
    if (gROOT->FindObject(label)) c[fiber-fibermin] = (TCanvas*)gROOT->FindObject(label);
    else   
      {c[fiber-fibermin]=new TCanvas (label,label,0,0,400,400);
	if (fiber<fibermax+1) c[fiber-fibermin]->Divide(3,2);
	if (fiber==fibermax+1) c[fiber-fibermin]->Divide(1,2);
	if (fiber==fibermax+2)c[fiber-fibermin]->Divide(2,2);
      }
  }
  for(fiber=fibermin;fiber<fibermax+1;fiber++){
    for(fiberch=0;fiberch<6;fiberch++){
      c[fiber-fibermin]->cd(fiberch+1);
      //   c[fiber-fibermin]->GetPad(fiberch+1)->SetLogy();
      if (whattodraw==0)
	{
	  int whattodrawfirst,TSmaxvalue=0,TSmax=0;
	  for(int TS=0;TS<20;TS++)
	    {
	      if  (avgPS[rbxnumber][fiberch][fiber][int(TS/10)]->GetBinContent(TS%10)>TSmax) {TSmax=TS;TSmaxvalue=avgPS[rbxnumber][fiberch][fiber][int(TS/10)]->GetBinContent(TS%10)>TSmax;}
	    }
	  avgPS[rbxnumber][fiberch][fiber][1]->SetLineColor(kRed);
	  avgPS[rbxnumber][fiberch][fiber][int(TSmax/10)]->Draw();
	  avgPS[rbxnumber][fiberch][fiber][1-int(TSmax/10)]->Draw("same");
	}
      if ((whattodraw==1)||(whattodraw==2))
	{  
	  c[fiber-fibermin]->GetPad(fiberch+1)->SetLogy();
	  Double_t par[9];
	  //   hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Rebin(4);
	  int binmax=   hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetMaximumBin();
	  int maximumvalue=  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetBinCenter(binmax);
	  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(0,maximumvalue-4*(maximumvalue>4));

	  while( (hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Integral())>((hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetEntries())*0.02))
	    {
	      if ((hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Integral())>=(hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetEntries())) break;

	      binmax=   hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetMaximumBin();
	      maximumvalue=  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetBinCenter(binmax);
	      if (maximumvalue<5) break;
	      hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(0,maximumvalue-4*(maximumvalue>4));
	    
	   
	      // printf("%d %d %d %f \n",fiberch,fiber,maximumvalue,hpulse[fiberch][fiber][whattodraw-1]->Integral());
	    }
	  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(maximumvalue-20,maximumvalue+20);
	  // hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Rebin(4);
	  binmax=   hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetMaximumBin();
	  maximumvalue=  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetBinCenter(binmax);
	  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(maximumvalue-20,maximumvalue+20);
	  TF1 *g1    = new TF1("g1","gaus",maximumvalue-10,maximumvalue+20);
	  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Fit(g1,"R");
	  g1->GetParameters(&par[0]);
	  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(maximumvalue+20,maximumvalue+100);
	  binmax=   hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetMaximumBin();
	  maximumvalue=  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetBinCenter(binmax);
	  TF1  *g2    = new TF1("g2","gaus",maximumvalue-15,maximumvalue+30);
	  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Fit(g2,"R+");
	  g2->GetParameters(&par[3]);
	  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(maximumvalue+20,maximumvalue+100);
	  binmax=   hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetMaximumBin();
	  maximumvalue=  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetBinCenter(binmax);
	  TF1  *g3    = new TF1("g3","gaus",maximumvalue-15,maximumvalue+30);
	  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Fit(g3,"R+");
	  g3->GetParameters(&par[6]);
	  // TF1 *total = new TF1("total","g1(0)+g2(3)",0,maximumvalue+20);
	  // total->SetParameters(par);
	  //hpulse[fiberch][fiber][whattodraw-1]->Fit(g1,"R");
	  // hpulse[fiberch][fiber][whattodraw-1]->Fit(g2,"R+");
	  // hpulse[fiberch][fiber][whattodraw-1]->Fit(total,"R+");
	  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(0,hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->GetMean()*3);
	  hpulse[rbxnumber][fiberch][fiber][whattodraw-1]->Draw();

	  //begin gain calculation
	 
	  //if (fiber<8) 
	      gain1d[rbxnumber][0]->Fill(par[4]-par[1]);

	  gain2d[rbxnumber][0]->SetBinContent(fiber-fibermin+1,fiberch+1,(par[4]-par[1]));
	  gain2d[rbxnumber][0]->GetZaxis()->SetRangeUser(40,60);
	  gain1d[rbxnumber][1]->Fill(par[7]-par[4]);

	  gain2d[rbxnumber][1]->SetBinContent(fiber-fibermin+1,fiberch+1,(par[7]-par[4]));
	  gain2d[rbxnumber][1]->GetZaxis()->SetRangeUser(40,60);
	  // end gain calculation

	}
    }
  }
  //c[8]->SetLogz();
 
  int maxrange[2]={0,0};
  for (int led=0;led<2;led++){
    for (int channelN=0;channelN<Numberofchannels_emap;channelN++){
      if (maxrange[led]<hpulse[rbxnumber][ch2eta[channelN]][ch2phi[channelN]][led]->GetMean()) 
	maxrange[led]=hpulse[rbxnumber][ch2eta[channelN]][ch2phi[channelN]][led]->GetMean();
    }
  }
  spectr[rbxnumber][0]->GetYaxis()->SetRangeUser(0,maxrange[0]*4);
  spectr[rbxnumber][1]->GetYaxis()->SetRangeUser(0,maxrange[1]*4);
  c[fibermax-fibermin+1]->cd(1);
  c[fibermax-fibermin+1]->GetPad(1)->SetLogz();
  spectr[rbxnumber][1]->Draw("colz");
  c[fibermax-fibermin+1]->cd(2);
  c[fibermax-fibermin+1]->GetPad(2)->SetLogz();
  spectr[rbxnumber][0]->Draw("colz");
  
  c[fibermax-fibermin+2]->cd(1);
  //meansignalmap[1]->GetXaxis()->CenterLabels();
  //meansignalmap[1]->GetYaxis()->CenterLabels();
  meansignalmap[rbxnumber][1]->SetStats(0);
  meansignalmap[rbxnumber][1]->Draw("colz text");
  c[fibermax-fibermin+2]->cd(2);
  //meansignalmap[0]->GetXaxis()->CenterLabels();
  //meansignalmap[0]->GetYaxis()->CenterLabels();
  meansignalmap[rbxnumber][0]->SetStats(0);
  meansignalmap[rbxnumber][0]->Draw("colz text");
  c[fibermax-fibermin+2]->cd(3);
  //rmssignalmap[1]->GetXaxis()->CenterLabels();
  //rmssignalmap[1]->GetYaxis()->CenterLabels();
  rmssignalmap[rbxnumber][1]->SetStats(0);
  rmssignalmap[rbxnumber][1]->Draw("colz text");
  c[fibermax-fibermin+2]->cd(4);
  //rmssignalmap[0]->GetXaxis()->CenterLabels();
  //rmssignalmap[0]->GetYaxis()->CenterLabels();
  rmssignalmap[rbxnumber][0]->SetStats(0);
  rmssignalmap[rbxnumber][0]->Draw("colz text");

} 
int maxts(float TSvalue[10],int number_TS_to_process){
  // return 2;
  int max[10]={0,0,0,0,0,0,0,0,0,0};
  int currentmax=0,currentmin=0,currentmaxcounter=0;
  for(TSn=0;TSn<10;TSn++){
    if (max[currentmax]<TSvalue[TSn]){currentmax=TSn;currentmin=TSn;max[currentmax]=TSvalue[currentmax];}
  }
  currentmaxcounter++;
  while (((currentmaxcounter<=number_TS_to_process)&&(currentmin>=0)&&(currentmax<=9))){
  
  if ((currentmaxcounter==number_TS_to_process)||(currentmin==0)) return currentmin;
  if (currentmax==9) return currentmax-number_TS_to_process;

  if (TSvalue[currentmin-1]>TSvalue[currentmax+1]){currentmin--;} else {currentmax++;}
  currentmaxcounter++;
 
  }

  return 0;
}
