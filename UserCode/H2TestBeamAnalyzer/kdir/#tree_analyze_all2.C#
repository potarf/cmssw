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
#include <limits.h>

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
TH1F *hsinglePS[6][32][10][2];//rbx fiberch fiber ts led
TH1F *TSmax[6][32][2];//rbx fiberch fiber led
TH1F *hpulse[6][32][2];

float pulse[Numberofchannels][10];
int TSn,fibermin=32,fibermax=0,rbxmin=0,rbxmax=0;
int numberofrbx;
//int phi;
//int channel[10][6];
TH1F *avgPS[6][32][2];
//string listn ;
int ch2fiberch[Numberofchannels];
int ch2fiber[Numberofchannels];
int ch2rbx[Numberofchannels];
TH2F *spectr[2],*meansignalmap[2],*rmssignalmap[2];
TH2F *meansignalmap_BV[2],*rmssignalmap_BV[2];
int fiber2BV[32][6];
TH1F *gain1d[5];//= new TH1F("gain1d","gain1d",100,0,100);
TH2F *gain2d[5],*gain2d_BV[5];//
TH2F *gain2derror[5],*gain2derror_BV[5];
int number_TS_to_process=2;
int maxts(float TSvalue[10],int number_TS_to_process);
void DrawPlots(int rbxnumber,int whattodraw);
void DrawPS(int rbxnumber);
void DrawPedSpectra(int rbxnumber);
void DrawLedSpectra(int rbxnumber);
void SaveRBX2file(int rbxnumber,int runnumber);
void FitPedSpectra(int rbxnumber);
void AutoScaleSpectra(int rbxnumber);
void ReadBVmap();
void makeBVmaps(int rbxnumber);
void tree_analyze(int runnumber,int uhtr,int rbxnumber);

void tree_analyze_all(int runnumber,int uhtr){
    ReadBVmap();
    for (int rbxN=rbxmin;rbxN<rbxmax+1;rbxN++){
        tree_analyze(runnumber,uhtr,rbxN);
        FitPedSpectra(rbxN);
        AutoScaleSpectra(rbxN);
        makeBVmaps(rbxN);
        SaveRBX2file(rbxN, runnumber);
    }
}

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
    sprintf(name,"tb_chanmap_EMAP-kalinin_HTR%d_rbx.py",uhtr);
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
        sscanf(label1,"chanmap[%d,%d,%d] = %d",&etaN,&phiN,&depth,&channelN);
        // channel[etaN][phiN]=channelN;
        //   printf("channel[%d][%d]=%d\n",etaN,phiN,channelN);
        ch2fiberch[channelN-1]=depth;
        ch2fiber[channelN-1]=phiN;
        ch2rbx[channelN-1]=etaN;


    }
    fclose(in);
    //getting fibermin fibermax
    for(int channel=0;channel<Numberofchannels_emap;channel++){
        if (ch2fiber[channel]>fibermax) fibermax = ch2fiber[channel];
        if (ch2fiber[channel]<fibermin) fibermin = ch2fiber[channel];
    }
    printf("fiber %d - %d\n",fibermin,fibermax);
    //getting rbxmin rbxmax
    for(int channel=0;channel<Numberofchannels_emap;channel++){
        if (ch2rbx[channel]>rbxmax) rbxmax = ch2rbx[channel];
        if (ch2rbx[channel]<rbxmin) rbxmin = ch2rbx[channel];
    }
    printf("rbx %d - %d\n",rbxmin,rbxmax);
    ///create histos

    for(int iFIBERCHn=0;iFIBERCHn<6;iFIBERCHn++){
        for(int FIBERn=fibermin;FIBERn<fibermax+1;FIBERn++){
            for(int led=0;led<2;led++){
                for(int TSnn=0;TSnn<10;TSnn++){
                    //printf("charge_rbx%d_FIBERCH%d_FIBER%d_TS%d_led=%d\n",rbxnumber,iFIBERCHn,FIBERn,TSnn,led);
                    sprintf(label,"charge_rbx%d_fiberch%d_fiber%d_TS%d_led=%d",rbxnumber,iFIBERCHn,FIBERn,TSnn,led);
                    if(gROOT->FindObject(label))  hsinglePS[iFIBERCHn][FIBERn][TSnn][led]=(TH1F*)gROOT->FindObject(label);
                    //	 else	 hsinglePS[iFIBERCHn][FIBERn][TSnn][led]= new TH1F(label,label,1000000,0,1000000);
                    else	 hsinglePS[iFIBERCHn][FIBERn][TSnn][led]= new TH1F(label,label,255,low_edge);
                    hsinglePS[iFIBERCHn][FIBERn][TSnn][led]->Reset();

                }


                int sipma=fiber2BV[FIBERn][iFIBERCHn];


                sprintf(label,"avgPS_rbx%d_f%d__fch%d_led=%d_RN=%d_sipm=%d",rbxnumber,FIBERn,iFIBERCHn,led,runnumber,sipma);
                if(gROOT->FindObject(label))  avgPS[iFIBERCHn][FIBERn][led] = (TH1F*)gROOT->FindObject(label);
                else   avgPS[iFIBERCHn][FIBERn][led]= new TH1F(label,label,10,0,10);
                avgPS[iFIBERCHn][FIBERn][led]->Reset();

                sprintf(label,"TSmax_rbx%d_f%d_fch%d_led=%d_RN=%d_sipm=%d",rbxnumber,FIBERn,iFIBERCHn,led,runnumber,sipma);
                if(gROOT->FindObject(label))  TSmax[iFIBERCHn][FIBERn][led] = (TH1F*)gROOT->FindObject(label);
                else   TSmax[iFIBERCHn][FIBERn][led]= new TH1F(label,label,10,0,10);
                TSmax[iFIBERCHn][FIBERn][led]->Reset();

                sprintf(label,"rbx%d_f%d_fch%d_led=%d_RN=%d_sipm=%d",rbxnumber,FIBERn,iFIBERCHn,led,runnumber,sipma);
                if(gROOT->FindObject(label)) hpulse[iFIBERCHn][FIBERn][led] = (TH1F*)gROOT->FindObject(label);
                else   hpulse[iFIBERCHn][FIBERn][led]=new TH1F(label,label,400,0,800); // orig binning: 1000000/4, 1000000
                hpulse[iFIBERCHn][FIBERn][led]->Reset();

            }

        }
    }

    printf("next step1\n");

    for (int led=0;led<2;led++){
        //printf("idepth=%d led=%d\n",rbxnumber,led);
        sprintf(label,"spectraVS_fiber__rbx%d_led=%d",rbxnumber,led);
        if(gROOT->FindObject(label))  spectr[led] =  (TH2F*)gROOT->FindObject(label);
        else  spectr[led] = new TH2F(label,label,10*(fibermax-fibermin+1),fibermin*10,(fibermax+1)*10,400,0,800); // orig binning: 1000000/4
        spectr[led]->Reset();
        //printf("idepth=%d led=%d finished1\n",rbxnumber,led);
        sprintf(label,"mean_rbx%d_led=%d",rbxnumber,led);
        if(gROOT->FindObject(label)) meansignalmap[led] =  (TH2F*)gROOT->FindObject(label);
        else   meansignalmap[led] = new TH2F(label,label,fibermax-fibermin+1,fibermin,fibermax+1,6,0,6);
        meansignalmap[led]->Reset();

        sprintf(label,"mean_rbx%d_led=%d_BV",rbxnumber,led);
        if(gROOT->FindObject(label)) meansignalmap_BV[led] =  (TH2F*)gROOT->FindObject(label);
        else   meansignalmap_BV[led] = new TH2F(label,label,fibermax-fibermin+1,fibermin,fibermax+1,6,0,6);
        meansignalmap_BV[led]->Reset();

        //printf("idepth=%d led=%d finished2\n",rbxnumber,led);
        sprintf(label,"rms_rbx%d_led=%d",rbxnumber,led);
        if(gROOT->FindObject(label)) rmssignalmap[led] = (TH2F*)gROOT->FindObject(label);
        else  rmssignalmap[led] = new TH2F(label,label,fibermax-fibermin+1,fibermin,fibermax+1,6,0,6);
        rmssignalmap[led]->Reset();

        sprintf(label,"rms_rbx%d_led=%d_BV",rbxnumber,led);
        if(gROOT->FindObject(label)) rmssignalmap_BV[led] = (TH2F*)gROOT->FindObject(label);
        else  rmssignalmap_BV[led] = new TH2F(label,label,fibermax-fibermin+1,fibermin,fibermax+1,6,0,6);
        rmssignalmap_BV[led]->Reset();

        //printf("idepth=%d led=%d finished3\n",rbxnumber,led);

    }

    printf("next step 2\n");

    for (int i=0;i<5;i++){
        sprintf(label,"rbx%d_gain1d[%d]",rbxnumber,i);
        if(gROOT->FindObject(label)) gain1d[i] = (TH1F*)gROOT->FindObject(label);
        else  gain1d[i] = new TH1F(label,label,100,0,100);
        gain1d[i]->Reset();


        sprintf(label,"rbx%d_gain2d[%d]",rbxnumber,i);
        if(gROOT->FindObject(label)) gain2d[i] = (TH2F*)gROOT->FindObject(label);
        else  gain2d[i] = new TH2F(label,label,fibermax-fibermin+1,fibermin,fibermax+1,6,0,6);
        gain2d[i]->Reset();

        sprintf(label,"rbx%d_gain2derror[%d]",rbxnumber,i);
        if(gROOT->FindObject(label)) gain2derror[i] = (TH2F*)gROOT->FindObject(label);
        else  gain2derror[i] = new TH2F(label,label,fibermax-fibermin+1,fibermin,fibermax+1,6,0,6);
        gain2derror[i]->Reset();

        sprintf(label,"rbx%d_gain2d[%d]_BV",rbxnumber,i);
        if(gROOT->FindObject(label)) gain2d_BV[i] = (TH2F*)gROOT->FindObject(label);
        else  gain2d_BV[i] = new TH2F(label,label,fibermax-fibermin+1,fibermin,fibermax+1,6,0,6);
        gain2d_BV[i]->Reset();

        sprintf(label,"rbx%d_gain2derror[%d]_BV",rbxnumber,i);
        if(gROOT->FindObject(label)) gain2derror_BV[i] = (TH2F*)gROOT->FindObject(label);
        else  gain2derror_BV[i] = new TH2F(label,label,fibermax-fibermin+1,fibermin,fibermax+1,6,0,6);
        gain2derror_BV[i]->Reset();

    }

    sprintf(name,"ana_h2_tb_run%.6d_EMAP-kalinin_HTR%d_rbx.root",runnumber,uhtr);
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
        ch2fiberch[channelN]=idep[channelN];
        ch2fiber[channelN]=phi[channelN];
        //ch2rbx[channelN]=idep[channelN];
        ch2rbx[channelN]=ieta[channelN];
        // printf("channelN= %d ieta=%d   phi=%d\n",channelN,ieta[channelN],phi[channelN]);
    }
    //end mapping



    for (Int_t i=0; i<nentries; i++) {
        if(i%(100+(900*(i>999)))==0)  printf("processing %dth event\n",i);
        Events->GetEntry(i);


        for (channelN=0;channelN<Numberofchannels_emap;channelN++){
            sumTS=0;
            sumall=0;
            //here we will skip unused RBX
            if (ieta[channelN]!=rbxnumber) {printf ("!!%d\n",ieta[channelN]);continue;}


            int maxtsn=maxts(pulse[channelN],number_TS_to_process);
            TSmax[ch2fiberch[channelN]][ch2fiber[channelN]][led_trigger]->Fill(maxtsn);

            for(int TS=0;TS<10;TS++){

                // if (led_trigger==1){
                // if (led_trigger==1 ) hsinglePS[ch2fiberch[channelN]][ch2fiber[channelN]][TS][led_trigger]->Fill(pulse[channelN][TS]);
                //       if ((led_trigger==0)&&(pulse[channelN][TS]<50))
                hsinglePS[ch2fiberch[channelN]][ch2fiber[channelN]][TS][led_trigger]->Fill(pulse[channelN][TS]);

                if((TS>=maxtsn)&&(TS<maxtsn+number_TS_to_process)) sumTS+=pulse[channelN][TS];
                sumall+=pulse[channelN][TS];
                // }
            }
            if (led_trigger==0) {
                spectr[led_trigger]->Fill(ch2fiber[channelN]*10+ch2fiberch[channelN],sumTS);
                hpulse[ch2fiberch[channelN]][ch2fiber[channelN]][led_trigger]->Fill(sumTS);
            }
            if (led_trigger==1) {
                spectr[led_trigger]->Fill(ch2fiber[channelN]*10+ch2fiberch[channelN],sumTS);
                hpulse[ch2fiberch[channelN]][ch2fiber[channelN]][led_trigger]->Fill(sumTS);
            }

        }
    }

    //Filling Pulse shape histo

    for(int iFIBERCHn=0;iFIBERCHn<6;iFIBERCHn++){
        for(int FIBERn=fibermin;FIBERn<fibermax+1;FIBERn++){
            for (int led=0;led<2;led++){
                for(int TSnn=0;TSnn<10;TSnn++){
                    avgPS[iFIBERCHn][FIBERn][led]->Fill(TSnn,hsinglePS[iFIBERCHn][FIBERn][TSnn][led]->GetMean());
                    avgPS[iFIBERCHn][FIBERn][led]->SetBinError(TSnn+1,hsinglePS[iFIBERCHn][FIBERn][TSnn][led]->GetRMS());
                }
                meansignalmap[led]->SetBinContent(FIBERn-fibermin+1,iFIBERCHn+1,hpulse[iFIBERCHn][FIBERn][led]->GetMean());
                rmssignalmap[led]->SetBinContent(FIBERn-fibermin+1,iFIBERCHn+1,hpulse[iFIBERCHn][FIBERn][led]->GetRMS());
            }
        }
    }

    fileinput->Close();
}

void DrawPS(int rbxnumber){DrawPlots(rbxnumber,0);}
void DrawPedSpectra(int rbxnumber){DrawPlots(rbxnumber,1);}
void DrawLedSpectra(int rbxnumber){DrawPlots(rbxnumber,2);}
void SaveRBX2file(int rbxnumber,int runnumber){
    // sprintf(name,"ana_averagePS_run%.6d_EMAP-kalinin_HTR%d_phi.root",runnumber,uhtr);
    //fileout = TFile::Open(name);
    char filename[30];
    sprintf(filename,"analysis_run_%d.root",runnumber);
    fileout = new TFile(filename,"update");
    fileout->cd();
    sprintf(name,"rbx%d",rbxnumber);
    TDirectory *rbxdir=fileout->mkdir(name);
    rbxdir->cd();
    TDirectory *pulsedir = rbxdir->mkdir("pulse");
    pulsedir->cd();

    for(fiber=fibermin;fiber<fibermax+1;fiber++){
        for(fiberch=0;fiberch<6;fiberch++){

            hpulse[fiberch][fiber][0]->Write();
            hpulse[fiberch][fiber][1]->Write();

        }
    }

    TDirectory *PSdir = rbxdir->mkdir("PS");
    PSdir->cd();
    for(fiber=fibermin;fiber<fibermax+1;fiber++){
        for(fiberch=0;fiberch<6;fiberch++){

            avgPS[fiberch][fiber][0]->Write();
            avgPS[fiberch][fiber][1]->Write();
        }
    }
    TDirectory *overviewdir = rbxdir->mkdir("Overview");
    overviewdir->cd();
    gain2d[0]->Write();
    gain2d[1]->Write();

    gain1d[0]->Write();
    gain1d[1]->Write();
    gain2derror[0]->Write();
    gain2derror[1]->Write();
    gain2derror[2]->Write();

    spectr[0]->Write();
    meansignalmap[0]->Write();
    rmssignalmap[0]->Write();

    spectr[1]->Write();
    meansignalmap[1]->Write();
    rmssignalmap[1]->Write();

    meansignalmap_BV[0]->Write();
    rmssignalmap_BV[0]->Write();
    meansignalmap_BV[1]->Write();
    rmssignalmap_BV[1]->Write();
    gain2d_BV[0]->Write();
    gain2d_BV[1]->Write();
    gain2derror_BV[0]-> Write();
    gain2derror_BV[1]->Write();
    gain2derror_BV[2]->Write();

    fileout->Close();
    /*
 char hostname[HOST_NAME_MAX];
 gethostname(hostname, HOST_NAME_MAX);
 char command[250];

 if(strcmp(hostname, "cmskam05") == 0)
 {
     printf("Copying file to JSroot... \n");
     sprintf(command, "cp ./%d /www/html/files/%d",name,name);
 } else
 {
      printf("Enter daqweb password if you wish to copy this file to JSroot: \n");
      sprintf(command, "scp ./%d daqweb@cmskam05.cern.ch:/www/html/files/%d",name,name);
 }

 printf("Executing command %d", command);
 system(command);*/

}


void FitPedSpectra(int rbxnumber=0){
    for(fiber=fibermin;fiber<fibermax+1;fiber++){
        for(fiberch=0;fiberch<6;fiberch++){
            Double_t par[9];
            int whattodraw=1;
            //hpulse[rbxnumber][fiberch][fiber][0]->Rebin(4);
            int binmax=   hpulse[fiberch][fiber][0]->GetMaximumBin();
            int maximumvalue=  hpulse[fiberch][fiber][0]->GetBinCenter(binmax);
            hpulse[fiberch][fiber][0]->GetXaxis()->SetRangeUser(0,maximumvalue-4*(maximumvalue>4));

            while( (hpulse[fiberch][fiber][whattodraw-1]->Integral())>((hpulse[fiberch][fiber][whattodraw-1]->GetEntries())*0.02))
            {
                if ((hpulse[fiberch][fiber][whattodraw-1]->Integral())>=(hpulse[fiberch][fiber][whattodraw-1]->GetEntries())) break;

                binmax=   hpulse[fiberch][fiber][whattodraw-1]->GetMaximumBin();
                maximumvalue=  hpulse[fiberch][fiber][whattodraw-1]->GetBinCenter(binmax);
                if (maximumvalue<5) break;
                hpulse[fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(0,maximumvalue-4*(maximumvalue>4));


                // printf("%d %d %d %f \n",fiberch,fiber,maximumvalue,hpulse[fiberch][fiber][whattodraw-1]->Integral());
            }
            hpulse[fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(maximumvalue-20,maximumvalue+20);
            TF1 *g1    = new TF1("g1","gaus",maximumvalue-20,maximumvalue+20);
            hpulse[fiberch][fiber][whattodraw-1]->Fit(g1,"RQ");
            g1->GetParameters(&par[0]);
            hpulse[fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(maximumvalue+20,maximumvalue+100);
            binmax=   hpulse[fiberch][fiber][whattodraw-1]->GetMaximumBin();
            maximumvalue=  hpulse[fiberch][fiber][whattodraw-1]->GetBinCenter(binmax);
            TF1  *g2    = new TF1("g2","gaus",maximumvalue-15,maximumvalue+30);
            hpulse[fiberch][fiber][whattodraw-1]->Fit(g2,"RQ+");
            g2->GetParameters(&par[3]);
            hpulse[fiberch][fiber][whattodraw-1]->GetXaxis()->SetRangeUser(maximumvalue+20,maximumvalue+100);
            binmax=   hpulse[fiberch][fiber][whattodraw-1]->GetMaximumBin();
            maximumvalue=  hpulse[fiberch][fiber][whattodraw-1]->GetBinCenter(binmax);
            TF1  *g3    = new TF1("g3","gaus",maximumvalue-15,maximumvalue+30);
            hpulse[fiberch][fiber][whattodraw-1]->Fit(g3,"RQ+");
            g3->GetParameters(&par[6]);

            //set info concerning peak diff

            Double_t peak_diff;

            peak_diff = par[4] - par[1];

            char * title;
            char title2[50] = "";
            char title3[50] = "";
            title = strdup(hpulse[fiberch][fiber][whattodraw-1]->GetTitle());

            strcpy(title2,title);
            sprintf(title3,"_peakdelta=%.2f",peak_diff);
            strcat(title2,title3);
            hpulse[fiberch][fiber][whattodraw-1]->SetTitle(title2);


            //gain fit functions

            binmax=   gain1d[0]->GetMaximumBin();
            maximumvalue=  gain1d[0]->GetBinCenter(binmax);
            TF1 *gainfit1 = new TF1("gainfit1","gaus",maximumvalue-20,maximumvalue+20);
            gain1d[0]->Fit(gainfit1,"RQ");

            binmax = gain1d[1]->GetMaximumBin();
            maximumvalue=  gain1d[1]->GetBinCenter(binmax);
            TF1  *gainfit2 = new TF1("gainfit2","gaus",maximumvalue-20,maximumvalue+20);
            gain1d[1]->Fit(gainfit2,"RQ");


            gain1d[0]->Fill(par[4]-par[1]);

            gain2d[0]->SetBinContent(fiber-fibermin+1,fiberch+1,(par[4]-par[1]));
            gain2d[0]->GetZaxis()->SetRangeUser(40,60);
            gain1d[1]->Fill(par[7]-par[4]);

            gain2d[1]->SetBinContent(fiber-fibermin+1,fiberch+1,(par[7]-par[4]));
            gain2d[1]->GetZaxis()->SetRangeUser(40,60);

            gain2derror[0]->SetBinContent(fiber-fibermin+1,fiberch+1,g1->GetParError(1));
            gain2derror[1]->SetBinContent(fiber-fibermin+1,fiberch+1,g2->GetParError(1));
            gain2derror[2]->SetBinContent(fiber-fibermin+1,fiberch+1,g3->GetParError(1));
        }
    }
}

void  AutoScaleSpectra(int rbxnumber){
    for(fiber=fibermin;fiber<fibermax+1;fiber++){
        for(fiberch=0;fiberch<6;fiberch++){
            hpulse[fiberch][fiber][0]->GetXaxis()->SetRangeUser(0,hpulse[fiberch][fiber][0]->GetMean()*3);
            hpulse[fiberch][fiber][1]->GetXaxis()->SetRangeUser(0,hpulse[fiberch][fiber][1]->GetMean()*3);
        }
    }
    int maxrange[2]={0,0};
    for (int led=0;led<2;led++){
        for (int channelN=0;channelN<Numberofchannels_emap;channelN++){
            if (maxrange[led]<hpulse[ch2fiberch[channelN]][ch2fiber[channelN]][led]->GetMean())
                maxrange[led]=hpulse[ch2fiberch[channelN]][ch2fiber[channelN]][led]->GetMean();
        }
    }
    spectr[0]->GetYaxis()->SetRangeUser(0,maxrange[0]*4);
    spectr[1]->GetYaxis()->SetRangeUser(0,maxrange[1]*4);

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
            printf("fiber-fibermin: %d, fiberch+1:%d \n",fiber-fibermin,fiberch+1);


            //   c[fiber-fibermin]->GetPad(fiberch+1)->SetLogy();
            if (whattodraw==0)
            {
                int whattodrawfirst,TSmaxvalue=0,TSmax=0;
                for(int TS=0;TS<20;TS++)
                {
                    if  (avgPS[fiberch][fiber][int(TS/10)]->GetBinContent(TS%10)>TSmax) {TSmax=TS;TSmaxvalue=avgPS[fiberch][fiber][int(TS/10)]->GetBinContent(TS%10)>TSmax;}
                }
                avgPS[fiberch][fiber][1]->SetLineColor(kRed);
                avgPS[fiberch][fiber][int(TSmax/10)]->Draw();
                avgPS[fiberch][fiber][1-int(TSmax/10)]->Draw("same");

                //set stats
                //gStyle->SetOptStat();

                /*TPaveStats *ps = (TPaveStats*)c[fiber-fibermin]->GetPrimitive("stats");
                ps->SetName("mystats");
                TList *listOfLines = ps->GetListOfLines();
                TLatex *myt = new TLatex(0,0,"Test = 10");
                myt ->SetTextFont(42);
                myt ->SetTextSize(0.04);
                myt ->SetTextColor(kRed);
                listOfLines->Add(myt);
                avgPS[fiberch][fiber][int(TSmax/10)]->SetStats(0);
                c[fiber-fibermin]->Modified();*/
            }
            if ((whattodraw==1)||(whattodraw==2))
            {
                c[fiber-fibermin]->GetPad(fiberch+1)->SetLogy();

                hpulse[fiberch][fiber][whattodraw-1]->Draw();
            }
        }
    }
    //c[8]->SetLogz();

    //set stats

    for(fiber=fibermin;fiber<fibermax+1;fiber++){


    }

    c[fibermax-fibermin+1]->cd(1);
    c[fibermax-fibermin+1]->GetPad(1)->SetLogz();
    spectr[1]->Draw("colz");
    c[fibermax-fibermin+1]->cd(2);
    c[fibermax-fibermin+1]->GetPad(2)->SetLogz();
    spectr[0]->Draw("colz");

    c[fibermax-fibermin+2]->cd(1);
    //meansignalmap[1]->GetXaxis()->CenterLabels();
    //meansignalmap[1]->GetYaxis()->CenterLabels();
    meansignalmap[1]->SetStats(0);
    meansignalmap[1]->Draw("colz text");
    c[fibermax-fibermin+2]->cd(2);
    //meansignalmap[0]->GetXaxis()->CenterLabels();
    //meansignalmap[0]->GetYaxis()->CenterLabels();
    meansignalmap[0]->SetStats(0);
    meansignalmap[0]->Draw("colz text");
    c[fibermax-fibermin+2]->cd(3);
    //rmssignalmap[1]->GetXaxis()->CenterLabels();
    //rmssignalmap[1]->GetYaxis()->CenterLabels();
    rmssignalmap[1]->SetStats(0);
    rmssignalmap[1]->Draw("colz text");
    c[fibermax-fibermin+2]->cd(4);
    //rmssignalmap[0]->GetXaxis()->CenterLabels();
    //rmssignalmap[0]->GetYaxis()->CenterLabels();
    rmssignalmap[0]->SetStats(0);
    rmssignalmap[0]->Draw("colz text");

}

void ReadBVmap(){
    FILE *fileintxt;
    int RM,card,fiber,fiberch,QIE,BV,sipma,sipmx,sipmy;
    sprintf(name,"fiber2BV.txt");
    fileintxt = fopen (name, "rt");
    rewind(fileintxt);
    for (int ifiber=0;ifiber<32;ifiber++){
        for(int ifiberch=0;ifiberch<6;ifiberch++){
            fgets(label,50,fileintxt);
            sscanf(label,"%d\t%d\t%d\t%d\t%d\t%d",&RM,&card,&fiber,&fiberch,&QIE,&BV);
            // printf("%d %d %d %d\n",RM,fiber,fiberch,BV);
            //fiber2BV[(RM-1)*8+fiber+((fiber)%2)-1*((fiber%2)==0)-1][fiberch]=BV-1;
            //fiber2BV[(RM-1)*8+fiber-1][3*(fiberch<3)+fiberch-3*(fiberch>=3)]=BV-1;
            fiber2BV[(RM-1)*8+fiber-1][fiberch]=BV-1;
        }
    }
    fclose(fileintxt);
}

void makeBVmaps(int rbxnumber){
    int RM,card,fiber,fiberch,QIE,BV,sipma,sipmx,sipmy;
    for (int ifiber=fibermin;ifiber<fibermax+1;ifiber++){
        for (int ifiberch=0;ifiberch<6;ifiberch++){
            RM=int(ifiber/8);
            sipma=fiber2BV[ifiber][ifiberch];
            sipmx=sipma%8;
            sipmy=int(sipma/8);
            meansignalmap_BV[0]->SetBinContent(RM*8+sipmx+1,sipmy+1,meansignalmap[0]->GetBinContent(ifiber+1,ifiberch+1));
            rmssignalmap_BV[0]->SetBinContent(RM*8+sipmx+1,sipmy+1,rmssignalmap[0]->GetBinContent(ifiber+1,ifiberch+1));
            meansignalmap_BV[1]->SetBinContent(RM*8+sipmx+1,sipmy+1,meansignalmap[1]->GetBinContent(ifiber+1,ifiberch+1));
            rmssignalmap_BV[1]->SetBinContent(RM*8+sipmx+1,sipmy+1,rmssignalmap[1]->GetBinContent(ifiber+1,ifiberch+1));
            gain2d_BV[0]->SetBinContent(RM*8+sipmx+1,sipmy+1,gain2d[0]->GetBinContent(ifiber+1,ifiberch+1));
            gain2d_BV[1]->SetBinContent(RM*8+sipmx+1,sipmy+1,gain2d[1]->GetBinContent(ifiber+1,ifiberch+1));
            gain2derror_BV[0]->SetBinContent(RM*8+sipmx+1,sipmy+1,gain2derror[0]->GetBinContent(ifiber+1,ifiberch+1));
            gain2derror_BV[1]->SetBinContent(RM*8+sipmx+1,sipmy+1,gain2derror[1]->GetBinContent(ifiber+1,ifiberch+1));
            gain2derror_BV[2]->SetBinContent(RM*8+sipmx+1,sipmy+1,gain2derror[2]->GetBinContent(ifiber+1,ifiberch+1));
            gain2d_BV[0]->GetZaxis()->SetRangeUser(0,100);
            gain2d_BV[1]->GetZaxis()->SetRangeUser(0,100);

        }

    }
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
