//include TMath;
const int fibermin=0,fibermax=31;
TFile *fileinput;
const int numberofruns=5;
TF1 *g1,*g2;
TH1F *hbreakdown;
TH2F *hbreakdown2d,*hbreakdown2dBV;
Double_t peak[2][6][fibermax-fibermin+1];
TGraphErrors *g[6][fibermax-fibermin+1];
int runnumber[numberofruns]={1000021569,1000021570,1000021571,1000021572,1000021573};
  Double_t BV[numberofruns]={69.,68.5,68.,67.5,67.};
char name[100],label[100];
void treecompare(){
  for (int ifiber=fibermin;ifiber<fibermax+1;ifiber++){
      for(int ifiberch=0;ifiberch<6;ifiberch++){
	sprintf(label,"pulse eta%d phi%d graph",ifiberch,ifiber);
	g[ifiberch][ifiber]=new TGraphErrors();

      }
  }
  for (int irun=0;irun<1;irun++){
    //  sprintf (name,"ana_averagePS_run%d_EMAP-kalinin_HTR1_phi.root",runnumber[irun]);
    sprintf (name,"test4.root");
    fileinput= TFile::Open(name);
    for (int ifiber=fibermin;ifiber<fibermax+1;ifiber++){
      for(int ifiberch=0;ifiberch<6;ifiberch++){
	Double_t par[9];
	sprintf(label,"rbx%d/pulse",0);
	fileinput->cd(label);
	sprintf(label,"pulse_rbx%d_eta%d_phi%d_led=0",8,ifiberch,ifiber);

	TH1F *hpulse=(TH1F*)gDirectory->Get(label);
	g1 = (TF1*)hpulse->FindObject("g1");
	g2 = (TF1*)hpulse->FindObject("g2");
	if (g1!=NULL)g1->GetParameters(&par[0]);
	if (g2!=NULL)g2->GetParameters(&par[3]);
	if((g1!=NULL)&&(g2!=NULL)){
	  if ((par[4]<0)||(par[1]<0)) continue;
	printf("%f\n",par[4]-par[1]);
	g[ifiberch][ifiber]->SetPoint(g[ifiberch][ifiber]->GetN(),BV[irun],par[4]-par[1]);
	g[ifiberch][ifiber]->SetPointError(g[ifiberch][ifiber]->GetN()-1,0.025,TMath::Sqrt(g1->GetParError(1)*g1->GetParError(1)+g2->GetParError(1)*g2->GetParError(1)));
	sprintf(label,"pulse eta%d phi%d graph",ifiberch,ifiber);

	g[ifiberch][ifiber]->SetName(label);
	g[ifiberch][ifiber]->SetTitle(label);
	g[ifiberch][ifiber]->SetMarkerStyle(20);

	}

      }
    }
    fileinput->Close();
  }
  for (int ifiber=fibermin;ifiber<fibermax+1;ifiber++){
    for(int ifiberch=0;ifiberch<6;ifiberch++){
      TF1 *gpol    = new TF1("gpol","pol1",60.,70.);
      g[ifiberch][ifiber]->Fit("gpol");
	Double_t par[2];
	gpol->GetParameters(&par[0]);
	peak[0][ifiberch][ifiber]=par[0];
	peak[1][ifiberch][ifiber]=par[1];
	printf("fiber%d fiberch%d pol %fx+%f\n",ifiber,ifiberch,par[1],par[0]);
	g[ifiberch][ifiber]->GetXaxis()->SetRangeUser(63,70);
	g[ifiberch][ifiber]->GetYaxis()->SetRangeUser(0,70);
	}
  }
  hbreakdown= new TH1F("breakdown","breakdown",50,45,70);
  hbreakdown2d= new TH2F("breakdown2d","breakdown2d",32,0,32,6,0,6);
  hbreakdown2dBV= new TH2F("breakdown2dBV","breakdown2dBV",32,0,32,6,0,6);
  hbreakdown2d->SetStats(0);
  hbreakdown2dBV->SetStats(0);
 
  gStyle->SetPaintTextFormat("2.1f");
    FILE *fileintxt;
  int fiber2BV[32][6];
  int RM,card,fiber,fiberch,QIE,BV,sipma,sipmx,sipmy;
 sprintf(name,"fiber2BV.txt");
 fileintxt = fopen (name, "rt");
rewind(fileintxt);
 for (int ifiber=fibermin;ifiber<fibermax+1;ifiber++){
    for(int ifiberch=0;ifiberch<6;ifiberch++){
      fgets(label,50,fileintxt);
      sscanf(label,"%d\t%d\t%d\t%d\t%d\t%d",&RM,&card,&fiber,&fiberch,&QIE,&BV);
      printf("%d %d %d %d\n",RM,fiber,fiberch,BV);
      //fiber2BV[(RM-1)*8+fiber+((fiber)%2)-1*((fiber%2)==0)-1][fiberch]=BV-1;
      //fiber2BV[(RM-1)*8+fiber-1][3*(fiberch<3)+fiberch-3*(fiberch>=3)]=BV-1;
      fiber2BV[(RM-1)*8+8-fiber][fiberch]=BV-1;
    }
    }
    fclose(fileintxt);
  char *fout="runcompare.txt";
  FILE *fc=fopen(fout,"w+");
  fprintf(fc,"fiber fiberch slope offset breakdown\n");
 for (int ifiber=fibermin;ifiber<fibermax+1;ifiber++){
    for(int ifiberch=0;ifiberch<6;ifiberch++){
      fprintf(fc,"%d %f\n",fiber2BV[ifiber][ifiberch],(50-peak[0][ifiberch][ifiber])/peak[1][ifiberch][ifiber]);
printf("%d %d %f %f %f\n",ifiber,ifiberch,peak[1][ifiberch][ifiber],peak[0][ifiberch][ifiber],-peak[0][ifiberch][ifiber]/peak[1][ifiberch][ifiber]);

 hbreakdown->Fill(-peak[0][ifiberch][ifiber]/peak[1][ifiberch][ifiber]);
 RM=int(ifiber/8);
 sipma=fiber2BV[ifiber][ifiberch];
 sipmx=sipma%8;
 sipmy=int(sipma/8);
 printf("RM %d %d %d\n",RM,sipmx,sipmy);
 hbreakdown2dBV->SetBinContent(RM*8+sipmx+1,sipmy+1,(50-peak[0][ifiberch][ifiber])/peak[1][ifiberch][ifiber]);
 hbreakdown2d->SetBinContent(ifiber+1,ifiberch+1,(50-peak[0][ifiberch][ifiber])/peak[1][ifiberch][ifiber]);
    }
 }
 fclose(fc);
}
