#include <sstream> 
#include <TGraph.h>
#include <TGraphErrors.h>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <TSystem.h> 
#include <TString.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TStyle.h>
#include <TFile.h>
#include <TLegend.h>
#include <TH1F.h>
#include <stdio.h>      /* printf */
#include <math.h>       /* floor */
#include <TChain.h>
#include <TTree.h>
#include <TF1.h>
#include <TNtuple.h>
#include <iostream>     //for using cout
#include <stdlib.h>     //for using the function sleep
using namespace std;

const 	int num_bins=20;
const   int ncuts=5;
		// run one deno run 2 is num
void X_bins_q(int run1=0, int run2=0,int quite =0,int quick=1) {
		time_t start = time(0) ;

	cout <<"\n";



	if(run1==0||run2==0){
	  cout << "Please input the run numbers of the disturbution you like to look at the ratios of (1st/2nd)." <<endl;
	 	cin >> run1 ; cin >> run2 ;cout <<endl;
	}

	time_t rawtime;
  	struct tm * timeinfo;
  	time (&rawtime);
  	timeinfo = localtime (&rawtime);
  	printf("Started on %s \n",asctime(timeinfo));

// Nameing the two runs/////
	double ones = ceil(((run1/100.0 - run1/100)-.000001)*100);
	int oness = ones; char* one;

	switch(oness){ 	case 2:  one =(char *)"D"; break;
					case 3:  one =(char *)"He3"; break;
					case 4:  one =(char *)"He4"; break;
					case 9:  one =(char *)"Be9"; break;
					case 12: one =(char *)"C12"; break;
					case 16: one =(char *)"O16"; break;
					case 40: one =(char *)"Ca40"; break;
					default: one =(char *)"No label "  ; }

	double twos = ceil(((run2/100.0 - run2/100)-.000001)*100);
	if(run2==4002){twos = floor((run2/100.0 - run2/100)*100);}
	int twoss = twos; char* two;

	switch(twoss){ 	case 2:  two =(char *)"D"; break;
					case 3:  two =(char *)"He3"; break;
					case 4:  two =(char *)"He4"; break;
					case 9:  two =(char *)"Be9"; break;
					case 12: two =(char *)"C12"; break;
					case 16: two =(char *)"O16"; break;
					case 40: two =(char *)"Ca40"; break;
					default: two =(char *)"No label"; }
///////////////////////////////////////////////////////////////////////

//Env vars...
  char* data_dir;
  	data_dir = getenv ("OUTPUT_DIR");
  char* root_dir;
  	root_dir = getenv ("OUT_DIR");
  char* img_dir;
  	img_dir = getenv ("IMG_DIR");
/////////////////////
	double min,max,count,count2;


//Double checking files
	ifstream run_1;
	run_1.open(Form("%s/IS_%d.root",root_dir,run1));
	ifstream run_2;
	run_2.open(Form("%s/IS_%d.root",root_dir,run2));
	if (!run_1.good()){cout << "The root file of run "<< run1 <<" does not exist."<<endl; return; }
	if (!run_2.good()){cout << "The root file of run "<< run2 <<" does not exist."<<endl; return; }
	run_1.close(); run_2.close();
///////////////////////////////////////////////////////////////////

//Opening ROOT files
	TFile *file_1 = new TFile(Form("%s/IS_%d.root",root_dir,run1),"read");
 	TFile *file_2 = new TFile(Form("%s/IS_%d.root",root_dir,run2),"read");
	TFile *f = new TFile(Form("%s/ratio_w_%d_%d.root",root_dir,run1,run2),"RECREATE");
///////////////////////////////////////////////////////////////

///tree vars
	double xb1,xb2,sig1,sig2, theta,theta2;
	//Float_t xb1,xb2,sig1,sig2, theta,theta2;
	float q1,q2;

///Make trees for two runs
 	TTree *T2 = (TTree*)file_2->Get("tree");
 		T2->SetBranchStatus("*",0);
		T2->SetBranchStatus("lambda",1);
		T2->SetBranchStatus("Beam",1);
		T2->SetBranchStatus("efinal_lab",1);
		T2->SetBranchStatus("qsquared",1);
		T2->SetBranchAddress("qsquared",&q1);
		T2->SetBranchStatus("invar_mass",1);
		T2->SetBranchStatus("theta_lab",1);
		T2->SetBranchStatus("theta_rest",1);
		T2->SetBranchAddress("theta_rest",&theta);
		//T2->SetBranchStatus("nuw_rest",1);
		T2->SetBranchStatus("Sigma_rest",1);
		T2->SetBranchAddress("Sigma_rest",&sig1);
		T2->SetBranchStatus("Xb_lab",1);
		T2->SetBranchAddress("Xb_lab",&xb1);
		T2->SetBranchStatus("Prot_momentum",1);


	TTree *tree = (TTree*)file_1->Get("tree");
		tree->SetBranchStatus("*",0);

		tree->SetBranchStatus("Beam",1);
		tree->SetBranchStatus("efinal_lab",1);
		tree->SetBranchStatus("lambda",1);
		tree->SetBranchStatus("qsquared",1);
		tree->SetBranchAddress("qsquared",&q2);
		tree->SetBranchStatus("invar_mass",1);
		tree->SetBranchStatus("theta_lab",1);
		tree->SetBranchStatus("theta_rest",1);
		tree->SetBranchAddress("theta_rest",&theta2);
		//tree->SetBranchStatus("nuw_rest",1);
		tree->SetBranchStatus("Sigma_rest",1);
		tree->SetBranchAddress("Sigma_rest",&sig2);
		tree->SetBranchStatus("Xb_lab",1);
		tree->SetBranchAddress("Xb_lab",&xb2);
		tree->SetBranchStatus("Prot_momentum",1);
//////////////////////////////////////////////////////////////////////////////////

//Cuts and varibles for cuts;
	double Qsq[5] ={0.7, 1.7 ,2.3 ,3.30 ,5.0};
	double DQsq[ncuts]={0.5,0.5,0.15,0.15,0.15};
	//double Qsq[ncuts] ={10.0, 20.0 ,30.0 ,35.0 ,45.0};
	//double DQsq[ncuts]={0.25,0.25,0.25,0.25,0.25};
	TCut Q[5];	
		for(int qi=0;qi<ncuts;qi++){Q[qi]= Form("(abs(qsquared-%f)/qsquared)<=%f",Qsq[qi],DQsq[qi]);}

//Cut constants
	double Qmax =6;	
	double Qmin =2;
	double Wmax = 500;
	double Wmin = 1.5;
	double theta_min = 35*3.14159/180;
	double theta_max = 50*3.14159/180;
	double deltaEmax = 5;
	double deltaEmin = 0;

	TCut QQ = Form("qsquared>=%f&&qsquared<=%f",Qmin,Qmax);
	TCut W  = Form("invar_mass>=%f&&invar_mass<=%f",Wmin,Wmax);
	TCut Theta = Form("theta_rest>=%f&&theta_rest<=%f",theta_min,theta_max);
	TCut deltaE= Form("(Beam-efinal_lab)>=%f&&(Beam-efinal_lab)<=%f",deltaEmin,deltaEmax);
	TCut Sigma_pos = "Sigma_rest>=0";
	TCut Pro_mo = "Prot_momentum>=0.0";
	TCut Lam ="1";//= "lambda <=3.14159";
	
	TCut Total = QQ&&Theta&&deltaE&&Sigma_pos&&W&&Pro_mo&&Lam;


//Histogram prep;
	TH1F *H1[5]; //Run 1 Full histos
	TH1F *H2[5]; //Run 2 Full histos
///
/*
	TCanvas *CC = new TCanvas("CC,");
	TH1F *R1 = new TH1F("R1","Xbins run1",30,0,1);
	TH1F *R2 = new TH1F("R2","Xbins run2",30,0,1);
	TCut Full= "1";
	T2->Draw("Xb_lab>>R2","Sigma_rest"*Total);
	double count1=R2->GetEntries();
	R2->Sumw2();
	R2->Scale(1/count1);
cout << "Denom has " << count1 <<" entries."<<endl;
	tree->Draw("Xb_lab>>R1","Sigma_rest"*Total,"same");
	double count22=R1->GetEntries();
	R1->Sumw2();
	R1->Scale(1/count22);

cout << "Numerator has " << count22 <<" entries."<<endl;


	TCanvas *CC1 = new TCanvas("CC1");
	TH1F *Ratio1 = new TH1F("Ratio1",Form("Xbins %d/%d",run1,run2),30,0,1);
	Ratio1->GetXaxis()->SetTitle("Xb");
	Ratio1->GetYaxis()->SetTitle(Form("%s/%s",two,one));
	Ratio1->Divide(R2,R1,1,1);
	Ratio1->Draw();
////
*/

	TCanvas *C[5];

////
	int N = T2->GetEntries();
	int N2= tree->GetEntries();
	int events=0;	
	if(N>=N2){events=N;}else{events=N2;}	

	
//Begining of xbin calc.
	for(int i=0;i<ncuts;i++){//Qsq loop- 5 values of qsq!!!!
//Histo creat...			
H1[i] = new TH1F(Form("Run1_%4.2f",Qsq[i]),Form("Xbins of Run%d with Q[%4.2f]",run1,Qsq[i]),1000,0 ,1.8);
H2[i] = new TH1F(Form("Run2_%4.2f",Qsq[i]),Form("Xbins of Run%d with Q[%4.2f]",run2,Qsq[i]),1000,0 , 1.8);	
///////////////////////////////////////////////
		}

///counters for skipping
	int skip=0,skip2=0;
	
	if(quite==0){cout <<"Starting the first event! Good luck! "<<endl;}
////Event by event loop for both trees
	for(int j=0;j<events;j++){
		tree->GetEntry(j);
		T2->GetEntry(j);
////////run1 qsq cuts	
		//if(theta>=35*3.1415/180.0&&theta<=50*3.1415/180.0){
		if(j<N){			
			if(q1> Qsq[4]+Qsq[4]*DQsq[4]){skip++;}
				else if(q1> Qsq[4]-Qsq[4]*DQsq[4]){H1[4]->Fill(xb1,sig1);}

				else if(q1> Qsq[3]-Qsq[3]*DQsq[3]&&q1< Qsq[3]+Qsq[3]*DQsq[3]){
					H1[3]->Fill(xb1,sig1);}
				else if(q1> Qsq[2]-Qsq[2]*DQsq[2]&&q1< Qsq[2]+Qsq[2]*DQsq[2]){
					H1[2]->Fill(xb1,sig1);}		
				else if(q1> Qsq[1]-Qsq[1]*DQsq[1]&&q1< Qsq[1]+Qsq[1]*DQsq[1]){
					H1[1]->Fill(xb1,sig1);}
				else if(q1> Qsq[0]-Qsq[0]*DQsq[0]&&q1< Qsq[0]+Qsq[0]*DQsq[0]){
					H1[0]->Fill(xb1,sig1);}
				else{skip++;}
			}//}
///////////////////
////////run2 qsq cuts	
	//	if(theta2>=35*3.1415/180.0&&theta2<=50*3.1415/180.0){
		if(j<N2){
			if(q2> Qsq[4]+Qsq[4]*DQsq[4]){skip2++;}
				else if(q2> Qsq[4]-Qsq[4]*DQsq[4]){H2[4]->Fill(xb2,sig2);}

				else if(q2> Qsq[3]-Qsq[3]*DQsq[3]&&q2< Qsq[3]+Qsq[3]*DQsq[3]){
					H2[3]->Fill(xb2,sig2);}
				else if(q2> Qsq[2]-Qsq[2]*DQsq[2]&&q2< Qsq[2]+Qsq[2]*DQsq[2]){
					H2[2]->Fill(xb2,sig2);}		
				else if(q2> Qsq[1]-Qsq[1]*DQsq[1]&&q2< Qsq[1]+Qsq[1]*DQsq[1]){
					H2[1]->Fill(xb2,sig2);}
				else if(q2> Qsq[0]-Qsq[0]*DQsq[0]&&q2< Qsq[0]+Qsq[0]*DQsq[0]){
					H2[0]->Fill(xb2,sig2);}
				else{skip2++;}
			}//}
///Print statment for santity of running long runs;
	if(quick==1){if(quite==0){if(floor(j/5000000.0) == ceil(j/5000000.0)){cout << j <<" "<<"\n";}}}
	if(quick!=1){if(floor(j/50000000.0) == ceil(j/50000000.0)){cout << j <<" "<<"\n";}}
		
//Debbugging forcing the tree to stop at some # 
		if(quick==1){if(j>20000000){break;}}
	}

	C[0] = new TCanvas("C0","histos");
	C[0]->Divide(1,5);
	//TPad *pad[5];
	if(quite==0){cout <<"Starting to Draw! So close don't mess up yet!"<<endl;}
///Draw the 5 Qsq values of x bins
	double padtop[5]={1.0,0.8,0.6,0.4,0.2};
	double padbot[5]={0.8,0.6,0.4,0.2,0.05};

	for(int i=0;i<5;i++){
	C[0]->cd(i+1);
	//pad[i] = new TPad(Form("pad%d",i),Form("pad%d",i),0,padbot[i],1,padtop[i]);
	//pad[i]->SetTopMargin(0);
	//pad[i]->SetBottomMargin(0);
	//pad[i]->Draw();
	//pad[i]->cd();	
		H1[i]->SetStats(0);
		H2[i]->SetStats(0);
		H1[i]->Draw("same");
		H2[i]->Draw("same");
		H2[i]->SetLineColor(2);
		H1[i]->SetLineWidth(2);
		count = H1[i]->GetEntries();
		H1[i]->Scale(1/count);
		count2 = H2[i]->GetEntries();
		H2[i]->Scale(1/count2);
	
	}
///////////////////////////

////////////Find the values for selected bins of x.
	//Aray varibles for TGraghs.
	double Xvalue1[5][num_bins];
	double Xvalue2[5][num_bins];
	double DX1[5][num_bins]; double E1;
	double DX2[5][num_bins]; double E2;

	
	double Xpos[num_bins]; for(int xi=0;xi<num_bins;xi++){Xpos[xi]=xi*0.1+0.1;}
	int bin1, bin2;

	Double_t x[num_bins];
	Double_t y[num_bins];
	Double_t dx[num_bins];
	Double_t dy[num_bins];
	Double_t y2[num_bins];
	Double_t dy2[num_bins];
	Double_t Rdy[num_bins];
	Double_t ratio[num_bins];
	
	TGraphErrors *G1[5];
	TGraphErrors *G2[5];
	TGraphErrors *R[5];

//File for ratio results
	char xfile[100];
	int xf = sprintf(xfile, "%s/Xbins_%d_%d.txt",data_dir,run2,run1);
	ofstream xbins;
	xbins.open(xfile);
	if (!xbins.good()){cout << "error.. error.. Will Robinson.... "<<endl;return;}
/////////////////////////

	xbins<<setw(7)<<"Qsq ";
	for(int m=0;m<num_bins;m++){xbins<<setprecision(3)<<setw(7)<<m*0.1+0.1<<" ";}
	xbins<<endl;
		 
	for(int i=0;i<5;i++){	
		xbins<<setw(7)<<setprecision(3)<<Qsq[i]<<" ";
		for(int k=0;k<num_bins;k++){
			bin1 = H1[i]->FindBin(Xpos[k]-Xpos[k]*0.15,0,0);
			bin2 = H1[i]->FindBin(Xpos[k]+Xpos[k]*0.15,0,0);
			Xvalue1[i][k]=H1[i]->IntegralAndError(bin1,bin2, E1,"width");
			Xvalue2[i][k]=H2[i]->IntegralAndError(bin1,bin2, E2,"width");
			DX1[i][k]=E1;DX2[i][k]=E2;
		//Graph run1 
			y[k] = Xvalue1[i][k]; 		dy[k]=DX1[i][k];	
			x[k] = Xpos[k]; 		dx[k]= Xpos[k]*0.15;
				
		//Graph run2
			y2[k] = Xvalue2[i][k]; 		dy2[k]=DX2[i][k];	
			x[k] = Xpos[k]; 		dx[k]= Xpos[k]*0.15;

		//ratio run2/run1
			ratio[k]=y[k]/y2[k];  Rdy[k]= sqrt(dy[k]*dy[k]+dy2[k]*dy2[k]);
			if(y[k]==0||y2[k]==0){ratio[k]=0;}

	if(quite==0){cout<<setprecision(3)<<setw(7) <<ratio[k]<<" ";}
	xbins<<setprecision(3)<<setw(7) <<ratio[k]<<" ";
	
		}//end of k or x points
	xbins<<endl;
	cout <<endl;
		G1[i]=new TGraphErrors(num_bins, x,y,dx,dy);
		G2[i]=new TGraphErrors(num_bins, x,y2,dx,dy2);
		R[i] =new TGraphErrors(num_bins, x,ratio,dx,Rdy);
	
	}//end of i or Qsq
	xbins.close();

// Draw qsq cuts for run1	
	C[1] = new TCanvas("C1","run1");
	C[1]->Divide(1,5);
	for(int j=0;j<5;j++){
		C[1]->cd(j+1);
		G1[j]->SetName(Form("Run%d_qsq[%4.2f]",run2,Qsq[j]));
		G1[j]->Write();
		G1[j]->Draw("ap");}
/////////////////////////////////////
//Draw qsq cuts for run2
	C[2] = new TCanvas("C2","run2");
	C[2]->Divide(1,5);
	for(int j=0;j<5;j++){
		C[2]->cd(j+1);
		G2[j]->SetName(Form("Run%d_qsq[%4.2f]",run1,Qsq[j]));
		G2[j]->Write();
		G2[j]->Draw("ap");}
//////////////////////////////////////
//Draw qsq cuts for ratio 2/1
	C[3] = new TCanvas("C3","Ratio");
	C[3]->Divide(1,5);
	for(int j=0;j<5;j++){
		C[3]->cd(j+1);
		R[j]->SetName(Form("Ratio(%d/%d)_qsq[%4.2f]",run2,run1,Qsq[j]));
		R[j]->GetYaxis()->SetRangeUser(0.90,1.05);
		R[j]->Write();		
		R[j]->Draw("ap");}
//////////////////////////////////////


	//

     f->Write();
	 time_t finish = time(0) ;

	if(quite==0){if(finish -start >= 60){
 		cout << "This program took " << floor((finish-start)/60.0) << " minutes and ";
		double secs = (finish-start)/60.0 - floor((finish-start)/60.0);
		cout << secs*60 << " seconds to run" <<"\n"; }
	else{ cout << "This program took " << finish-start << " seconds to run"<<"\n";}}

	/*delete C[0];
	delete C[1];	
	delete C[2];
	delete C[3];
	*/

	TCanvas *Rat[5];
	for(int i=0;i<5;i++){
		Rat[i]=new TCanvas(Form("Ra%d",i),Form("Ratio(%d/%d) for qsq[%4.2f]",run2,run1,Qsq[i]));
		R[i]->GetXaxis()->SetRangeUser(0,1.2);
		R[i]->SetTitle(Form("%s to %s ratio - Qsq[%4.2f]",two,one,Qsq[i]));
		R[i]->GetXaxis()->SetTitle("Xb");
		R[i]->GetYaxis()->SetTitle(" A/D ");
		R[i]->Draw("ap");
		//delete R[i];
		}
	
	}	





/*	

//Histogram prep;
	TH1F *H1[5][num_bins]; //Run 1 Full histos
	TH1F *H2[5][num_bins]; //Run 2 Full histos
//Begining of xbin calc.
	for(int i=0;i<5;i++){//Qsq loop- 5 values of qsq!!!!
	//New Tcanvas for each Q
		C[i] = new TCanvas(Form("C%d",i));

		for(int j=0;j<num_bins;j++){//Xbin loop, looking to find num_bins points, between 0.1 and 0.13!!!

//Histo creat...			
H1[i][j] = new TH1F(Form("Run1_%f_%f",Qsq[i],Xpos[j]),Form("Xbins of Run%d with Q[%f] and X[%f]",run1,Qsq[i],Xpos[j]),10, Xpos[j]-0.04, Xpos[j]-0.05);			
H2[i][j] = new TH1F(Form("Run2_%f_%f",Qsq[i],Xpos[j]),Form("Xbins of Run%d with Q[%f] and X[%f]",run2,Qsq[i],Xpos[j]),10, Xpos[j]-0.04, Xpos[j]-0.05);	
///////////////////////////////////////////////


		T->Draw(Form("Xb_lab>>Run1_%f_%f",Qsq[i],Xpos[j]),Q[i],"same");
		tree->Draw(Form("Xb_lab>>Run2_%f_%f",Qsq[i],Xpos[j]),Q[i],"same");
	
		H2[i][j]->SetLineColor(4);

					
				


		} //Closing xbins
	}//Closing qsq


*/

	



