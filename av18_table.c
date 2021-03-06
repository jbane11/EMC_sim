#include <sstream> 
#include <string>
#include <cmath>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <TSystem.h> 
#include <TGraph.h>
#include <TString.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TChain.h>
#include <TTree.h>
#include <TF1.h>
#include <TNtuple.h>
#include <iostream>     //for using cout
#include <stdlib.h>     //for using the function sleep
const Int_t steps =20000;
const int normalizer=100;
const double pi = 3.14159265359; 
using namespace std;


void av18_table(char* table,int run, int Des,int NoP=1){
		cout<<endl<<endl;	
// NoP -1 for proton, 0 for Neutron(Ns only on He3 and Be9 (may 5th 2017))!
	std::ostringstream str;
    str << table<<"_table.txt"; 
   	TString Table = str.str().c_str();
	//int Des =0;

  	cin.clear();
	char* data_dir; data_dir = getenv ("OUTPUT_DIR");
  	char* img_dir;  img_dir = getenv ("IMG_DIR");
  	char* root_dir; root_dir = getenv ("OUT_DIR");

/*	ofstream runlist;
	runlist.open(Form("%s/Runlist.txt",img_dir),std::fstream::in | std::fstream::out |std::fstream::app);

	if (!runlist.good()){ return; }
cout <<"check";
*/

	TString Titles_T[10]={"K","RHOK","RHOKP","RHOKN","DRHOKP","DRHOKN","junk","junk","junk","junk"} ;



//Opening and determine the number of columns for the av18 mom dist.
	//table = "Be9";
	int columns=0;
	TString Titles[10][2];
	int Type[10]={9};
	ifstream in_f;	in_f.open(Form("%s",Table.Data()));
	if (!in_f.good()) {cout << "That table is not here. :(" <<endl<<endl;return ;}
	if(Des==1){
		cout<< Table.Data();
		cout<<" Has the following columns,"<<endl;
		cout<<"and they are stored in the following columns:"<<endl;}
	for(int jki =0;jki<10;jki++){
		in_f >> Titles[jki][0];if(Titles[jki][0]=="****"||Titles[jki][0]=="*****"){columns=jki;break;}}
//Print the column headers;
	if(Des==1){		
		for(int jki =0;jki<columns;jki++){cout <<Titles[jki][0] <<" ";if(jki==columns){cout <<endl;}}
	}
// Determining what the columns mean
	//Titles[0][jki]="int"; 0 =momentum, 1 = Rho(k),2 = Rho(k) for the P,3 = Rho(k) for the N
	// 4 = Delta Rho(k) for the P,5 = Delta Rho(k) for the N,
	// 10 = Other
	int Selector=100;
	for(int jki =0;jki<columns;jki++){
		if(Titles[jki][0]=="K"){Type[jki]=0;}
		else if(Titles[jki][0]=="RHOK") {Type[jki]=1;Selector=1;}
		else if(Titles[jki][0]=="RHOKP"||Titles[jki][0]=="RHOKP_V18"){	
				Type[jki]=2;if(Selector==3){Selector=4;}else{Selector=2;}}
		else if(Titles[jki][0]=="RHOKN"){Type[jki]=3;if(Selector==2){Selector=4;}else{Selector=3;}}
		else if(Titles[jki][0]=="DRHOKP"||Titles[jki][0]=="DRHOKP_V18"){Type[jki]=4;}
		else if(Titles[jki][0]=="DRHOKN"){Type[jki]=5;}
		else{Type[jki]=9;}}

//Print the column =assignments;
	if(Des==1){
		cout<<endl;
		for(int jki =0;jki<columns;jki++){cout <<Type[jki] <<"    ";if(jki==columns){cout <<endl;}}
		}	
	cout <<endl<<endl;


//Get pass the header seperator.
	char star;int st_1=1;int n=0;
	while(st_1){
		in_f >> star;bool A ; A = (star=='*');
		if(A != 1){break;}n++;
		if(n>=100){break;}}
///////////////////////////////////////////////////////////////////////////////

//Begin inputing data;
	Double_t K[steps], Rho[steps],Scaled_Rho[steps];
	double Table_in[steps][10];
	Table_in[0][Type[0]]=star - '0';
	int o=0,p=0,DD=1;
//Varibles for intergration;
	double Cumaltive_int_Rho[steps], int_Rho=0;
	double int_SR_Rho=0,Cumaltive_int_SR_Rho[steps];
	int sle=Selector;
/////////////////////////////////////////////////////////////////

//Neutron or Proton;
	 //1 for proton, 0 for Neutron(Ns only on He3 and Be9 (may 5th 2017))!
	if(sle==4){if(NoP==0){sle=3;}}
	cout<<endl<<"What type of Rho are we dealing with: "<<sle<<endl;

//Print statment to title columns
	if(Des==1&&DD==1){for(int jki =0;jki<9;jki++){
		cout<<setiosflags(ios::left)<<setw(8)<<Titles_T[jki]<<"  ";}}
	cout<<endl;

//Begin the loop to insert data.
	while(1){if (!in_f.good()){cout<<endl << "Input Done" <<endl;break;}
		for(o=0;o<columns;o++){
		if(p==0&&o==0){o++;}

		in_f >>Table_in[p][Type[o]];
//Debuggin O16
	//cout << Table_in[p][Type[o]] << "  "<< p << " " << Type[o] <<endl;
		}
// Print statment to check the correct inputs
		if(p<=10){if(Des==1&&DD==1){	for(o=0;o<9;o++){
			cout<<setiosflags(ios::left) << setw(8)<<Table_in[p][o]<<"  ";}cout<<endl;}}
//Calculate the intergral of Rho(k);				
		if(sle==4){sle=2;} //Use RhoP if both are provided
		if(p>0){int_Rho += (Table_in[p][0]-Table_in[p-1][0])*((Table_in[p][sle]+Table_in[p-1][sle])/2.0);
		Cumaltive_int_Rho[p]=int_Rho;}
			else{Cumaltive_int_Rho[0]=0;}
//Finding the Sum rule intergral K^2*P(k)		 
		if(p>0){double K1=Table_in[p][0],K2=Table_in[p-1][0];
		int_SR_Rho += (K1-K2)*((K1*K1*Table_in[p][sle]+K2*K2*Table_in[p-1][sle])/2.0);}
		Cumaltive_int_SR_Rho[p]=int_SR_Rho;		
		p++;}//End of filling the input 2D Array.


///////////////////////////////////////////////////////////////////////////////

//Informative print statments:
	cout <<endl<<"There are "<< p<< " data points in " << Table.Data()  <<endl;
	cout << "The intergral of Rho is "<<int_Rho<<endl;
	cout << "The intergral of K^2*Rho is "<<int_SR_Rho<<" Normalization Constant = ";
	double Norm_cost = int_SR_Rho*4*pi/(pow(2*pi,3));
	cout << Norm_cost <<endl;
	cout << table<<" Has been normalized to :: "<<int_SR_Rho/Norm_cost<<endl<<endl;

//In order to use the cumlative sum method of randomly selecting:
//The cumaltive sum needs to max out at 1.
//Need to Normalize Rho by the intergal of Rho. Also fill one D array for K and Rho
 	Double_t K_1,Rho_1,Cumaltive_int_Rho_1;
	double K_2[steps], Cumaltive_int_Rho_2[steps],Cumaltive_Scaled_Rho[steps];
	int N;
	for(N=0;N<=p-1;N++){ 
		Cumaltive_int_Rho[N] = Cumaltive_int_Rho[N]/int_Rho/1.00000018651;
		K[N]=Table_in[N][0];Rho[N]=Table_in[N][sle];Scaled_Rho[N]=Rho[N]*K[N]*K[N];
		Cumaltive_Scaled_Rho[N]=Cumaltive_int_SR_Rho[N]/int_SR_Rho;
	}

	cout << "Now, the Cumalitive sum as been normalized to max out at : ";	
	cout << Cumaltive_Scaled_Rho[N-1]<<endl;
////////////////////////////////////////////////////////////////////////////
//Graphing out a few graphs.
	TCanvas *C1 = new TCanvas("C1","C1",0,0,1000,900);
	C1->Divide(1,2);
	C1->cd(1);
	TGraph *Rho_K = new TGraph(p-1,K,Rho);
	Rho_K->Draw();
	gPad->SetLogy();
	gPad->SetGridx();
	gPad->SetGridy();
	C1->cd(2);
	TGraph *Cumal_inv= new TGraph(p-1,Cumaltive_Scaled_Rho,K);
	Cumal_inv->Draw();
	gPad->SetGridx();
	TGraph *KK_Rho_K= new TGraph(p-1,K,Scaled_Rho);
	//Determining the scale for differnet k cuts
	
	int min = KK_Rho_K->GetXaxis()->FindBin(0.25*5.068);
	int max = KK_Rho_K->GetXaxis()->FindBin(10);
	//double int_12 = KK_Rho_K->Integral(min,max)/19.7392;

	 min = KK_Rho_K->GetXaxis()->FindBin(0.3*5.068);
	//double int_22 = KK_Rho_K->Integral(min,max)/19.7392;

	  min = KK_Rho_K->GetXaxis()->FindBin(0.35*5.068);
	//double int_32 = KK_Rho_K->Integral(min,max)/19.7392;
	
	double min1=0.25*5.068,min2=0.30*5.068,min3=0.35*5.068,min_check=2.0;
        //double  min_check = KK_Rho_K->GetXaxis()->FindBin(2.0);
	 double int_check=0,int_1=0,int_2=0,int_3=0; 
	 N=0;
	double K1,K2;
	for(N=0;N<=p-1;N++){ 
	  if(N>0){ K1=Table_in[N][0]; K2=Table_in[N-1][0];}
		else{K1=Table_in[N][0];K2=0.0;}
	  //	cout <<K[N]<<endl;
	        if(K[N]>=min1){
		  int_1 += (K1-K2)*((K1*K1*Table_in[N][sle]+(K2*K2*Table_in[N-1][sle]))/2.0);}
	        if(K[N]>=min2){
		  int_2 += (K1-K2)*((K1*K1*Table_in[N][sle]+K2*K2*Table_in[N-1][sle])/2.0);}
	        if(K[N]>=min3){
		  int_3 += (K1-K2)*((K1*K1*Table_in[N][sle]+K2*K2*Table_in[N-1][sle])/2.0);}
	        if(K[N]>=min_check){
		  int_check += (K1-K2)*((K1*K1*Table_in[N][sle]+K2*K2*Table_in[N-1][sle])/2.0);}
	}


	cout <<endl;
	cout<< "The kk*Rho(k) integral for > than 2.0 (1/fm) = "<< int_check/19.739<<endl;
	cout<< "The kk*Rho(k) integral for > than 0.25 GeV = "<< int_1/19.739<<endl;
	cout<< "The kk*Rho(k) integral for > than 0.30 GeV = "<< int_2/19.739<<endl;
	cout<< "The kk*Rho(k) integral for > than 0.35 GeV = "<< int_3/19.739<<endl;


//////////////////////////////////////////////////////////////////////////////

	TFile *f = new TFile(Form("%s/Av18_%d.root",root_dir,run),"recreate");
	TTree *tree = new TTree("tree",Form("Tree for run %d", run));
		tree->Branch("K",&K_1,"K_1/D");		
		tree->Branch("Cumaltive_Sum",&Cumaltive_int_Rho_1,"Cumaltive_int_Rho_1/D");
		tree->Branch("Rho",&Rho_1,"Rho_1/D");

//File to use in Sim_IS_6. In contains the cumaltive sum and k in units of (1/fm)!!!!
//And Fill the Tree
	ofstream lists;
	lists.open(Form("%s/mo_chart_%d.txt",data_dir,run));
	for(int i= 0;i<p-1;i++){Cumaltive_int_Rho_1=Cumaltive_int_Rho[i];
		Cumaltive_int_Rho_2[i]=Cumaltive_int_Rho[i]; K_2[i]=K[i];
		K_1=K[i];Rho_1=Rho[i];
		lists << setprecision(20)<<Cumaltive_int_Rho_2[i] <<" "<<setprecision(8)<< K_2[i]<<endl;
		tree->Fill();}
	lists.close();
///////////////////////////////////////////
	Rho_K->SetName("Rho_K");
	Cumal_inv->SetName("Cumal_inv");
	KK_Rho_K->SetName("Scaled_Rho");

	Rho_K->Write();
	Cumal_inv->Write();
	KK_Rho_K->Write();
	f->Write();
	


}
