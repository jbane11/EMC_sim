#include <sstream> 
#include <cmath>
#include <fstream>
#include <ctime>
#include <sys/stat.h>
#include <unistd.h>
#include <iomanip>
#include <TSystem.h> 
#include <TString.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TChain.h>
#include <TTree.h>
#include <fstream>
#include <TF1.h>
#include <TNtuple.h>
#include <iostream>     //for using cout
#include <stdlib.h>     //for using the function sleep
using namespace std;


double Sigma_L_Up(double e_final, double theta, double e_in,double F_two,double xb,string a);

double look_up(string a, double x);
double Ftwo_mod(double Ftwo, double x, int A);
double	cross_section(double e_final, double theta, double e_in,double F_two,double xb,int nth);

void Sim_root_6() {

		time_t start = time(0) ;
		cout << "\n" << "\n";

  char* data_dir;
  data_dir = getenv ("OUTPUT_DIR");
  char* root_dir;
  root_dir = getenv ("OUT_DIR");

int run =0;
cout << "Please input the run number of the disturbution used." <<endl;
cin >> run;
//run = 260;
	ifstream in;
	in.open(Form("%s/ISroot_run_%d.txt",data_dir,run));
	if(!in){cout << "Run " << run << " does not exsit." <<endl;return;}


int yes = 1;
cout << "Would you like to modify F_two? 1 for yes 0 for no"<<endl;
cin >> yes;

string nucli ="He4";
if(yes ==1){
cout << "Which nucleus do want to use for the modification? \n"; cout <<endl;
   	getline (cin, nucli);
    getline (cin, nucli);
	nucli += ".csv";
	cout << "file name = " <<nucli <<endl;

	ifstream nucleus;
	nucleus.open (nucli.c_str());

	if (!nucleus){cout << "No txt file exist for this nucleous. Please try again.\n";
		getline (cin, nucli);nucli += ".csv";
		cout << "file name = " <<nucli <<endl;					
		nucleus.open (nucli.c_str());
			if(!nucleus){cout <<" Sorry, try again later" <<endl; return;}}}


double  DeltaELab,DeltaERest,W;
//int bad=0;
//int k=0;
//int j=0;
int i=0;
//double mp = .938;



	char output_file[100];
	int n = sprintf(output_file, "%s/IS_%d.root",root_dir,run);
	ifstream ifile(output_file);
	if (ifile) {
		cout << "Run " << run << " already exists. Would you like to replace it? 1 for yes 0 for no"<<endl;
		int replace; 
		cin >> replace ;
		cout <<endl;
	if(replace == 1){cout << "Replacing run "<< run<<endl<<endl;}
		else{return;}	
	 }
	
  	  time_t rawtime;
	  struct tm * timeinfo;
	  time (&rawtime);
	  timeinfo = localtime (&rawtime);
	  printf("Started run %d on %s",run, asctime(timeinfo));


		std::fstream hist; 
		hist.open("./Simulation_Root_history.txt", std::fstream::app);
		hist << "Started run "<<run <<"  on "<<asctime(timeinfo)<<"!"<<endl;
		hist.close();



	Float_t lambda,qsquared,theta_lab,Xb_lab,efinal_lab, Prot_momentum, Diffcross_lab;
	Float_t Xb_rest, theta_rest,eEri,ISE,Eloss,Beam,F_two,Sigma_mod[7],F_two_mod[7];
	Float_t Xa_rest, Xa_lab, Ftwo_mod_new,Sigma_mod_new;
	double R1,R2,R3,R4;
	int event=0;
		
//	double ones = ceil((run/10.0 - run/10)*10);
//	if(ones == 0){ones ==1};
//	int A = ones;

	Int_t nlines = 0;

	remove(Form("%s/IS_%d.root",root_dir,run));
	TFile *f = new TFile(Form("%s/IS_%d.root",root_dir,run),"RECREATE");

	TTree *tree = new TTree("tree",Form("Tree for run %d", run));
		tree->Branch("Run",&run,"Run/I");
		tree->Branch("Prot_momentum",&Prot_momentum,"Prot_momentum/F");
	 	tree->Branch("lambda",&lambda,"lambda/F");
	 	tree->Branch("theta_rest",&theta_rest,"theta_rest/F");
	 	tree->Branch("E_rest",&eEri,"E_rest/F");
	 	tree->Branch("E_prime_rest",&ISE,"E_prime_rest/F");
	 	tree->Branch("Random_E_loss",&Eloss,"Random_E_loss/F");
	 	tree->Branch("qsquared",&qsquared,"qsquared/F");
	 	tree->Branch("Xb_rest",&Xb_rest,"Xb_rest/F");
	 	tree->Branch("Xa_rest",&Xa_rest,"Xa_rest/F");
	 	tree->Branch("Sigma_rest",&Diffcross_lab,"Sigma_rest/F");
	 	tree->Branch("F_two",&F_two,"F_two/F");
	 	tree->Branch("theta_lab",&theta_lab,"theta_lab/F");
	 	tree->Branch("Beam",&Beam,"Beam/F");
	 	tree->Branch("efinal_lab",&efinal_lab,"efinal_lab/F");
	 	tree->Branch("Xb_lab",&Xb_lab,"Xb_lab/F");
	 	tree->Branch("invar_mass",&W,"invar_mass/D");
		tree->Branch("Sigma_mod_minus",&Sigma_mod[5],"Sigma_mod_minus/F");
		tree->Branch("Sigma_mod_plus",&Sigma_mod[6],"Sigma_mod_plus/F");
		tree->Branch("Ftwo_mod_minus",&F_two_mod[5],"Ftwo_mod_minus/F");
		tree->Branch("Ftwo_mod_plus",&F_two_mod[6],"Ftwo_mod_plus/F");
		tree->Branch("event",&event,"event/I");
		tree->Branch("Ftwo_mod_new",&Ftwo_mod_new,"Ftwo_mod_new/F");
		tree->Branch("Sigma_mod_new",&Sigma_mod_new,"Sigma_mod_new/F");

//	TNtuple *Trun = new TNtuple("T2","data from ascii file","lambda:qsquared:theta_lab:theta_rest:xb_lab:xb_rest:efinal_lab: Proton_mom:Sigma_rest:F_two:e_rest_i:e_rest_f:random_loss:Beam:Invar_mass");

TNtuple *T = new TNtuple("T","data from ascii file","lambda:qsquared:theta_lab:theta_rest:xb_lab:xb_rest:efinal_lab:Proton_mom:Sigma_rest:F_two:e_rest_i:e_rest_f:random_loss:Beam:invar_mass");

TNtuple *Randoms = new TNtuple("Randoms","Random_nums","Prot_mom:random_lam:scatt_theta:random_loss");

TH1F *Xb_lab_W_low = new TH1F("Xb_lab_Weighted_low",Form("Counts of Lab frame Xb weighted with cross section for dist %d, xb less than one",run),70,0,1);
Xb_lab_W_low->GetXaxis()->SetTitle("Xb - scaling variable");
Xb_lab_W_low->GetYaxis()->SetTitle("Counts");


TH1F *Xb_rest_W_low = new TH1F("Xb_rest_Weighted_low",Form("Counts of rest frame Xb weighted with cross section for dist %d, xb less than one",run),70,0,1);
Xb_rest_W_low->GetXaxis()->SetTitle("Xb - scaling variable");
Xb_rest_W_low->GetYaxis()->SetTitle("Counts");

TH1F *Xb_lab_W_high = new TH1F("Xb_lab_Weighted_high",Form("Counts of Lab frame Xb weighted with cross section for dist %d, xb greater than one",run),70,1.0,1.5);
Xb_lab_W_high->GetXaxis()->SetTitle("Xb - scaling variable");
Xb_lab_W_high->GetYaxis()->SetTitle("Counts");


while (1) { //Xb-> lab frame, x-> rest fram

 in >> Prot_momentum>> lambda>>theta_rest>>eEri>>ISE>>Eloss>>qsquared >>Xb_rest>>Diffcross_lab >>F_two >>theta_lab>>Beam>>efinal_lab >>Xb_lab>>R1>>R2>>R3>>R4;

//DeltaELab= ,DeltaERest
//if(Diffcross_lab < 0.0){continue;}
if (!in.good()){  cout << i <<endl;break; }


 Xb_lab_W_low  ->Fill(Xb_lab ,Diffcross_lab);
 Xb_rest_W_low ->Fill(Xb_rest,Diffcross_lab); 
 Xb_lab_W_high ->Fill(Xb_lab ,Diffcross_lab);


W=.938*.938 + 2*.938*(eEri-ISE)-qsquared;
DeltaELab =Beam-efinal_lab;
DeltaERest =eEri-ISE;

Sigma_mod[4] = Diffcross_lab;
Sigma_mod[5] = Diffcross_lab;
Sigma_mod[6] = Diffcross_lab;
F_two_mod[6]  = F_two;
Ftwo_mod_new = F_two;
Sigma_mod_new = Diffcross_lab;



//Xb_rest = 0.60-0.05*i;
if(yes==1){
if(Xb_rest >= 0.325 && Xb_rest <= 0.700){
	for(int ii=5; ii<7;ii++){	Sigma_mod[ii]=cross_section(ISE, theta_rest, eEri, F_two, Xb_rest,ii);
		F_two_mod[ii] = Ftwo_mod(F_two,Xb_rest,ii); }}
	if(yes==1){
		double factor = look_up( nucli, Xb_rest);
		Ftwo_mod_new = F_two - F_two*factor;
		Sigma_mod_new = Sigma_L_Up(ISE, theta_rest, eEri,Ftwo_mod_new,Xb_rest,nucli);
}
}

//Trun->Fill(lambda,qsquared,theta_lab,theta_rest,Xb_lab,Xb_rest,efinal_lab , Prot_momentum,Diffcross_lab,F_two,eEri,ISE,Eloss,Beam,W);

	tree->Fill();
	T->Fill(lambda,qsquared,theta_lab,theta_rest,Xb_lab,Xb_rest,efinal_lab, Prot_momentum,Diffcross_lab,F_two,eEri,ISE,Eloss,Beam,W);
	Randoms->Fill(R2,R1,R3,R4);
nlines++;
i++;
event++;


//Counter
	if(i/250000 == i/250000.0){ if(i >= 1000000){cout << "Event count  " << i/1000000.0<< " "<<"mil" << "\n";}  
		else{ cout << "Event count  " << i << "\n";	}}

//Time per event calculator, avg. over 10 mil events
		if(i/10000000 == i/10000000.0){
  	 		 time_t rawtime1;
	  		struct tm * timeinfo1;
	  		time (&rawtime1);
	  		timeinfo1 = localtime (&rawtime1);
	  		printf("At %s", asctime(timeinfo1));
			time_t now = time(0) ;
				if(now -start >= 60){double runtime = (now - start)/60.0;
				cout <<"Run " << run <<" has been running for "<< runtime << " mins."<<endl;}
				else{ double runtime = now - start;
				cout <<"Run " << run <<" has been running for "<< runtime << " seconds."<<endl;}
			double eps = i/(now-start);
			cout <<"Currently running at "<< eps << " events per second."<<endl; 
			}
//if(i>=5){break;}

}printf(" found %d points\n",nlines);
	
in.close();f->Write();

time_t finish = time(0) ;

if(finish -start >= 60){
 	cout << "This program took " << floor((finish-start)/60) << " minutes and ";
	double secs = (finish-start)/60.0 - floor((finish-start)/60);
	cout << secs*60 << " seconds to run" <<endl; }
else{ cout << "This program took " << finish-start << " seconds to run"<<endl;}

	cout << "You just completed run " <<run <<" . :)"<<endl<<endl; 
// 	double	eps = i/(finish-start);
//	cout <<"Currently running at "<< eps << " events per second."<<endl; 

}
//////////////////////Start of Sub ros.!!!


double Ftwo_mod(double Ftwo, double x, int A){
	double	F_two[7];
	double mod[7]={0};
	double fourth[5]={ -49.49489, 89.39285 , -58.43405 , 17.0892   , -1.8123 };
	double fivth[6] ={-189.68558, 425.24736, -376.73102, 165.78331 , -35.69145,	3.041914};
	double sixth[7] ={-611.70385, 1646.367 , -1835.1431, 1083.17658, -355.6528,	61.91305 ,	-4.42943} ;	

	mod[4] = fourth[0]*pow(x,4) + fourth[1]*pow(x,3) + fourth[2]*x*x + fourth[3]*x + fourth[4];
	mod[5] = fivth[0]*pow(x,5) + fivth[1]*pow(x,4) + fivth[2]*pow(x,3) + fivth[3]*x*x + fivth[4]*x + fivth[5];
	mod[6] = sixth[0]*pow(x,6) + sixth[1]*pow(x,5) + sixth[2]*pow(x,4) + sixth[3]*pow(x,3) + sixth[4]*x*x + sixth[5]*x + sixth[6];
	F_two[5] = Ftwo - mod[A]*Ftwo ;
	F_two[6] = Ftwo + mod[A]*Ftwo ;
	
return F_two[A];}

double	cross_section(double e_final, double theta, double e_in,double F_two,double xb,int nth){
// This calculation use equations from the jlab pac30-petratos.pdf.

//Inelastic Cross Section
//	int nth = 4; //4, 5, or 6
	double nu = e_in - e_final;
	double alpha = (1.0/137.0);
	double mp = .938;
	F_two = Ftwo_mod(F_two,xb,nth);
	double A = alpha*alpha/(4*e_in*e_in*pow(sin(theta/2.0),4));
	double B = (F_two/nu)*pow(cos(theta/2.0),2);
	double C = (F_two/(mp*xb)) *pow(sin(theta/2),2);
	double diffcross = A*(B+C);
	// Convert 1/Gev^3 	to nb/(sr*GeV)
		   diffcross*=0.01973*0.01973*1e9;
return diffcross;}




double Sigma_L_Up(double e_final, double theta, double e_in,double F_two,double xb,string a){
// This calculation use equations from the jlab pac30-petratos.pdf.

//Inelastic Cross Section
//	int nth = 4; //4, 5, or 6
	double nu = e_in - e_final;
	double alpha = (1.0/137.0);
	double mp = .938;
	double A = alpha*alpha/(4*e_in*e_in*pow(sin(theta/2.0),4));
	double B = (F_two/nu)*pow(cos(theta/2.0),2);
	double C = (F_two/(mp*xb)) *pow(sin(theta/2),2);
	double diffcross = A*(B+C);
	// Convert 1/Gev^3 	to nb/(sr*GeV)
		   diffcross*=0.01973*0.01973*1e9;
return diffcross;}


double look_up(string a, double x){
	double factor, xb[20], mod1[20],blank[2][20]; 
	int i =0;
	int j =0;
	int k =0;
	string b;
	ifstream mod;
	mod.open (a.c_str());
	

			if(!mod.good()){factor = 0;}
			while (mod.good()){
				//getline(mod, b);

				if(j>=20){break;}
				mod >> xb[i] >> blank[0][i] >>blank[1][i] >> mod1[i];
//				cout<< x << " \t"<< xb[i]<< "\t" << mod1[i] <<endl;
				bool ijk = xb[i] == x;
				double kk = abs(xb[i]-x);
//				cout << ijk << "     "<< kk  <<endl;
				j++;
				
				if(xb[i] >= x){double x1 = xb[i];
					if(x1 == x){factor = mod1[i];}
					else{factor = (mod1[i] + mod1[i-1])/2;}
					k = i;break;}
				i++;

				}


//	cout<<endl<< "results " << xb[i]<< "\t" << factor <<endl<<endl;
return factor;}









//////////////////////////////////JUNK///////////////////////////////////////

/*
					//cout << xb[i] - x<< endl;
					if(i>=20){break;}
					if( xb[i] < x){continue;i++;}
					 else if(xb[i] == x){factor=mod1[i];break;}
					 else{factor = (mod1[i]+mod1[i-1])/2.0;
						cout<<endl<<i<< " " << j<< " "<< factor<<endl;break;}
*/
