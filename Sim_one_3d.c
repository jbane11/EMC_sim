#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <time.h>
#include <ctime>
#include <sys/stat.h>
#include <unistd.h>
#include <TSystem.h> 
#include <TString.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TFile.h>
#include <TH1F.h>
#include <TChain.h>
#include <TTree.h>
#include <TF1.h>
#include <TNtuple.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <math.h>

using namespace std;
const 	int number_of_electrons =500000000; 
//#include "nr3.h"
//#include "ran.h"

//Function protype.
double distribution(double A,double number[500],double F_x[500],double help[50],double rndm,int k, int run);
//uses mo_chart_%d.txt to randomize the proton momentum from a distribution.
double rotationZ(double A[4], double angle,int i);
// Rotates about the Z axis
double rotationX(double A[], double angle, int i);
// Rotates about the x axis
double rotationY(double A[], double angle, int i);
// Rotates about the x axis
void printvectors(double A[4], double B[4],int f, char* string);
//Prints the 4 vecotr of e. Debugging!
double Ftwo(double Beam, double Eprime, double theta);
double Boostx(double A[4], double gamma, double beta, int i);

double Ftwo_mod(double Ftwo, double x);
double	cross_section(double e_final, double theta, double e_in,double F_two,double xb);


////	run = name -> run1 = mo_dist
void Sim_one_3d(int run =0,int run1 =0, int pri=0){
	int no_O=1;	

	time_t start = time(0) ;
	cout << "\n" << "\n";

// Enviromental varibles for 
  char* data_dir;  data_dir = getenv ("OUTPUT_DIR");
  char* root_dir;  root_dir = getenv ("OUT_DIR");
	if(run1==0){
cout << "Please input the number of the momentum disturbution you would like to use." <<"\n";;
		cin >> run1;}
	if(run==0){
cout << "Please input the run number, you would like to call this run." <<"\n";;
		cin >> run;}

	double spread = 359.9998;
	double bottom = 0.0001;

	char output_file[100];
	int n = sprintf(output_file,"%s/IS_3d_%d.root",root_dir,run);
	ifstream ifile(output_file);


	if (ifile) {
cout << "Run " << run << " 3D, already exists. Would you like to replace it? 1 for yes 0 for no"<<endl;
		int replace=1; 
		cin >> replace ;
		cout <<endl;
	if(replace == 0){return;}	
	 }
	
//Output root file!
	remove(Form("%s/IS_3d_%d.root",root_dir,run));
	TFile *f = new TFile(Form("%s/IS_3d_%d.root",root_dir,run),"RECREATE");








//Input for momentun distribution
	char chart[100];
	int Tt = sprintf(chart, "%s/mo_chart_%d.txt",data_dir,run1);
	ifstream lists;
	lists.open(chart);
	if (!lists.good()) {cout << "No mo_chart_"<<run1<<" file: error" <<endl;return;}
	cout << "\n" << "\n";

	double number[500]={0};
	double F_x[500]={0};
	for(int T=0;T<=500;T++){if (!lists.good()){ break;}
		//		lists >>setprecision(14)>> number[T] >>setprecision(8)>> F_x[T];
		lists >> number[T] >> F_x[T];
		//		cout<<setprecision(14) << number[T] <<" "<<endl;
		if (!lists.good()){ break;}}

////////////////////////////////////////////////////////////////////////////////////
	//Root random number gen.	
	TRandom *R = new TRandom3();

//cheeky << "Lambda" << "\t"<< "Q^2"<<"\t"<<  "theta"<<"\t"<< "xb-lab" <<"\t"<<"e_f"<<"\t"<<"F_two"<< "\t"<< "P_mo"<< "\t"<<"d_cross"<<"\t"<<"dcross_mod"<<"\t" <<"deltaE_rest"<<"\t"<< "QQ" << "\t"<<"x_RF"<< "\t"<<"W_RF"<< "\t"<< "W_LF" << "\t"<<"dcross_lab"<<"\t"<<"dcross_m_lab"<<" "<<"\n";

// Constants
	double Fit[16];
	double intersect[4];
	double help[50] ={0.0};	
	double random_max;

    double mp = 0.938;            // Mass of Proton in GeV/c^2
	double me = .0005;            // Mass of Electron in GeC/c^2
	double c  = 3.0 * pow(10.0,8.0);   // Speed of light in m/s
	double pi = 3.14159265359;     // Pi
   	double rad = pi/180;           // Convert from degrees to radians
   	double e[4];                   //Electron vector array
    double p[4];                   //Proton vector array
    double e_prime[4];             //ELectron vector array in prime frame
    double p_prime[4];             //Proton vector array in prime frame
    double e_prime2[4];            //ELectron vector array in prime frame
    double p_prime2[4];            //Proton vector array in prime frame
	double e_final[4],p_final[4];	//Final vectors for the electron and proton
    float elec_Beam_momentum = 5.76;//Beam momentum in GeV
    double Prot_momentum = 0.250     ;  //Proton momentum in GeV 
	float Lambda; 					//Proton angle
	float phi; 					//Proton angle
	double eEri;					//Electron energy before scattering
	double scatter_theta;          //Code - angle of scattered electron in rest frame
	double theta_rest;				//tree varible ....
	double min_angle = 0;// By defaults scattered angle covers 360 degrees
	double angle_range = 360;  
	int check =0;	
	int swag=0;
	// counters for misses. 
	int N = 0;int M = 0;int O = 0;int OO = 0;	
	int kk =0 ;int MOK =0 ;char trash[10];
	int print = pri*10;// 0; //If you want to check the 4 vectors -> 10;
	int scatt=1;int jj = 0;	int i;	int Low_P=0;

	double scattered_p_theta; 	double scattered_p_momentum;
	double diff_cross; 	double diffcross_mod; 	double deltaE_lab ;
	double deltaE_rest;
	float Qsquared; 	double xb; double Xb_rest; double Xb_lab;
	float invar_mass ;	double invar_mass_RF;
	double F_two_lab ; 	double F_two_mod_lab;
	double diff_cross_lab;	double diffcross_mod_lab;
	float theta; 	double F_two; 	double QQ ; double x ;
	double F_two_mod ;
	double eErf;	 // Electron energy after scattering in rest frame!
	double ISE ; 	float nuw; 	double Eloss;
	double efinal_lab;
	cout << "\n" << "\n";

//Possible floats to save space.
	




  	  time_t rawtime;
	  struct tm * timeinfo;
	  time (&rawtime);
	  timeinfo = localtime (&rawtime);
	  printf("Started run %d on %s",run, asctime(timeinfo));



	TTree *tree = new TTree("tree",Form("Tree for run %d", run));
		tree->Branch("Run",&run,"Run/I");
		tree->Branch("Prot_momentum",&Prot_momentum,"Prot_momentum/D");
	 	tree->Branch("lambda",&Lambda,"lambda/f");
	 	tree->Branch("Phi",&phi,"phi/f");
	 	tree->Branch("theta_rest",&theta_rest,"theta_rest/D");
	 	tree->Branch("E_rest",&eEri,"E_rest/D");
	 	tree->Branch("E_prime_rest",&ISE,"E_prime_rest/D");
	 	tree->Branch("Random_E_loss",&Eloss,"Random_E_loss/D");
	 	tree->Branch("qsquared",&Qsquared,"qsquared/F");
	 	tree->Branch("Xb_rest",&Xb_rest,"Xb_rest/D");
	 	//tree->Branch("Xa_rest",&Xa_rest,"Xa_rest/F");
	 	tree->Branch("Sigma_rest",&diff_cross,"Sigma_rest/D");
	 	tree->Branch("F_two",&F_two,"F_two/D");
	 	tree->Branch("theta_lab",&theta,"theta_lab/F");
	 	tree->Branch("Beam",&elec_Beam_momentum,"Beam/F");
	 	tree->Branch("efinal_lab",&efinal_lab,"efinal_lab/D");
	 	tree->Branch("Xb_lab",&Xb_lab,"Xb_lab/D");
		tree->Branch("nuw",&nuw,"nuw/F");	 	
		tree->Branch("invar_mass",&invar_mass,"invar_mass/F");
		//tree->Branch("Sigma_mod_minus",&Sigma_mod[5],"Sigma_mod_minus/F");
		//tree->Branch("Sigma_mod_plus",&Sigma_mod[6],"Sigma_mod_plus/F");
		//tree->Branch("Ftwo_mod_minus",&F_two_mod[5],"Ftwo_mod_minus/F");
		//tree->Branch("Ftwo_mod_plus",&F_two_mod[6],"Ftwo_mod_plus/F");
		tree->Branch("event",&i,"event/I");
		//tree->Branch("Ftwo_mod_new",&Ftwo_mod_new,"Ftwo_mod_new/F");
		//tree->Branch("Sigma_mod_new",&Sigma_mod_new,"Sigma_mod_new/F");





//////////////////////////////////////////////////////////////////////////////////////////
// Loop for the number of electrons
	for(i=1; i <= number_of_electrons; i++){

/////////////////////////////////////////////	
	// Random proton momentum from a disturbution
		long double rndm= R->Rndm();	
		Prot_momentum = distribution(i, number,F_x,help,rndm,kk,run)*0.1973;



	//if(Prot_momentum <= 0.20){Low_P++;continue;}
//	if(i/1000000==i/1000000.0){cout << " mo = "<<Prot_momentum<<"  "<<endl;}
// Random angle in degrees used for Lambda the incoming proton angle
		double random_lam =R->Rndm();
		Lambda = ((random_lam)*360)*rad;
		double random_phi =R->Rndm();
		phi =((random_phi*180)-90)*rad;		

// incomeing electron and proton 4 vectors.
    e[0] = elec_Beam_momentum;      //electron energy
    e[1] = elec_Beam_momentum;      //Electron Beam Momentum
    e[2] = 0 ;
    e[3] = 0 ;

    p[1] = Prot_momentum*cos(Lambda+pi)*sin(pi/2.0-phi);  //Proton x momentum
    p[2] = Prot_momentum*sin(Lambda+pi)*sin(pi/2.0-phi);  //Proton y momentum
    p[3] = Prot_momentum*cos(pi/2.0-phi);//Proton Z momentum
    p[0] = sqrt( mp*mp+Prot_momentum*Prot_momentum); //Proton energy

if(print == 10){cout << setw(25)<< "Electron 4 vector" <<"\t" << "\t     "<<  "Proton 4 vector " <<"\n";}
	printvectors(e, p, print," Initial vectors ");
////////////////////////////////



// Rotating by Lambda,  clockwise!! This sets Py to 0.
for(int k = 0; k < 4; k++){ 
    e_prime[k] = rotationZ(e , -Lambda, k);
    p_prime[k] = rotationZ(p , -Lambda, k); }
	printvectors(e_prime, p_prime, print," rotate around z ");


  	rotationY( p_prime, phi,1);
 	rotationY( e_prime, phi,1);
	printvectors(e_prime, p_prime, print," rotate around y ");


// Lorentz factors!
	double gamma = p[0]/mp;
	double beta  = p_prime[1]/p[0];

//Lorentz boost in the x direction.
for(int k = 0; k < 4; k++){ 
    e_prime2[k] = Boostx(e_prime , gamma, beta, k);
    p_prime2[k] = Boostx(p_prime , gamma, beta, k); }
	printvectors(e_prime2, p_prime2, print," boost ");

// Calculate new angle Delta, angle between e and the x axis. 
    double Delta = atan(e_prime2[2]/e_prime2[1]);

// Rotate the electron to the x axis. 
for(int k = 0; k < 4; k++){ 
    e_prime[k] = rotationZ(e_prime2 , -Delta, k);
    p_prime[k] = rotationZ(p_prime2 , -Delta, k); }
	printvectors(e_prime, p_prime, print," rotate around z ");

//Rotate by phi-
	double phi2 = -atan(e_prime[3]/e_prime[1]);
 	rotationY( p_prime, phi2,1);
 	rotationY( e_prime, phi2,1);
	printvectors(e_prime, p_prime, print," rotate around y ");


	eEri = e_prime[0]; //Electron energy in rest frame before

// Due to limitations of atan, forceing the e[1](x momentum) to be postive. This makes the beam always head to the right, in my internal drawling. This will also let me know if we need to flip back.
	int flip =0;
	if(e_prime[1] < 0){flip = 10; e_prime[1] = abs(e_prime[1]);}  

	printvectors(e_prime, p_prime, print," flip if needed ");

//WE ARE NOW IN THE REST FRAME OF THE TARGET IN THE ORENTATION OF A FIXED TARGET SCATTERING
//Randomly select scattered angle
	double R2=R->Rndm();
	scatter_theta = R2*spread*rad + bottom*rad;
	theta_rest=scatter_theta;
// Scatter!! Scatter the electron to a random angle, using conservation of momentum to calulate all varibles.
    e_prime2[0] = -1*(e_prime[0]*mp/(cos(scatter_theta)*e_prime[0] - e_prime[0]-mp));

// Inelastic scattering! Some energy will be lost!
	//double E_loss_fun =1- exp(-1.2*RANDOM(i*run*2));
	double R3=R->Rndm();
	double E_loss_fun =R3;
	 Eloss = e_prime2[0]*E_loss_fun;

//Recalculate energy for the e`
	 ISE = e_prime2[0]-Eloss;
	 nuw = e_prime[0] - ISE;

	e_prime2[0] =  ISE ;
    e_prime2[1] =  e_prime2[0]*cos(scatter_theta);
    e_prime2[2] =  e_prime2[0]*sin(scatter_theta);
    e_prime2[3] =  0.0 ;   

//Calculate x,Q,F_2. Call a function, sending
	 F_two = Ftwo(e_prime[0], e_prime2[0],scatter_theta);
	 QQ = 4*e_prime[0]*e_prime2[0]*sin(scatter_theta/2)*sin(scatter_theta/2);
	 x = QQ/(2*0.938*((e_prime[0]-e_prime2[0])));
	 F_two_mod = F_two;
	if(x >= 0.3 && x <= 0.7){F_two_mod = Ftwo_mod(F_two,x);}
	Xb_rest = x;



 	diff_cross=cross_section( e_prime2[0], scatter_theta,e_prime[0] ,F_two, x);
 	diffcross_mod=cross_section( e_prime2[0], scatter_theta, e_prime[0],F_two_mod, x);

//Using results for e, calculate p's info. 

    scattered_p_momentum = sqrt(e_prime[0]*e_prime[0] + e_prime2[0]*e_prime2[0] + 2*mp*e_prime[0] - 2*e_prime[0]*e_prime2[0] - 2*mp*e_prime2[0]);
    p_prime2[0] = sqrt(scattered_p_momentum*scattered_p_momentum +mp*mp );
   scattered_p_theta = asin(-e_prime2[0]*sin(scatter_theta)/scattered_p_momentum);
    p_prime2[1] = scattered_p_momentum*cos(scattered_p_theta);
    p_prime2[2] = scattered_p_momentum*sin(scattered_p_theta);
    p_prime2[3] = 0.0;  

printvectors(e_prime2, p_prime2, print," scatter ");

	 eErf = e_prime2[0]; // Electron energy after scattering in rest frame!

// We begin transforming back to the lab frame!
// If we needed to flip before, we unflip.
if(flip == 10){for(int k = 0; k < 4; k++){ 
    e_prime[k] = rotationZ(e_prime2 , 180*rad, k);
    p_prime[k] = rotationZ(p_prime2 , 180*rad, k); }   
    for(int k = 0; k < 4; k++) {
        e_prime2[k]=e_prime[k];
        p_prime2[k]=p_prime[k];}}

 	rotationY( p_prime2, phi2,1);
 	rotationY( e_prime2, phi2,1);
	printvectors(e_prime2, p_prime2, print," rotate around y ");

// Rotated back by delta.
for(int k = 0; k < 4; k++){ 
    e_prime[k] = rotationZ(e_prime2 , Delta, k);
    p_prime[k] = rotationZ(p_prime2 , Delta, k); }

printvectors(e_prime, p_prime, print," rotate around Z ");

//Boost back!   
for(int k = 0; k < 4; k++){ 
    e_prime2[k] = Boostx(e_prime , gamma, -beta, k);
    p_prime2[k] = Boostx(p_prime , gamma, -beta, k); }

printvectors(e_prime2, p_prime2, print," boost back");

 	rotationY( p_prime2, phi,1);
 	rotationY( e_prime2, phi,1);
	printvectors(e_prime2, p_prime2, print," rotate around y ");

//Rotate back by Lambda
for(int k = 0; k < 4; k++){ 
    e_final[k] = rotationZ(e_prime2 , Lambda, k);
    p_final[k] = rotationZ(p_prime2 , Lambda, k); }

	 theta = atan(e_final[2]/e_final[1]); // Final lab frame scattered angle
// To determine the right quadrent. 
if(e_final[1] < 0 && e_final[2] < 0){theta = pi + theta;}	
if(e_final[1] < 0 && e_final[2] > 0){theta =  (pi+theta);}
if(e_final[1] > 0 && e_final[2] < 0){theta = 2*pi + theta;}
   
printvectors(e_final, p_final,  print, " Final vecotors");
	 efinal_lab= e_final[0];
	 deltaE_lab  = elec_Beam_momentum -e_final[0];
	 deltaE_rest = eEri-eErf;
	 Qsquared = 4 * elec_Beam_momentum * e_final[0] * sin(theta/2) *sin(theta/2);
	 xb	= Qsquared/(2*mp*((elec_Beam_momentum-e_final[0])));
	 invar_mass = sqrt(mp*mp - Qsquared + 2*mp*deltaE_lab);
	 invar_mass_RF = sqrt(mp*mp - Qsquared + 2*mp*deltaE_rest);
	 Xb_lab = xb;
//##################################################

	 F_two_lab = Ftwo(elec_Beam_momentum , e_final[0],theta);
	 F_two_mod_lab = F_two_lab;
	if(x > 0.3&& xb <= 0.7){F_two_mod_lab = Ftwo_mod(F_two,xb);}

	 diff_cross_lab=cross_section(e_final[0],theta, elec_Beam_momentum,F_two_lab, xb);
	 diffcross_mod_lab=cross_section(e_final[0],theta, elec_Beam_momentum,F_two_mod_lab, xb);

//Cross section
double	Cross_0 =  (1.0/137.0)*(1.0/137.0)*pow(cos(theta/2.0),2)/(4.0*elec_Beam_momentum*elec_Beam_momentum*pow(sin(theta/2.0),4));
double Tau = Qsquared/(4*mp*mp);
double Ge = pow((1+Qsquared)/.71,-2);
double Gm = 2.79*Ge;
double frec = 1 + (2*elec_Beam_momentum/mp)*pow(sin(theta/2),2);
double diffcross_elastic_lab = Cross_0*(e_final[0]/elec_Beam_momentum)*((Ge*Ge+Tau*Gm*Gm)/(1+Tau) + 2*Tau*Gm*Gm*pow(tan(theta/2),2));

//if(xb>=.850){diff_cross_lab=diffcross + diff_cross_lab;}
// if(xb>=1){diff_cross_lab= diffcross;}




// If the scattered angle is really close to the beam, I get NaN fo P(0).
	double L = p_final[0];



//if(Prot_momentum > 1){kk++;  continue;}
if(L != L){ N++;  continue;   } //Forces the code to not record any NaN
if(Qsquared < 0.6){M++;  continue;} //Does not record any unwanted q^2 values.
//if(invar_mass != invar_mass){OO++; continue;}// invar_mass=0.0; }//
if(F_two != F_two){O++;   continue;} 
if(F_two/F_two != F_two/F_two){O++;  continue;} 
if(diff_cross_lab != diff_cross_lab || diff_cross_lab/diff_cross_lab != diff_cross_lab/diff_cross_lab ){MOK++; diff_cross_lab=0; continue;}
if(diffcross_mod_lab != diffcross_mod_lab || diffcross_mod_lab/diffcross_mod_lab != diffcross_mod_lab/diffcross_mod_lab ){MOK++; diffcross_mod_lab =0; continue;} 



if(floor(i/(number_of_electrons/50.0)) == ceil(i/(number_of_electrons/50.0))){cout << i <<" "<<"\n";}
tree->Fill();
scatt++;

//Time per event calculator, avg. over 10 mil events
if(floor(i/(number_of_electrons/4.0)) == ceil(i/(number_of_electrons/4.0))){
  	 		 time_t rawtime1;
	  		struct tm * timeinfo1;
	  		time (&rawtime1);
	  		timeinfo1 = localtime (&rawtime1);
	  		printf("At %s", asctime(timeinfo1));
			time_t now = time(0) ;
				if(now -start >= 60){double runtime = (now - start)/60.0;
				cout <<"Run " << run <<" has been running for "<< runtime << " mins."<<"\n";}
				else{ double runtime = now - start;
				cout <<"Run " << run <<" has been running for "<< runtime << " seconds."<<"\n";}
			double eps = i/(now-start);
			cout <<"Currently running at "<< eps << " events per second."<<"\n"; 
			double time_left = (number_of_electrons-i)/eps;
		cout <<"This simulation should be done in approximently " << time_left/60.0<<" minutes!"<<"\n";
			}




//Break from loop for debugging
//if(i >= 1000000){break;}

}// End of the beam loop! ////////////////////////////////////////

int missed = M + N + kk+OO + O+Low_P+ MOK;

lists.close();f->Write();

//Output statment for general information
 time_t finish = time(0) ;
 cout <<"\n" << i-1 << " " << "electrons were sent." <<"\n";
 cout <<missed << " events have been thrown out due to q^2 Value or unphysical results." << "\n";


if(finish -start >= 60){
 	cout << "This program took " << floor((finish-start)/60.0) << " minutes and ";
	double secs = (finish-start)/60.0 - floor((finish-start)/60.0);
	cout << secs*60 << " seconds to run" <<"\n"; }
else{ cout << "This program took " << finish-start << " seconds to run"<<"\n";}

cout <<scatt<< "  Electrons were scattered."<<endl<<endl<<endl;
cout << "Q^2" <<"\t "<< "P_f" <<"\t "<<  "F2"<< "\t"<<"Low_P"<< "\t "<<"Sigma"<<"\t  "<<"W"<<endl;
cout << M <<"\t "<< N<< "\t " <<  O<< "\t "<<Low_P<< "\t "<<MOK<<"\t  "<<OO<<endl;
}//End of main program
///////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////
////////////////////////////////////
///////////////////////
///////////

////////////////////////////////////////////////////////////////////////////////////////
////Functions!!!

// Rotation about the Z axis! Rotate CCW from the postive X axis.
double rotationZ(double A[], double angle,int i) { 
    double Prime[4] = {0.0};
    double rot[4][4] = {0.0};
    rot[1][1] = rot[2][2] = cos(angle);
    rot[1][2] = - sin(angle);
    rot[2][1] = sin(angle);
    rot[0][0] = rot[3][3] =1.0;
        for (int j = 0; j < 4; j++) {
          Prime[i] = Prime[i]+ A[j]*rot[i][j]  ; }
return Prime[i];
}

double rotationX(double A[], double angle, int i){
	angle =  -angle;    
	double Prime[4] = {0.0};
    double rot[4][4] = {0.0};
    rot[3][3] = rot[2][2] = cos(angle);
    rot[2][3] = - sin(angle);
    rot[3][2] = sin(angle);
    rot[0][0] = rot[1][1] =1.0;
	
	for(int k=0; k<4;k++){		
		for (int j = 0; j < 4; j++) {
          Prime[k] = Prime[k] + A[j]*rot[k][j]  ; }
		}
	for(int k=0;k<4;k++){A[k]=Prime[k];}

}

double rotationY(double A[], double angle, int i){
	angle =  -angle;    
	double Prime[4] = {0.0};
    double rot[4][4] = {0.0};
    rot[1][1] = rot[3][3] = cos(angle);
    rot[3][1] = - sin(angle);
    rot[1][3] = sin(angle);
    rot[0][0] = rot[2][2] =1.0;
	
	for(int k=0; k<4;k++){		
		for (int j = 0; j < 4; j++) {
          Prime[k] = Prime[k] + A[j]*rot[k][j]  ; }
		}
	for(int k=0;k<4;k++){A[k]=Prime[k];}

}


// Function to print out two four vectors E|x|y|z
void printvectors(double A[], double B[],int f ,char *string){
if(f==10){
for(int i=0; i < 4; i++){ if(abs(A[i]) <= 0.00001){A[i] = 0.0;}cout <<setprecision(4)<< setw(10)<< A[i] << "|"  ;}
cout << "       ";
for(int j=0; j < 4; j++){ if(abs(B[j]) <= 0.00001){B[j] = 0.0;}cout <<setprecision(4)<< setw(10)<< B[j] << "|  "  ;}
cout << string <<"\n";
}}

// Boost in the x axis!
double Boostx(double A[], double gamma, double beta, int i){
    double Prime[4] = {0.0};
    double B[4][4] = {0.0};
    B[0][1] = B[1][0] = -beta*gamma;
    B[0][0] = B[1][1] =gamma;
    B[2][2]=B[3][3] = 1.0;
        for (int j = 0; j < 4; j++) {
          Prime[i] = Prime[i]+ A[j]*B[i][j]  ; }
return Prime[i];
}

//Calculate Ftwo-Need Beam energy,scattered electron energy and angle
double Ftwo(double Beam, double Eprime, double theta){
	double QQ = 4*Beam*Eprime*sin(theta/2)*sin(theta/2);
	double xb = QQ/(2*0.938*((Beam-Eprime)));  
	double A = 1.22*exp(3.2*xb);
	double Ftwo_thr = 0;
	double lambda[2] = {0};
	double C[12] = {0.948,-0.115,1.861,-4.733,2.348,-0.065,-0.224,1.085,0.213,-1.687,3.409,-3.255 };
	double beta =1;
	for(int i=0; i < 4; i++){lambda[0] = lambda[0] + C[i+8]*pow(xb,i); }
	if(QQ < A){lambda[1] = C[5] + C[6]*xb + C[7]*xb*xb;}
		else{lambda[1] = 0;} 
	for(int i=1; i < 6; i++){Ftwo_thr = Ftwo_thr + C[i-1]*pow((1-xb),i+2);}
	double Ftwo = beta*Ftwo_thr*(1+lambda[0]*log(QQ/A) + lambda[1]*pow(log(QQ/A),2));
return Ftwo; }

double Ftwo_mod(double Ftwo, double x){
	double A[6] ={ -189.69, 425.25, -376.76, 165.78,-35.69,3.05};
	double mod = A[0]*pow(x,5)+A[1]*pow(x,4)+A[2]*pow(x,3)+A[3]*x*x +A[4]*x+ A[5];
	double F_two = mod + Ftwo;
return F_two;}

//Produces a momemtum disturbution to roughly match av18. WIP
double distribution(double A,double number[],double F_x[],double help[], double rndm,int k, int run){
/* This function brings; A, a seed for the random number generator. It is currently the electron number.*/
// A is electorn number, number[] is Rho(k), F_x[] is K, help is a divider for large arrays
// max,k,run are not used!	F and num have 500 elements, only about 100 to 200 are filled.
	long double B;long double x;
// Generate a random number
	x= rndm;


	int count=0;
	int x_floor = 0;//48;
//Loop through the help array to speed up the final look up. help[50]
//	if(A>0){for(int h=0;h<=50;h++){if(x >= help[h]){x_floor=10*h;if(x<help[h+1]){break;}}}}
	double last_o,last_2;								
	int small = x_floor;	
	double diff = abs(x-number[small]);
	for(int i = 0;i<=500;i++){
				diff = abs(x-number[small]);

		if(	abs(x-number[i]) < diff){small = i;diff=abs(x-number[i]);}
			//if(i/10 == i/10.0){int k = floor(i/10);help[k] =number[small];} 

		if(	abs(x-number[i+1]) > diff){break;}
	count=i;
	if(i>2){
	if(F_x[i]==last_2&&F_x[i]==last_o){break;}
	last_o = F_x[i];
	last_2 = F_x[i-1];}
		}

	B=(F_x[small] + F_x[small+1])/2.0 ;
	double p = abs(B);

//cout << x << " "<< p<< " "<<small<<endl;
return p; }

double	cross_section(double e_final, double theta, double e_in,double F_two,double xb){
// This calculation use equations from the jlab pac30-petratos.pdf.

//Inelastic Cross Section
	double nu = e_in - e_final;
	double alpha = (1.0/137.0);
	double mp = .938;
	double A = alpha*alpha/(4*e_in*e_in*pow(sin(theta/2.0),4));
	double  B = (F_two/nu)*pow(cos(theta/2.0),2);
	double C =  (F_two/(mp*xb)) *pow(sin(theta/2),2);
	double diffcross = A*(B+C);
	// Convert 1/Gev^3 	to nb/(sr*GeV)
		   diffcross*=0.01973*0.01973*1e9;

return diffcross;}

