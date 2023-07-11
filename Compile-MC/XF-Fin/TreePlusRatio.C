#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <TH1F.h>
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <TMath.h>
#include <TLegend.h>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TPad.h>
#include <TF1.h>


void TreePlusRatio() {

Float_t e_p, mc_e_p, e_theta, mc_e_theta, p1_p, mc_p1_p, p1_theta, mc_p1_theta, Q2, mc_Q2, W, mc_W, Mx, mc_Mx, x, mc_x, y, mc_y, z, mc_z, xF, mc_xF, pT, mc_pT, zeta, mc_zeta, eta, mc_eta, trento_phi, mc_trento_phi, matching_e_pid, matching_p1_pid, mc_p1_parent, dummy, mc_dummy, delta_p, ratio_pT, ratio_eta, ratio_Q2, ratio_x, ratio_xF, ratio_zeta;

// Read in file
FILE *fp = fopen("/work/clas12/fatiha/MCfiles/MC_Eloss-Long.txt","r");//Open Eloss file

 TFile *hfile = hfile = TFile::Open("MCEloss.root","RECREATE");
 
TTree *outTree = new TTree("T","kinematics tree");
outTree->Branch("e_p",&e_p,"e_p/F");//Electron momentum
outTree->Branch("mc_e_p",&mc_e_p,"mc_e_p/F");//MonteCarlo Electron Momentum
outTree->Branch("e_theta",&e_theta,"e_theta/F");//Electron theta
outTree->Branch("mc_e_theta",&mc_e_theta,"mc_e_theta/F");//MonteCarlo Electron theta
outTree->Branch("p1_p",&p1_p,"p1_p/F");//p1 Momentum
outTree->Branch("mc_p1_p",&mc_p1_p,"mc_p1_p/F");//MonteCarlo p1 momentum
outTree->Branch("p1_theta",&p1_theta,"p1_theta/F");//p1 Theta
outTree->Branch("mc_p1_theta",&mc_p1_theta,"mc_p1_theta/F");//MonteCarlo p1 theta
outTree->Branch("Q2",&Q2,"Q2/F");//Q2
outTree->Branch("mc_Q2",&mc_Q2,"mc_Q2/F");//MonteCarlo Q2nk
outTree->Branch("W",&W,"W/F");//W
outTree->Branch("mc_W",&mc_W,"mc_W/F");//MonteCarlo W
outTree->Branch("Mx",&Mx,"Mx/F");//Mx
outTree->Branch("mc_Mx",&mc_Mx,"mc_Mx/F");//MonteCarlo Mx
outTree->Branch("x",&x,"x/F");//x
outTree->Branch("mc_x",&mc_x,"mc_x/F");//MonteCarlo x
outTree->Branch("y",&y,"y/F");//y
outTree->Branch("mc_y",&mc_y,"mc_y/F");//MonteCarlo y
outTree->Branch("z",&z,"z/F");//z
outTree->Branch("mc_z",&mc_z,"mc_z/F");//MonteCarlo z
outTree->Branch("xF",&xF,"xF/F");//xF
outTree->Branch("mc_xF",&mc_xF,"mc_xF/F");//MonteCarlo xF
outTree->Branch("pT",&pT,"pT/F");//pT
outTree->Branch("mc_pT",&mc_pT,"mc_pT/F");//MonteCarlo pT
outTree->Branch("zeta",&zeta,"zeta/F");//zeta
outTree->Branch("mc_zeta",&mc_zeta,"mc_zeta/F");//MonteCarlo zeta
outTree->Branch("eta",&eta,"eta/F");//eta
outTree->Branch("mc_eta",&mc_eta,"mc_eta/F");//MonteCarlo eta
outTree->Branch("trento_phi",&trento_phi,"trento_phi/F");//Trento phi?
outTree->Branch("mc_trento_phi",&mc_trento_phi,"mc_trento_phi/F");//MonteCarlo Trento Phi
outTree->Branch("matching_e_pid",&matching_e_pid,"matching_e_pid/F");//Electron PID
outTree->Branch("matching_p1_pid",&matching_p1_pid,"matching_p1_pid/F");//p1 PID
outTree->Branch("mc_p1_parent",&mc_p1_parent,"mc_p1_parent/F");//MonteCarlo PID?
outTree->Branch("delta_p",&delta_p,"delta_p/F");
outTree->Branch("ratio_pT",&ratio_pT,"ratio_pT/F");
outTree->Branch("ratio_eta",&ratio_eta,"ratio_eta/F");
outTree->Branch("ratio_Q2",&ratio_Q2,"ratio_Q2/F");
outTree->Branch("ratio_x",&ratio_x,"ratio_x/F");
outTree->Branch("ratio_xF",&ratio_xF,"ratio_xF/F");
outTree->Branch("ratio_zeta",&ratio_zeta,"ratio_zeta/F");

//"mc_ " are generated and the other are reconstructed
for (int k=0; k<5000000; k++) {//2441322 total events for ProtonFallinbending_MC_Eloss.txt
//for(int k=0; k<100000; k++){

	fscanf(fp, " %f %f %f %f %f %f", &e_p, &mc_e_p, &e_theta, &mc_e_theta, &p1_p, &mc_p1_p);

	fscanf(fp, " %f %f %f %f %f %f", &dummy, &mc_dummy, &Q2, &mc_Q2, &W, &mc_W);

	fscanf(fp, " %f %f %f %f %f %f", &Mx, &mc_Mx, &x, &mc_x, &y, &mc_y);

	fscanf(fp, " %f %f %f %f %f %f", &z, &mc_z, &xF, &mc_xF, &pT, &mc_pT);

	fscanf(fp, " %f %f %f %f %f %f", &zeta, &mc_zeta, &eta, &mc_eta, &trento_phi, &mc_trento_phi);

	fscanf(fp, " %f %f %f", &matching_e_pid, &matching_p1_pid, &mc_p1_parent);

	delta_p=mc_p1_p - p1_p;//take the difference in generated and reconstructed momentum
	p1_theta = dummy*(180/(22/7));
	mc_p1_theta = mc_dummy*(180/(22/7));

	ratio_eta=mc_eta / eta;
	ratio_pT=mc_pT / pT;
	ratio_Q2=mc_Q2 / Q2;
	ratio_x=mc_x / x;
	ratio_xF=mc_xF / xF;
	ratio_zeta=mc_zeta / zeta;

	if(mc_xF !=0) {
	outTree->Fill();
	}

}//end file read

cout << "ELoss File Read Successful" << endl;

outTree->Write();//write Tree to file MC.root
fclose(fp);//close txt file
hfile->Write();//write the f 

// CUTS & LOOPING -> CREATES TEXT FILE
///////////////////////////////////////////////////////////////
TCut NegHeli = ("helicity == -1");  // Negative helicity cut
TCut PosHeli = ("helicity == 1");  // Positive helicity cut
TCut Mxg1 = ("");            // Cut applied to all cuts, change depending on parameters

char cut[20];
strcpy(cut, "-1.0<mc_xF&&mc_xF<-0.9");   // First cut applied on variable of focus, change for different variable

float start = -0.9;      // Larger value of first cut made on target variable
float end = 1.1;        // Number one increment greater than the larger value of the final cut (getting around the for-loop being not inclusive)
float increment = 0.1;  // Size of cut

float mean1; // Create variables for printing
float mean2;
FILE *outputfile;
float ratio_mean;
float ratio_err;
float std1;
float std2;
float prop_error;
float N1;
float N2;
	
	outputfile = fopen("xF-Tree.txt" , "wb");  // Text file to print to

	for (float j = start; j < end; j += increment) {  // Loop to run through the target variable cuts
		
        	TCut VARcut = (cut);

        	TH1F*h1 = new TH1F("h1","h1",500, -1.0, 1.0);  // Histogram for mean of GENERATED xF
        	outTree->Draw("mc_xF>>h1", VARcut);
		mean1 = h1->GetMean();
		std1 = h1->GetStdDev();
		N1 = h1->GetEntries();

        	TH1F*h2 = new TH1F("h2","h2",500, -1.0, 1.0);  // Histogram for mean of RECONSTRUCTED xF
        	outTree->Draw("xF>>h2", VARcut);
		mean2 = h2->GetMean();
		std2 = h2->GetStdDev();
		N2 = h2->GetEntries();

		ratio_mean = mean1 / mean2; // generated divided by reconstructed

		prop_error = ratio_mean * sqrt( pow((std1/sqrt(N1))/mean1, 2) + pow((std2/sqrt(N2))/mean2, 2));

		// Print into a file for plotting
		fprintf(outputfile, "%f \t", mean1); 
		fprintf(outputfile, "%f \t", ratio_mean); //ratio_mean
		fprintf(outputfile, "%f \n", prop_error);
	
		sprintf(cut, "%f<mc_xF&&mc_xF<%f", j, j + increment);
		printf("%s \n", cut);
		
	}
fclose(outputfile);

// FORMAT & PLOTTING -> CREATES CANVAS
//////////////////////////////////////////////
double x_plot, y_plot, erry_plot;

FILE * fp1 = fopen("xF-Tree.txt", "r");

int n=0;

gStyle->SetPadTopMargin(0.1);
gStyle->SetPadBottomMargin(0.2);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadRightMargin(0.05);

TCanvas*c01 = new TCanvas("c01", "Asymmetry", 1000,700);
c01->Divide(1,1);
TGraphErrors*h3 = new TGraphErrors();

for(int k=0; k<25; k++){ 
fscanf(fp1, "%lf %lf %lf ", &x_plot, &y_plot, &erry_plot);

n = h3->GetN();

h3->SetPoint(n, x_plot, y_plot);
h3->SetPointError(n ,0, erry_plot);

}
fclose(fp1);
c01->cd(1);

h3->GetXaxis()->SetTitle("<x_{F}>");
h3->GetXaxis()->SetTitleSize(0.09);
h3->GetXaxis()->SetLabelSize(0.06);
h3->GetXaxis()->SetTitleOffset(0.8);
h3->GetYaxis()->SetTitle("<Gen>/<Recon>");
h3->GetYaxis()->SetTitleSize(0.09);
h3->GetYaxis()->SetLabelSize(0.06);
h3->GetYaxis()->SetTitleOffset(0.7);
h3->SetLineWidth(1.5);
h3->SetMarkerColor(1);
h3->SetMarkerSize(2.5);
h3->SetMarkerStyle(20);

h3->Draw("AP");
h3->GetXaxis()->SetRangeUser(-1.0,1.0);
h3->GetYaxis()->SetRangeUser(0.8, 1.2);
gPad->Modified();
gPad->Update();
c01->Print("xF-Tree.png");

}




