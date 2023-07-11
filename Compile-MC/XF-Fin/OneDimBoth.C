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

void OneDimBoth(){

	Float_t e_p, mc_e_p, e_theta, mc_e_theta, p1_p, mc_p1_p, p1_theta, mc_p1_theta, Q2, mc_Q2, W, mc_W, Mx, mc_Mx, x, mc_x, y, mc_y, z, mc_z, xF, mc_xF, pT, mc_pT, zeta, mc_zeta, eta, mc_eta, trento_phi, mc_trento_phi, matching_e_pid, matching_p1_pid, mc_p1_parent, dummy, mc_dummy, delta_p, delta_xF, delta_Q2;
// number of events 4075759,
// Read in file
	FILE *fp = fopen("/work/clas12/fatiha/MCfiles/MC_Eloss-Long.txt","r");

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
	outTree->Branch("mc_Q2",&mc_Q2,"mc_Q2/F");//MonteCarlo Q2
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
	outTree->Branch("delta_xF",&delta_xF,"delta_xF/F");
	outTree->Branch("delta_Q2",&delta_Q2,"delta_Q2/F");

// Determine number of lines to read in the imported text file
//number if events is 4075759 total events, but will avoid the end, in case of corrupted file. 
//for (int k=0; k<4075000; k++) {
for (int k=0; k<5000000; k++) {

	fscanf(fp, " %f %f %f %f %f %f", &e_p, &mc_e_p, &e_theta, &mc_e_theta, &p1_p, &mc_p1_p);

	fscanf(fp, " %f %f %f %f %f %f", &dummy, &mc_dummy, &Q2, &mc_Q2, &W, &mc_W);

	fscanf(fp, " %f %f %f %f %f %f", &Mx, &mc_Mx, &x, &mc_x, &y, &mc_y);

	fscanf(fp, " %f %f %f %f %f %f", &z, &mc_z, &xF, &mc_xF, &pT, &mc_pT);

	fscanf(fp, " %f %f %f %f %f %f", &zeta, &mc_zeta, &eta, &mc_eta, &trento_phi, &mc_trento_phi);

	fscanf(fp, " %f %f %f", &matching_e_pid, &matching_p1_pid, &mc_p1_parent);

	delta_p=mc_p1_p - p1_p;//take the difference in generated and reconstructed momentum
	p1_theta = dummy*(180/(22/7));
	mc_p1_theta = mc_dummy*(180/(22/7));

	delta_xF = xF - mc_xF; //reconstructed - generated

	delta_Q2 = Q2 - mc_Q2;

	if(mc_xF !=0){
	outTree->Fill();
	}

}
// If variables need to be created, they must be added to the tree and the loop
// Example of delta_pT on lines 60 and 82

cout << "ELoss File Read Successful" << endl;

outTree->Write();//write Tree to file MC.root
fclose(fp);//close txt file
hfile->Write();//write the f 

gStyle->SetPadTopMargin(0.1);
gStyle->SetPadBottomMargin(0.15);
gStyle->SetPadLeftMargin(0.12);
gStyle->SetPadRightMargin(0.12);
TGaxis::SetMaxDigits(3);
////////////////////////////////////////////////////////////////
// ONE-DIMENSIONAL HISTOGRAM
//TCanvas*c01 = new TCanvas("c01", "", 1000,700);
TCanvas*c00 = new TCanvas("c00", "", 1000,700);
gStyle->SetPalette(1);
c00->Divide(1,1);
TH1D*h1 = new TH1D("h1","", 500, -1.0, 1.0);
TH1D*h2 = new TH1D("h2","", 500, -1.0, 1.0);
TH1D*h3 = new TH1D("h3","", 400, -0.1, 0.1);

c00->cd(1);
h1->GetXaxis()->SetTitle("mc_x_{F}");  // The "#" format will print the symbol
h1->GetXaxis()->SetTitleSize(0.1);
h1->GetXaxis()->SetLabelSize(0.05);
h1->GetXaxis()->SetTitleOffset(0.6);
h1->GetYaxis()->SetTitle("Counts");
h1->GetYaxis()->SetTitleSize(0.1);
h1->GetYaxis()->SetLabelSize(0.05);
h1->GetYaxis()->SetTitleOffset(0.6);

h1->SetLineWidth(5);
h1->SetStats(1);
h1->SetLineColor(4);
h1->SetName("h1");
outTree->Draw("mc_xF>>h1");
h1->Draw();

h2->GetXaxis()->SetTitle("x_{F}");  // The "#" format will print the symbol
h2->GetXaxis()->SetTitleSize(0.1);
h2->GetXaxis()->SetLabelSize(0.05);
h2->GetXaxis()->SetTitleOffset(0.6);
h2->GetYaxis()->SetTitle("Counts");
h2->GetYaxis()->SetTitleSize(0.1);
h2->GetYaxis()->SetLabelSize(0.05);
h2->GetYaxis()->SetTitleOffset(0.6);


h2->SetLineWidth(5);
h2->SetStats(1);
h2->SetLineColor(2);
h2->SetName("h2");
outTree->Draw("xF>>h2");
h2->Draw();

h2->GetXaxis()->SetTitle("x_{F}");

h1->Draw("Same");

//legend to show what is blue and what is red 
auto legend = new TLegend(0.78,0.65,0.98,0.75);
legend->AddEntry("h1","Generated","l"); // mc_ = generated = blue
legend->AddEntry("h2","Reconstructed","l"); // no mc_ = reconstructed = red



legend->Draw(); 
c00->Update();
gPad->Modified();
gPad->Update();

///////////////////////////////
//PLOT THE DIFFERENCE HERE 

TCanvas*c01 = new TCanvas("c01", "", 1000,700);
c01->cd(1);
h3->GetXaxis()->SetTitle("#delta x_{F}");  // The "#" format will print the symbol
h3->GetXaxis()->SetTitleSize(0.1);
h3->GetXaxis()->SetLabelSize(0.05);
h3->GetXaxis()->SetTitleOffset(0.6);
h3->GetYaxis()->SetTitle("Counts");
h3->GetYaxis()->SetTitleSize(0.1);
h3->GetYaxis()->SetLabelSize(0.05);
h3->GetYaxis()->SetTitleOffset(0.6);

h3->SetLineWidth(5);
h3->SetStats(1);
h3->SetLineColor(1);
outTree->Draw("delta_xF>>h3");


c00->Print("xF-Gen-Recon.png");
c01->Print("Delta-xF.png");
}






