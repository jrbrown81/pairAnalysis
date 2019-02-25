/**********************************************************************
Created by Alex Cherlin, 19/07/2012
**********************************************************************/

#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TSystem.h"
#include <TLatex.h>
#include <TApplication.h>
#include <math.h>
#include <TSystemDirectory.h>
#include <TChain.h>
#include <TLegend.h>
#include <fstream>
#include <TGraph.h>
#include <list>
#include <vector>
#include <limits>
#include <TRandom3.h>
#include <dirent.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TTree.h>
#include <TString.h>

using std::setw;
using namespace std;
ifstream in_file;

//
// Running parameters
//

TString pathRoot = "./";

const TString ver = "v4";
const TString setup_file = "./analysis_setup_pairAnalysis_" + ver + ".txt";

const TString fname = "data";

const Float_t G_FWHM = 2.35482;
const Float_t screen_width_def = 1600;
const Float_t screen_height_def = 900;

////////////////////////////////////////////////////////////////////////////////////////
// Running parameters, set here some default values.
// These parameters are re-set with values read from the analysis setup file.
////////////////////////////////////////////////////////////////////////////////////////

Bool_t printOutSetupFileParameters = kFALSE;
TString spectrumName = "";
TString pathRoot_in_calib = "";
TString TcalibFileName = "";
TString fname_pixel_mapping_file = "";
TString num_str;
Int_t prntOutF = 100000;

Int_t firstPixel = 1;
Int_t lastPixel = 968;
Int_t nPixXY = 22;
Int_t nAMs = 2;

Float_t minE = 0;
Float_t maxE = 600; // keV
Int_t nBinsE = 1200;
Float_t minEcath = 0; // keV
Float_t maxEcath = 700; // keV
Int_t nBinsEcath = 1400;
Float_t minPhotopeakE4PeakSearch = 460;
Float_t maxPhotopeakE4PeakSearch = 550;

Int_t minNTrigs4H = 0;
Int_t maxNTrigs4H = 10;
Int_t minNTrigs4H_2D = 0;
Int_t maxNTrigs4H_2D = 6;
Int_t minClusterSize4H = 0;
Int_t maxClusterSize4H = 8;
Int_t nBinsTheta4H = 18;

Float_t Na22Energy = 511;
Float_t angleOfSecondHead = 90;
Float_t relativePhiAngle = 90;
Float_t EwindowFor4PhiAnalyis_min = 530;
Float_t EwindowFor4PhiAnalyis_max = 530;
Float_t ThetaWindowFor4PhiAnalyis_min = 60;
Float_t ThetaWindowFor4PhiAnalyis_max = 90;

Bool_t makeTimingStuff = kTRUE;

// 0 - all, 1 - fully internal (disable two most external rows), 2 - internal (disable one external row), 3 - NOT fully internal, 4 - NOT internal
Int_t typeOfPixelsToUse = 0;

Int_t nBins_dPhi = 18;
Float_t maxDistanceBetweenClusters4ComptonPair = 7;

// Cluster reconstruction stuff
Int_t maxClusterSize = 10;
Float_t minPixelEnergyThr4ClusterReconstruction = 15;
Int_t neighbourSearchWindow = 1;
Int_t maxClusterSizeToSkipFrame = 25;
Bool_t useOnlyTriggeredPixelsInCluster = kFALSE;

////////////////////////////////////////////////////////////////////////////////////////
// End of running parameters
////////////////////////////////////////////////////////////////////////////////////////

Int_t TempChannelBase = 62000;
Int_t Cath1ChannelBase = 60000;

TString str;
string line;
Bool_t readAnalysisSetupFile(TString);
Int_t getCath1Channel(const Int_t, const Int_t, const Int_t);
Float_t getPhiAngleDeg(const Float_t, const Float_t);

const Int_t nCath = 1100;
const Int_t cathChannelNumberOffset = 60000;

Float_t **tCalibAnodes;
Float_t tCalibCathodes[nCath][2] = {0};
const Float_t cathodeSlopeADCperKEV = 3.104;
const Float_t cathodeOffsetADC = 500.33;
const Int_t newNeighbourFlag = 5000;
const Int_t buffSizeSingleHits = 50;
const Int_t nVarsSingleHits = 4; // energy, x, y, z
const Int_t buffSizeComptonChains = 50;
const Int_t nLegsComptonChains = 10;
const Int_t nVarsComptonChains = 4; // energy, x, y, z

Bool_t ***pixelStatus;
TH2F *image;
TH2F *imageNEG;
TH2F *imageBG;
TH2F *imageBGneg;
TH2I *imageTiming;
TH2I *imagePos;
TH2I *imageTrig;
TH2I **pixelPattern;
TH2I **pixelTypePattern;
TH2I **usedCentrePixelPattern;
Float_t ***singleHits;
Float_t ****comptonChains;
Int_t *comptonChainsLength;
Int_t *singleHitsLength;
Int_t **comptonLegsLength;

std::vector<std::vector<Int_t>> matrix_trigs;
std::vector<std::vector<Int_t>> matrix_flags;
std::vector<std::vector<Float_t>> matrix_E;
std::vector<Int_t> buffX;
std::vector<Int_t> buffY;
std::vector<Int_t> buffXnew;
std::vector<Int_t> buffYnew;

std::vector<std::vector<Float_t>> buffClusterX;
std::vector<std::vector<Float_t>> buffClusterY;
std::vector<std::vector<Float_t>> buffClusterE;
std::vector<std::vector<Float_t>> buffClusterArea;
std::vector<std::vector<Float_t>> buffClusterFlag;
std::vector<std::vector<Float_t>> buffClusterTrigs;
std::vector<Float_t> buffClusterXloc;
std::vector<Float_t> buffClusterYloc;
std::vector<Float_t> buffClusterEloc;
std::vector<Float_t> buffClusterArealoc;
std::vector<Float_t> buffClusterFlagloc;
std::vector<Float_t> buffClusterTrigsloc;

TChain* chain_events;
TChain** chain1p;
Long64_t *event;
Int_t *nTrigPixels;
Int_t *AM;
Int_t *GM;
Int_t *nAMsInEvent;
Int_t *pixel;
Float_t *E;
Float_t *E_neg;
Int_t *timeStamp;
Int_t *triggerFlag;
Int_t *nloop;
Int_t *timeDetect;
Int_t *timeDetectPos;
Float_t *Temp;
Int_t *pos;
Long64_t event_eev;
Int_t *AM_flag_eev;
Int_t *cath_flag_eev;
Int_t nAMsInEvent_eev;
Int_t *nTrigPixels_eev;
Int_t *cntEntryTree;
Long64_t *Nev;

TH1F **singleHitPixelSpec;
TH1F **allClustersSpec;
TH1F **singleClusterEventsSpec;
TH1F **summedClusterEventsSpec;
TH2F *singleHitPixelSpecCorr;
TH1F **nTrigsInEvent;
TH2F **nTrigsInEvent_maxE;
TH2F **nTrigsInEvent_maxEcluster;
TH1F **comptonSummedSpecClusters;
TH1F **comptonHitsSpecClusters;
TH2F **comptonHitsSpecClustersCorr;
TH2F *comptonSummedSpecClustersCorr;
TH1F **comptonPixelSpec;
TH1F **phiAngle;
TH1F **phiAnglePairs;
TH1F *dPhiAngle;
TH1F *dPhiAngle_summed;
TH1F *dPhiAngle_Ewin;
TH2F *dPhiAngle_2ClusterE;
TH2F *dPhiAngle_totalE;
TH2F *phiAngleCorr;
TH2F **dPhiPlaneAM;
TH2F **firstComptonCluster;
TH2F **secondComptonCluster;
TH2F **clusterE_vsSize;
TH2F *singleClusterSpecCorr;
TH2F *summedClusterSpecCorr;
TH1F **thetaFromFirstClusterE;
TH2F **thetaFromFirstClusterE_dPhi;
TH1F **clusterSize;
TH2F *theta0_theta1_FromFirstClusterE;
TH2F **nTrigs_clusterSize;
TH2F **clusterE_nTrigs;
TH2F **clusterMaxE_vsSize;
TH2F **clusterSecE_vsSize;
TH2F **clusterSizeMaxE_SecE;

TRandom3 *rand3;
TStopwatch localTimer;
Float_t totalTimeElapced = 0;

Float_t getFWHM(const TH1F*, const Float_t, const Float_t, const Float_t, const Float_t, const Float_t, const Float_t, Float_t&, Float_t&);
void getHistoPeak(const TH1F*, Float_t, Float_t, Float_t&, Float_t&);
Float_t getDistanceBetweenPixels(const Int_t, const Int_t, const Int_t);
Bool_t getPixel2Dcoordinates(const Int_t, const Int_t, Int_t&, Int_t&);
Bool_t analyseNextEvent(const Int_t);
Bool_t analyseAM(const Int_t, const Int_t);
Float_t getPhiAngleDeg(const Float_t, const Float_t);
void getNewNeighbours(Int_t, Int_t, Int_t);
Float_t getThetaFromEnergy(const Float_t, const Float_t);
Float_t getScatteredEnergyFromAngle(const Float_t, const Float_t);

int main()
{
	if (!readAnalysisSetupFile(setup_file))
	{
		cerr << "ERROR reading analysis setup file " << pathRoot + "/" + setup_file << ". Exiting." << endl;
		return 0;
	}
	
	pixelPattern = new TH2I*[nAMs];
	pixelTypePattern = new TH2I*[nAMs];
	usedCentrePixelPattern = new TH2I*[nAMs];
	pixelStatus = new Bool_t**[nAMs];
	for (Int_t i = 0; i < nAMs; i++)
	{
		pixelPattern[i] = new TH2I(Form("pixelPattern_AM%d",i),Form("Pixel pattern, AM%d",i),nPixXY,0,nPixXY,nPixXY,0,nPixXY);
		pixelPattern[i]->GetXaxis()->SetTitle("Pixel number");
		pixelPattern[i]->GetXaxis()->SetTitleOffset(1.2);
		pixelPattern[i]->GetYaxis()->SetTitle("Pixel number");
		pixelPattern[i]->GetYaxis()->SetTitleOffset(1.2);
		pixelTypePattern[i] = new TH2I(Form("pixelTypePattern_AM%d",i),Form("Pixel type pattern, AM%d",i),nPixXY,0,nPixXY,nPixXY,0,nPixXY);
		pixelTypePattern[i]->SetTitle("Pattern of pixels used in the analysis");
		pixelTypePattern[i]->GetXaxis()->SetTitle("Pixel number");
		pixelTypePattern[i]->GetXaxis()->SetTitleOffset(1.2);
		pixelTypePattern[i]->GetYaxis()->SetTitle("Pixel number");
		pixelTypePattern[i]->GetYaxis()->SetTitleOffset(1.2);
		usedCentrePixelPattern[i] = new TH2I(Form("usedCentrePixelPattern_AM%d",i),Form("Used centre pixel pattern, AM%d",i),nPixXY,0,nPixXY,nPixXY,0,nPixXY);
		usedCentrePixelPattern[i]->SetTitle("Pattern of pixels used in the analysis");
		usedCentrePixelPattern[i]->GetXaxis()->SetTitle("Pixel number");
		usedCentrePixelPattern[i]->GetXaxis()->SetTitleOffset(1.2);
		usedCentrePixelPattern[i]->GetYaxis()->SetTitle("Pixel number");
		usedCentrePixelPattern[i]->GetYaxis()->SetTitleOffset(1.2);

		pixelStatus[i] = new Bool_t*[nPixXY];
		for (Int_t n = 0; n < nPixXY; n++)
		{
			pixelStatus[i][n] = new Bool_t[nPixXY];
			for (Int_t m = 0; m < nPixXY; m++)
			{
				pixelStatus[i][n][m] = 0;
			}
		}
	}
	
	while (pathRoot_in_calib.BeginsWith(" ")) pathRoot_in_calib.Remove(0,1);
	if (!((TString) pathRoot_in_calib[pathRoot_in_calib.Length()]).IsAlnum()) pathRoot_in_calib = pathRoot_in_calib.Chop();
	
	while (fname_pixel_mapping_file.BeginsWith(" ")) fname_pixel_mapping_file.Remove(0,1);
	if (!((TString) fname_pixel_mapping_file[fname_pixel_mapping_file.Length()]).IsAlnum()) fname_pixel_mapping_file = fname_pixel_mapping_file.Chop();
	for (Int_t im = 0; im < nAMs; im++)
	{
		TString fname_pxl_cfg = Form("%s/%s_AM%d.txt",pathRoot_in_calib.Data(),fname_pixel_mapping_file.Data(),im);
		cout << "fname_pxl_cfg  = " << fname_pxl_cfg << endl;
		if (gSystem->AccessPathName(fname_pxl_cfg))
		{
			cerr << "ERROR: Pixel config file " << fname_pxl_cfg << " doesn't exist. Exiting." << endl;
			return 0;
		}
		in_file.open(fname_pxl_cfg, ios::in);
		for (Int_t i = nPixXY; i >= 1 && !in_file.eof(); i--)
		{
			for (Int_t j = 0; j < nPixXY && !in_file.eof(); j++)
			{
				Int_t pix = -1;
				in_file >> pix;
				if (pix < firstPixel || pix > lastPixel) continue;
				pixelPattern[im]->SetBinContent(j+1,i,pix);
			}
		}
		in_file.close();
		in_file.clear();
	}
	
	gStyle->SetOptStat(0);
	//gStyle->SetOptStat(10);
	//gStyle->SetOptFit(111111);
	TH1::AddDirectory(0);
	gErrorIgnoreLevel = 1001;
	gStyle->SetPalette(1);
	rand3 = new TRandom3();
	rand3->SetSeed(0);
	
	if (makeTimingStuff)
	{
		tCalibAnodes = new Float_t*[lastPixel+1];
		for (Int_t i = 0; i <= lastPixel; i++)
		{
			tCalibAnodes[i] = new Float_t[2];
			tCalibAnodes[i][0] = 0;
			tCalibAnodes[i][1] = 0;
		}
		
		while (TcalibFileName.BeginsWith(" ")) TcalibFileName.Remove(0,1);
		if (!((TString) TcalibFileName[TcalibFileName.Length()]).IsAlnum()) TcalibFileName = TcalibFileName.Chop();
		TString fname2 = Form("%s/%s.txt",pathRoot_in_calib.Data(),TcalibFileName.Data());
		if (gSystem->AccessPathName(fname2))  // Strange convention - this function return 0;s 1 (true) if path name doesn't exist !!!!
		{
			cerr << "ERROR: Timing calibration data file " << fname2 << " doesn't exist. Exiting." << endl;
			return 0;
		}
		in_file.open(fname2, ios::in);
		while (!in_file.eof())
		{
			Int_t pix = -1;
			Float_t a1 = -1, a2 = -1;
			in_file >> pix >> a2 >> a1;
			if (pix > 0 && pix <= lastPixel)
			{
				tCalibAnodes[pix][0] = a2; // slope
				tCalibAnodes[pix][1] = a1; // intercept
			}
			if (pix >= cathChannelNumberOffset)
			{
				tCalibCathodes[pix - cathChannelNumberOffset][0] = a2;
				tCalibCathodes[pix - cathChannelNumberOffset][1] = a1;
			}
		}
		in_file.close();
		in_file.clear();
	}
	
	for (Int_t im = 0; im < nAMs; im++)
	{
		if (typeOfPixelsToUse <= 2)
		{
			for (Int_t i = 0; i <= lastPixel; i++)
			{
				Int_t xx1 = -1, yy1 = -1;
				if (getPixel2Dcoordinates(i,im,xx1,yy1))
				{
					pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (xx1 == 1 && (typeOfPixelsToUse == 1 || typeOfPixelsToUse == 2)) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (xx1 == 2 && typeOfPixelsToUse == 1) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (xx1 == 10 && typeOfPixelsToUse == 1) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (xx1 == 11 && (typeOfPixelsToUse == 1 || typeOfPixelsToUse == 2)) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (xx1 == 12 && (typeOfPixelsToUse == 1 || typeOfPixelsToUse == 2)) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (xx1 == 13 && typeOfPixelsToUse == 1) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (xx1 == 21 && typeOfPixelsToUse == 1) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (xx1 == 22 && (typeOfPixelsToUse == 1 || typeOfPixelsToUse == 2)) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (yy1 == 1 && (typeOfPixelsToUse == 1 || typeOfPixelsToUse == 2)) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (yy1 == 2 && typeOfPixelsToUse == 1) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (yy1 == 10 && typeOfPixelsToUse == 1) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (yy1 == 11&& (typeOfPixelsToUse == 1 || typeOfPixelsToUse == 2)) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (yy1 == 12 && (typeOfPixelsToUse == 1 || typeOfPixelsToUse == 2)) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (yy1 == 13 && typeOfPixelsToUse == 1) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (yy1 == 21 && typeOfPixelsToUse == 1) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (yy1 == 22 && (typeOfPixelsToUse == 1 || typeOfPixelsToUse == 2)) pixelStatus[im][xx1-1][yy1-1] = kFALSE;
				}
			}
		}
		if (typeOfPixelsToUse == 3 || typeOfPixelsToUse == 4)
		{
			for (Int_t i = 0; i <= lastPixel; i++)
			{
				Int_t xx1 = -1, yy1 = -1;
				if (getPixel2Dcoordinates(i,im,xx1,yy1))
				{
					pixelStatus[im][xx1-1][yy1-1] = kFALSE;
					if (xx1 == 1 && (typeOfPixelsToUse == 3 || typeOfPixelsToUse == 4)) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (xx1 == 2 && typeOfPixelsToUse == 3) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (xx1 == 10 && typeOfPixelsToUse == 3) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (xx1 == 11 && (typeOfPixelsToUse == 3 || typeOfPixelsToUse == 4)) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (xx1 == 12 && (typeOfPixelsToUse == 3 || typeOfPixelsToUse == 4)) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (xx1 == 13 && typeOfPixelsToUse == 3) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (xx1 == 21 && typeOfPixelsToUse == 3) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (xx1 == 22 && (typeOfPixelsToUse == 3 || typeOfPixelsToUse == 4)) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (yy1 == 1 && (typeOfPixelsToUse == 3 || typeOfPixelsToUse == 4)) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (yy1 == 2 && typeOfPixelsToUse == 3) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (yy1 == 10 && typeOfPixelsToUse == 3) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (yy1 == 11&& (typeOfPixelsToUse == 3 || typeOfPixelsToUse == 4)) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (yy1 == 12 && (typeOfPixelsToUse == 3 || typeOfPixelsToUse == 4)) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (yy1 == 13 && typeOfPixelsToUse == 3) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (yy1 == 21 && typeOfPixelsToUse == 3) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
					if (yy1 == 22 && (typeOfPixelsToUse == 3 || typeOfPixelsToUse == 4)) pixelStatus[im][xx1-1][yy1-1] = kTRUE;
				}
			}
		}
		
		for (Int_t im = 0; im < nAMs; im++)
		{
			for (Int_t i = 0; i < nPixXY; i++)
			{
				for (Int_t j = 0; j < nPixXY; j++)
				{
					pixelTypePattern[im]->SetBinContent(i+1,j+1,pixelStatus[im][i][j]);
				}
			}
		}
	}
	
	TString subdir = "Pair analysis " + ver;
	subdir = Form("%s pixType%d",subdir.Data(),typeOfPixelsToUse);
	TString outputPathPairs = Form("%s%s/pairAnalysis",pathRoot.Data(),subdir.Data());
	if (subdir.Length() > 1)
	{
		gSystem->MakeDirectory(pathRoot+subdir);
		gSystem->MakeDirectory(outputPathPairs);
	}
	
	TString *outputPathAM = new TString[nAMs];
	for (Int_t im = 0; im < nAMs; im++)
	{
		outputPathAM[im] = Form("%s%s/AM%d",pathRoot.Data(),subdir.Data(),im);
		gSystem->MakeDirectory(outputPathAM[im]);
	}
	TFile *hfile = new TFile(Form("%s%s/histos.root",pathRoot.Data(),subdir.Data()),"recreate");
	
	image = new TH2F("image","Event image",nPixXY,0,nPixXY,nPixXY,0,nPixXY);
	image->GetXaxis()->SetTitle("Pixel number");
	image->GetXaxis()->SetTitleOffset(1.2);
	image->GetYaxis()->SetTitle("Pixel number");
	image->GetYaxis()->SetTitleOffset(1.2);
	
	imageNEG = new TH2F("imageNEG","Event image, negative energy",nPixXY,0,nPixXY,nPixXY,0,nPixXY);
	imageNEG->GetXaxis()->SetTitle("Pixel number");
	imageNEG->GetXaxis()->SetTitleOffset(1.2);
	imageNEG->GetYaxis()->SetTitle("Pixel number");
	imageNEG->GetYaxis()->SetTitleOffset(1.2);
	
	imageBG = new TH2F("imageBG","Event image BG",nPixXY,0,nPixXY,nPixXY,0,nPixXY);
	imageBG->GetXaxis()->SetTitle("Pixel number");
	imageBG->GetXaxis()->SetTitleOffset(1.2);
	imageBG->GetYaxis()->SetTitle("Pixel number");
	imageBG->GetYaxis()->SetTitleOffset(1.2);
	
	imageBGneg = new TH2F("imageBGneg","Event image BG negative",nPixXY,0,nPixXY,nPixXY,0,nPixXY);
	imageBGneg->GetXaxis()->SetTitle("Pixel number");
	imageBGneg->GetXaxis()->SetTitleOffset(1.2);
	imageBGneg->GetYaxis()->SetTitle("Pixel number");
	imageBGneg->GetYaxis()->SetTitleOffset(1.2);
	
	imagePos = new TH2I("imagePos","Image of pixel energy ordering in event",nPixXY,0,nPixXY,nPixXY,0,nPixXY);
	imagePos->GetXaxis()->SetTitle("Pixel number");
	imagePos->GetXaxis()->SetTitleOffset(1.2);
	imagePos->GetYaxis()->SetTitle("Pixel number");
	imagePos->GetYaxis()->SetTitleOffset(1.2);
	
	imageTrig = new TH2I("imageTrig","Triggered pixels in event",nPixXY,0,nPixXY,nPixXY,0,nPixXY);
	imageTrig->GetXaxis()->SetTitle("Pixel number");
	imageTrig->GetXaxis()->SetTitleOffset(1.2);
	imageTrig->GetYaxis()->SetTitle("Pixel number");
	imageTrig->GetYaxis()->SetTitleOffset(1.2);
	
	imageTiming = new TH2I("imageTiming","Anode timing event image",nPixXY,0,nPixXY,nPixXY,0,nPixXY);
	imageTiming->GetXaxis()->SetTitle("Pixel number");
	imageTiming->GetXaxis()->SetTitleOffset(1.2);
	imageTiming->GetYaxis()->SetTitle("Pixel number");
	imageTiming->GetYaxis()->SetTitleOffset(1.2);
	
	chain_events = new TChain("tree2_events");
	chain1p = new TChain*[nAMs];
	for (Int_t im = 0; im < nAMs; im++) chain1p[im] = new TChain(Form("tree2_AM%d",im));
	
	TString fname2do;
	Int_t Nfiles = 0;
	TSystemDirectory dir("./","./");
	TList *files = dir.GetListOfFiles();
	if (files)
	{
		TSystemFile *file;
		TIter next(files);
		cout << "Counting valid root files in " << pathRoot << " ... ";
		while (file=(TSystemFile*)next())
		{
			fname2do = file->GetName();
			if (file->IsDirectory()) files->Remove(file);
			if (!fname2do.EndsWith(".root")) files->Remove(file);
			if (!fname2do.Contains("Ecal")) files->Remove(file);
		}
		next.Reset();
		cout << "List of files to read:" << endl;
		while (file=(TSystemFile*)next())
		{
			fname2do = file->GetName();
			Nfiles++;
		}
		next.Reset();
		if (Nfiles < 1)
		{
			cerr << "No valid ROOT files to read. Exiting." << endl;
			return kFALSE;
		}
		cout << "Total of " << Nfiles << " files found." << endl;
		while (file=(TSystemFile*)next())
		{
			fname2do = pathRoot + file->GetName();
			cout << "Reading file:" << endl;
			cout << fname2do << endl;
			for (Int_t im = 0; im < nAMs; im++) chain1p[im]->Add(fname2do);
			chain_events->Add(fname2do);
		}
	}
	
	event = new Long64_t[nAMs];
	nTrigPixels = new Int_t[nAMs];
	AM = new Int_t[nAMs];
	GM = new Int_t[nAMs];
	nAMsInEvent = new Int_t[nAMs];
	pixel = new Int_t[nAMs];
	E = new Float_t[nAMs];
	E_neg = new Float_t[nAMs];
	timeStamp = new Int_t[nAMs];
	triggerFlag = new Int_t[nAMs];
	nloop = new Int_t[nAMs];
	timeDetect = new Int_t[nAMs];
	timeDetectPos = new Int_t[nAMs];
	Temp = new Float_t[nAMs];
	pos = new Int_t[nAMs];
	
	AM_flag_eev = new Int_t[nAMs];
	cath_flag_eev = new Int_t[nAMs];
	nTrigPixels_eev = new Int_t[nAMs];
	cntEntryTree = new Int_t[nAMs];
	chain_events->SetBranchAddress("event",&event_eev);
	chain_events->SetBranchAddress("nAMsInEvent",&nAMsInEvent_eev);
	chain_events->SetBranchAddress("nTrigPixels",nTrigPixels_eev);
	chain_events->SetBranchAddress("AM_flag",AM_flag_eev);
	chain_events->SetBranchAddress("cath_flag",cath_flag_eev);
	chain_events->SetBranchStatus("*",1);
	Long_t nEvents2Analyse = chain_events->GetEntries();
	TString nev_str = Form("%ld",nEvents2Analyse);
	for (Int_t is = nev_str.Length(); is > 1; is--)
		if ((is-nev_str.Length()+2)%4 == 0) nev_str.Insert(is-1,",",1);
	cout << "nEvents2Analyse = " << nev_str << endl;

	Nev = new Long64_t[nAMs];
	singleHits = new Float_t**[nAMs];
	//singleHitsLength = (Int_t *) malloc(sizeof(Int_t)*nAMs);
	singleHitsLength = new Int_t[nAMs];
	//comptonChainsLength = (Int_t *) malloc(sizeof(Int_t)*nAMs);
	comptonChainsLength = new Int_t[nAMs];
	//comptonLegsLength = (Int_t **) malloc(sizeof(Int_t*)*nAMs);
	comptonLegsLength = new Int_t*[nAMs];
	comptonChains = new Float_t***[nAMs];
	for (Int_t im = 0; im < nAMs; im++)
	{
		num_str = Form("%d",chain1p[im]->GetNbranches());
		for (Int_t is = num_str.Length(); is > 1; is--)
			if ((is-num_str.Length()+2)%4 == 0) num_str.Insert(is-1,",",1);
		cout << "Nvar = " << num_str << endl;
		
		Nev[im] = chain1p[im]->GetEntries();
		num_str = Form("%lld",Nev[im]);
		for (Int_t is = num_str.Length(); is > 1; is--)
			if ((is-num_str.Length()+2)%4 == 0) num_str.Insert(is-1,",",1);
		cout << "Nev = " << num_str << endl;
		
		// Assign local variable addresses to variable pointers in the ntuple
		chain1p[im]->SetBranchAddress("event",&event[im]);
		chain1p[im]->SetBranchAddress("nTrigPixels",&nTrigPixels[im]);
		chain1p[im]->SetBranchAddress("AM",&AM[im]);
		chain1p[im]->SetBranchAddress("GM",&GM[im]);
		chain1p[im]->SetBranchAddress("nAMsInEvent",&nAMsInEvent[im]);
		chain1p[im]->SetBranchAddress("pixel",&pixel[im]);
		chain1p[im]->SetBranchAddress("E",&E[im]);
		chain1p[im]->SetBranchAddress("E_neg",&E_neg[im]);
		chain1p[im]->SetBranchAddress("pos",&pos[im]);
		chain1p[im]->SetBranchAddress("time",&timeStamp[im]);
		chain1p[im]->SetBranchAddress("triggerFlag",&triggerFlag[im]);
		chain1p[im]->SetBranchAddress("nloop",&nloop[im]);
		chain1p[im]->SetBranchAddress("timeDetect",&timeDetect[im]);
		chain1p[im]->SetBranchAddress("timeDetectPos",&timeDetectPos[im]);
		chain1p[im]->SetBranchAddress("Temp",&Temp[im]);
		chain1p[im]->SetBranchStatus("*",1);
		chain1p[im]->GetEntry(0);
		cntEntryTree[im] = 1;
		
		singleHits[im] = new Float_t*[buffSizeSingleHits];
		for (Int_t b = 0; b < buffSizeSingleHits; b++)
		{
			singleHits[im][b] = new Float_t[nVarsSingleHits];
			for (Int_t b2 = 0; b2 < nVarsSingleHits; b2++) singleHits[im][b][b2] = 0;
		}
		
		comptonLegsLength[im] = (Int_t*) malloc(sizeof(Int_t)*buffSizeComptonChains);
		comptonLegsLength[im] = new Int_t[buffSizeComptonChains];
		comptonChains[im] = new Float_t**[buffSizeComptonChains];
		for (Int_t b = 0; b < buffSizeComptonChains; b++)
		{
			comptonChains[im][b] = new Float_t*[nLegsComptonChains];
			for (Int_t b2 = 0; b2 < nLegsComptonChains; b2++) comptonChains[im][b][b2] = new Float_t[nVarsComptonChains];
		}
	}
	
	nTrigsInEvent = new TH1F*[nAMs];
	nTrigsInEvent_maxE = new TH2F*[nAMs];
	clusterE_vsSize = new TH2F*[nAMs];
	clusterMaxE_vsSize = new TH2F*[nAMs];
	clusterSecE_vsSize = new TH2F*[nAMs];
	clusterSizeMaxE_SecE = new TH2F*[nAMs];
	nTrigsInEvent_maxEcluster = new TH2F*[nAMs];
	singleHitPixelSpec = new TH1F*[nAMs];
	allClustersSpec = new TH1F*[nAMs];
	singleClusterEventsSpec = new TH1F*[nAMs];
	summedClusterEventsSpec = new TH1F*[nAMs];
	comptonSummedSpecClusters = new TH1F*[nAMs];
	comptonHitsSpecClusters = new TH1F*[nAMs];
	comptonHitsSpecClustersCorr = new TH2F*[nAMs];
	comptonPixelSpec = new TH1F*[nAMs];
	phiAngle = new TH1F*[nAMs];
	phiAnglePairs = new TH1F*[nAMs];
	dPhiPlaneAM = new TH2F*[nAMs];
	firstComptonCluster = new TH2F*[nAMs];
	secondComptonCluster = new TH2F*[nAMs];
	thetaFromFirstClusterE = new TH1F*[nAMs];
	thetaFromFirstClusterE_dPhi = new TH2F*[nAMs];
	clusterSize = new TH1F*[nAMs];
	nTrigs_clusterSize = new TH2F*[nAMs];
	clusterE_nTrigs = new TH2F*[nAMs];
	
	for (Int_t im = 0; im < nAMs; im++)
	{
		thetaFromFirstClusterE[im] = new TH1F(Form("thetaFromFirstClusterE_AM%d",im),Form("%s, #theta calculated from first interaction energy, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180);
		thetaFromFirstClusterE[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromFirstClusterE[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));

		thetaFromFirstClusterE_dPhi[im] = new TH2F(Form("thetaFromFirstClusterE_dPhi_AM%d",im),Form("%s, #theta calculated from first interaction energy, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180,nBins_dPhi,-180,180);
		thetaFromFirstClusterE_dPhi[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromFirstClusterE_dPhi[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));
		thetaFromFirstClusterE_dPhi[im]->GetYaxis()->SetTitleOffset(1.2);
		thetaFromFirstClusterE_dPhi[im]->GetYaxis()->SetTitle(Form("#Delta#varphi AM%d, degrees",im));

		nTrigsInEvent[im] = new TH1F(Form("nTrigsInEvent_AM%d",im),Form("%s, number of triggers per event, AM%d",spectrumName.Data(),im),maxNTrigs4H-minNTrigs4H,minNTrigs4H,maxNTrigs4H);
		nTrigsInEvent[im]->GetXaxis()->SetTitleOffset(1.2);
		nTrigsInEvent[im]->GetXaxis()->SetTitle(Form("Number of triggers in AM%d",im));
		
		nTrigsInEvent_maxE[im] = new TH2F(Form("nTrigsInEvent_maxE_AM%d",im),Form("%s, number of triggers per event vs pixel with maximum energy, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,maxNTrigs4H_2D-minNTrigs4H_2D,minNTrigs4H_2D,maxNTrigs4H_2D);
		nTrigsInEvent_maxE[im]->GetXaxis()->SetTitleOffset(1.2);
		nTrigsInEvent_maxE[im]->GetXaxis()->SetTitle(Form("AM%d pixel energy, keV",im));		
		nTrigsInEvent_maxE[im]->GetYaxis()->SetTitleOffset(1.4);
		nTrigsInEvent_maxE[im]->GetYaxis()->SetTitle(Form("Number of triggers in AM%d",im));		
		
		nTrigs_clusterSize[im] = new TH2F(Form("nTrigs_clusterSize_AM%d",im),Form("%s, number of triggers in cluster vs cluster size, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H,maxNTrigs4H_2D-minNTrigs4H_2D,minNTrigs4H_2D,maxNTrigs4H_2D);
		nTrigs_clusterSize[im]->GetXaxis()->SetTitleOffset(1.2);
		nTrigs_clusterSize[im]->GetXaxis()->SetTitle(Form("AM%d cluster size, pixels",im));	
		nTrigs_clusterSize[im]->GetYaxis()->SetTitleOffset(1.4);
		nTrigs_clusterSize[im]->GetYaxis()->SetTitle("Number of triggers in cluster");		
		
		clusterSize[im] = new TH1F(Form("clusterSize_AM%d",im),Form("%s, number of triggers per event vs pixel with maximum energy, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		clusterSize[im]->GetXaxis()->SetTitleOffset(1.4);
		clusterSize[im]->GetXaxis()->SetTitle(Form("AM%d cluster size, pixels",im));		
		
		clusterE_vsSize[im] = new TH2F(Form("clusterE_vsSize_AM%d",im),Form("%s, cluster size vs cluster energy, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		clusterE_vsSize[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterE_vsSize[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		
		clusterE_vsSize[im]->GetYaxis()->SetTitleOffset(1.4);
		clusterE_vsSize[im]->GetYaxis()->SetTitle(Form("AM%d cluster size, pixels",im));		
		
		clusterMaxE_vsSize[im] = new TH2F(Form("clusterMaxE_vsSize_AM%d",im),Form("%s, cluster size with maximum energy vs its size, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		clusterMaxE_vsSize[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterMaxE_vsSize[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		
		clusterMaxE_vsSize[im]->GetYaxis()->SetTitleOffset(1.4);
		clusterMaxE_vsSize[im]->GetYaxis()->SetTitle(Form("AM%d cluster size, pixels",im));		
		
		clusterSecE_vsSize[im] = new TH2F(Form("clusterSecE_vsSize_AM%d",im),Form("%s, cluster size second in energy vs its size, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		clusterSecE_vsSize[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterSecE_vsSize[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		
		clusterSecE_vsSize[im]->GetYaxis()->SetTitleOffset(1.4);
		clusterSecE_vsSize[im]->GetYaxis()->SetTitle(Form("AM%d cluster size, pixels",im));		
	
		clusterSizeMaxE_SecE[im] = new TH2F(Form("clusterSizeMaxE_SecE_AM%d",im),Form("%s, first vs second cluster size, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H,maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		clusterSizeMaxE_SecE[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterSizeMaxE_SecE[im]->GetXaxis()->SetTitle(Form("AM%d first cluster energy, pixels",im));		
		clusterSizeMaxE_SecE[im]->GetYaxis()->SetTitleOffset(1.4);
		clusterSizeMaxE_SecE[im]->GetYaxis()->SetTitle(Form("AM%d second cluster size, pixels",im));		
	
		clusterE_nTrigs[im] = new TH2F(Form("clusterE_nTrigs_AM%d",im),Form("%s, number of triggers in cluster vs cluster energy, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		clusterE_nTrigs[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterE_nTrigs[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		
		clusterE_nTrigs[im]->GetYaxis()->SetTitleOffset(1.4);
		clusterE_nTrigs[im]->GetYaxis()->SetTitle("Number of triggers in cluster");		
		
		nTrigsInEvent_maxEcluster[im] = new TH2F(Form("nTrigsInEvent_maxEcluster_AM%d",im),Form("%s, number of triggers per event vs cluster with maximum energy, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,maxNTrigs4H_2D-minNTrigs4H_2D,minNTrigs4H_2D,maxNTrigs4H_2D);
		nTrigsInEvent_maxEcluster[im]->GetXaxis()->SetTitleOffset(1.2);
		nTrigsInEvent_maxEcluster[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		
		nTrigsInEvent_maxEcluster[im]->GetYaxis()->SetTitleOffset(1.4);
		nTrigsInEvent_maxEcluster[im]->GetYaxis()->SetTitle(Form("Number of trigger in AM%d",im));		

		singleHitPixelSpec[im] = new TH1F(Form("singleHitPixelSpec_AM%d",im),Form("%s, single hit events, triggered pixel spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		singleHitPixelSpec[im]->GetXaxis()->SetTitleOffset(1.2);
		singleHitPixelSpec[im]->GetXaxis()->SetTitle(Form("AM%d pixel energy, keV",im));		

		allClustersSpec[im] = new TH1F(Form("allClustersSpec_AM%d",im),Form("%s, all cluster spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		allClustersSpec[im]->GetXaxis()->SetTitleOffset(1.2);
		allClustersSpec[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		

		singleClusterEventsSpec[im] = new TH1F(Form("singleClusterEventsSpec_AM%d",im),Form("%s, single hit events, cluster spectum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		singleClusterEventsSpec[im]->GetXaxis()->SetTitleOffset(1.2);
		singleClusterEventsSpec[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		

		summedClusterEventsSpec[im] = new TH1F(Form("summedClusterEventsSpec_AM%d",im),Form("%s, single hit events + summed Compton events, cluster spectum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		summedClusterEventsSpec[im]->GetXaxis()->SetTitleOffset(1.2);
		summedClusterEventsSpec[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		

		comptonHitsSpecClusters[im] = new TH1F(Form("comptonHitsSpecClusters_AM%d",im),Form("%s, unsummed Compton interaction cluster spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		comptonHitsSpecClusters[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonHitsSpecClusters[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		

		comptonHitsSpecClustersCorr[im] = new TH2F(Form("comptonHitsSpecClustersCorr_AM%d",im),Form("%s, unsummed Compton interactions, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,nBinsE,minE,maxE);
		comptonHitsSpecClustersCorr[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonHitsSpecClustersCorr[im]->GetXaxis()->SetTitle(Form("AM%d first cluster energy, keV",im));		
		comptonHitsSpecClustersCorr[im]->GetYaxis()->SetTitleOffset(1.4);
		comptonHitsSpecClustersCorr[im]->GetYaxis()->SetTitle(Form("AM%d second cluster energy, keV",im));		

		comptonSummedSpecClusters[im] = new TH1F(Form("comptonSummedSpecClusters_AM%d",im),Form("%s, two interaction summed Compton spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		comptonSummedSpecClusters[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonSummedSpecClusters[im]->GetXaxis()->SetTitle(Form("AM%d summed cluster energy, keV",im));		

		comptonPixelSpec[im] = new TH1F(Form("comptonPixelSpec_AM%d",im),Form("%s, triggered pixel spectrum of Compton events, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		comptonPixelSpec[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonPixelSpec[im]->GetXaxis()->SetTitle(Form("AM%d pixel energy, keV",im));		

		phiAngle[im] = new TH1F(Form("phiAngle_AM%d",im),Form("%s, cluster azimuthal angle, AM%d",spectrumName.Data(),im),360,0,360);
		phiAngle[im]->GetXaxis()->SetTitleOffset(1.2);
		phiAngle[im]->GetXaxis()->SetTitle(Form("#varphi AM%d, degrees",im));

		phiAnglePairs[im] = new TH1F(Form("phiAnglePairs_AM%d",im),Form("%s, cluster azimuthal angle from pairs, AM%d",spectrumName.Data(),im),360,0,360);
		phiAnglePairs[im]->GetXaxis()->SetTitleOffset(1.2);
		phiAnglePairs[im]->GetXaxis()->SetTitle(Form("#varphi AM%d, degrees",im));
		
		dPhiPlaneAM[im] = new TH2F(Form("dPhiPlaneAM_AM%d",im),"First - Second Compton scattering interaction distribution",nPixXY/2,0,nPixXY/2,nPixXY/2,0,nPixXY/2);
		dPhiPlaneAM[im]->GetXaxis()->SetTitle("Pixel coordinate");
		dPhiPlaneAM[im]->GetXaxis()->SetTitleOffset(1.2);
		dPhiPlaneAM[im]->GetYaxis()->SetTitle("Pixel coordinate");
		dPhiPlaneAM[im]->GetYaxis()->SetTitleOffset(1.4);
		
		firstComptonCluster[im] = new TH2F(Form("firstComptonCluster_AM%d",im),Form("First Compton scattering cluster positions, AM%d",im),nPixXY,0,nPixXY,nPixXY,0,nPixXY);
		firstComptonCluster[im]->GetXaxis()->SetTitle("Pixel coordinate");
		firstComptonCluster[im]->GetXaxis()->SetTitleOffset(1.2);
		firstComptonCluster[im]->GetYaxis()->SetTitle("Pixel coordinate");
		firstComptonCluster[im]->GetYaxis()->SetTitleOffset(1.4);
		
		secondComptonCluster[im] = new TH2F(Form("secondComptonCluster_AM%d",im),Form("Second Compton scattering cluster positions, AM%d",im),nPixXY,0,nPixXY,nPixXY,0,nPixXY);
		secondComptonCluster[im]->GetXaxis()->SetTitle("Pixel coordinate");
		secondComptonCluster[im]->GetXaxis()->SetTitleOffset(1.2);
		secondComptonCluster[im]->GetYaxis()->SetTitle("Pixel coordinate");
		secondComptonCluster[im]->GetYaxis()->SetTitleOffset(1.4);
	}
	
	singleHitPixelSpecCorr = new TH2F("singleHitPixelSpecCorr",Form("%s, energy spectra of single hit events in AMs",spectrumName.Data()),nBinsE,minE,maxE,nBinsE,minE,maxE);
	singleHitPixelSpecCorr->GetXaxis()->SetTitleOffset(1.2);
	singleHitPixelSpecCorr->GetXaxis()->SetTitle("AM0 pixel energy, keV");
	singleHitPixelSpecCorr->GetYaxis()->SetTitleOffset(1.4);
	singleHitPixelSpecCorr->GetYaxis()->SetTitle("AM1 pixel energy, keV");
	
	singleClusterSpecCorr = new TH2F("singleClusterSpecCorr",Form("%s, energy spectra of single cluster events in AMs",spectrumName.Data()),nBinsE,minE,maxE,nBinsE,minE,maxE);
	singleClusterSpecCorr->GetXaxis()->SetTitleOffset(1.2);
	singleClusterSpecCorr->GetXaxis()->SetTitle("AM0 cluster energy, keV");
	singleClusterSpecCorr->GetYaxis()->SetTitleOffset(1.4);
	singleClusterSpecCorr->GetYaxis()->SetTitle("AM1 cluster energy, keV");
	
	summedClusterSpecCorr = new TH2F("summedClusterSpecCorr",Form("%s, energy spectra of single cluster events and summed Compton events in AMs",spectrumName.Data()),nBinsE,minE,maxE,nBinsE,minE,maxE);
	summedClusterSpecCorr->GetXaxis()->SetTitleOffset(1.2);
	summedClusterSpecCorr->GetXaxis()->SetTitle("AM0 cluster energy, keV");
	summedClusterSpecCorr->GetYaxis()->SetTitleOffset(1.4);
	summedClusterSpecCorr->GetYaxis()->SetTitle("AM1 cluster energy, keV");
	
	comptonSummedSpecClustersCorr = new TH2F("comptonSummedSpecClustersCorr",Form("%s, two interaction summed Compton spectra in AMs",spectrumName.Data()),nBinsE,minE,maxE,nBinsE,minE,maxE);
	comptonSummedSpecClustersCorr->GetXaxis()->SetTitleOffset(1.2);
	comptonSummedSpecClustersCorr->GetXaxis()->SetTitle("AM0 summed cluster energy, keV");
	comptonSummedSpecClustersCorr->GetYaxis()->SetTitleOffset(1.4);
	comptonSummedSpecClustersCorr->GetYaxis()->SetTitle("AM1 summed cluster energy, keV");

	dPhiAngle = new TH1F("dPhiAngle",Form("%s, #Delta#varphi",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiAngle->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle_summed = new TH1F("dPhiAngle_summed",Form("%s, #Delta#varphi summed",spectrumName.Data()),nBins_dPhi,0,180);
	dPhiAngle_summed->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_summed->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle_Ewin = new TH1F("dPhiAngle_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV",spectrumName.Data(),EwindowFor4PhiAnalyis_min,EwindowFor4PhiAnalyis_max),nBins_dPhi,-180,180);
	dPhiAngle_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle_2ClusterE = new TH2F("dPhiAngle_2ClusterE",Form("%s, #Delta#varphi vs summed energy of two Compton clusters",spectrumName.Data()),nBinsE,minE,maxE,nBins_dPhi,-180,180);
	dPhiAngle_2ClusterE->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_2ClusterE->GetXaxis()->SetTitle("Total summed energy, keV");
	dPhiAngle_2ClusterE->GetYaxis()->SetTitleOffset(1.4);
	dPhiAngle_2ClusterE->GetYaxis()->SetTitle("#Delta#varphi, degrees");		
	dPhiAngle_2ClusterE->RebinX(10);
	
	dPhiAngle_totalE = new TH2F("dPhiAngle_totalE",Form("%s, #Delta#varphi vs total energy of Compton event",spectrumName.Data()),nBinsE,minE,maxE,nBins_dPhi,-180,180);
	dPhiAngle_totalE->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_totalE->GetXaxis()->SetTitle("Total energy, keV");
	dPhiAngle_totalE->GetYaxis()->SetTitleOffset(1.4);
	dPhiAngle_totalE->GetYaxis()->SetTitle("#Delta#varphi, degrees");		
	dPhiAngle_totalE->RebinX(10);
	
	phiAngleCorr = new TH2F("phiAngleCorr",Form("%s, #Delta#varphi",spectrumName.Data()),36,0,360,36,0,360);
	phiAngleCorr->GetXaxis()->SetTitleOffset(1.2);
	phiAngleCorr->GetXaxis()->SetTitle("#varphi AM0, degrees");	
	phiAngleCorr->GetYaxis()->SetTitleOffset(1.2);
	phiAngleCorr->GetYaxis()->SetTitle("#varphi AM1, degrees");	

	theta0_theta1_FromFirstClusterE = new TH2F("theta0_theta1_FromFirstClusterE",Form("%s, #theta_{0} - #theta_{1} calculated from first interaction energy",spectrumName.Data()),nBinsTheta4H,0,180,nBinsTheta4H,0,180);
	theta0_theta1_FromFirstClusterE->GetXaxis()->SetTitleOffset(1.2);
	theta0_theta1_FromFirstClusterE->GetXaxis()->SetTitle("#theta in AM0, degrees");
	theta0_theta1_FromFirstClusterE->GetYaxis()->SetTitleOffset(1.2);
	theta0_theta1_FromFirstClusterE->GetYaxis()->SetTitle("#theta in AM1, degrees");
	
	localTimer.Start();

	matrix_trigs.clear();
	matrix_trigs.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_flags.clear();
	matrix_flags.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_E.clear();
	matrix_E.resize(nPixXY, vector<Float_t>(nPixXY));
	
	for (Int_t ie = 0; ie < nEvents2Analyse; ie++)
	//for (Int_t ie = 0; ie < 100; ie++)
	{
		chain_events->GetEntry(ie);
		//if (ie == 3) continue;
		if (ie%prntOutF == 0)
		{
			num_str = Form("%d",ie);
			for (Int_t is = num_str.Length(); is > 1; is--)
			{
				if ((is-num_str.Length()+2)%4 == 0) num_str.Insert(is-1,",",1);
			}
			cout << Form("%s, Nev = %s, RealTime: %.3fs, CPUTime: %.3fs",num_str.Data(),nev_str.Data(),localTimer.RealTime(), localTimer.CpuTime()) << endl;
			totalTimeElapced += localTimer.RealTime();
			localTimer.ResetRealTime();
			localTimer.ResetCpuTime();
			localTimer.Start();
		}
		//cout << Form("Event %lld, nTrigPixels_eev[0] = %d, nTrigPixels_eev[1] = %d, AM_flag_eev[0] = %d, AM_flag_eev[1] = %d",event_eev,nTrigPixels_eev[0],nTrigPixels_eev[1],AM_flag_eev[0],AM_flag_eev[1]) << endl;
		if (!analyseNextEvent(event_eev))
		{
			cerr << "ERROR analysing event " << ie << ". Exiting." << endl;
			return 0;
		}
		//cout << "------------------------" << endl;
	}
	
	UInt_t screen_width = screen_width_def;
	UInt_t screen_height = screen_height_def;
	Float_t size_factor = 0.95;
	Float_t size_ratio = 1.333;
	const Int_t ix = 3;
	const Int_t iy = 3;
	TCanvas *c0 = new TCanvas("c0","",10,10,Int_t(screen_height*size_factor*size_ratio),Int_t(screen_height*size_factor));
	Int_t curr_Pad = 0;
	c0->Divide(ix,iy);
	for (curr_Pad = 0; curr_Pad <= ix*iy; curr_Pad++)
	{
		c0->cd(curr_Pad);
		c0->GetPad(curr_Pad)->SetGridx();
		c0->GetPad(curr_Pad)->SetGridy();
		c0->GetPad(curr_Pad)->SetLeftMargin(0.12);
		c0->GetPad(curr_Pad)->SetRightMargin(0.12);
		c0->GetPad(curr_Pad)->SetTopMargin(0.11);
		c0->GetPad(curr_Pad)->SetBottomMargin(0.13);
	}
	TLatex *txt = new TLatex();
	txt->SetTextSize(.035);
	txt->SetTextColor(2);
	txt->SetNDC();
	TLine *line = new TLine();
	line->SetLineColor(1);
	line->SetLineWidth(2);
	TLine *line2 = new TLine();
	line2->SetLineColor(3);
	line2->SetLineWidth(2);

	TCanvas *c1 = new TCanvas("c1","",10,10,Int_t(screen_height*size_factor*size_ratio),Int_t(screen_height*size_factor));
	curr_Pad = 0;
	c1->cd(curr_Pad);
	c1->GetPad(curr_Pad)->SetGridx();
	c1->GetPad(curr_Pad)->SetGridy();
	c1->GetPad(curr_Pad)->SetLeftMargin(0.12);
	c1->GetPad(curr_Pad)->SetRightMargin(0.12);
	c1->GetPad(curr_Pad)->SetTopMargin(0.11);
	c1->GetPad(curr_Pad)->SetBottomMargin(0.13);
	
	size_ratio = 1;
	TCanvas *c2 = new TCanvas("c2","",10,10,Int_t(screen_height*size_factor*size_ratio),Int_t(screen_height*size_factor));
	curr_Pad = 0;
	c2->Divide(2,1);
	for (curr_Pad = 0; curr_Pad <= 2; curr_Pad++)
	{
		c2->cd(curr_Pad);
		c2->GetPad(curr_Pad)->SetGridx();
		c2->GetPad(curr_Pad)->SetGridy();
		c2->GetPad(curr_Pad)->SetLeftMargin(0.12);
		c2->GetPad(curr_Pad)->SetRightMargin(0.12);
		c2->GetPad(curr_Pad)->SetTopMargin(0.11);
		c2->GetPad(curr_Pad)->SetBottomMargin(0.13);
	}
	
	Float_t X_text = 0.21;
	//Float_t X_text = 0.68;
	Float_t Y_text_top = 0.84;
	Float_t Y_text_step = 0.05;

	Float_t Peak_FWHM_raw = 0, Peak_FWHM_raw_leftEdge = 0, Peak_FWHM_raw_rightEdge = 0, locMaxPeak = 0, locMaxPeakHeight = 0;

	TH1F **nTrigsInEvent_norm = new TH1F*[nAMs];
	TH1F **clusterSize_norm = new TH1F*[nAMs];
	for (Int_t im = 0; im < nAMs; im++)
	{
		c2->cd(0);
		pixelPattern[im]->Draw("text");
		line->DrawLine(pixelPattern[im]->GetXaxis()->GetXmin(),pixelPattern[im]->GetYaxis()->GetXmax()/2,pixelPattern[im]->GetXaxis()->GetXmax(),pixelPattern[im]->GetYaxis()->GetXmax()/2);
		line->DrawLine(pixelPattern[im]->GetXaxis()->GetXmax()/2,pixelPattern[im]->GetYaxis()->GetXmin(),pixelPattern[im]->GetXaxis()->GetXmax()/2,pixelPattern[im]->GetYaxis()->GetXmax());
		c2->SaveAs(Form("%s/pixelPattern.gif",outputPathAM[im].Data()));

		comptonHitsSpecClustersCorr[im]->Draw("colz");
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(comptonHitsSpecClustersCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(comptonHitsSpecClustersCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
	
		dPhiPlaneAM[im]->Draw("colz");
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dPhiPlaneAM[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dPhiPlaneAM[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
	
		clusterSizeMaxE_SecE[im]->Draw("colz");
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(clusterSizeMaxE_SecE[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(clusterSizeMaxE_SecE[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
		
		c1->cd(0);
		singleHitPixelSpec[im]->Draw();
		Y_text_top = 0.82;
		getHistoPeak(singleHitPixelSpec[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(singleHitPixelSpec[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/singleHitPixelSpec[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(singleHitPixelSpec[im]->GetName()).Data()));
		
		singleClusterEventsSpec[im]->Draw();
		Y_text_top = 0.82;
		getHistoPeak(singleClusterEventsSpec[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(singleClusterEventsSpec[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/singleClusterEventsSpec[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(singleClusterEventsSpec[im]->GetName()).Data()));
		
		summedClusterEventsSpec[im]->Draw();
		Y_text_top = 0.82;
		getHistoPeak(summedClusterEventsSpec[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(summedClusterEventsSpec[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/summedClusterEventsSpec[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(summedClusterEventsSpec[im]->GetName()).Data()));
		
		allClustersSpec[im]->Draw();
		Y_text_top = 0.82;
		getHistoPeak(allClustersSpec[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(allClustersSpec[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/allClustersSpec[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(allClustersSpec[im]->GetName()).Data()));

		clusterE_nTrigs[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(clusterE_nTrigs[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(clusterE_nTrigs[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		comptonHitsSpecClusters[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonHitsSpecClusters[im]->GetName()).Data()));
				
		nTrigsInEvent[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(nTrigsInEvent[im]->GetName()).Data()));
		nTrigsInEvent_norm[im] = (TH1F *) nTrigsInEvent[im]->Clone();
		nTrigsInEvent_norm[im]->SetName(Form("nTrigsInEvent_norm_AM%d",im));
		nTrigsInEvent_norm[im]->SetNormFactor(1);
		nTrigsInEvent_norm[im]->GetYaxis()->SetTitle("Normalised counts");
		nTrigsInEvent_norm[im]->GetYaxis()->SetTitleOffset(1.3);
		nTrigsInEvent_norm[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(nTrigsInEvent_norm[im]->GetName()).Data()));
		
		nTrigsInEvent_maxE[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(nTrigsInEvent_maxE[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(nTrigsInEvent_maxE[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		nTrigs_clusterSize[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(nTrigs_clusterSize[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(nTrigs_clusterSize[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		clusterE_vsSize[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(clusterE_vsSize[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(clusterE_vsSize[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
	
		clusterMaxE_vsSize[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(clusterMaxE_vsSize[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(clusterMaxE_vsSize[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
	
		clusterSecE_vsSize[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(clusterSecE_vsSize[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(clusterSecE_vsSize[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
	
		nTrigsInEvent_maxEcluster[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(nTrigsInEvent_maxEcluster[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(nTrigsInEvent_maxEcluster[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		c1->cd(0);
		comptonSummedSpecClusters[im]->Draw();
		Y_text_top = 0.82;
		getHistoPeak(comptonSummedSpecClusters[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(comptonSummedSpecClusters[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpecClusters[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSummedSpecClusters[im]->GetName()).Data()));
		
		comptonPixelSpec[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonPixelSpec[im]->GetName()).Data()));
		
		phiAngle[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(phiAngle[im]->GetName()).Data()));
		
		phiAnglePairs[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(phiAnglePairs[im]->GetName()).Data()));
				
		firstComptonCluster[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(firstComptonCluster[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(firstComptonCluster[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		secondComptonCluster[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(secondComptonCluster[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(secondComptonCluster[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		thetaFromFirstClusterE_dPhi[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_dPhi[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_dPhi[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		thetaFromFirstClusterE[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(thetaFromFirstClusterE[im]->GetName()).Data()));
		
		clusterSize[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(clusterSize[im]->GetName()).Data()));
		clusterSize_norm[im] = (TH1F *) clusterSize[im]->Clone();
		clusterSize_norm[im]->SetName(Form("clusterSize_norm_AM%d",im));
		clusterSize_norm[im]->SetNormFactor(1);
		clusterSize_norm[im]->GetYaxis()->SetTitle("Normalised counts");
		clusterSize_norm[im]->GetYaxis()->SetTitleOffset(1.3);
		clusterSize_norm[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(clusterSize_norm[im]->GetName()).Data()));
	}

	c1->cd(0);
	dPhiAngle->SetMinimum(0);
	dPhiAngle->Draw();
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle->GetName()).Data()));

	dPhiAngle_summed->SetMinimum(0);
	dPhiAngle_summed->Draw();
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle_summed->GetName()).Data()));

	dPhiAngle_Ewin->SetMinimum(0);
	dPhiAngle_Ewin->Draw();
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle_Ewin->GetName()).Data()));	
	
	c2->cd(0);
	singleHitPixelSpecCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(singleHitPixelSpecCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(singleHitPixelSpecCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
	
	singleClusterSpecCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(singleClusterSpecCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(singleClusterSpecCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
	
	comptonSummedSpecClustersCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(comptonSummedSpecClustersCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(comptonSummedSpecClustersCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
	
	summedClusterSpecCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(summedClusterSpecCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(summedClusterSpecCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
	
	phiAngleCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(phiAngleCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(phiAngleCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
		
	dPhiAngle_2ClusterE->Draw("colz");
	dPhiAngle_2ClusterE->GetXaxis()->SetRangeUser(300,dPhiAngle_2ClusterE->GetXaxis()->GetXmax());
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle_2ClusterE->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(dPhiAngle_2ClusterE->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
	
	dPhiAngle_totalE->Draw("colz");
	dPhiAngle_totalE->GetXaxis()->SetRangeUser(300,dPhiAngle_totalE->GetXaxis()->GetXmax());
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle_totalE->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(dPhiAngle_totalE->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);

	theta0_theta1_FromFirstClusterE->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(theta0_theta1_FromFirstClusterE->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(theta0_theta1_FromFirstClusterE->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
	
	delete c0;
	delete c1;
	delete c2;
	
	if (totalTimeElapced < 60) cout << Form("Total running time is %.2f sec",totalTimeElapced) << endl;
	if (totalTimeElapced >= 60 && totalTimeElapced < 3600) cout << Form("Total running time is %dm:%.0f sec",Int_t(totalTimeElapced/60),totalTimeElapced - Int_t(totalTimeElapced/60)*60) << endl;
	if (totalTimeElapced >= 3600) cout << Form("Total running time is %dh:%dm:%.0fs",Int_t(totalTimeElapced/3600),Int_t((totalTimeElapced - Int_t(totalTimeElapced/3600)*3600)/60),
		totalTimeElapced - Int_t(totalTimeElapced/3600)*3600 - Int_t((totalTimeElapced - Int_t(totalTimeElapced/3600)*3600)/60)*60) << endl;
	
	hfile->cd();
	for (Int_t im = 0; im < nAMs; im++)
	{
		singleHitPixelSpec[im]->Write();
		allClustersSpec[im]->Write();
		singleClusterEventsSpec[im]->Write();
		summedClusterEventsSpec[im]->Write();
		nTrigsInEvent[im]->Write();
		nTrigsInEvent_maxE[im]->Write();
		nTrigs_clusterSize[im]->Write();
		nTrigsInEvent_maxEcluster[im]->Write();
		comptonSummedSpecClusters[im]->Write();
		comptonHitsSpecClusters[im]->Write();
		comptonHitsSpecClustersCorr[im]->Write();
		comptonPixelSpec[im]->Write();
		phiAngle[im]->Write();
		phiAnglePairs[im]->Write();
		dPhiPlaneAM[im]->Write();
		firstComptonCluster[im]->Write();
		secondComptonCluster[im]->Write();
		clusterE_vsSize[im]->Write();
		clusterMaxE_vsSize[im]->Write();
		thetaFromFirstClusterE[im]->Write();
		thetaFromFirstClusterE_dPhi[im]->Write();
		nTrigsInEvent_norm[im]->Write();
		clusterSize[im]->Write();
		clusterE_nTrigs[im]->Write();
		clusterSecE_vsSize[im]->Write();
		clusterSizeMaxE_SecE[im]->Write();
	}
	singleHitPixelSpecCorr->Write();
	singleClusterSpecCorr->Write();
	summedClusterSpecCorr->Write();
	dPhiAngle->Write();
	dPhiAngle_summed->Write();
	dPhiAngle_Ewin->Write();
	dPhiAngle_2ClusterE->Write();
	dPhiAngle_totalE->Write();
	phiAngleCorr->Write();
	comptonSummedSpecClustersCorr->Write();
	theta0_theta1_FromFirstClusterE->Write();
	hfile->Close();
}

Bool_t analyseNextEvent(const Int_t ievent)
{
	for (Int_t im = 0; im < nAMs; im++)
	{
		singleHitsLength[im] = 0;
		comptonChainsLength[im] = 0;
		for (Int_t b = 0; b < buffSizeSingleHits; b++)
		{
			for (Int_t b2 = 0; b2 < nVarsSingleHits; b2++) singleHits[im][b][b2] = 0;
		}
		for (Int_t b = 0; b < buffSizeComptonChains; b++)
		{
			comptonLegsLength[im][b] = 0;
			for (Int_t b2 = 0; b2 < nLegsComptonChains; b2++)
			{
				for (Int_t b3 = 0; b3 < nVarsComptonChains; b3++) comptonChains[im][b][b2][b3] = 0;
			}
		}
	}
	buffClusterX.clear();
	buffClusterY.clear();
	buffClusterE.clear();
	buffClusterArea.clear();
	buffClusterFlag.clear();
	buffClusterTrigs.clear();
	for (Int_t im = 0; im < nAMs; im++)	if (!analyseAM(ievent, im)) return kFALSE;

	// Just for self-check
	for (Int_t im = 0; im < nAMs; im++)
	{
		for (Int_t ic = 0; ic < comptonChainsLength[im]; ic++)
			if (comptonLegsLength[im][ic] != nTrigPixels[im])
					cout << Form("ERROR!!! event = %d, comptonLegsLength[%d][%d] = %d, nTrigPixels[%d] = %d",ievent,im,ic,comptonLegsLength[im][ic],im,nTrigPixels[im]) << endl;
	}
	
	// So far undeveloped
	/*
	Float_t summedE0 = -999;
	Float_t summedE1 = -999;
	for (Int_t im = 0; im < nAMs; im++)
	{
		if (singleHitsLength[im] > 0)
		{	
			singleHitSpec[im]->Fill(singleHits[im][0][0]);
		}
		
		for (Int_t ic = 0; ic < comptonChainsLength[im]; ic++)
		{
			if (comptonLegsLength[im][ic] >= 2)
			{
				Float_t xc2 = comptonChains[im][ic][0][1]-1;
				Float_t yc2 = comptonChains[im][ic][0][2]-1;
				Float_t xc1 = comptonChains[im][ic][1][1]-1;
				Float_t yc1 = comptonChains[im][ic][1][2]-1;
				Float_t dist = TMath::Hypot(xc1-xc2,yc1-yc2);
				Float_t phiAng = getPhiAngleDeg(xc2-xc1,yc2-yc1);
				phiAngle[im]->Fill(phiAng);
			}
			Float_t summedE = 0;
			for (Int_t p = 0; p < comptonLegsLength[im][ic]&&ic < buffSizeComptonChains&&p < nLegsComptonChains; p++) summedE += comptonChains[im][ic][p][0];
			if (summedE > 1) comptonSummedSpec[im]->Fill(summedE);
			if (summedE > 1 && im == 1) summedE0 = summedE;
			if (summedE > 1 && im == 0) summedE1 = summedE;
		}
	}
	if (summedE0 > 1 && summedE1 > 1) comptonSummedSpecCorr->Fill(summedE0,summedE1);
	*/
	
	if (singleHitsLength[0] > 0 && singleHitsLength[1] > 0) singleHitPixelSpecCorr->Fill(singleHits[0][0][0],singleHits[1][0][0]);

	// For now, make pairs just from two first clusters in cluster array sorted by energy
	// First Compton interaction - second cluster (second in energy), second interaction - first cluster (maximum energy), assuming forward Compton scattering

	//cout << "=============================" << endl;

	for (Int_t im = 0; im < nAMs; im++)
	{
		if (buffClusterX[im].size() >= 2)
		{
			Float_t xc2 = buffClusterX[im][0]-1;
			Float_t yc2 = buffClusterY[im][0]-1;
			Float_t xc1 = buffClusterX[im][1]-1;
			Float_t yc1 = buffClusterY[im][1]-1;
			
			//cout << im << " xc2: " << xc2 << " xc1: " << xc1 << " yc2: " << yc2 << " yc1: " << yc1 << endl;
			
			Float_t dist = TMath::Hypot(xc1-xc2,yc1-yc2);
			Float_t phiAng = getPhiAngleDeg(xc2-xc1,yc2-yc1);
			phiAngle[im]->Fill(phiAng);
		}
		if (nTrigPixels[im] > 1)
		{
			for (Int_t p = 0; p < nLegsComptonChains; p++)
				if (comptonChains[im][0][p][0] > 1) comptonPixelSpec[im]->Fill(comptonChains[im][0][p][0]);
		}
		if (nTrigPixels[im] == 1) singleHitPixelSpec[im]->Fill(singleHits[im][0][0]);
	}

	if (buffClusterX[0].size() == 1 && buffClusterX[1].size() == 1)
	{
		singleClusterSpecCorr->Fill(buffClusterE[0][0],buffClusterE[1][0]);
		summedClusterSpecCorr->Fill(buffClusterE[0][0],buffClusterE[1][0]);
	}

	if (buffClusterX[0].size() >= 2 && buffClusterX[1].size() >= 2)
	{
		Float_t xcentre = 6, ycentre = 17;
		Float_t phiAng0 = -9999;
		Float_t xc2 = buffClusterX[0][0]-1;
		Float_t yc2 = buffClusterY[0][0]-1;
		Float_t xc1 = buffClusterX[0][1]-1;
		Float_t yc1 = buffClusterY[0][1]-1;
		
		//cout << "Pairs: 0 " << "xc2: " << xc2 << " xc1: " << xc1 << " yc2: " << yc2 << " yc1: " << yc1 << endl;
		
		Float_t dist0 = TMath::Hypot(xc1-xc2,yc1-yc2);
		Bool_t phi0_valid = kFALSE;
		if (fabs(xc2-xc1) < maxDistanceBetweenClusters4ComptonPair && dist0 > 0.1 && buffClusterArea[0][1] <= 4)
		{
			phiAng0 = getPhiAngleDeg(xc2-xc1,yc2-yc1);
			phiAnglePairs[0]->Fill(phiAng0);
			firstComptonCluster[0]->Fill(xc1+1,yc1+1);
			secondComptonCluster[0]->Fill(xc2+1,yc2+1);
			dPhiPlaneAM[0]->Fill(xcentre,ycentre-11);
			dPhiPlaneAM[0]->Fill(xc2+xcentre-xc1,yc2+ycentre-yc1-11);
			phi0_valid = kTRUE;
		}
		
		Float_t phiAng1 = -9999;
		Float_t phiAng1mirrored = -9999;
		Float_t phiAng1mirrored_shifted = -9999;
		xc2 = buffClusterX[1][0]-1;
		yc2 = buffClusterY[1][0]-1;
		xc1 = buffClusterX[1][1]-1;
		yc1 = buffClusterY[1][1]-1;
		
		//cout << "Pairs: 1 " << "xc2: " << xc2 << " xc1: " << xc1 << " yc2: " << yc2 << " yc1: " << yc1 << endl;

		Float_t dist1 = TMath::Hypot(xc1-xc2,yc1-yc2);
		Bool_t phi1_valid = kFALSE;
		if (fabs(xc2-xc1) < maxDistanceBetweenClusters4ComptonPair && dist0 > 0.1 && buffClusterArea[0][1] <= 4)
		{
			phiAng1 = getPhiAngleDeg(xc2-xc1,yc2-yc1);
			phiAng1mirrored = getPhiAngleDeg(xc1-xc2,yc2-yc1);
			phiAng1mirrored_shifted = phiAng1mirrored + relativePhiAngle;
			if (phiAng1mirrored_shifted > 360) phiAng1mirrored_shifted -= 360;
			phiAnglePairs[1]->Fill(phiAng1mirrored_shifted);
			firstComptonCluster[1]->Fill(xc1+1,yc1+1);
			secondComptonCluster[1]->Fill(xc2+1,yc2+1);
			dPhiPlaneAM[1]->Fill(xcentre,ycentre-11);
			dPhiPlaneAM[1]->Fill(xc2+xcentre-xc1,yc2+ycentre-yc1-11);
			phi1_valid = kTRUE;
		}
		if (phi0_valid && phi1_valid)
		{
			Float_t initialE1 = getScatteredEnergyFromAngle(Na22Energy, angleOfSecondHead);
			Float_t total2E0 = buffClusterE[0][0]+buffClusterE[0][1];
			Float_t total2E1 = buffClusterE[1][0]+buffClusterE[1][1];
			Float_t theta0 = getThetaFromEnergy(Na22Energy, buffClusterE[0][1]);
			Float_t theta1 = getThetaFromEnergy(initialE1, buffClusterE[1][1]);
			
			phiAngleCorr->Fill(phiAng0,phiAng1mirrored_shifted);
			Float_t dPh = phiAng0-phiAng1mirrored_shifted;
			if (dPh > 180) dPh = 360 - dPh;
			if (dPh < -180) dPh = -360 - dPh;
			if (theta0 > ThetaWindowFor4PhiAnalyis_min && theta0 < ThetaWindowFor4PhiAnalyis_max && theta1 > ThetaWindowFor4PhiAnalyis_min && theta1 < ThetaWindowFor4PhiAnalyis_max) 
			{
				dPhiAngle->Fill(dPh);
				if(dPh>=0) dPhiAngle_summed->Fill(dPh);
				else  dPhiAngle_summed->Fill(-dPh);
				if (total2E0 >= EwindowFor4PhiAnalyis_min && total2E0 <= EwindowFor4PhiAnalyis_max && total2E1 >= EwindowFor4PhiAnalyis_min && total2E1 <= EwindowFor4PhiAnalyis_max) dPhiAngle_Ewin->Fill(dPh);
			}
			dPhiAngle_2ClusterE->Fill(total2E0,dPh);
			dPhiAngle_2ClusterE->Fill(total2E1,dPh);
			thetaFromFirstClusterE[0]->Fill(theta0);
			
			thetaFromFirstClusterE[1]->Fill(theta1);
			thetaFromFirstClusterE_dPhi[0]->Fill(theta0,dPh);
			thetaFromFirstClusterE_dPhi[1]->Fill(theta1,dPh);
			theta0_theta1_FromFirstClusterE->Fill(theta0,theta1);
			Float_t totalE = 0;
			for (Int_t p = 0; p < buffClusterX[0].size(); p++) totalE += buffClusterE[0][p];
			dPhiAngle_totalE->Fill(totalE,dPh);
			totalE = 0;
			for (Int_t p = 0; p < buffClusterX[1].size(); p++) totalE += buffClusterE[1][p];
			dPhiAngle_totalE->Fill(totalE,dPh);
		}
		comptonSummedSpecClustersCorr->Fill(buffClusterE[0][0]+buffClusterE[0][1],buffClusterE[1][0]+buffClusterE[1][1]);
		summedClusterSpecCorr->Fill(buffClusterE[0][0]+buffClusterE[0][1],buffClusterE[1][0]+buffClusterE[1][1]);
	}
	for (Int_t im = 0; im < nAMs; im++)
	{
		if (buffClusterX[im].size() >= 2)
		{
			for (Int_t p = 0; p < 2; p++) comptonHitsSpecClusters[im]->Fill(buffClusterE[im][p]);
			comptonSummedSpecClusters[im]->Fill(buffClusterE[im][0] + buffClusterE[im][1]);
			summedClusterEventsSpec[im]->Fill(buffClusterE[im][0] + buffClusterE[im][1]);
			comptonHitsSpecClustersCorr[im]->Fill(buffClusterE[im][0],buffClusterE[im][1]);
			clusterMaxE_vsSize[im]->Fill(buffClusterE[im][0],buffClusterArea[im][0]);
			clusterSecE_vsSize[im]->Fill(buffClusterE[im][1],buffClusterArea[im][1]);
			clusterSizeMaxE_SecE[im]->Fill(buffClusterArea[im][0],buffClusterArea[im][1]);
		}
		if (buffClusterX[im].size() == 1)
		{
			singleClusterEventsSpec[im]->Fill(buffClusterE[im][0]);
			summedClusterEventsSpec[im]->Fill(buffClusterE[im][0]);
		}
		for (Int_t p = 0; p < buffClusterX[im].size(); p++)
		{
			clusterE_vsSize[im]->Fill(buffClusterE[im][p],buffClusterArea[im][p]);
			clusterSize[im]->Fill(buffClusterArea[im][p]);
			allClustersSpec[im]->Fill(buffClusterE[im][p]);
			nTrigs_clusterSize[im]->Fill(buffClusterArea[im][p],buffClusterTrigs[im][p]);
			clusterE_nTrigs[im]->Fill(buffClusterE[im][p],buffClusterTrigs[im][p]);
		}
	}

	// So far undeveloped
	/*
	if (nTrigPixels[0] == 2 && nTrigPixels[1] == 2)
	{
		if (summedE0 > 1 && summedE1 > 1) comptonSummedSpecCorr->Fill(summedE0,summedE1);
		
		Float_t xcentre = 6, ycentre = 17;
		Float_t phiAng0 = -9999;
		Float_t xc2 = comptonChains[0][0][0][1]-1;
		Float_t yc2 = comptonChains[0][0][0][2]-1;
		Float_t xc1 = comptonChains[0][0][1][1]-1;
		Float_t yc1 = comptonChains[0][0][1][2]-1;
		Float_t dist0 = TMath::Hypot(xc1-xc2,yc1-yc2);
		if (fabs(xc2-xc1) < 7 && dist0 > 0.1)
		{
			phiAng0 = getPhiAngleDeg(xc2-xc1,yc2-yc1);
			phiAnglePairs[0]->Fill(phiAng0);
			dPhiPlaneAM[0]->Fill(xcentre,ycentre-11);
			dPhiPlaneAM[0]->Fill(xc2+xcentre-xc1,yc2+ycentre-yc1-11);
			firstComptonHit[0]->Fill(xc1,yc1);
			secondComptonHit[0]->Fill(xc2,yc1);
		}
		
		Float_t phiAng1 = -9999;
		Float_t phiAng1mirrored = -9999;
		Float_t phiAng1mirrored_shifted = -9999;
		xc2 = comptonChains[1][0][0][1]-1;
		yc2 = comptonChains[1][0][0][2]-1;
		xc1 = comptonChains[1][0][1][1]-1;
		yc1 = comptonChains[1][0][1][2]-1;
		Float_t dist1 = TMath::Hypot(xc1-xc2,yc1-yc2);
		if (fabs(xc2-xc1) < 5 && dist1 > 0.1)
		{
			phiAng1 = getPhiAngleDeg(xc2-xc1,yc2-yc1);
			phiAng1mirrored = getPhiAngleDeg(xc1-xc2,yc2-yc1);
			phiAng1mirrored_shifted = phiAng1mirrored + relativePhiAngle;
			if (phiAng1mirrored_shifted > 360) phiAng1mirrored_shifted -= 360;
			phiAnglePairs[1]->Fill(phiAng1);
			dPhiPlaneAM[1]->Fill(xcentre,ycentre-11);
			dPhiPlaneAM[1]->Fill(xc2+xcentre-xc1,yc2+ycentre-yc1-11);
			firstComptonHit[1]->Fill(xc1,yc1);
			secondComptonHit[1]->Fill(xc2,yc1);
		}
		if (dist0 > 2 && dist1 > 2)
		{
			phiAngleCorr->Fill(phiAng0,phiAng1);
			Float_t dPh = fabs(phiAng0-phiAng1mirrored_shifted);
			if (dPh > 180) dPh = 360 - dPh;
			if (dPh > 0.01 && dPh < 179.99) dPhiAngle->Fill(dPh);
		}
	}
	*/
	return kTRUE;
}

Bool_t analyseAM(const Int_t evt, const Int_t AM2do)
{
	Int_t cntComptonChain = 0;
	Bool_t comptonFound = kFALSE;
	Int_t centrePixelX = -1;
	Int_t centrePixelY = -1;
	buffClusterXloc.clear();
	buffClusterYloc.clear();
	buffClusterEloc.clear();
	buffClusterArealoc.clear();
	buffClusterFlagloc.clear();
	buffClusterTrigsloc.clear();
	buffClusterXloc.resize(0,0);
	buffClusterYloc.resize(0,0);
	buffClusterEloc.resize(0,0);
	buffClusterArealoc.resize(0,0);
	buffClusterFlagloc.resize(0,0);
	buffClusterTrigsloc.resize(0,0);
	matrix_trigs.clear();
	matrix_trigs.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_flags.clear();
	matrix_flags.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_E.clear();
	matrix_E.resize(nPixXY, vector<Float_t>(nPixXY));
	for (Int_t ii = 0; ii < Nev[AM2do] && event[AM2do] != evt; ii++)
	{
		cntEntryTree[AM2do]++;
		chain1p[AM2do]->GetEntry(cntEntryTree[AM2do]);
		if (event[AM2do] == evt) break;
	}
	
	while (event[AM2do] == evt && cntEntryTree[AM2do] <= Nev[AM2do])
	{
		if (singleHitsLength[AM2do] == 0 && pos[AM2do] == 0)
		{
			nTrigsInEvent[AM2do]->Fill(nTrigPixels[AM2do]);
			nTrigsInEvent_maxE[AM2do]->Fill(E[AM2do],nTrigPixels[AM2do]);
		}
		if (pixel[AM2do] >= firstPixel && pixel[AM2do] <= lastPixel)
		{
			Int_t xx1 = -1, yy1 = -1;
			if (getPixel2Dcoordinates(pixel[AM2do],AM2do,xx1,yy1))
			{
				if (E[AM2do] > minPixelEnergyThr4ClusterReconstruction)
				{
				if ((useOnlyTriggeredPixelsInCluster && triggerFlag[AM2do] == 1) || !useOnlyTriggeredPixelsInCluster)
					{
						matrix_E[xx1-1][yy1-1] = E[AM2do];
						matrix_trigs[xx1-1][yy1-1] = triggerFlag[AM2do];
						matrix_flags[xx1-1][yy1-1] = newNeighbourFlag;
					}
				}
				if (triggerFlag[AM2do] == 1 && pos[AM2do] == 0)
				{
					centrePixelX = xx1;
					centrePixelY = yy1;
				}
			}
			if (nTrigPixels[AM2do] == 1)
			{
				if (triggerFlag[AM2do] == 1)
				{
					singleHits[AM2do][singleHitsLength[AM2do]][0] = E[AM2do];
					singleHits[AM2do][singleHitsLength[AM2do]][1] = xx1;
					singleHits[AM2do][singleHitsLength[AM2do]][2] = yy1;
					singleHitsLength[AM2do]++;
				}
			}
			if (nTrigPixels[AM2do] > 1) 
			{
				comptonFound = kTRUE;
				if (triggerFlag[AM2do] == 1)
				{
					if (comptonLegsLength[AM2do][comptonChainsLength[AM2do]] < nLegsComptonChains)
					{
						comptonChains[AM2do][comptonChainsLength[AM2do]][comptonLegsLength[AM2do][comptonChainsLength[AM2do]]][0] = E[AM2do];
						comptonChains[AM2do][comptonChainsLength[AM2do]][comptonLegsLength[AM2do][comptonChainsLength[AM2do]]][1] = xx1;
						comptonChains[AM2do][comptonChainsLength[AM2do]][comptonLegsLength[AM2do][comptonChainsLength[AM2do]]][2] = yy1;
						comptonChains[AM2do][comptonChainsLength[AM2do]][comptonLegsLength[AM2do][comptonChainsLength[AM2do]]][3] = pixel[AM2do];  // temporary
					}
					comptonLegsLength[AM2do][comptonChainsLength[AM2do]]++;
				}
			}
		}
		cntEntryTree[AM2do]++;
		chain1p[AM2do]->GetEntry(cntEntryTree[AM2do]);
	}
	if (comptonFound) comptonChainsLength[AM2do]++; // assuming one chain for now
	cntEntryTree[AM2do]--;
	chain1p[AM2do]->GetEntry(cntEntryTree[AM2do]);
	
	Int_t clusterCounter = 0;
	Bool_t skipFrame = kFALSE;
	Int_t nPrimaryCluster = 0;
	Float_t primaryClusterE = 0;
	Float_t E1 = -1, E2 = -1;
	Float_t xc1 = -1, yc1 = -1;
	Float_t xc2 = -1, yc2 = -1;
	for (Int_t ii = 0; ii < nPixXY&&(!skipFrame); ii++)
	{
		for (Int_t jj = 0; jj < nPixXY&&(!skipFrame); jj++)
		{
			if (matrix_flags[ii][jj] != newNeighbourFlag) continue;
			Float_t CoGx = 0;
			Float_t CoGy = 0;
			Int_t area = 0;
			Float_t totAmp = 0;
			Float_t totAmpEcal = 0;
			Float_t CoGxEcal = 0;
			Float_t CoGyEcal = 0;
			Bool_t isPrimaryCluster = kFALSE;
			Int_t nTrigsInCluster = 0;
			//Int_t nPrimaryCluster = 0;
			
			buffX.clear();
			buffY.clear();
			
			buffX.resize(0,0);
			buffY.resize(0,0);
			
			getNewNeighbours(ii,jj,neighbourSearchWindow);
			Int_t nNb = buffX.size();
			clusterCounter++;
			
			matrix_flags[ii][jj] = clusterCounter;
			CoGx += (ii+0.5)*matrix_E[ii][jj];
			CoGy += (jj+0.5)*matrix_E[ii][jj];
			if (ii == centrePixelX-1 && jj == centrePixelY-1) isPrimaryCluster = kTRUE;
			if (matrix_trigs[ii][jj] == 1) nTrigsInCluster++;
			area++;
			
			totAmp += matrix_E[ii][jj];
			while (nNb > 0 && !skipFrame)
			{
				nNb = 0;
				buffXnew.swap(buffX);
				buffYnew.swap(buffY);
				buffX.clear();
				buffY.clear();
				buffX.resize(0,0);
				buffY.resize(0,0);
				for (Int_t ic = 0; ic < buffXnew.size(); ic++)
				{
					if (matrix_flags[buffXnew[ic]][buffYnew[ic]] == newNeighbourFlag)
					{
						CoGx += (buffXnew[ic]+0.5)*matrix_E[buffXnew[ic]][buffYnew[ic]];
						CoGy += (buffYnew[ic]+0.5)*matrix_E[buffXnew[ic]][buffYnew[ic]];
						if (buffXnew[ic] == centrePixelX-1 && buffYnew[ic] == centrePixelY-1) isPrimaryCluster = kTRUE;
						if (matrix_trigs[buffXnew[ic]][buffYnew[ic]] == 1) nTrigsInCluster++;
						area++;
						if (area > maxClusterSizeToSkipFrame)
						{
							skipFrame = kTRUE;
							cout << "Skipping frame, max cluster size exceeds " << maxClusterSizeToSkipFrame << endl;
							break;
						}
						totAmp += matrix_E[buffXnew[ic]][buffYnew[ic]];
					}
					matrix_flags[buffXnew[ic]][buffYnew[ic]] = clusterCounter;
				}
				for (Int_t ic = 0; ic < buffXnew.size(); ic++)
				{
					getNewNeighbours(buffXnew[ic],buffYnew[ic],neighbourSearchWindow);
					nNb += buffX.size();
				}
				buffXnew.clear();
				buffYnew.clear();
				buffXnew.resize(0,0);
				buffYnew.resize(0,0);
			}
			if (skipFrame) break;
			if (E1 < 0 && E2 < 0)
			{
				E1 = totAmp;
				xc1 = CoGx/totAmp;
				yc1 = CoGy/totAmp;
			}
			else if (E1 > 0 && E2 < 0)
			{
				E2 = totAmp;
				xc2 = CoGx/totAmp;
				yc2 = CoGy/totAmp;
			}
									
			if (isPrimaryCluster) nPrimaryCluster++;
			buffClusterXloc.push_back(CoGx/totAmp);
			buffClusterYloc.push_back(CoGy/totAmp);
			buffClusterEloc.push_back(totAmp);
			buffClusterArealoc.push_back(area);
			if (isPrimaryCluster) buffClusterFlagloc.push_back(1);
			else buffClusterFlagloc.push_back(0);
			buffClusterTrigsloc.push_back(nTrigsInCluster);
		}
	}

	if (clusterCounter > 0)
	{
		Int_t niter = 0;
		while (niter < buffClusterXloc.size())
		{
			Float_t maxE = -9999;
			Int_t idx = -1, maxArea, maxFlag, maxTrigs;
			Float_t maxX, maxY;
			for (Int_t ic = niter; ic < buffClusterXloc.size(); ic++)
			{
				if (buffClusterEloc[ic] > maxE)
				{
					maxE = buffClusterEloc[ic];
					idx = ic;
					maxX = buffClusterXloc[ic];
					maxY = buffClusterYloc[ic];
					maxArea = buffClusterArealoc[ic];
					maxFlag = buffClusterFlagloc[ic];
					maxTrigs = buffClusterTrigsloc[ic];
				}
			}
			buffClusterEloc.erase(buffClusterEloc.begin()+idx);
			buffClusterXloc.erase(buffClusterXloc.begin()+idx);
			buffClusterYloc.erase(buffClusterYloc.begin()+idx);
			buffClusterArealoc.erase(buffClusterArealoc.begin()+idx);
			buffClusterFlagloc.erase(buffClusterFlagloc.begin()+idx);
			buffClusterTrigsloc.erase(buffClusterTrigsloc.begin()+idx);
			buffClusterEloc.insert(buffClusterEloc.begin()+niter,maxE);
			buffClusterXloc.insert(buffClusterXloc.begin()+niter,maxX);
			buffClusterYloc.insert(buffClusterYloc.begin()+niter,maxY);
			buffClusterArealoc.insert(buffClusterArealoc.begin()+niter,maxArea);
			buffClusterFlagloc.insert(buffClusterFlagloc.begin()+niter,maxFlag);
			buffClusterTrigsloc.insert(buffClusterTrigsloc.begin()+niter,maxTrigs);
			niter++;
		}
		nTrigsInEvent_maxEcluster[AM2do]->Fill(buffClusterEloc[0],nTrigPixels[AM2do]);
		buffClusterX.push_back(buffClusterXloc);
		buffClusterY.push_back(buffClusterYloc);
		buffClusterE.push_back(buffClusterEloc);
		buffClusterArea.push_back(buffClusterArealoc);
		buffClusterFlag.push_back(buffClusterFlagloc);
		buffClusterTrigs.push_back(buffClusterTrigsloc);
	}
	if (nPrimaryCluster > 1) cout << "ERROR: more than one primary cluster detected, nPrimaryCluster = " << nPrimaryCluster << endl;
	return kTRUE;
}

Float_t getFWHM(const TH1F *histo, const Float_t peak_pos, const Float_t peak_H, const Float_t minADC_4search_sigma, const Float_t maxADC_4search_sigma, const Float_t minADC_4search,
							const Float_t maxADC_4search, Float_t &leftEdge, Float_t &rightEdge)
{
	leftEdge = 0;
	rightEdge = 0;
	Int_t peakBin = histo->GetXaxis()->FindBin(peak_pos);
	Int_t minBin = histo->GetXaxis()->FindBin(peak_pos - minADC_4search);
	if (histo->GetXaxis()->FindBin(peak_pos - minADC_4search_sigma) > minBin) minBin = histo->GetXaxis()->FindBin(peak_pos - minADC_4search_sigma);
	if (minBin < 2) minBin = 2;
	Int_t maxBin = histo->GetXaxis()->FindBin(maxADC_4search);
	if (histo->GetXaxis()->FindBin(peak_pos + minADC_4search_sigma) < maxBin) maxBin = histo->GetXaxis()->FindBin(peak_pos + minADC_4search_sigma);
	if (maxBin >= histo->GetNbinsX()) maxBin = histo->GetNbinsX() - 1;
	for (Int_t n = peakBin; n >= minBin; n--)
	{
		if (peak_H/2 >= histo->GetBinContent(n-1) && peak_H/2 <= histo->GetBinContent(n))
		{
			leftEdge = histo->GetXaxis()->GetBinCenter(n)-(histo->GetXaxis()->GetBinCenter(n) - histo->GetXaxis()->GetBinCenter(n-1))*
				(histo->GetBinContent(n) - peak_H/2)/(histo->GetBinContent(n) - histo->GetBinContent(n-1));
			break;
		}
	}
	for (Int_t n = peakBin; n <= maxBin; n++)
	{
		if (histo->GetXaxis()->GetBinCenter(n) > peak_pos && peak_H/2 <= histo->GetBinContent(n) && peak_H/2 >= histo->GetBinContent(n+1))
		{
			rightEdge = histo->GetXaxis()->GetBinCenter(n)+(histo->GetXaxis()->GetBinCenter(n+1) - histo->GetXaxis()->GetBinCenter(n))*
				(histo->GetBinContent(n) - peak_H/2)/(histo->GetBinContent(n) - histo->GetBinContent(n+1));
				break;
		}
	}
	if (rightEdge >= maxADC_4search) return 0;
	if (leftEdge <= minADC_4search) return 0;
	if (rightEdge - leftEdge > 0) return rightEdge - leftEdge;
	else return 0;
}

void getHistoPeak(const TH1F *histo, Float_t minX, Float_t maxX, Float_t &peak_pos1, Float_t &peak_H)
{ 
	peak_H = -9999;
	for (Int_t n = 1; n	< histo->GetNbinsX(); n++)
	{
		if (histo->GetXaxis()->GetBinCenter(n) < minX) continue;
		if (histo->GetXaxis()->GetBinCenter(n) > maxX) break;
		if (histo->GetBinContent(n) > peak_H)
		{
			peak_pos1 = histo->GetXaxis()->GetBinCenter(n);
			peak_H = histo->GetBinContent(n);
		}
	}
}

Bool_t readAnalysisSetupFile(TString fname)
{
	cout << "Reading analysis setup file:" << endl;
	cout << fname << endl;
	if (gSystem->AccessPathName(fname))  // Strange convention - this function return 0;s 1 (true) if path name doesn't exist !!!!
	{
		cerr << "ERROR: Analysis setup file " << fname << " doesn't exist. Exiting." << endl;
		return kFALSE;
	}
	
	Int_t buf = 0;
	in_file.open(fname, ios::in);
	while (!in_file.eof())
	{
		getline(in_file,line);
		TString sline = TString(line);
		if (sline.BeginsWith("/")) continue;
		
		if (sline.Contains("printOutSetupFileParameters ="))
		{
			sline.ReplaceAll("printOutSetupFileParameters =","");
			buf = sline.Atoi();
			printOutSetupFileParameters = kTRUE;
			if (buf == 0) printOutSetupFileParameters = kFALSE;
			if (printOutSetupFileParameters) cout << "printOutSetupFileParameters = " << printOutSetupFileParameters << endl;
		}

		if (sline.Contains("spectrumName ="))
		{
			sline.ReplaceAll("spectrumName =","");
			spectrumName = sline;
			if (printOutSetupFileParameters) cout << "spectrumName = " << spectrumName << endl;
		}
				
		if (sline.Contains("pathRoot_in_calib ="))
		{
			sline.ReplaceAll("pathRoot_in_calib =","");
			pathRoot_in_calib = sline;
			if (printOutSetupFileParameters) cout << "pathRoot_in_calib = " << pathRoot_in_calib << endl;
		}

		if (sline.Contains("TcalibFileName ="))
		{
			sline.ReplaceAll("TcalibFileName =","");
			TcalibFileName = sline;
			if (printOutSetupFileParameters) cout << "TcalibFileName = " << TcalibFileName << endl;
		}
			
		if (sline.Contains("fname_pixel_mapping_file ="))
		{
			sline.ReplaceAll("fname_pixel_mapping_file =","");
			fname_pixel_mapping_file = sline;
			if (printOutSetupFileParameters) cout << "fname_pixel_mapping_file = " << fname_pixel_mapping_file << endl;
		}

		if (sline.Contains("firstPixel ="))
		{
			sline.ReplaceAll("firstPixel =","");
			firstPixel = sline.Atoi();
			if (printOutSetupFileParameters) cout << "firstPixel = " << firstPixel << endl;
		}
		if (sline.Contains("lastPixel ="))
		{
			sline.ReplaceAll("lastPixel =","");
			lastPixel = sline.Atoi();
			if (printOutSetupFileParameters) cout << "lastPixel = " << lastPixel << endl;
		}
		if (sline.Contains("nPixXY ="))
		{
			sline.ReplaceAll("nPixXY =","");
			nPixXY = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nPixXY = " << nPixXY << endl;
		}
		if (sline.Contains("nAMs ="))
		{
			sline.ReplaceAll("nAMs =","");
			nAMs = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nAMs = " << nAMs << endl;
		}
		
		if (sline.Contains("minE ="))
		{
			sline.ReplaceAll("minE =","");
			minE = sline.Atof();
			if (printOutSetupFileParameters) cout << "minE = " << minE << endl;
		}
		
		if (sline.Contains("maxE ="))
		{
			sline.ReplaceAll("maxE =","");
			maxE = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxE = " << maxE << endl;
		}
		
		if (sline.Contains("nBinsE ="))
		{
			sline.ReplaceAll("nBinsE =","");
			nBinsE = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nBinsE = " << nBinsE << endl;
		}

		if (sline.Contains("minEcath ="))
		{
			sline.ReplaceAll("minEcath =","");
			minEcath = sline.Atof();
			if (printOutSetupFileParameters) cout << "minEcath = " << minEcath << endl;
		}
				
		if (sline.Contains("maxEcath ="))
		{
			sline.ReplaceAll("maxEcath =","");
			maxEcath = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxEcath = " << maxEcath << endl;
		}
				
		if (sline.Contains("nBinsEcath ="))
		{
			sline.ReplaceAll("nBinsEcath =","");
			nBinsEcath = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nBinsEcath = " << nBinsEcath << endl;
		}
		
		if (sline.Contains("minPhotopeakE4PeakSearch ="))
		{
			sline.ReplaceAll("minPhotopeakE4PeakSearch =","");
			minPhotopeakE4PeakSearch = sline.Atof();
			if (printOutSetupFileParameters) cout << "minPhotopeakE4PeakSearch = " << minPhotopeakE4PeakSearch << endl;
		}
		
		if (sline.Contains("maxPhotopeakE4PeakSearch ="))
		{
			sline.ReplaceAll("maxPhotopeakE4PeakSearch =","");
			maxPhotopeakE4PeakSearch = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxPhotopeakE4PeakSearch = " << maxPhotopeakE4PeakSearch << endl;
		}		
		
		if (sline.Contains("minClusterSize4H ="))
		{
			sline.ReplaceAll("minClusterSize4H =","");
			minClusterSize4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "minClusterSize4H = " << minClusterSize4H << endl;
		}
		
		if (sline.Contains("minNTrigs4H ="))
		{
			sline.ReplaceAll("minNTrigs4H =","");
			minNTrigs4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "minNTrigs4H = " << minNTrigs4H << endl;
		}
		
		if (sline.Contains("maxNTrigs4H ="))
		{
			sline.ReplaceAll("maxNTrigs4H =","");
			maxNTrigs4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "maxNTrigs4H = " << maxNTrigs4H << endl;
		}
		
		if (sline.Contains("minNTrigs4H_2D ="))
		{
			sline.ReplaceAll("minNTrigs4H_2D =","");
			minNTrigs4H_2D = sline.Atoi();
			if (printOutSetupFileParameters) cout << "minNTrigs4H_2D = " << minNTrigs4H_2D << endl;
		}
		
		if (sline.Contains("maxNTrigs4H_2D ="))
		{
			sline.ReplaceAll("maxNTrigs4H_2D =","");
			maxNTrigs4H_2D = sline.Atoi();
			if (printOutSetupFileParameters) cout << "maxNTrigs4H_2D = " << maxNTrigs4H_2D << endl;
		}
		
		if (sline.Contains("maxClusterSize4H ="))
		{
			sline.ReplaceAll("maxClusterSize4H =","");
			maxClusterSize4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "maxClusterSize4H = " << maxClusterSize4H << endl;
		}
		
		if (sline.Contains("nBinsTheta4H ="))
		{
			sline.ReplaceAll("nBinsTheta4H =","");
			nBinsTheta4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nBinsTheta4H = " << nBinsTheta4H << endl;
		}

		if (sline.Contains("Na22Energy ="))
		{
			sline.ReplaceAll("Na22Energy =","");
			Na22Energy = sline.Atof();
			if (printOutSetupFileParameters) cout << "Na22Energy = " << Na22Energy << endl;
		}
		
		if (sline.Contains("angleOfSecondHead ="))
		{
			sline.ReplaceAll("angleOfSecondHead =","");
			angleOfSecondHead = sline.Atof();
			if (printOutSetupFileParameters) cout << "angleOfSecondHead = " << angleOfSecondHead << endl;
		}
		
		if (sline.Contains("relativePhiAngle ="))
		{
			sline.ReplaceAll("relativePhiAngle =","");
			relativePhiAngle = sline.Atof();
			if (printOutSetupFileParameters) cout << "relativePhiAngle = " << relativePhiAngle << endl;
		}
		
		if (sline.Contains("EwindowFor4PhiAnalyis_min ="))
		{
			sline.ReplaceAll("EwindowFor4PhiAnalyis_min =","");
			EwindowFor4PhiAnalyis_min = sline.Atof();
			if (printOutSetupFileParameters) cout << "EwindowFor4PhiAnalyis_min = " << EwindowFor4PhiAnalyis_min << endl;
		}
		
		if (sline.Contains("EwindowFor4PhiAnalyis_max ="))
		{
			sline.ReplaceAll("EwindowFor4PhiAnalyis_max =","");
			EwindowFor4PhiAnalyis_max = sline.Atof();
			if (printOutSetupFileParameters) cout << "EwindowFor4PhiAnalyis_max = " << EwindowFor4PhiAnalyis_max << endl;
		}
		
		if (sline.Contains("ThetaWindowFor4PhiAnalyis_min ="))
		{
			sline.ReplaceAll("ThetaWindowFor4PhiAnalyis_min =","");
			ThetaWindowFor4PhiAnalyis_min = sline.Atof();
			if (printOutSetupFileParameters) cout << "ThetaWindowFor4PhiAnalyis_min = " << ThetaWindowFor4PhiAnalyis_min << endl;
		}
		
		if (sline.Contains("ThetaWindowFor4PhiAnalyis_max ="))
		{
			sline.ReplaceAll("ThetaWindowFor4PhiAnalyis_max =","");
			ThetaWindowFor4PhiAnalyis_max = sline.Atof();
			if (printOutSetupFileParameters) cout << "ThetaWindowFor4PhiAnalyis_max = " << ThetaWindowFor4PhiAnalyis_max << endl;
		}
		
		if (sline.Contains("nBins_dPhi ="))
		{
			sline.ReplaceAll("nBins_dPhi =","");
			nBins_dPhi = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nBins_dPhi = " << nBins_dPhi << endl;
		}
		
		if (sline.Contains("maxDistanceBetweenClusters4ComptonPair ="))
		{
			sline.ReplaceAll("maxDistanceBetweenClusters4ComptonPair =","");
			maxDistanceBetweenClusters4ComptonPair = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxDistanceBetweenClusters4ComptonPair = " << maxDistanceBetweenClusters4ComptonPair << endl;
		}

		if (sline.Contains("maxClusterSize ="))
		{
			sline.ReplaceAll("maxClusterSize =","");
			maxClusterSize = sline.Atoi();
			if (printOutSetupFileParameters) cout << "maxClusterSize = " << maxClusterSize << endl;
		}
		
		if (sline.Contains("neighbourSearchWindow ="))
		{
			sline.ReplaceAll("neighbourSearchWindow =","");
			neighbourSearchWindow = sline.Atoi();
			if (printOutSetupFileParameters) cout << "neighbourSearchWindow = " << neighbourSearchWindow << endl;
		}
		
		if (sline.Contains("minPixelEnergyThr4ClusterReconstruction ="))
		{
			sline.ReplaceAll("minPixelEnergyThr4ClusterReconstruction =","");
			minPixelEnergyThr4ClusterReconstruction = sline.Atof();
			if (printOutSetupFileParameters) cout << "minPixelEnergyThr4ClusterReconstruction = " << minPixelEnergyThr4ClusterReconstruction << endl;
		}
		
		if (sline.Contains("maxClusterSizeToSkipFrame ="))
		{
			sline.ReplaceAll("maxClusterSizeToSkipFrame =","");
			maxClusterSizeToSkipFrame = sline.Atoi();
			if (printOutSetupFileParameters) cout << "maxClusterSizeToSkipFrame = " << maxClusterSizeToSkipFrame << endl;
		}
		
		if (sline.Contains("makeTimingStuff ="))
		{
			sline.ReplaceAll("makeTimingStuff =","");
			buf = sline.Atoi();
			makeTimingStuff = kTRUE;
			if (buf == 0) makeTimingStuff = kFALSE;
			if (printOutSetupFileParameters) cout << "makeTimingStuff = " << makeTimingStuff << endl;
		}
		
		if (sline.Contains("typeOfPixelsToUse ="))
		{
			sline.ReplaceAll("typeOfPixelsToUse =","");
			typeOfPixelsToUse = sline.Atoi();
			if (printOutSetupFileParameters) cout << "typeOfPixelsToUse = " << typeOfPixelsToUse << endl;
		}
		
		if (sline.Contains("prntOutF ="))
		{
			sline.ReplaceAll("prntOutF =","");
			prntOutF = sline.Atoi();
			if (printOutSetupFileParameters) cout << "prntOutF = " << prntOutF << endl;
		}
		
		if (sline.Contains("useOnlyTriggeredPixelsInCluster ="))
		{
			sline.ReplaceAll("useOnlyTriggeredPixelsInCluster =","");
			buf = sline.Atoi();
			useOnlyTriggeredPixelsInCluster = kTRUE;
			if (buf == 0) useOnlyTriggeredPixelsInCluster = kFALSE;
			if (printOutSetupFileParameters) cout << "useOnlyTriggeredPixelsInCluster = " << useOnlyTriggeredPixelsInCluster << endl;
		}
	}
	in_file.close();
	in_file.clear();
	return kTRUE;
}

Int_t getCath1Channel(const Int_t AM_number, const Int_t ASIC_number, const Int_t serial)
{
	return Cath1ChannelBase + AM_number*16 + ASIC_number*4 + serial;
}

Float_t getDistanceBetweenPixels(const Int_t pix1, const Int_t pix2, const Int_t am)
{
	Int_t xx1 = -1, yy1 = -1;
	if (getPixel2Dcoordinates(pix1,am,xx1,yy1))
	{
		//cout << "pix1 = " << pix1 << ": xx1 = " << xx1 << ", yy1 = " << yy1 << endl;
		Int_t xx2 = -1, yy2 = -1;
		if (getPixel2Dcoordinates(pix2,am,xx2,yy2))
		{
			//cout << "pix2 = " << pix2 << ": xx2 = " << xx2 << ", yy2 = " << yy2 << endl;
			return TMath::Hypot(Float_t(xx1-xx2),Float_t(yy1-yy2));
		}
	}
	return -100;
}

Bool_t getPixel2Dcoordinates(const Int_t pixx, const Int_t am, Int_t &x1, Int_t &y1)
{
	for (Int_t ix = 1; ix <= pixelPattern[am]->GetNbinsX(); ix++)
	{
		for (Int_t iy = 1; iy <= pixelPattern[am]->GetNbinsY(); iy++)
		{
			if (pixelPattern[am]->GetBinContent(ix,iy) == pixx)
			{
				x1 = ix;
				y1 = iy;
				return kTRUE;
			}
		}
	}
	return kFALSE;
}

Float_t getPhiAngleDeg(const Float_t x, const Float_t y)
{
	Float_t ang = -9999;
	if (x >= 0 && y >= 0) ang = TMath::ATan(y/x)*TMath::RadToDeg();
	if (x < 0 && y >= 0) ang = 180 - TMath::ATan(y/fabs(x))*TMath::RadToDeg();
	if (x < 0 && y < 0) ang = 180 + TMath::ATan(fabs(y/x))*TMath::RadToDeg();
	if (x >= 0 && y < 0) ang = 360 - TMath::ATan(fabs(y)/x)*TMath::RadToDeg();
	return ang;
}

void getNewNeighbours(Int_t ix, Int_t iy, Int_t win)
{
	for (Int_t k = ix-win; k <= ix+win; k++)
	{
		if (k >= lastPixel || k < 0) continue;
		for (Int_t m = iy-win; m <= iy+win; m++)
		{
			if (m >= lastPixel || m < 0) continue;
			if (matrix_flags[k][m] == newNeighbourFlag)
			{
				buffX.push_back(k);
				buffY.push_back(m);				
			}
		}		
	}
}

Float_t getThetaFromEnergy(const Float_t eneRef, const Float_t ene)
{
	return TMath::ACos(1.-ene/(eneRef-ene))*TMath::RadToDeg();
}

Float_t getScatteredEnergyFromAngle(const Float_t eneRef, const Float_t ang)
{
	return eneRef - eneRef*(1. - TMath::Cos(ang*TMath::DegToRad()))/(2. - TMath::Cos(ang*TMath::DegToRad()));
}
