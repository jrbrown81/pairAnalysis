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
ofstream logfile;

//
// Running parameters
//

TString pathRoot = "./";

const TString ver = "v8";
const TString setup_file = "analysis_setup_pairAnalysis_" + ver + ".txt";

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
Int_t nEvents2Do = -1;

Int_t firstPixel = 1;
Int_t lastPixel = 968;
Int_t nPixXY = 22;
Int_t nAMs = 2;
Bool_t usePixelListToDisable = kFALSE;
TString pixelListFileName = "pixel_list";
Bool_t saveSelectedComptonRootTree = kFALSE;
Int_t savingPointLocation = 1;

Bool_t enableEventMonitor = kFALSE;
Int_t nEvents2Display = 10;
Int_t eventNumber2Display = 0;

Float_t minE = 0;
Float_t maxE = 600; // keV
Int_t nBinsE = 1200;
Float_t minEcath = 0; // keV
Float_t maxEcath = 700; // keV
Int_t nBinsEcath = 1400;
Float_t minPhotopeakE4PeakSearch = 460;
Float_t maxPhotopeakE4PeakSearch = 550;
Bool_t plotPixelClusterSpectra = kFALSE;

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
Bool_t makeCosineFit = kFALSE;
Float_t factor2ConvertAnodeTime2Distance = 120;

Bool_t makeRawTimingStuff = kTRUE;
Int_t minRawTiming4H = 0;
Int_t maxRawTiming4H = 1200;
Int_t nBinsRawTiming4H = 100;
Int_t minDeltaRawTiming4H = 0;
Int_t maxDeltaRawTiming4H = 500;
Int_t nBinsDeltaRawTiming4H = 100;

Bool_t makeTimingCalibrationStuff = kTRUE;
Float_t minAnodeTiming4H = 0;
Float_t maxAnodeTiming4H = 1200;
Int_t nBinsAnodeTiming4H = 100;
Float_t minDeltaAnodeTiming4H = 0;
Float_t maxDeltaAnodeTiming4H = 500;
Int_t nBinsDeltaAnodeTiming4H = 100;

// 0 - all, 1 - fully internal (disable two most external rows), 2 - internal (disable one external row), 3 - NOT fully internal, 4 - NOT internal
Int_t typeOfPixelsToUse = 0;

Int_t nBins_dPhi = 18;
Float_t maxDistanceBetweenClusters4ComptonPair = 7;

// Cluster reconstruction stuff
Int_t maxClusterSize = 10;
Float_t minPixelEnergyThr4ClusterReconstruction = 15;
Int_t neighbourSearchWindow = 1;
Bool_t doNotUseCornerTouchingPixels4Clustering = kTRUE;
Int_t maxClusterSizeToSkipFrame = 25;
Bool_t useOnlyTriggeredPixelsInCluster = kFALSE;
Int_t minClusterSize4PairAnalysis = -1;
Int_t maxClusterSize4PairAnalysis = -1;
Bool_t updateSinglePixelClusterCOGEneg = kTRUE;
Bool_t updateTwoPixelClusterCOGEneg = kTRUE;
Bool_t updateTwoPixelClusterCOG_1Donly = kTRUE;
Float_t maxDtBetweenPixel4Clustering = 200;
Bool_t applyMaxDtCut4Clustering = kFALSE;
Bool_t doNotUseCornerPixelsInPixelClusterCOGEneg = kFALSE;
Int_t minNClusters4PairAnalysis = 1;
Int_t maxNClusters4PairAnalysis = 2;

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
Int_t ***pixelStatusGlobal;
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

std::vector<std::vector<Int_t>> matrix_trigs;
std::vector<std::vector<Int_t>> matrix_flags;
std::vector<std::vector<Int_t>> matrix_time;
std::vector<std::vector<Float_t>> matrix_timeCalib;
std::vector<std::vector<Float_t>> matrix_E;
std::vector<std::vector<Float_t>> matrix_Eneg;
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
std::vector<std::vector<Int_t>> buffClusterIsSplit;
std::vector<std::vector<Float_t>> buffClusterAnodeTiming;
std::vector<Float_t> buffClusterXloc;
std::vector<Float_t> buffClusterYloc;
std::vector<Float_t> buffClusterEloc;
std::vector<Float_t> buffClusterArealoc;
std::vector<Float_t> buffClusterFlagloc;
std::vector<Float_t> buffClusterTrigsloc;
std::vector<Int_t> buffClusterXarr;
std::vector<Int_t> buffClusterYarr;
std::vector<Int_t> buffClusterTimearr;
std::vector<Int_t> buffClusterAnodeTimearr;
std::vector<Int_t> buffClusterTrigsarr;
std::vector<Int_t> buffClusterPixelarr;
std::vector<Int_t> buffClusterIsSplitloc;
std::vector<Float_t> buffClusterAnodeTimingloc;

TChain* chain_events;
TChain** chain1p;
Long64_t *event;
Int_t *nTrigPixels;
Float_t *cathodeE;
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
TLine *line2;
Bool_t *disabledPixels;
Int_t nDisabledPixels = 0;

TH1F **cathodeSpecAllEvents;
TH1F **cathodeSpecPair;
TH1F **cathodeSpecSelPairs;
TH1F **allClustersSpec;
TH1F **singleClusterEventsSpec;
TH1F **summedClusterEventsSpec;
TH1F **nTrigsInEvent;
TH2F **nTrigsInEvent_maxEcluster;
TH1F **comptonSummedSpec2Clusters;
TH1F **comptonSummedSpecEventClusters;
TH1F **comptonSpecClusters;
TH2F **comptonSpec2ClustersCorr;
TH2F *comptonSummedSpec2ClustersCorr;
TH2F *comptonSummedSpecEventCorr;
TH1F **phiAngle;
TH1F **phiAngleSelEvents;
TH1F *dPhiAngle;
TH1F *dPhiAngle_Ewin;
TH1F *dPhiAngleNorm;
TH1F *dPhiAngleNorm_Ewin;
TH1F *dPhiAngle1;
TH1F *dPhiAngle1_Ewin;
TH2F *dPhiAngle_2ClusterE;
TH2F *dPhiAngle_totalE;
TH2F *phiAngleCorr;
TH2F *phiAngleSelEventsCorr;
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
TH2F **imageE;
TH2F **imageT;
TH2F **imageC;
TH2I **imageN;
TH2F **comptonHitsSpecClustersSizeCorr;
TH2F **comptonHitsClusters_nTrigsCorr;
TH1F **nClustersInEvent;
TH2F *nClustersInEventCorr;
TH1F **nClustersInSelEvent;
TH2F *nClustersInSelEventCorr;
TH1F **firstComptonClusterSize;
TH1F **secondComptonClusterSize;
TH1F **firstComptonClusterSizeSelEvents;
TH1F **secondComptonClusterSizeSelEvents;
TH1F **thetaFromFirstClusterE_w;
TH1F **thetaFromFirstClusterE_1pixOnly;
TH1F **thetaFromFirstClusterE_1pixOnly_w;
TH1F **thetaFromFirstClusterE_2pixOnly;
TH1F **thetaFromFirstClusterE_1_2pixOnly;
TH2I **rawTimingTwoPixelClustersCorr;
TH1I **rawTimingDiff2PixClusters;
TH1I **rawTimingDiff2PixDiagClusters;
TH1I **rawTimingDiffAllClusters;
TH1I **rawTimingDiff3PixClusters;
TH2F **anodeTimingTwoPixelClustersCorr;
TH1F **anodeTimingDiff2PixDiagClusters;
TH1F **anodeTimingDiff2PixClusters;
TH1F **anodeTimingDiffAllClusters;
TH1F **anodeTimingDiff3PixClusters;
TH2F **anodeTiming2PixDiagClustersCorr;
TH2F **anodeTiming3PixClustersCorr;
TH2F **anodeTiming2PixClustersCorr;
TH2F **allClustersCOGImage;
TH2F **allPixelsImage;
TH2F **allClustersFineCOGImage;
TH2F **allPixelsFreqImage;
TH2F **allClustersCOGFreqImage;
TH2F **allClustersFineCOGFreqImage;
TH1F ***clusterSpec;
TH2F **anodeTimingTwoClustersAMCorr;
TH2F **dAnodeTiming_vs_theta;
TH2F **dAnodeTiming_vs_E1;
TH2F **dAnodeTiming_vs_E2;
TH2F **dAnodeTiming_vs_E12;
TH2F **dAnodeTimingAbs_vs_E12;
TH2F **dist12_vs_E1;
TH2F **dist12_vs_E2;
TH2F **dist12_vs_E12;
TH2F **dist12_over_E12sq;

TRandom3 *rand3;
TStopwatch localTimer;
Float_t totalTimeElapced = 0;
TString evMonDir;
TCanvas *c3;
Int_t nEventsDisplayed = 0;
Int_t firstClusterIdx[2] = {0};
Int_t secondClusterIdx[2] = {0};
Bool_t eventDisplayFlag = kFALSE;
Long64_t *selEvt;
Long64_t **selEvtAM;
Int_t *firstAMTreeEventIdx;
Int_t *lastAMTreeEventIdx;
Long64_t cntEv2BSaved = 0;
Long64_t evtTreeIdx;

Double_t fitfun(Double_t *, Double_t *);
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
Float_t updateSinglePixelClusterCOGwithEneg(const Int_t, const Int_t, Float_t &, Float_t &);
void updateTwoPixelClusterCOGwithEneg(const Int_t, const Float_t, const Float_t, const Int_t, const Int_t, const Int_t, const Int_t, Float_t &, Float_t &);
void sortClusters(const Int_t);

int main()
{
	if (!readAnalysisSetupFile("./" + setup_file))
	{
		cerr << "ERROR reading analysis setup file " << pathRoot + "/" + setup_file << ". Exiting." << endl;
		return 0;
	}
	
	if (minClusterSize4PairAnalysis > maxClusterSize4PairAnalysis && minClusterSize4PairAnalysis >= 0 && maxClusterSize4PairAnalysis >= 0)
	{
		cerr << "ERROR: Wrong settings for the cluster size in the analysis: minClusterSize4PairAnalysis = " << minClusterSize4PairAnalysis
				<< ",  maxClusterSize4PairAnalysis = " << maxClusterSize4PairAnalysis << ". Exiting." << endl;
		return 0;
	}
	
	if (minNClusters4PairAnalysis > maxNClusters4PairAnalysis && minNClusters4PairAnalysis >= 0 && maxNClusters4PairAnalysis >= 0)
	{
		cerr << "ERROR: Wrong settings for the number of clusters in the analysis: minNClusters4PairAnalysis = " << minNClusters4PairAnalysis
				<< ",  maxNClusters4PairAnalysis = " << maxNClusters4PairAnalysis << ". Exiting." << endl;
		return 0;
	}
	
	if (enableEventMonitor)
	{
		cout << endl;
		cout << "***********************************************" << endl;
		cout << "***** Careful!!! Event Monitor is enabled *****" << endl;
		cout << "***********************************************" << endl;
		cout << endl;
	}
	
	pixelPattern = new TH2I*[nAMs];
	pixelTypePattern = new TH2I*[nAMs];
	usedCentrePixelPattern = new TH2I*[nAMs];
	pixelStatus = new Bool_t**[nAMs];
	pixelStatusGlobal = new Int_t**[nAMs];
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
		pixelStatusGlobal[i] = new Int_t*[nPixXY];
		for (Int_t n = 0; n < nPixXY; n++)
		{
			pixelStatus[i][n] = new Bool_t[nPixXY];
			pixelStatusGlobal[i][n] = new Int_t[nPixXY];
			for (Int_t m = 0; m < nPixXY; m++)
			{
				pixelStatus[i][n][m] = 0;
				pixelStatusGlobal[i][n][m] = 0;
			}
		}
	}
	
	while (pathRoot_in_calib.BeginsWith(" ")) pathRoot_in_calib.Remove(0,1);
	if (!((TString) pathRoot_in_calib[pathRoot_in_calib.Length()]).IsAlnum()) pathRoot_in_calib = pathRoot_in_calib.Chop();
	
	while (pixelListFileName.BeginsWith(" ")) pixelListFileName.Remove(0,1);
	if (!((TString) pixelListFileName[pixelListFileName.Length()]).IsAlnum()) pixelListFileName = pixelListFileName.Chop();
	
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
	
	for (Int_t im = 0; im < nAMs; im++)
	{
		for (Int_t i = firstPixel; i <= lastPixel; i++)
		{
			Int_t xx1 = -1, yy1 = -1;
			if (getPixel2Dcoordinates(i,im,xx1,yy1))
			{
				pixelStatusGlobal[im][xx1-1][yy1-1] = 1;
				if (xx1 == 1) pixelStatusGlobal[im][xx1-1][yy1-1] = 0;
				if (xx1 == 11) pixelStatusGlobal[im][xx1-1][yy1-1] = 0;
				if (xx1 == 12) pixelStatusGlobal[im][xx1-1][yy1-1] = 0;
				if (xx1 == 22) pixelStatusGlobal[im][xx1-1][yy1-1] = 0;
				if (yy1 == 1) pixelStatusGlobal[im][xx1-1][yy1-1] = 0;
				if (yy1 == 11) pixelStatusGlobal[im][xx1-1][yy1-1] = 0;
				if (yy1 == 12) pixelStatusGlobal[im][xx1-1][yy1-1] = 0;
				if (yy1 == 22) pixelStatusGlobal[im][xx1-1][yy1-1] = 0;
				if (xx1 == 1 && yy1 == 1) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 1 && yy1 == 11) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 1 && yy1 == 12) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 1 && yy1 == 22) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 11 && yy1 == 1) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 11 && yy1 == 11) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 11 && yy1 == 12) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 11 && yy1 == 22) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 12 && yy1 == 1) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 12 && yy1 == 11) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 12 && yy1 == 12) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 12 && yy1 == 22) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 22 && yy1 == 1) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 22 && yy1 == 11) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 22 && yy1 == 12) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
				if (xx1 == 22 && yy1 == 22) pixelStatusGlobal[im][xx1-1][yy1-1] = -1;
			}
		}
	}

	gStyle->SetOptStat(0);
	//gStyle->SetOptStat(10);
	//gStyle->SetOptFit(111111);
	TH1::AddDirectory(0);
	gErrorIgnoreLevel = 1001;
	gStyle->SetPalette(1);
	rand3 = new TRandom3();
	rand3->SetSeed(0);
	
	if (makeTimingCalibrationStuff)
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
			if (pix >= firstPixel && pix <= lastPixel)
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
			for (Int_t i = firstPixel; i <= lastPixel; i++)
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
			for (Int_t i = firstPixel; i <= lastPixel; i++)
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
		
		for (Int_t i = 0; i < nPixXY; i++)
		{
			for (Int_t j = 0; j < nPixXY; j++)
			{
				pixelTypePattern[im]->SetBinContent(i+1,j+1,pixelStatus[im][i][j]);
			}
		}
	}
	
	TString subdir = "Pair analysis " + ver;
	subdir = Form("%s pixType%d",subdir.Data(),typeOfPixelsToUse);
	if (useOnlyTriggeredPixelsInCluster) subdir += " onlyTrigPixels";
	if (minClusterSize4PairAnalysis > 0 && maxClusterSize4PairAnalysis > 0)
	{
		if (minClusterSize4PairAnalysis == maxClusterSize4PairAnalysis) subdir += Form(" %d-pixCls",minClusterSize4PairAnalysis);
		if (minClusterSize4PairAnalysis < maxClusterSize4PairAnalysis) subdir += Form(" %d-%d-pixCls",minClusterSize4PairAnalysis,maxClusterSize4PairAnalysis);
	}
	if (minNClusters4PairAnalysis > 0 && maxNClusters4PairAnalysis > 0)
	{
		if (minNClusters4PairAnalysis == maxNClusters4PairAnalysis) subdir += Form(" %d-nCls",minNClusters4PairAnalysis);
		if (minNClusters4PairAnalysis < maxNClusters4PairAnalysis) subdir += Form(" %d-%d-nCls",minNClusters4PairAnalysis,maxNClusters4PairAnalysis);
	}
	TString outputPathPairs = Form("%s%s/pairAnalysis",pathRoot.Data(),subdir.Data());
	if (subdir.Length() > 1)
	{
		gSystem->MakeDirectory(pathRoot+subdir);
		gSystem->MakeDirectory(outputPathPairs);
		if (enableEventMonitor)
		{
			evMonDir = pathRoot + subdir + "/eventMonitor";
			gSystem->MakeDirectory(evMonDir);
		}
	}
	TString setup_file_copy = pathRoot + subdir + "/" + setup_file;
	setup_file_copy.ReplaceAll(".txt","_USED.txt");
	gSystem->CopyFile(setup_file, setup_file_copy, kTRUE);

	if (usePixelListToDisable)
	{
		TString pixelList_fname = "./" + pixelListFileName;
		// Check whether pixel list file exists and source the master file is defined so.
		if (gSystem->AccessPathName(pixelList_fname))  // Strange convention - this function return 0;s 1 (true) if path name doesn't exist !!!!
		{
			cerr << "ERROR: Pixel list file (" << pixelList_fname << ") doesn't exist. Exiting." << endl;
			return kFALSE;
		}
		TString pixelList_fname_copy = pathRoot + subdir + "/" + pixelListFileName;
		pixelList_fname_copy.ReplaceAll(".txt","_USED.txt");
		gSystem->CopyFile(pixelList_fname, pixelList_fname_copy, kTRUE);
		disabledPixels = new Bool_t[lastPixel+1];
		in_file.open(pixelList_fname, ios::in);
		Int_t cnnt = 0;
		for (Int_t nn = 0; nn < lastPixel+1; nn++) disabledPixels[nn] = kFALSE;
		while (!in_file.eof())
		{
			Int_t px = -1;
			in_file >> px;
			if (px >= firstPixel && px <= lastPixel)
			{
				disabledPixels[px] = kTRUE;
				nDisabledPixels++;
			}
		}
		in_file.close();
		in_file.clear();
		
		cout << "nDisabledPixels = " << nDisabledPixels << ": " << endl;
		for (Int_t nn = 0; nn < lastPixel+1; nn++) if (disabledPixels[nn]) cout << nn << endl;
	}
	
	TString *outputPathAM = new TString[nAMs];
	TString *outputPixelPathAM = new TString[nAMs];
	for (Int_t im = 0; im < nAMs; im++)
	{
		outputPathAM[im] = Form("%s%s/AM%d",pathRoot.Data(),subdir.Data(),im);
		gSystem->MakeDirectory(outputPathAM[im]);
		if (plotPixelClusterSpectra)
		{
			outputPixelPathAM[im] = Form("%s/pixelClusterSpectra",outputPathAM[im].Data());
			gSystem->MakeDirectory(outputPixelPathAM[im]);
		}
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
	if (saveSelectedComptonRootTree)
	{	selEvt = new Long64_t[nEvents2Analyse];
		selEvtAM = new Long64_t*[nAMs];
	}
	firstAMTreeEventIdx = new Int_t[nAMs];
	lastAMTreeEventIdx = new Int_t[nAMs];

	Nev = new Long64_t[nAMs];
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
		if (saveSelectedComptonRootTree)
		{
			selEvtAM[im] = new Long64_t[Nev[im]];
			for (Long64_t pr = 0; pr < Nev[im]; pr++) selEvtAM[im][pr] = 0;
		}
		
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
	}
	
	cathodeE = new Float_t[nAMs];
	nTrigsInEvent = new TH1F*[nAMs];
	clusterE_vsSize = new TH2F*[nAMs];
	clusterMaxE_vsSize = new TH2F*[nAMs];
	clusterSecE_vsSize = new TH2F*[nAMs];
	clusterSizeMaxE_SecE = new TH2F*[nAMs];
	nTrigsInEvent_maxEcluster = new TH2F*[nAMs];
	cathodeSpecAllEvents = new TH1F*[nAMs];
	cathodeSpecPair = new TH1F*[nAMs];
	cathodeSpecSelPairs = new TH1F*[nAMs];
	allClustersSpec = new TH1F*[nAMs];
	singleClusterEventsSpec = new TH1F*[nAMs];
	summedClusterEventsSpec = new TH1F*[nAMs];
	comptonSummedSpec2Clusters = new TH1F*[nAMs];
	comptonSummedSpecEventClusters = new TH1F*[nAMs];
	comptonSpecClusters = new TH1F*[nAMs];
	comptonSpec2ClustersCorr = new TH2F*[nAMs];
	comptonHitsSpecClustersSizeCorr = new TH2F*[nAMs];
	comptonHitsClusters_nTrigsCorr = new TH2F*[nAMs];
	phiAngle = new TH1F*[nAMs];
	phiAngleSelEvents = new TH1F*[nAMs];
	dPhiPlaneAM = new TH2F*[nAMs];
	firstComptonCluster = new TH2F*[nAMs];
	allClustersCOGImage = new TH2F*[nAMs];
	allClustersFineCOGImage = new TH2F*[nAMs];
	allPixelsImage = new TH2F*[nAMs];
	secondComptonCluster = new TH2F*[nAMs];
	thetaFromFirstClusterE = new TH1F*[nAMs];
	thetaFromFirstClusterE_dPhi = new TH2F*[nAMs];
	clusterSize = new TH1F*[nAMs];
	nTrigs_clusterSize = new TH2F*[nAMs];
	clusterE_nTrigs = new TH2F*[nAMs];
	nClustersInEvent = new TH1F*[nAMs];
	nClustersInSelEvent = new TH1F*[nAMs];	
	firstComptonClusterSize = new TH1F*[nAMs];
	secondComptonClusterSize = new TH1F*[nAMs];
	firstComptonClusterSizeSelEvents = new TH1F*[nAMs];
	secondComptonClusterSizeSelEvents = new TH1F*[nAMs];
	thetaFromFirstClusterE_w = new TH1F*[nAMs];
	thetaFromFirstClusterE_1pixOnly = new TH1F*[nAMs];
	thetaFromFirstClusterE_1pixOnly_w = new TH1F*[nAMs];
	thetaFromFirstClusterE_2pixOnly = new TH1F*[nAMs];
	thetaFromFirstClusterE_1_2pixOnly = new TH1F*[nAMs];
	allPixelsFreqImage = new TH2F*[nAMs];
	allClustersCOGFreqImage = new TH2F*[nAMs];
	allClustersFineCOGFreqImage = new TH2F*[nAMs];
	if (plotPixelClusterSpectra) clusterSpec = new TH1F**[nAMs];
	
	if (makeRawTimingStuff)
	{
		rawTimingTwoPixelClustersCorr = new TH2I*[nAMs];
		rawTimingDiff2PixClusters = new TH1I*[nAMs];
		rawTimingDiffAllClusters = new TH1I*[nAMs];
		if (!doNotUseCornerPixelsInPixelClusterCOGEneg) rawTimingDiff2PixDiagClusters = new TH1I*[nAMs];
		rawTimingDiff3PixClusters = new TH1I*[nAMs];
	}
	if (makeTimingCalibrationStuff)
	{
		anodeTimingTwoPixelClustersCorr = new TH2F*[nAMs];
		anodeTimingDiff2PixClusters = new TH1F*[nAMs];
		anodeTimingDiffAllClusters = new TH1F*[nAMs];
		anodeTimingDiff3PixClusters = new TH1F*[nAMs];
		if (!doNotUseCornerPixelsInPixelClusterCOGEneg)
		{
			anodeTiming2PixDiagClustersCorr = new TH2F*[nAMs];
			anodeTimingDiff2PixDiagClusters = new TH1F*[nAMs];
		}
		anodeTiming3PixClustersCorr = new TH2F*[nAMs];
		anodeTiming2PixClustersCorr = new TH2F*[nAMs];
		anodeTimingTwoClustersAMCorr = new TH2F*[nAMs];
		dAnodeTiming_vs_theta = new TH2F*[nAMs];
		dAnodeTiming_vs_E1 = new TH2F*[nAMs];
		dAnodeTiming_vs_E2 = new TH2F*[nAMs];
		dAnodeTiming_vs_E12 = new TH2F*[nAMs];
		dAnodeTimingAbs_vs_E12 = new TH2F*[nAMs];
		dist12_vs_E1 = new TH2F*[nAMs];
		dist12_vs_E2 = new TH2F*[nAMs];
		dist12_vs_E12 = new TH2F*[nAMs];
		dist12_over_E12sq = new TH2F*[nAMs];
	}

	if (enableEventMonitor)
	{
		imageE = new TH2F*[nAMs];
		imageT = new TH2F*[nAMs];
		imageC = new TH2F*[nAMs];
		imageN = new TH2I*[nAMs];
	}
	
	for (Int_t im = 0; im < nAMs; im++)
	{
		if (plotPixelClusterSpectra)
		{
			clusterSpec[im] = new TH1F*[lastPixel+1];
			for (Int_t b = 0; b < lastPixel+1; b++)
			{
				if (b >= firstPixel)
				{
					clusterSpec[im][b] = new TH1F(Form("clusterSpec_AM%d_pix%d",im,b),Form("%s, pixel %d cluster spectrum, AM%d",spectrumName.Data(),b,im),nBinsE,minE,maxE);
					clusterSpec[im][b]->GetXaxis()->SetTitleOffset(1.2);
					clusterSpec[im][b]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));
				}
				if (b < firstPixel) clusterSpec[im][b] = new TH1F();
			}
		}
		
		cathodeE[im] = 0;
		if (enableEventMonitor)
		{
			imageE[im] = new TH2F(Form("imageE_AM%d",im),"Event image",nPixXY,0,nPixXY,nPixXY,0,nPixXY);
			imageE[im]->GetXaxis()->SetTitle("Pixel number");
			imageE[im]->GetXaxis()->SetTitleOffset(1.2);
			imageE[im]->GetYaxis()->SetTitle("Pixel number");
			imageE[im]->GetYaxis()->SetTitleOffset(1.3);
			
			imageT[im] = new TH2F(Form("imageT_AM%d",im),"Event image",nPixXY*20,0,nPixXY,nPixXY*20,0,nPixXY);
			imageT[im]->GetXaxis()->SetTitle("Pixel number");
			imageT[im]->GetXaxis()->SetTitleOffset(1.2);
			imageT[im]->GetYaxis()->SetTitle("Pixel number");
			imageT[im]->GetYaxis()->SetTitleOffset(1.3);
			imageT[im]->SetMarkerStyle(20);
			imageT[im]->SetMarkerSize(0.8);
			
			imageC[im] = new TH2F(Form("imageC_AM%d",im),"Event image",nPixXY*20,0,nPixXY,nPixXY*20,0,nPixXY);
			imageC[im]->GetXaxis()->SetTitle("Pixel number");
			imageC[im]->GetXaxis()->SetTitleOffset(1.2);
			imageC[im]->GetYaxis()->SetTitle("Pixel number");
			imageC[im]->GetYaxis()->SetTitleOffset(1.3);
			imageC[im]->SetMarkerStyle(5);
			imageC[im]->SetMarkerSize(1.6);
			
			imageN[im] = new TH2I(Form("imageN_AM%d",im),"Event image",nPixXY,0,nPixXY,nPixXY,0,nPixXY);
			imageN[im]->GetXaxis()->SetTitle("Pixel number");
			imageN[im]->GetXaxis()->SetTitleOffset(1.2);
			imageN[im]->GetYaxis()->SetTitle("Pixel number");
			imageN[im]->GetYaxis()->SetTitleOffset(1.3);
		}
		thetaFromFirstClusterE[im] = new TH1F(Form("thetaFromFirstClusterE_AM%d",im),Form("%s, #theta calculated from first interaction energy, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180);
		thetaFromFirstClusterE[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromFirstClusterE[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));

		thetaFromFirstClusterE_w[im] = new TH1F(Form("thetaFromFirstClusterE_w_AM%d",im),Form("%s, #theta calculated from first interaction energy, weighted, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180);
		thetaFromFirstClusterE_w[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromFirstClusterE_w[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));

		thetaFromFirstClusterE_1pixOnly[im] = new TH1F(Form("thetaFromFirstClusterE_1pixOnly_AM%d",im),Form("%s, #theta calculated from first interaction energy, 1-pixel clusters only, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180);
		thetaFromFirstClusterE_1pixOnly[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromFirstClusterE_1pixOnly[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));

		thetaFromFirstClusterE_1pixOnly_w[im] = new TH1F(Form("thetaFromFirstClusterE_1pixOnly_w_AM%d",im),Form("%s, #theta calculated from first interaction energy, 1-pixel clusters only, weighted, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180);
		thetaFromFirstClusterE_1pixOnly_w[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromFirstClusterE_1pixOnly_w[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));

		thetaFromFirstClusterE_2pixOnly[im] = new TH1F(Form("thetaFromFirstClusterE_2pixOnly_AM%d",im),Form("%s, #theta calculated from first interaction energy, 2-pixel clusters only, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180);
		thetaFromFirstClusterE_2pixOnly[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromFirstClusterE_2pixOnly[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));

		thetaFromFirstClusterE_1_2pixOnly[im] = new TH1F(Form("thetaFromFirstClusterE_1_2pixOnly_AM%d",im),Form("%s, #theta calculated from first interaction energy, 1-2 pixel clusters combinations, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180);
		thetaFromFirstClusterE_1_2pixOnly[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromFirstClusterE_1_2pixOnly[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));

		thetaFromFirstClusterE_dPhi[im] = new TH2F(Form("thetaFromFirstClusterE_dPhi_AM%d",im),Form("%s, #theta calculated from first interaction energy, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180,nBins_dPhi,-180,180);
		thetaFromFirstClusterE_dPhi[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromFirstClusterE_dPhi[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));
		thetaFromFirstClusterE_dPhi[im]->GetYaxis()->SetTitleOffset(1.4);
		thetaFromFirstClusterE_dPhi[im]->GetYaxis()->SetTitle(Form("#Delta#varphi AM%d, degrees",im));

		nTrigsInEvent[im] = new TH1F(Form("nTrigsInEvent_AM%d",im),Form("%s, number of triggers per event, AM%d",spectrumName.Data(),im),maxNTrigs4H-minNTrigs4H,minNTrigs4H,maxNTrigs4H);
		nTrigsInEvent[im]->GetXaxis()->SetTitleOffset(1.2);
		nTrigsInEvent[im]->GetXaxis()->SetTitle(Form("Number of triggers in AM%d",im));

		nClustersInEvent[im] = new TH1F(Form("nClustersInEvent_AM%d",im),Form("%s, number of clusters per event, AM%d",spectrumName.Data(),im),maxNTrigs4H-minNTrigs4H,minNTrigs4H,maxNTrigs4H);
		nClustersInEvent[im]->GetXaxis()->SetTitleOffset(1.2);
		nClustersInEvent[im]->GetXaxis()->SetTitle(Form("Number of clusters in AM%d",im));

		nClustersInSelEvent[im] = new TH1F(Form("nClustersInSelEvent_AM%d",im),Form("%s, number of clusters/event in selected events, AM%d",spectrumName.Data(),im),maxNTrigs4H-minNTrigs4H,minNTrigs4H,maxNTrigs4H);
		nClustersInSelEvent[im]->GetXaxis()->SetTitleOffset(1.2);
		nClustersInSelEvent[im]->GetXaxis()->SetTitle(Form("Number of clusters in AM%d",im));
		
		nTrigs_clusterSize[im] = new TH2F(Form("nTrigs_clusterSize_AM%d",im),Form("%s, number of triggers in cluster vs cluster size, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H,maxNTrigs4H_2D-minNTrigs4H_2D,minNTrigs4H_2D,maxNTrigs4H_2D);
		nTrigs_clusterSize[im]->GetXaxis()->SetTitleOffset(1.2);
		nTrigs_clusterSize[im]->GetXaxis()->SetTitle(Form("AM%d cluster size, pixels",im));	
		nTrigs_clusterSize[im]->GetYaxis()->SetTitleOffset(1.4);
		nTrigs_clusterSize[im]->GetYaxis()->SetTitle("Number of triggers in cluster");		
		
		if (makeRawTimingStuff) 
		{
			rawTimingTwoPixelClustersCorr[im] = new TH2I(Form("rawTimingTwoPixelClustersCorr_AM%d",im),Form("%s, Correlation between triggered pixel raw timing, two triggered pixel clusters only, AM%d",spectrumName.Data(),im),nBinsRawTiming4H,minRawTiming4H,maxRawTiming4H,nBinsRawTiming4H,minRawTiming4H,maxRawTiming4H);
			rawTimingTwoPixelClustersCorr[im]->GetXaxis()->SetTitleOffset(1.2);
			rawTimingTwoPixelClustersCorr[im]->GetXaxis()->SetTitle(Form("Raw timing, first pixel, AM%d",im));	
			rawTimingTwoPixelClustersCorr[im]->GetYaxis()->SetTitleOffset(1.4);
			rawTimingTwoPixelClustersCorr[im]->GetYaxis()->SetTitle(Form("Raw timing, second pixel, AM%d",im));	

			rawTimingDiff2PixClusters[im] = new TH1I(Form("rawTimingDiff2PixClusters_AM%d",im),Form("%s, Difference in triggered pixel raw timing, two triggered pixel clusters only, AM%d",spectrumName.Data(),im),nBinsDeltaRawTiming4H,minDeltaRawTiming4H,maxDeltaRawTiming4H);
			rawTimingDiff2PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			rawTimingDiff2PixClusters[im]->GetXaxis()->SetTitle(Form("Raw timing difference between pixels, AM%d",im));	

			rawTimingDiff3PixClusters[im] = new TH1I(Form("rawTimingDiff3PixClusters_AM%d",im),Form("%s, Difference in triggered pixel raw timing, three and more triggered pixel clusters only, AM%d",spectrumName.Data(),im),nBinsDeltaRawTiming4H,minDeltaRawTiming4H,maxDeltaRawTiming4H);
			rawTimingDiff3PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			rawTimingDiff3PixClusters[im]->GetXaxis()->SetTitle(Form("Raw timing difference between pixels, AM%d",im));	
			
			rawTimingDiffAllClusters[im] = new TH1I(Form("rawTimingDiffAllClusters_AM%d",im),Form("%s, Difference in triggered pixel raw timing, all clusters, AM%d",spectrumName.Data(),im),nBinsDeltaRawTiming4H,minDeltaRawTiming4H,maxDeltaRawTiming4H);
			rawTimingDiffAllClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			rawTimingDiffAllClusters[im]->GetXaxis()->SetTitle(Form("Raw timing difference between pixels, AM%d",im));
			
			if (!doNotUseCornerPixelsInPixelClusterCOGEneg) 
			{
				rawTimingDiff2PixDiagClusters[im] = new TH1I(Form("rawTimingDiff2PixDiagClusters_AM%d",im),Form("%s, Difference in triggered pixel raw timing, two triggered diagonal pixel clusters only, AM%d",spectrumName.Data(),im),nBinsDeltaRawTiming4H,minDeltaRawTiming4H,maxDeltaRawTiming4H);
				rawTimingDiff2PixDiagClusters[im]->GetXaxis()->SetTitleOffset(1.2);
				rawTimingDiff2PixDiagClusters[im]->GetXaxis()->SetTitle(Form("Raw timing difference between pixels, AM%d",im));
			}
		}
		
		if (makeTimingCalibrationStuff)
		{
			anodeTimingTwoPixelClustersCorr[im] = new TH2F(Form("anodeTimingTwoPixelClustersCorr_AM%d",im),Form("%s, Correlation between triggered pixel anode timing, two triggered pixel clusters only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H);
			anodeTimingTwoPixelClustersCorr[im]->GetXaxis()->SetTitleOffset(1.2);
			anodeTimingTwoPixelClustersCorr[im]->GetXaxis()->SetTitle(Form("Anode timing, first pixel, ns, AM%d",im));	
			anodeTimingTwoPixelClustersCorr[im]->GetYaxis()->SetTitleOffset(1.4);
			anodeTimingTwoPixelClustersCorr[im]->GetYaxis()->SetTitle(Form("Anode timing, second pixel, ns, AM%d",im));	

			if (!doNotUseCornerPixelsInPixelClusterCOGEneg) 
			{
				anodeTimingDiff2PixDiagClusters[im] = new TH1F(Form("anodeTimingDiff2PixDiagClusters_AM%d",im),Form("%s, Difference in triggered pixel anode timing, two triggered diagonal pixel clusters only, AM%d",spectrumName.Data(),im),nBinsDeltaRawTiming4H,minDeltaRawTiming4H,maxDeltaRawTiming4H);
				anodeTimingDiff2PixDiagClusters[im]->GetXaxis()->SetTitleOffset(1.2);
				anodeTimingDiff2PixDiagClusters[im]->GetXaxis()->SetTitle(Form("Anode timing difference between pixels, ns, AM%d",im));	
			
				anodeTiming2PixDiagClustersCorr[im] = new TH2F(Form("anodeTiming2PixDiagClustersCorr_AM%d",im),Form("%s, Correlation between triggered pixel anode timing, two triggered diagonal pixel clusters only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H);
				anodeTiming2PixDiagClustersCorr[im]->GetXaxis()->SetTitleOffset(1.2);
				anodeTiming2PixDiagClustersCorr[im]->GetXaxis()->SetTitle(Form("Anode timing, first pixel, ns, AM%d",im));	
				anodeTiming2PixDiagClustersCorr[im]->GetYaxis()->SetTitleOffset(1.4);
				anodeTiming2PixDiagClustersCorr[im]->GetYaxis()->SetTitle(Form("Anode timing, second pixel, ns, AM%d",im));	
			}
			
			anodeTimingDiff2PixClusters[im] = new TH1F(Form("anodeTimingDiff2PixClusters_AM%d",im),Form("%s, Difference in triggered pixel anode timing, two triggered pixel clusters only, AM%d",spectrumName.Data(),im),nBinsDeltaAnodeTiming4H,minDeltaAnodeTiming4H,maxDeltaAnodeTiming4H);
			anodeTimingDiff2PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			anodeTimingDiff2PixClusters[im]->GetXaxis()->SetTitle(Form("Anode timing difference between pixels, ns, AM%d",im));	

			anodeTimingDiff3PixClusters[im] = new TH1F(Form("anodeTimingDiff3PixClusters_AM%d",im),Form("%s, Difference in triggered pixel anode timing, three and more triggered pixel clusters only, AM%d",spectrumName.Data(),im),nBinsDeltaAnodeTiming4H,minDeltaAnodeTiming4H,maxDeltaAnodeTiming4H);
			anodeTimingDiff3PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			anodeTimingDiff3PixClusters[im]->GetXaxis()->SetTitle(Form("Anode timing difference between pixels, ns, AM%d",im));	
			
			anodeTimingDiffAllClusters[im] = new TH1F(Form("anodeTimingDiffAllClusters_AM%d",im),Form("%s, Difference in triggered pixel anode timing, all clusters, AM%d",spectrumName.Data(),im),nBinsDeltaAnodeTiming4H,minDeltaAnodeTiming4H,maxDeltaAnodeTiming4H);
			anodeTimingDiffAllClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			anodeTimingDiffAllClusters[im]->GetXaxis()->SetTitle(Form("Anode timing difference between pixels, ns, AM%d",im));	

			anodeTimingTwoPixelClustersCorr[im] = new TH2F(Form("anodeTimingTwoPixelClustersCorr_AM%d",im),Form("%s, Correlation between triggered pixel anode timing, two triggered pixel clusters only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H);
			anodeTimingTwoPixelClustersCorr[im]->GetXaxis()->SetTitleOffset(1.2);
			anodeTimingTwoPixelClustersCorr[im]->GetXaxis()->SetTitle(Form("Anode timing, first pixel, ns, AM%d",im));	
			anodeTimingTwoPixelClustersCorr[im]->GetYaxis()->SetTitleOffset(1.4);
			anodeTimingTwoPixelClustersCorr[im]->GetYaxis()->SetTitle(Form("Anode timing, second pixel, ns, AM%d",im));	

			anodeTiming2PixClustersCorr[im] = new TH2F(Form("anodeTiming2PixClustersCorr_AM%d",im),Form("%s, Correlation between triggered pixel anode timing, two triggered pixel clusters only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H);
			anodeTiming2PixClustersCorr[im]->GetXaxis()->SetTitleOffset(1.2);
			anodeTiming2PixClustersCorr[im]->GetXaxis()->SetTitle(Form("Anode timing, first pixel, ns, AM%d",im));	
			anodeTiming2PixClustersCorr[im]->GetYaxis()->SetTitleOffset(1.4);
			anodeTiming2PixClustersCorr[im]->GetYaxis()->SetTitle(Form("Anode timing, second pixel, ns, AM%d",im));	

			anodeTiming3PixClustersCorr[im] = new TH2F(Form("anodeTiming3PixClustersCorr_AM%d",im),Form("%s, Correlation between triggered pixel anode timing, two triggered pixel clusters only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H);
			anodeTiming3PixClustersCorr[im]->GetXaxis()->SetTitleOffset(1.2);
			anodeTiming3PixClustersCorr[im]->GetXaxis()->SetTitle(Form("Anode timing, first pixel, ns, AM%d",im));	
			anodeTiming3PixClustersCorr[im]->GetYaxis()->SetTitleOffset(1.4);
			anodeTiming3PixClustersCorr[im]->GetYaxis()->SetTitle(Form("Anode timing, second pixel, ns, AM%d",im));
			
			anodeTimingTwoClustersAMCorr[im] = new TH2F(Form("anodeTimingTwoClustersAMCorr_AM%d",im),Form("%s, Correlation between cluster anode timing, two cluster events only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H);
			anodeTimingTwoClustersAMCorr[im]->GetXaxis()->SetTitleOffset(1.2);
			anodeTimingTwoClustersAMCorr[im]->GetXaxis()->SetTitle(Form("Anode timing, first cluster, ns, AM%d",im));	
			anodeTimingTwoClustersAMCorr[im]->GetYaxis()->SetTitleOffset(1.4);
			anodeTimingTwoClustersAMCorr[im]->GetYaxis()->SetTitle(Form("Anode timing, second cluster, ns, AM%d",im));

			dAnodeTiming_vs_theta[im] = new TH2F(Form("dAnodeTiming_vs_theta_AM%d",im),Form("%s, #Delta(two cluster anode timings) vs scattering angle, two cluster events only, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180,100,-500,500);
			dAnodeTiming_vs_theta[im]->GetXaxis()->SetTitleOffset(1.2);
			dAnodeTiming_vs_theta[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));
			dAnodeTiming_vs_theta[im]->GetYaxis()->SetTitleOffset(1.4);
			dAnodeTiming_vs_theta[im]->GetYaxis()->SetTitle(Form("#Delta(Anode timing), ns, AM%d",im));
			
			dAnodeTiming_vs_E1[im] = new TH2F(Form("dAnodeTiming_vs_E1_AM%d",im),Form("%s, #Delta(Anode timing) = t_{2} - t_{1} vs first cluster energy, two cluster events only, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,100,-500,500);
			dAnodeTiming_vs_E1[im]->GetXaxis()->SetTitleOffset(1.2);
			dAnodeTiming_vs_E1[im]->GetXaxis()->SetTitle(Form("First cluster energy, AM%d, keV",im));
			dAnodeTiming_vs_E1[im]->GetYaxis()->SetTitleOffset(1.5);
			dAnodeTiming_vs_E1[im]->GetYaxis()->SetTitle(Form("#Delta(Anode timing) = t_{2} - t_{1}, ns, AM%d",im));
			
			dAnodeTiming_vs_E2[im] = new TH2F(Form("dAnodeTiming_vs_E2_AM%d",im),Form("%s, #Delta(Anode timing) = t_{1} - t_{2} vs second cluster energy, two cluster events only, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,100,-500,500);
			dAnodeTiming_vs_E2[im]->GetXaxis()->SetTitleOffset(1.2);
			dAnodeTiming_vs_E2[im]->GetXaxis()->SetTitle(Form("Second cluster energy, AM%d, keV",im));
			dAnodeTiming_vs_E2[im]->GetYaxis()->SetTitleOffset(1.5);
			dAnodeTiming_vs_E2[im]->GetYaxis()->SetTitle(Form("#Delta(Anode timing) = t_{1} - t_{2}, ns, AM%d",im));
			
			dAnodeTiming_vs_E12[im] = new TH2F(Form("dAnodeTiming_vs_E12_AM%d",im),Form("%s, #Delta(Anode timing) vs 'second' cluster energy, two cluster events only, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,100,-500,500);
			dAnodeTiming_vs_E12[im]->GetXaxis()->SetTitleOffset(1.2);
			dAnodeTiming_vs_E12[im]->GetXaxis()->SetTitle(Form("'Second' cluster energy, AM%d, keV",im));
			dAnodeTiming_vs_E12[im]->GetYaxis()->SetTitleOffset(1.5);
			dAnodeTiming_vs_E12[im]->GetYaxis()->SetTitle(Form("#Delta(Anode timing), ns, AM%d",im));
			
			dAnodeTimingAbs_vs_E12[im] = new TH2F(Form("dAnodeTimingAbs_vs_E12_AM%d",im),Form("%s, |#Delta(Anode timing)| vs 'second' cluster energy, two cluster events only, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,100,0,500);
			dAnodeTimingAbs_vs_E12[im]->GetXaxis()->SetTitleOffset(1.2);
			dAnodeTimingAbs_vs_E12[im]->GetXaxis()->SetTitle(Form("'Second' cluster energy, AM%d, keV",im));
			dAnodeTimingAbs_vs_E12[im]->GetYaxis()->SetTitleOffset(1.5);
			dAnodeTimingAbs_vs_E12[im]->GetYaxis()->SetTitle(Form("|#Delta(Anode timing)|, ns, AM%d",im));
			
			dist12_vs_E1[im] = new TH2F(Form("dist12_vs_E1_AM%d",im),Form("%s, 'Distance 1-2' vs first cluster energy, two cluster events only, AM%d",spectrumName.Data(),im),nBinsE,minE,2e5,100,0,11);
			dist12_vs_E1[im]->GetXaxis()->SetTitleOffset(1.2);
			dist12_vs_E1[im]->GetXaxis()->SetTitle(Form("First cluster energy, AM%d, keV",im));
			dist12_vs_E1[im]->GetYaxis()->SetTitleOffset(1.5);
			dist12_vs_E1[im]->GetYaxis()->SetTitle(Form("'Distance 1-2', 'mm', AM%d",im));
			
			dist12_vs_E2[im] = new TH2F(Form("dist12_vs_E2_AM%d",im),Form("%s, 'Distance 1-2' vs second cluster energy, two cluster events only, AM%d",spectrumName.Data(),im),nBinsE,minE,5e5,100,0,11);
			dist12_vs_E2[im]->GetXaxis()->SetTitleOffset(1.2);
			dist12_vs_E2[im]->GetXaxis()->SetTitle(Form("Second cluster energy, AM%d, keV",im));
			dist12_vs_E2[im]->GetYaxis()->SetTitleOffset(1.5);
			dist12_vs_E2[im]->GetYaxis()->SetTitle(Form("'Distance 1-2', 'mm', AM%d",im));
			
			dist12_vs_E12[im] = new TH2F(Form("dist12_vs_E12_AM%d",im),Form("%s, 'Distance 1-2' vs 'second' cluster energy, two cluster events only, AM%d",spectrumName.Data(),im),nBinsE,minE,5e5,100,0,11);
			dist12_vs_E12[im]->GetXaxis()->SetTitleOffset(1.2);
			dist12_vs_E12[im]->GetXaxis()->SetTitle(Form("'Second' cluster energy, AM%d, keV",im));
			dist12_vs_E12[im]->GetYaxis()->SetTitleOffset(1.5);
			dist12_vs_E12[im]->GetYaxis()->SetTitle(Form("'Distance 1-2', 'mm', AM%d",im));

			dist12_over_E12sq[im] = new TH2F(Form("dist12_over_E12sq_AM%d",im),Form("%s, 'Distance 1-2' in #mum over E_{cluster}^{2.2}, two cluster events only, AM%d",spectrumName.Data(),im),100,0,0.08,100,0,0.5);
			dist12_over_E12sq[im]->GetXaxis()->SetTitleOffset(1.2);
			dist12_over_E12sq[im]->GetXaxis()->SetTitle(Form("Ratio with E_{1}, AM%d",im));
			dist12_over_E12sq[im]->GetYaxis()->SetTitleOffset(1.5);
			dist12_over_E12sq[im]->GetYaxis()->SetTitle(Form("Ratio with E_{2}, AM%d",im));
		}
		
		clusterSize[im] = new TH1F(Form("clusterSize_AM%d",im),Form("%s, cluster size, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
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

		allClustersSpec[im] = new TH1F(Form("allClustersSpec_AM%d",im),Form("%s, all cluster spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		allClustersSpec[im]->GetXaxis()->SetTitleOffset(1.2);
		allClustersSpec[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		

		cathodeSpecAllEvents[im] = new TH1F(Form("cathodeSpecAllEvents_AM%d",im),Form("%s, cathode spectrum, all evetns, AM%d",spectrumName.Data(),im),200,0,600);
		cathodeSpecAllEvents[im]->GetXaxis()->SetTitleOffset(1.2);
		cathodeSpecAllEvents[im]->GetXaxis()->SetTitle(Form("AM%d cathode energy, keV",im));		

		cathodeSpecPair[im] = new TH1F(Form("cathodeSpecPair_AM%d",im),Form("%s, cathode spectrum, Compton pair events, AM%d",spectrumName.Data(),im),200,0,600);
		cathodeSpecPair[im]->GetXaxis()->SetTitleOffset(1.2);
		cathodeSpecPair[im]->GetXaxis()->SetTitle(Form("AM%d cathode energy, keV",im));		

		cathodeSpecSelPairs[im] = new TH1F(Form("cathodeSpecSelPairs_AM%d",im),Form("%s, cathode spectrum, selected Compton pair events, AM%d",spectrumName.Data(),im),200,0,600);
		cathodeSpecSelPairs[im]->GetXaxis()->SetTitleOffset(1.2);
		cathodeSpecSelPairs[im]->GetXaxis()->SetTitle(Form("AM%d cathode energy, keV",im));		

		singleClusterEventsSpec[im] = new TH1F(Form("singleClusterEventsSpec_AM%d",im),Form("%s, single cluster events, cluster spectum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		singleClusterEventsSpec[im]->GetXaxis()->SetTitleOffset(1.2);
		singleClusterEventsSpec[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		

		summedClusterEventsSpec[im] = new TH1F(Form("summedClusterEventsSpec_AM%d",im),Form("%s, single cluster events + summed Compton events, cluster spectum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		summedClusterEventsSpec[im]->GetXaxis()->SetTitleOffset(1.2);
		summedClusterEventsSpec[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		

		comptonSpecClusters[im] = new TH1F(Form("comptonSpecClusters_AM%d",im),Form("%s, unsummed Compton interaction cluster spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		comptonSpecClusters[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonSpecClusters[im]->GetXaxis()->SetTitle(Form("AM%d cluster energy, keV",im));		

		comptonSpec2ClustersCorr[im] = new TH2F(Form("comptonSpec2ClustersCorr_AM%d",im),Form("%s, unsummed Compton interactions, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,nBinsE,minE,maxE);
		comptonSpec2ClustersCorr[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonSpec2ClustersCorr[im]->GetXaxis()->SetTitle(Form("AM%d first cluster energy, keV",im));		
		comptonSpec2ClustersCorr[im]->GetYaxis()->SetTitleOffset(1.4);
		comptonSpec2ClustersCorr[im]->GetYaxis()->SetTitle(Form("AM%d second cluster energy, keV",im));		

		comptonHitsSpecClustersSizeCorr[im] = new TH2F(Form("comptonHitsSpecClustersSizeCorr_AM%d",im),Form("%s, unsummed Compton interactions, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H,maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		comptonHitsSpecClustersSizeCorr[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonHitsSpecClustersSizeCorr[im]->GetXaxis()->SetTitle(Form("AM%d first cluster size, pixels",im));	
		comptonHitsSpecClustersSizeCorr[im]->GetYaxis()->SetTitleOffset(1.4);
		comptonHitsSpecClustersSizeCorr[im]->GetYaxis()->SetTitle(Form("AM%d second cluster size, pixels",im));				

		firstComptonClusterSize[im] = new TH1F(Form("firstComptonClusterSize_AM%d",im),Form("%s, cluster size of first Compton cluster, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		firstComptonClusterSize[im]->GetXaxis()->SetTitleOffset(1.4);
		firstComptonClusterSize[im]->GetXaxis()->SetTitle(Form("AM%d cluster size, pixels",im));		

		secondComptonClusterSize[im] = new TH1F(Form("secondComptonClusterSize_AM%d",im),Form("%s, cluster size of second Compton cluster, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		secondComptonClusterSize[im]->GetXaxis()->SetTitleOffset(1.4);
		secondComptonClusterSize[im]->GetXaxis()->SetTitle(Form("AM%d cluster size, pixels",im));		

		firstComptonClusterSizeSelEvents[im] = new TH1F(Form("firstComptonClusterSizeSelEvents_AM%d",im),Form("%s, cluster size of first Compton cluster, selected events, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		firstComptonClusterSizeSelEvents[im]->GetXaxis()->SetTitleOffset(1.4);
		firstComptonClusterSizeSelEvents[im]->GetXaxis()->SetTitle(Form("AM%d cluster size, pixels",im));		

		secondComptonClusterSizeSelEvents[im] = new TH1F(Form("secondComptonClusterSizeSelEvents_AM%d",im),Form("%s, cluster size of second Compton cluster, selected events, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		secondComptonClusterSizeSelEvents[im]->GetXaxis()->SetTitleOffset(1.4);
		secondComptonClusterSizeSelEvents[im]->GetXaxis()->SetTitle(Form("AM%d cluster size, pixels",im));		
		
		comptonHitsClusters_nTrigsCorr[im] = new TH2F(Form("comptonHitsClusters_nTrigsCorr_AM%d",im),Form("%s, unsummed Compton interactions, AM%d",spectrumName.Data(),im),maxNTrigs4H_2D-minNTrigs4H_2D,minNTrigs4H_2D,maxNTrigs4H_2D,maxNTrigs4H_2D-minNTrigs4H_2D,minNTrigs4H_2D,maxNTrigs4H_2D);
		comptonHitsClusters_nTrigsCorr[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonHitsClusters_nTrigsCorr[im]->GetXaxis()->SetTitle(Form("Number of triggers in first cluster, AM%d",im));	
		comptonHitsClusters_nTrigsCorr[im]->GetYaxis()->SetTitleOffset(1.4);
		comptonHitsClusters_nTrigsCorr[im]->GetYaxis()->SetTitle(Form("Number of triggers in second cluster, AM%d",im));		
		
		comptonSummedSpec2Clusters[im] = new TH1F(Form("comptonSummedSpec2Clusters_AM%d",im),Form("%s, two interaction summed Compton spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		comptonSummedSpec2Clusters[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonSummedSpec2Clusters[im]->GetXaxis()->SetTitle(Form("AM%d summed cluster energy, keV",im));		
		
		comptonSummedSpecEventClusters[im] = new TH1F(Form("comptonSummedSpecEventClusters_AM%d",im),Form("%s, summed Compton spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		comptonSummedSpecEventClusters[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonSummedSpecEventClusters[im]->GetXaxis()->SetTitle(Form("AM%d summed cluster energy, keV",im));		
		
		phiAngle[im] = new TH1F(Form("phiAngle_AM%d",im),Form("%s, cluster azimuthal angle, AM%d",spectrumName.Data(),im),360,0,360);
		phiAngle[im]->GetXaxis()->SetTitleOffset(1.2);
		phiAngle[im]->GetXaxis()->SetTitle(Form("#varphi AM%d, degrees",im));

		phiAngleSelEvents[im] = new TH1F(Form("phiAngleSelEvents_AM%d",im),Form("%s, cluster azimuthal angle, selected events, AM%d",spectrumName.Data(),im),360,0,360);
		phiAngleSelEvents[im]->GetXaxis()->SetTitleOffset(1.2);
		phiAngleSelEvents[im]->GetXaxis()->SetTitle(Form("#varphi AM%d, degrees",im));
		
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
		
		allClustersCOGImage[im] = new TH2F(Form("allClustersCOGImage_AM%d",im),Form("All clusters COG positions, integral energy image, AM%d",im),nPixXY,0,nPixXY,nPixXY,0,nPixXY);
		allClustersCOGImage[im]->GetXaxis()->SetTitle("Pixel coordinate");
		allClustersCOGImage[im]->GetXaxis()->SetTitleOffset(1.2);
		allClustersCOGImage[im]->GetYaxis()->SetTitle("Pixel coordinate");
		allClustersCOGImage[im]->GetYaxis()->SetTitleOffset(1.4);
		
		allClustersFineCOGImage[im] = new TH2F(Form("allClustersFineCOGImage_AM%d",im),Form("All clusters fine COG positions, integral energy image, AM%d",im),nPixXY*10,0,nPixXY,nPixXY*10,0,nPixXY);
		allClustersFineCOGImage[im]->GetXaxis()->SetTitle("Pixel coordinate");
		allClustersFineCOGImage[im]->GetXaxis()->SetTitleOffset(1.2);
		allClustersFineCOGImage[im]->GetYaxis()->SetTitle("Pixel coordinate");
		allClustersFineCOGImage[im]->GetYaxis()->SetTitleOffset(1.4);
		
		allPixelsImage[im] = new TH2F(Form("allPixelsImage_AM%d",im),Form("Integral energy image of all single pixel signals, AM%d",im),nPixXY,0,nPixXY,nPixXY,0,nPixXY);
		allPixelsImage[im]->GetXaxis()->SetTitle("Pixel coordinate");
		allPixelsImage[im]->GetXaxis()->SetTitleOffset(1.2);
		allPixelsImage[im]->GetYaxis()->SetTitle("Pixel coordinate");
		allPixelsImage[im]->GetYaxis()->SetTitleOffset(1.4);
		
		allPixelsFreqImage[im] = new TH2F(Form("allPixelsFreqImage_AM%d",im),Form("Image of all single pixel frequencies, AM%d",im),nPixXY,0,nPixXY,nPixXY,0,nPixXY);
		allPixelsFreqImage[im]->GetXaxis()->SetTitle("Pixel coordinate");
		allPixelsFreqImage[im]->GetXaxis()->SetTitleOffset(1.2);
		allPixelsFreqImage[im]->GetYaxis()->SetTitle("Pixel coordinate");
		allPixelsFreqImage[im]->GetYaxis()->SetTitleOffset(1.4);
		
		allClustersCOGFreqImage[im] = new TH2F(Form("allClustersCOGFreqImage_AM%d",im),Form("All clusters COG position frequencies, AM%d",im),nPixXY,0,nPixXY,nPixXY,0,nPixXY);
		allClustersCOGFreqImage[im]->GetXaxis()->SetTitle("Pixel coordinate");
		allClustersCOGFreqImage[im]->GetXaxis()->SetTitleOffset(1.2);
		allClustersCOGFreqImage[im]->GetYaxis()->SetTitle("Pixel coordinate");
		allClustersCOGFreqImage[im]->GetYaxis()->SetTitleOffset(1.4);
		
		allClustersFineCOGFreqImage[im] = new TH2F(Form("allClustersFineCOGFreqImage_AM%d",im),Form("All clusters fine COG position frequencies, AM%d",im),nPixXY*10,0,nPixXY,nPixXY*10,0,nPixXY);
		allClustersFineCOGFreqImage[im]->GetXaxis()->SetTitle("Pixel coordinate");
		allClustersFineCOGFreqImage[im]->GetXaxis()->SetTitleOffset(1.2);
		allClustersFineCOGFreqImage[im]->GetYaxis()->SetTitle("Pixel coordinate");
		allClustersFineCOGFreqImage[im]->GetYaxis()->SetTitleOffset(1.4);
		
		secondComptonCluster[im] = new TH2F(Form("secondComptonCluster_AM%d",im),Form("Second Compton scattering cluster positions, AM%d",im),nPixXY,0,nPixXY,nPixXY,0,nPixXY);
		secondComptonCluster[im]->GetXaxis()->SetTitle("Pixel coordinate");
		secondComptonCluster[im]->GetXaxis()->SetTitleOffset(1.2);
		secondComptonCluster[im]->GetYaxis()->SetTitle("Pixel coordinate");
		secondComptonCluster[im]->GetYaxis()->SetTitleOffset(1.4);
	}
	
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
	
	comptonSummedSpec2ClustersCorr = new TH2F("comptonSummedSpec2ClustersCorr",Form("%s, two interaction summed Compton spectra in AMs",spectrumName.Data()),nBinsE,minE,maxE,nBinsE,minE,maxE);
	comptonSummedSpec2ClustersCorr->GetXaxis()->SetTitleOffset(1.2);
	comptonSummedSpec2ClustersCorr->GetXaxis()->SetTitle("AM0 summed cluster energy, keV");
	comptonSummedSpec2ClustersCorr->GetYaxis()->SetTitleOffset(1.4);
	comptonSummedSpec2ClustersCorr->GetYaxis()->SetTitle("AM1 summed cluster energy, keV");
	
	comptonSummedSpecEventCorr = new TH2F("comptonSummedSpecEventCorr",Form("%s, summed Compton spectra in AMs",spectrumName.Data()),nBinsE,minE,maxE,nBinsE,minE,maxE);
	comptonSummedSpecEventCorr->GetXaxis()->SetTitleOffset(1.2);
	comptonSummedSpecEventCorr->GetXaxis()->SetTitle("AM0 summed cluster energy, keV");
	comptonSummedSpecEventCorr->GetYaxis()->SetTitleOffset(1.4);
	comptonSummedSpecEventCorr->GetYaxis()->SetTitle("AM1 summed cluster energy, keV");
	
	nClustersInEventCorr = new TH2F("nClustersInEventCorr",Form("%s, number of clusters/event in AMs",spectrumName.Data()),maxNTrigs4H_2D-minNTrigs4H_2D,minNTrigs4H_2D,maxNTrigs4H_2D,maxNTrigs4H_2D-minNTrigs4H_2D,minNTrigs4H_2D,maxNTrigs4H_2D);
	nClustersInEventCorr->GetXaxis()->SetTitleOffset(1.2);
	nClustersInEventCorr->GetXaxis()->SetTitle("Number of clusters in AM0");
	nClustersInEventCorr->GetYaxis()->SetTitleOffset(1.4);
	nClustersInEventCorr->GetYaxis()->SetTitle("Number of clusters in AM1");

	nClustersInSelEventCorr = new TH2F("nClustersInSelEventCorr",Form("%s, number of clusters/event in selected events in AMs",spectrumName.Data()),maxNTrigs4H_2D-minNTrigs4H_2D,minNTrigs4H_2D,maxNTrigs4H_2D,maxNTrigs4H_2D-minNTrigs4H_2D,minNTrigs4H_2D,maxNTrigs4H_2D);
	nClustersInSelEventCorr->GetXaxis()->SetTitleOffset(1.2);
	nClustersInSelEventCorr->GetXaxis()->SetTitle("Number of clusters in AM0");
	nClustersInSelEventCorr->GetYaxis()->SetTitleOffset(1.4);
	nClustersInSelEventCorr->GetYaxis()->SetTitle("Number of clusters in AM1");
	
	dPhiAngle = new TH1F("dPhiAngle",Form("%s, #Delta#varphi",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiAngle->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngleNorm = new TH1F("dPhiAngleNorm",Form("%s, #Delta#varphi, normalised by #Delta#varphi = 0",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiAngleNorm->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngleNorm->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle1 = new TH1F("dPhiAngle1",Form("%s, |#Delta#varphi|",spectrumName.Data()),nBins_dPhi/2,0,180);
	dPhiAngle1->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle1->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle_Ewin = new TH1F("dPhiAngle_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV",spectrumName.Data(),EwindowFor4PhiAnalyis_min,EwindowFor4PhiAnalyis_max),nBins_dPhi,-180,180);
	dPhiAngle_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngleNorm_Ewin = new TH1F("dPhiAngleNorm_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, normalised by #Delta#varphi = 0",spectrumName.Data(),EwindowFor4PhiAnalyis_min,EwindowFor4PhiAnalyis_max),nBins_dPhi,-180,180);
	dPhiAngleNorm_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngleNorm_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle1_Ewin = new TH1F("dPhiAngle1_Ewin",Form("%s, |#Delta#varphi| within  %.0f keV < E < %.0f keV",spectrumName.Data(),EwindowFor4PhiAnalyis_min,EwindowFor4PhiAnalyis_max),nBins_dPhi/2,0,180);
	dPhiAngle1_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle1_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle_2ClusterE = new TH2F("dPhiAngle_2ClusterE",Form("%s, #Delta#varphi vs summed energy of two first Compton clusters",spectrumName.Data()),nBinsE,minE,maxE,nBins_dPhi,-180,180);
	dPhiAngle_2ClusterE->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_2ClusterE->GetXaxis()->SetTitle("Total summed energy, keV");
	dPhiAngle_2ClusterE->GetYaxis()->SetTitleOffset(1.4);
	dPhiAngle_2ClusterE->GetYaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle_totalE = new TH2F("dPhiAngle_totalE",Form("%s, #Delta#varphi vs total energy of Compton events",spectrumName.Data()),nBinsE,minE,maxE,nBins_dPhi,-180,180);
	dPhiAngle_totalE->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_totalE->GetXaxis()->SetTitle("Total energy, keV");
	dPhiAngle_totalE->GetYaxis()->SetTitleOffset(1.4);
	dPhiAngle_totalE->GetYaxis()->SetTitle("#Delta#varphi, degrees");		
	
	phiAngleCorr = new TH2F("phiAngleCorr",Form("%s, #varphi angle",spectrumName.Data()),36,0,360,36,0,360);
	phiAngleCorr->GetXaxis()->SetTitleOffset(1.2);
	phiAngleCorr->GetXaxis()->SetTitle("#varphi AM0, degrees");	
	phiAngleCorr->GetYaxis()->SetTitleOffset(1.4);
	phiAngleCorr->GetYaxis()->SetTitle("#varphi AM1, degrees");	
	
	phiAngleSelEventsCorr = new TH2F("phiAngleSelEventsCorr",Form("%s, #varphi angle, selected events",spectrumName.Data()),36,0,360,36,0,360);
	phiAngleSelEventsCorr->GetXaxis()->SetTitleOffset(1.2);
	phiAngleSelEventsCorr->GetXaxis()->SetTitle("#varphi AM0, degrees");	
	phiAngleSelEventsCorr->GetYaxis()->SetTitleOffset(1.4);
	phiAngleSelEventsCorr->GetYaxis()->SetTitle("#varphi AM1, degrees");	
	
	theta0_theta1_FromFirstClusterE = new TH2F("theta0_theta1_FromFirstClusterE",Form("%s, #theta_{0} - #theta_{1} calculated from first interaction energy",spectrumName.Data()),nBinsTheta4H,0,180,nBinsTheta4H,0,180);
	theta0_theta1_FromFirstClusterE->GetXaxis()->SetTitleOffset(1.2);
	theta0_theta1_FromFirstClusterE->GetXaxis()->SetTitle("#theta in AM0, degrees");
	theta0_theta1_FromFirstClusterE->GetYaxis()->SetTitleOffset(1.4);
	theta0_theta1_FromFirstClusterE->GetYaxis()->SetTitle("#theta in AM1, degrees");

	matrix_trigs.clear();
	matrix_trigs.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_flags.clear();
	matrix_flags.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_E.clear();
	matrix_E.resize(nPixXY, vector<Float_t>(nPixXY));
	matrix_Eneg.clear();
	matrix_Eneg.resize(nPixXY, vector<Float_t>(nPixXY));
	matrix_time.clear();
	matrix_time.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_timeCalib.clear();
	matrix_timeCalib.resize(nPixXY, vector<Float_t>(nPixXY));

	UInt_t screen_width = screen_width_def;
	UInt_t screen_height = screen_height_def;
	Float_t size_factor = 0.95;
	Float_t size_ratio = 1.333;
	const Int_t ix = 3;
	const Int_t iy = 3;
	Int_t curr_Pad = 0;
	
	if (enableEventMonitor)
	{
		c3 = new TCanvas("c3","",10,10,Int_t(screen_height*size_factor*2),Int_t(screen_height*size_factor));
		curr_Pad = 0;
		c3->Divide(2,1);
		for (curr_Pad = 0; curr_Pad <= 2; curr_Pad++)
		{
			c3->cd(curr_Pad);
			c3->GetPad(curr_Pad)->SetGridx();
			c3->GetPad(curr_Pad)->SetGridy();
			c3->GetPad(curr_Pad)->SetLeftMargin(0.12);
			c3->GetPad(curr_Pad)->SetRightMargin(0.12);
			c3->GetPad(curr_Pad)->SetTopMargin(0.11);
			c3->GetPad(curr_Pad)->SetBottomMargin(0.13);
		}
	}	
	line2 = new TLine();
	line2->SetLineColor(1);
	line2->SetLineWidth(3);
	
	if (nEvents2Do > 0) nEvents2Analyse = nEvents2Do;
	localTimer.Start();
	for (Int_t ie = 0; ie < nEvents2Analyse; ie++)
	{
		if (saveSelectedComptonRootTree) selEvt[ie] = 0;
		chain_events->GetEntry(ie);
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
		evtTreeIdx = ie;
		if (!analyseNextEvent(event_eev))
		{
			//cerr << "ERROR analysing event " << ie << ". Exiting." << endl;
			//return 0;
		}
	}
	
	if (saveSelectedComptonRootTree)
	{
		TFile *hfileOutRoot = new TFile(Form("%s%s/newTrees_pt%d_Ecal.root",pathRoot.Data(),subdir.Data(),savingPointLocation),"recreate");
		TTree* chain_events_new = chain_events->CloneTree(0);
		TTree **chain1p_new = new TTree*[nAMs];
		
		for (Long_t ie = 0; ie < nEvents2Analyse; ie++)
		{
			chain_events->GetEntry(ie);
			if (selEvt[ie] < 1) continue;
			event_eev = selEvt[ie];
			chain_events_new->Fill();
		}
		chain_events_new->Write();
		cout << "Saved " << cntEv2BSaved << " events in the new selected event Tree." << endl;
		cout << Form("RealTime: %.3fs, CPUTime: %.3fs",localTimer.RealTime(), localTimer.CpuTime()) << endl;
		
		for (Int_t im = 0; im < nAMs; im++)
		{
			chain1p_new[im] = chain1p[im]->CloneTree(0);
			for (Long_t ie = 0; ie < Nev[im]; ie++)
			{
				if (selEvtAM[im][ie] < 1) continue;
				chain1p[im]->GetEntry(ie);
				event[im] = selEvtAM[im][ie];
				chain1p_new[im]->Fill();
			}
			chain1p_new[im]->Write();
			cout << Form("Saved AM%d data, RealTime: %.3fs, CPUTime: %.3fs",im,localTimer.RealTime(), localTimer.CpuTime()) << endl;
		}
		hfileOutRoot->Close();
	}
	
	TCanvas *c0 = new TCanvas("c0","",10,10,Int_t(screen_height*size_factor*size_ratio),Int_t(screen_height*size_factor));
	curr_Pad = 0;
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
	TH1F **nClustersInSelEvent_norm = new TH1F*[nAMs];
	TH1F **nClustersInEvent_norm = new TH1F*[nAMs];
	TH1F **firstComptonClusterSize_norm = new TH1F*[nAMs];
	TH1F **secondComptonClusterSize_norm = new TH1F*[nAMs];
	TH1F **firstComptonClusterSizeSelEvents_norm = new TH1F*[nAMs];
	TH1F **secondComptonClusterSizeSelEvents_norm = new TH1F*[nAMs];

	for (Int_t im = 0; im < nAMs; im++)
	{
		if (enableEventMonitor)
		{
			c3->cd(im+1);
			pixelPattern[im]->Draw("text");
			line2->DrawLine(pixelPattern[im]->GetXaxis()->GetXmin(),pixelPattern[im]->GetYaxis()->GetXmax()/2,pixelPattern[im]->GetXaxis()->GetXmax(),pixelPattern[im]->GetYaxis()->GetXmax()/2);
			line2->DrawLine(pixelPattern[im]->GetXaxis()->GetXmax()/2,pixelPattern[im]->GetYaxis()->GetXmin(),pixelPattern[im]->GetXaxis()->GetXmax()/2,pixelPattern[im]->GetYaxis()->GetXmax());
		}
		
		c2->cd(0);
		pixelPattern[im]->Draw("text");
		line2->DrawLine(pixelPattern[im]->GetXaxis()->GetXmin(),pixelPattern[im]->GetYaxis()->GetXmax()/2,pixelPattern[im]->GetXaxis()->GetXmax(),pixelPattern[im]->GetYaxis()->GetXmax()/2);
		line2->DrawLine(pixelPattern[im]->GetXaxis()->GetXmax()/2,pixelPattern[im]->GetYaxis()->GetXmin(),pixelPattern[im]->GetXaxis()->GetXmax()/2,pixelPattern[im]->GetYaxis()->GetXmax());
		c2->SaveAs(Form("%s/pixelPattern.gif",outputPathAM[im].Data()));
					
		comptonSpec2ClustersCorr[im]->Draw("colz");
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(comptonSpec2ClustersCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(comptonSpec2ClustersCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
	
		comptonHitsSpecClustersSizeCorr[im]->Draw("colz");
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(comptonHitsSpecClustersSizeCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(comptonHitsSpecClustersSizeCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
		
		comptonHitsClusters_nTrigsCorr[im]->Draw("colz");
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(comptonHitsClusters_nTrigsCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(comptonHitsClusters_nTrigsCorr[im]->GetName()).Data()));
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
					
		firstComptonCluster[im]->Draw("colz");
		line2->DrawLine(firstComptonCluster[im]->GetXaxis()->GetXmin(),firstComptonCluster[im]->GetYaxis()->GetXmax()/2,firstComptonCluster[im]->GetXaxis()->GetXmax(),firstComptonCluster[im]->GetYaxis()->GetXmax()/2);
		line2->DrawLine(firstComptonCluster[im]->GetXaxis()->GetXmax()/2,firstComptonCluster[im]->GetYaxis()->GetXmin(),firstComptonCluster[im]->GetXaxis()->GetXmax()/2,firstComptonCluster[im]->GetYaxis()->GetXmax());
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(firstComptonCluster[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(firstComptonCluster[im]->GetName()).Data()));
		pixelPattern[im]->Draw("sametext");
		c2->SaveAs(Form("%s/%s_logz_pix.png",outputPathAM[im].Data(),TString(firstComptonCluster[im]->GetName()).Data()));		
		c2->GetPad(0)->SetLogz(0);
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(firstComptonCluster[im]->GetName()).Data()));		
		
		secondComptonCluster[im]->Draw("colz");
		line2->DrawLine(secondComptonCluster[im]->GetXaxis()->GetXmin(),secondComptonCluster[im]->GetYaxis()->GetXmax()/2,secondComptonCluster[im]->GetXaxis()->GetXmax(),secondComptonCluster[im]->GetYaxis()->GetXmax()/2);
		line2->DrawLine(secondComptonCluster[im]->GetXaxis()->GetXmax()/2,secondComptonCluster[im]->GetYaxis()->GetXmin(),secondComptonCluster[im]->GetXaxis()->GetXmax()/2,secondComptonCluster[im]->GetYaxis()->GetXmax());
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(secondComptonCluster[im]->GetName()).Data()));
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(secondComptonCluster[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(secondComptonCluster[im]->GetName()).Data()));
		pixelPattern[im]->Draw("sametext");
		c2->SaveAs(Form("%s/%s_logz_pix.png",outputPathAM[im].Data(),TString(secondComptonCluster[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(secondComptonCluster[im]->GetName()).Data()));
		
		allClustersCOGImage[im]->Draw("colz");
		line2->DrawLine(allClustersCOGImage[im]->GetXaxis()->GetXmin(),allClustersCOGImage[im]->GetYaxis()->GetXmax()/2,allClustersCOGImage[im]->GetXaxis()->GetXmax(),allClustersCOGImage[im]->GetYaxis()->GetXmax()/2);
		line2->DrawLine(allClustersCOGImage[im]->GetXaxis()->GetXmax()/2,allClustersCOGImage[im]->GetYaxis()->GetXmin(),allClustersCOGImage[im]->GetXaxis()->GetXmax()/2,allClustersCOGImage[im]->GetYaxis()->GetXmax());
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(allClustersCOGImage[im]->GetName()).Data()));
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allClustersCOGImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(allClustersCOGImage[im]->GetName()).Data()));
		pixelPattern[im]->Draw("sametext");
		c2->SaveAs(Form("%s/%s_logz_pix.png",outputPathAM[im].Data(),TString(allClustersCOGImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allClustersCOGImage[im]->GetName()).Data()));
		
		allPixelsImage[im]->Draw("colz");
		line2->DrawLine(allPixelsImage[im]->GetXaxis()->GetXmin(),allPixelsImage[im]->GetYaxis()->GetXmax()/2,allPixelsImage[im]->GetXaxis()->GetXmax(),allPixelsImage[im]->GetYaxis()->GetXmax()/2);
		line2->DrawLine(allPixelsImage[im]->GetXaxis()->GetXmax()/2,allPixelsImage[im]->GetYaxis()->GetXmin(),allPixelsImage[im]->GetXaxis()->GetXmax()/2,allPixelsImage[im]->GetYaxis()->GetXmax());
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(allPixelsImage[im]->GetName()).Data()));
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allPixelsImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(allPixelsImage[im]->GetName()).Data()));
		pixelPattern[im]->Draw("sametext");
		c2->SaveAs(Form("%s/%s_logz_pix.png",outputPathAM[im].Data(),TString(allPixelsImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allPixelsImage[im]->GetName()).Data()));
		
		allClustersFineCOGImage[im]->Draw("colz");
		line2->DrawLine(allClustersCOGImage[im]->GetXaxis()->GetXmin(),allClustersCOGImage[im]->GetYaxis()->GetXmax()/2,allClustersCOGImage[im]->GetXaxis()->GetXmax(),allClustersCOGImage[im]->GetYaxis()->GetXmax()/2);
		line2->DrawLine(allClustersCOGImage[im]->GetXaxis()->GetXmax()/2,allClustersCOGImage[im]->GetYaxis()->GetXmin(),allClustersCOGImage[im]->GetXaxis()->GetXmax()/2,allClustersCOGImage[im]->GetYaxis()->GetXmax());
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(allClustersFineCOGImage[im]->GetName()).Data()));
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allClustersFineCOGImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(allClustersFineCOGImage[im]->GetName()).Data()));
		pixelPattern[im]->Draw("sametext");
		c2->SaveAs(Form("%s/%s_logz_pix.png",outputPathAM[im].Data(),TString(allClustersFineCOGImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allClustersFineCOGImage[im]->GetName()).Data()));
		
		allPixelsFreqImage[im]->Draw("colz");
		line2->DrawLine(allPixelsImage[im]->GetXaxis()->GetXmin(),allPixelsImage[im]->GetYaxis()->GetXmax()/2,allPixelsImage[im]->GetXaxis()->GetXmax(),allPixelsImage[im]->GetYaxis()->GetXmax()/2);
		line2->DrawLine(allPixelsImage[im]->GetXaxis()->GetXmax()/2,allPixelsImage[im]->GetYaxis()->GetXmin(),allPixelsImage[im]->GetXaxis()->GetXmax()/2,allPixelsImage[im]->GetYaxis()->GetXmax());
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(allPixelsFreqImage[im]->GetName()).Data()));
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allPixelsFreqImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(allPixelsFreqImage[im]->GetName()).Data()));
		pixelPattern[im]->Draw("sametext");
		c2->SaveAs(Form("%s/%s_logz_pix.png",outputPathAM[im].Data(),TString(allPixelsFreqImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allPixelsFreqImage[im]->GetName()).Data()));
		
		allClustersFineCOGFreqImage[im]->Draw("colz");
		line2->DrawLine(allClustersCOGImage[im]->GetXaxis()->GetXmin(),allClustersCOGImage[im]->GetYaxis()->GetXmax()/2,allClustersCOGImage[im]->GetXaxis()->GetXmax(),allClustersCOGImage[im]->GetYaxis()->GetXmax()/2);
		line2->DrawLine(allClustersCOGImage[im]->GetXaxis()->GetXmax()/2,allClustersCOGImage[im]->GetYaxis()->GetXmin(),allClustersCOGImage[im]->GetXaxis()->GetXmax()/2,allClustersCOGImage[im]->GetYaxis()->GetXmax());
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(allClustersFineCOGFreqImage[im]->GetName()).Data()));
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allClustersFineCOGFreqImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(allClustersFineCOGFreqImage[im]->GetName()).Data()));
		pixelPattern[im]->Draw("sametext");
		c2->SaveAs(Form("%s/%s_logz_pix.png",outputPathAM[im].Data(),TString(allClustersFineCOGFreqImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allClustersFineCOGFreqImage[im]->GetName()).Data()));
		
		allClustersCOGFreqImage[im]->Draw("colz");
		line2->DrawLine(allClustersCOGImage[im]->GetXaxis()->GetXmin(),allClustersCOGImage[im]->GetYaxis()->GetXmax()/2,allClustersCOGImage[im]->GetXaxis()->GetXmax(),allClustersCOGImage[im]->GetYaxis()->GetXmax()/2);
		line2->DrawLine(allClustersCOGImage[im]->GetXaxis()->GetXmax()/2,allClustersCOGImage[im]->GetYaxis()->GetXmin(),allClustersCOGImage[im]->GetXaxis()->GetXmax()/2,allClustersCOGImage[im]->GetYaxis()->GetXmax());
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(allClustersCOGFreqImage[im]->GetName()).Data()));
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allClustersCOGFreqImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(allClustersCOGFreqImage[im]->GetName()).Data()));
		pixelPattern[im]->Draw("sametext");
		c2->SaveAs(Form("%s/%s_logz_pix.png",outputPathAM[im].Data(),TString(allClustersCOGFreqImage[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
		c2->SaveAs(Form("%s/%s_pix.png",outputPathAM[im].Data(),TString(allClustersCOGFreqImage[im]->GetName()).Data()));
		
		c1->cd(0);
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
		
		cathodeSpecAllEvents[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(cathodeSpecAllEvents[im]->GetName()).Data()));
		
		cathodeSpecPair[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(cathodeSpecPair[im]->GetName()).Data()));
		
		cathodeSpecSelPairs[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(cathodeSpecSelPairs[im]->GetName()).Data()));
		
		comptonSpecClusters[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSpecClusters[im]->GetName()).Data()));
		
		nClustersInEvent[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(nClustersInEvent[im]->GetName()).Data()));
		nClustersInEvent_norm[im] = (TH1F *) nClustersInEvent[im]->Clone();
		nClustersInEvent_norm[im]->SetName(Form("nClustersInEvent_norm_AM%d",im));
		nClustersInEvent_norm[im]->SetNormFactor(1);
		nClustersInEvent_norm[im]->GetYaxis()->SetTitle("Normalised number of clusters");
		nClustersInEvent_norm[im]->GetYaxis()->SetTitleOffset(1.3);
		nClustersInEvent_norm[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(nClustersInEvent_norm[im]->GetName()).Data()));
		
		nClustersInSelEvent[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(nClustersInSelEvent[im]->GetName()).Data()));
		nClustersInSelEvent_norm[im] = (TH1F *) nClustersInSelEvent[im]->Clone();
		nClustersInSelEvent_norm[im]->SetName(Form("nClustersInSelEvent_norm_AM%d",im));
		nClustersInSelEvent_norm[im]->SetNormFactor(1);
		nClustersInSelEvent_norm[im]->GetYaxis()->SetTitle("Normalised number of clusters");
		nClustersInSelEvent_norm[im]->GetYaxis()->SetTitleOffset(1.3);
		nClustersInSelEvent_norm[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(nClustersInSelEvent_norm[im]->GetName()).Data()));
		
		nTrigsInEvent[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(nTrigsInEvent[im]->GetName()).Data()));
		nTrigsInEvent_norm[im] = (TH1F *) nTrigsInEvent[im]->Clone();
		nTrigsInEvent_norm[im]->SetName(Form("nTrigsInEvent_norm_AM%d",im));
		nTrigsInEvent_norm[im]->SetNormFactor(1);
		nTrigsInEvent_norm[im]->GetYaxis()->SetTitle("Normalised counts");
		nTrigsInEvent_norm[im]->GetYaxis()->SetTitleOffset(1.3);
		nTrigsInEvent_norm[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(nTrigsInEvent_norm[im]->GetName()).Data()));
		
		nTrigs_clusterSize[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(nTrigs_clusterSize[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(nTrigs_clusterSize[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		if (makeRawTimingStuff) 
		{
			c2->cd(0);
			rawTimingTwoPixelClustersCorr[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(rawTimingTwoPixelClustersCorr[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(rawTimingTwoPixelClustersCorr[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			c1->cd(0);
			rawTimingDiff2PixClusters[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(rawTimingDiff2PixClusters[im]->GetName()).Data()));

			rawTimingDiff3PixClusters[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(rawTimingDiff3PixClusters[im]->GetName()).Data()));
			
			rawTimingDiffAllClusters[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(rawTimingDiffAllClusters[im]->GetName()).Data()));
			
			if (!doNotUseCornerPixelsInPixelClusterCOGEneg) 
			{
				rawTimingDiff2PixDiagClusters[im]->Draw();
				c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(rawTimingDiff2PixDiagClusters[im]->GetName()).Data()));
			}
		}
		if (makeTimingCalibrationStuff)
		{			
			c2->cd(0);
			anodeTimingTwoPixelClustersCorr[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(anodeTimingTwoPixelClustersCorr[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(anodeTimingTwoPixelClustersCorr[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			if (!doNotUseCornerPixelsInPixelClusterCOGEneg) 
			{
				anodeTiming2PixDiagClustersCorr[im]->Draw("colz");
				c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(anodeTiming2PixDiagClustersCorr[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(1);
				c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(anodeTiming2PixDiagClustersCorr[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(0);
				
				anodeTimingDiff2PixDiagClusters[im]->Draw();
				c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(anodeTimingDiff2PixDiagClusters[im]->GetName()).Data()));
			}
			
			anodeTiming2PixClustersCorr[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(anodeTiming2PixClustersCorr[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(anodeTiming2PixClustersCorr[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			anodeTiming3PixClustersCorr[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(anodeTiming3PixClustersCorr[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(anodeTiming3PixClustersCorr[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			anodeTimingTwoClustersAMCorr[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(anodeTimingTwoClustersAMCorr[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(anodeTimingTwoClustersAMCorr[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			dAnodeTiming_vs_theta[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_theta[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_theta[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			dAnodeTiming_vs_E2[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_E2[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_E2[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			dAnodeTiming_vs_E12[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_E12[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_E12[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			dAnodeTimingAbs_vs_E12[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dAnodeTimingAbs_vs_E12[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dAnodeTimingAbs_vs_E12[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			dist12_vs_E1[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dist12_vs_E1[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dist12_vs_E1[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			dist12_vs_E2[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dist12_vs_E2[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dist12_vs_E2[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			dist12_vs_E12[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dist12_vs_E12[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dist12_vs_E12[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			dist12_over_E12sq[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dist12_over_E12sq[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dist12_over_E12sq[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
			
			dAnodeTiming_vs_E1[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_E1[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_E1[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
						
			c1->cd(0);
			anodeTimingDiff2PixClusters[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(anodeTimingDiff2PixClusters[im]->GetName()).Data()));

			anodeTimingDiff3PixClusters[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(anodeTimingDiff3PixClusters[im]->GetName()).Data()));
			
			anodeTimingDiffAllClusters[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(anodeTimingDiffAllClusters[im]->GetName()).Data()));
		}
		
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
		comptonSummedSpec2Clusters[im]->Draw();
		Y_text_top = 0.82;
		getHistoPeak(comptonSummedSpec2Clusters[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(comptonSummedSpec2Clusters[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpec2Clusters[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSummedSpec2Clusters[im]->GetName()).Data()));
		
		comptonSummedSpecEventClusters[im]->Draw();
		Y_text_top = 0.82;
		getHistoPeak(comptonSummedSpecEventClusters[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(comptonSummedSpecEventClusters[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpecEventClusters[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSummedSpecEventClusters[im]->GetName()).Data()));
		
		phiAngle[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(phiAngle[im]->GetName()).Data()));
		
		phiAngleSelEvents[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(phiAngleSelEvents[im]->GetName()).Data()));
		
		thetaFromFirstClusterE_dPhi[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_dPhi[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_dPhi[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		thetaFromFirstClusterE[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(thetaFromFirstClusterE[im]->GetName()).Data()));
		
		thetaFromFirstClusterE_w[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_w[im]->GetName()).Data()));
		
		thetaFromFirstClusterE_1pixOnly[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_1pixOnly[im]->GetName()).Data()));
		
		thetaFromFirstClusterE_1pixOnly_w[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_1pixOnly_w[im]->GetName()).Data()));
		
		thetaFromFirstClusterE_2pixOnly[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_2pixOnly[im]->GetName()).Data()));
		
		thetaFromFirstClusterE_1_2pixOnly[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_1_2pixOnly[im]->GetName()).Data()));
		
		clusterSize[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(clusterSize[im]->GetName()).Data()));
		clusterSize_norm[im] = (TH1F *) clusterSize[im]->Clone();
		clusterSize_norm[im]->SetName(Form("clusterSize_norm_AM%d",im));
		clusterSize_norm[im]->SetNormFactor(1);
		clusterSize_norm[im]->GetYaxis()->SetTitle("Normalised counts");
		clusterSize_norm[im]->GetYaxis()->SetTitleOffset(1.3);
		clusterSize_norm[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(clusterSize_norm[im]->GetName()).Data()));
		
		firstComptonClusterSize[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(firstComptonClusterSize[im]->GetName()).Data()));
		firstComptonClusterSize_norm[im] = (TH1F *) firstComptonClusterSize[im]->Clone();
		firstComptonClusterSize_norm[im]->SetName(Form("firstComptonClusterSize_norm_AM%d",im));
		firstComptonClusterSize_norm[im]->SetNormFactor(1);
		firstComptonClusterSize_norm[im]->GetYaxis()->SetTitle("Normalised counts");
		firstComptonClusterSize_norm[im]->GetYaxis()->SetTitleOffset(1.3);
		firstComptonClusterSize_norm[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(firstComptonClusterSize_norm[im]->GetName()).Data()));
		
		secondComptonClusterSize[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(secondComptonClusterSize[im]->GetName()).Data()));
		secondComptonClusterSize_norm[im] = (TH1F *) secondComptonClusterSize[im]->Clone();
		secondComptonClusterSize_norm[im]->SetName(Form("secondComptonClusterSize_norm_AM%d",im));
		secondComptonClusterSize_norm[im]->SetNormFactor(1);
		secondComptonClusterSize_norm[im]->GetYaxis()->SetTitle("Normalised counts");
		secondComptonClusterSize_norm[im]->GetYaxis()->SetTitleOffset(1.3);
		secondComptonClusterSize_norm[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(secondComptonClusterSize_norm[im]->GetName()).Data()));
		
		firstComptonClusterSizeSelEvents[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(firstComptonClusterSizeSelEvents[im]->GetName()).Data()));
		firstComptonClusterSizeSelEvents_norm[im] = (TH1F *) firstComptonClusterSizeSelEvents[im]->Clone();
		firstComptonClusterSizeSelEvents_norm[im]->SetName(Form("firstComptonClusterSizeSelEvents_norm_AM%d",im));
		firstComptonClusterSizeSelEvents_norm[im]->SetNormFactor(1);
		firstComptonClusterSizeSelEvents_norm[im]->GetYaxis()->SetTitle("Normalised counts");
		firstComptonClusterSizeSelEvents_norm[im]->GetYaxis()->SetTitleOffset(1.3);
		firstComptonClusterSizeSelEvents_norm[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(firstComptonClusterSizeSelEvents_norm[im]->GetName()).Data()));
		
		secondComptonClusterSizeSelEvents[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(secondComptonClusterSizeSelEvents[im]->GetName()).Data()));
		secondComptonClusterSizeSelEvents_norm[im] = (TH1F *) secondComptonClusterSizeSelEvents[im]->Clone();
		secondComptonClusterSizeSelEvents_norm[im]->SetName(Form("secondComptonClusterSizeSelEvents_norm_AM%d",im));
		secondComptonClusterSizeSelEvents_norm[im]->SetNormFactor(1);
		secondComptonClusterSizeSelEvents_norm[im]->GetYaxis()->SetTitle("Normalised counts");
		secondComptonClusterSizeSelEvents_norm[im]->GetYaxis()->SetTitleOffset(1.3);
		secondComptonClusterSizeSelEvents_norm[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(secondComptonClusterSizeSelEvents_norm[im]->GetName()).Data()));
	}
	if (enableEventMonitor) c3->SaveAs(Form("%s/pixelPatterns.gif",evMonDir.Data()));	
			
	if (plotPixelClusterSpectra)
	{
		for (Int_t im = 0; im < nAMs; im++)
		{
			c1->cd(0);
			for (Int_t b = firstPixel; b < lastPixel+1; b++)
			{
				if (clusterSpec[im][b]->GetEntries() < 10) continue; 
				clusterSpec[im][b]->Draw();
				Y_text_top = 0.82;
				getHistoPeak(clusterSpec[im][b], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
				Peak_FWHM_raw = getFWHM(clusterSpec[im][b], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
				txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
				Y_text_top -= Y_text_step;
				txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
				Y_text_top -= Y_text_step;
				txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
				Y_text_top -= Y_text_step;
				txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/clusterSpec[im][b]->GetBinWidth(1)));
				line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
				c1->SaveAs(Form("%s/pixel_%d.gif",outputPixelPathAM[im].Data(),b));
			}
		}
	}
	
	c1->cd(0);
	dPhiAngle->SetMinimum(0);
	dPhiAngle->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc = new TF1("fit",fitfun,dPhiAngle->GetXaxis()->GetXmin(),dPhiAngle->GetXaxis()->GetXmax(),2);
		dPhiAngle->Fit(ffunc,"RQ0");
		ffunc->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle->GetName()).Data()));
	
	dPhiAngleNorm->SetMinimum(0);
	Int_t b1 = dPhiAngleNorm->GetXaxis()->FindBin(-0.01);
	Int_t b2 = dPhiAngleNorm->GetXaxis()->FindBin(0.01);
	Float_t nf = (dPhiAngleNorm->GetBinContent(b1) + dPhiAngleNorm->GetBinContent(b2))/2;
	for (Int_t ib = 1; ib <= dPhiAngleNorm->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiAngleNorm->GetBinContent(ib)/nf;
		dPhiAngleNorm->SetBinContent(ib,bc);
	}
	dPhiAngleNorm->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc = new TF1("fit",fitfun,dPhiAngleNorm->GetXaxis()->GetXmin(),dPhiAngleNorm->GetXaxis()->GetXmax(),2);
		dPhiAngleNorm->Fit(ffunc,"RQ0");
		ffunc->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngleNorm->GetName()).Data()));

	dPhiAngle_Ewin->SetMinimum(0);
	dPhiAngle_Ewin->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc_Ewin = new TF1("fit_Ewin",fitfun,dPhiAngle_Ewin->GetXaxis()->GetXmin(),dPhiAngle_Ewin->GetXaxis()->GetXmax(),2);
		dPhiAngle_Ewin->Fit(ffunc_Ewin,"RQ0");
		ffunc_Ewin->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle_Ewin->GetName()).Data()));

	dPhiAngleNorm_Ewin->SetMinimum(0);
	Int_t b1w = dPhiAngleNorm_Ewin->GetXaxis()->FindBin(-0.01);
	Int_t b2w = dPhiAngleNorm_Ewin->GetXaxis()->FindBin(0.01);
	Float_t nf_w = (dPhiAngleNorm_Ewin->GetBinContent(b1w) + dPhiAngleNorm_Ewin->GetBinContent(b2w))/2;
	for (Int_t ib = 1; ib <= dPhiAngleNorm_Ewin->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiAngleNorm_Ewin->GetBinContent(ib)/nf_w;
		dPhiAngleNorm_Ewin->SetBinContent(ib,bc);
	}
	dPhiAngleNorm_Ewin->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc_Ewin = new TF1("fit_Ewin",fitfun,dPhiAngleNorm_Ewin->GetXaxis()->GetXmin(),dPhiAngleNorm_Ewin->GetXaxis()->GetXmax(),2);
		dPhiAngleNorm_Ewin->Fit(ffunc_Ewin,"RQ0");
		ffunc_Ewin->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngleNorm_Ewin->GetName()).Data()));
	
	dPhiAngle1->SetMinimum(0);
	dPhiAngle1->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc = new TF1("fit",fitfun,dPhiAngle1->GetXaxis()->GetXmin(),dPhiAngle1->GetXaxis()->GetXmax(),2);
		dPhiAngle1->Fit(ffunc,"RQ0");
		ffunc->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle1->GetName()).Data()));

	dPhiAngle1_Ewin->SetMinimum(0);
	dPhiAngle1_Ewin->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc_Ewin = new TF1("fit_Ewin",fitfun,dPhiAngle1_Ewin->GetXaxis()->GetXmin(),dPhiAngle1_Ewin->GetXaxis()->GetXmax(),2);
		dPhiAngle1_Ewin->Fit(ffunc_Ewin,"RQ0");
		ffunc_Ewin->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle1_Ewin->GetName()).Data()));	
	
	c2->cd(0);
	singleClusterSpecCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(singleClusterSpecCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(singleClusterSpecCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
	
	comptonSummedSpec2ClustersCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(comptonSummedSpec2ClustersCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(comptonSummedSpec2ClustersCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
	
	comptonSummedSpecEventCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(comptonSummedSpecEventCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(comptonSummedSpecEventCorr->GetName()).Data()));
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
	
	phiAngleSelEventsCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(phiAngleSelEventsCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(phiAngleSelEventsCorr->GetName()).Data()));
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
	
	nClustersInEventCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(nClustersInEventCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(nClustersInEventCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
	
	nClustersInSelEventCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(nClustersInSelEventCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(nClustersInSelEventCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
	
	delete c0;
	delete c1;
	delete c2;
	//logfile.close();
	
	if (totalTimeElapced < 60) cout << Form("Total running time is %.2f sec",totalTimeElapced) << endl;
	if (totalTimeElapced >= 60 && totalTimeElapced < 3600) cout << Form("Total running time is %dm:%.0f sec",Int_t(totalTimeElapced/60),totalTimeElapced - Int_t(totalTimeElapced/60)*60) << endl;
	if (totalTimeElapced >= 3600) cout << Form("Total running time is %dh:%dm:%.0fs",Int_t(totalTimeElapced/3600),Int_t((totalTimeElapced - Int_t(totalTimeElapced/3600)*3600)/60),
		totalTimeElapced - Int_t(totalTimeElapced/3600)*3600 - Int_t((totalTimeElapced - Int_t(totalTimeElapced/3600)*3600)/60)*60) << endl;
	
	hfile->cd();
	for (Int_t im = 0; im < nAMs; im++)
	{
		cathodeSpecAllEvents[im]->Write();
		cathodeSpecPair[im]->Write();
		cathodeSpecSelPairs[im]->Write();
		allClustersSpec[im]->Write();
		singleClusterEventsSpec[im]->Write();
		summedClusterEventsSpec[im]->Write();
		nTrigsInEvent[im]->Write();
		nTrigs_clusterSize[im]->Write();
		nTrigsInEvent_maxEcluster[im]->Write();
		comptonSummedSpec2Clusters[im]->Write();
		comptonSummedSpecEventClusters[im]->Write();
		comptonSpecClusters[im]->Write();
		comptonSpec2ClustersCorr[im]->Write();
		comptonHitsSpecClustersSizeCorr[im]->Write();
		comptonHitsClusters_nTrigsCorr[im]->Write();
		phiAngle[im]->Write();
		phiAngleSelEvents[im]->Write();
		dPhiPlaneAM[im]->Write();
		firstComptonCluster[im]->Write();
		secondComptonCluster[im]->Write();
		allClustersCOGImage[im]->Write();
		allClustersFineCOGImage[im]->Write();
		allPixelsFreqImage[im]->Write();
		allClustersCOGFreqImage[im]->Write();
		allClustersFineCOGFreqImage[im]->Write();
		allPixelsImage[im]->Write();
		clusterE_vsSize[im]->Write();
		clusterMaxE_vsSize[im]->Write();
		thetaFromFirstClusterE[im]->Write();
		thetaFromFirstClusterE_dPhi[im]->Write();
		nTrigsInEvent_norm[im]->Write();
		clusterSize[im]->Write();
		clusterE_nTrigs[im]->Write();
		clusterSecE_vsSize[im]->Write();
		clusterSizeMaxE_SecE[im]->Write();
		nClustersInEvent[im]->Write();
		nClustersInSelEvent[im]->Write();
		firstComptonClusterSize[im]->Write();
		secondComptonClusterSize[im]->Write();
		firstComptonClusterSizeSelEvents[im]->Write();
		secondComptonClusterSizeSelEvents[im]->Write();
		thetaFromFirstClusterE_w[im]->Write();
		thetaFromFirstClusterE_1pixOnly[im]->Write();
		thetaFromFirstClusterE_1pixOnly_w[im]->Write();
		thetaFromFirstClusterE_2pixOnly[im]->Write();
		thetaFromFirstClusterE_1_2pixOnly[im]->Write();
		if (makeRawTimingStuff)
		{
			rawTimingTwoPixelClustersCorr[im]->Write();
			rawTimingDiff2PixClusters[im]->Write();
			rawTimingDiff3PixClusters[im]->Write();
			if (!doNotUseCornerPixelsInPixelClusterCOGEneg) rawTimingDiff2PixDiagClusters[im]->Write();
			rawTimingDiffAllClusters[im]->Write();			
		}
		if (makeTimingCalibrationStuff)
		{
			anodeTimingTwoPixelClustersCorr[im]->Write();
			if (!doNotUseCornerPixelsInPixelClusterCOGEneg)
			{
				anodeTimingDiff2PixDiagClusters[im]->Write();
				anodeTiming2PixDiagClustersCorr[im]->Write();
			}
			anodeTimingDiff2PixClusters[im]->Write();
			anodeTimingDiff3PixClusters[im]->Write();
			anodeTiming2PixClustersCorr[im]->Write();
			anodeTiming3PixClustersCorr[im]->Write();
			anodeTimingDiffAllClusters[im]->Write();
			anodeTimingTwoClustersAMCorr[im]->Write();
			dAnodeTiming_vs_theta[im]->Write();
			dAnodeTiming_vs_E1[im]->Write();
			dAnodeTiming_vs_E2[im]->Write();
			dAnodeTiming_vs_E12[im]->Write();
			dAnodeTimingAbs_vs_E12[im]->Write();
			dist12_vs_E1[im]->Write();
			dist12_vs_E2[im]->Write();
			dist12_vs_E12[im]->Write();
			dist12_over_E12sq[im]->Write();
		}
		if (plotPixelClusterSpectra)
		{
			for (Int_t b = firstPixel; b < lastPixel+1; b++)
			{
				if (clusterSpec[im][b]->GetEntries() < 10) continue; 
				clusterSpec[im][b]->Write();
			}
		}
	}
	singleClusterSpecCorr->Write();
	summedClusterSpecCorr->Write();
	dPhiAngle->Write();
	dPhiAngleNorm->Write();
	dPhiAngle_Ewin->Write();
	dPhiAngleNorm_Ewin->Write();
	dPhiAngle1->Write();
	dPhiAngle1_Ewin->Write();
	dPhiAngle_2ClusterE->Write();
	dPhiAngle_totalE->Write();
	phiAngleCorr->Write();
	phiAngleSelEventsCorr->Write();
	comptonSummedSpec2ClustersCorr->Write();
	comptonSummedSpecEventCorr->Write();
	theta0_theta1_FromFirstClusterE->Write();
	nClustersInEventCorr->Write();
	nClustersInSelEventCorr->Write();
	hfile->Close();
}

Bool_t analyseNextEvent(const Int_t ievent)
{
	buffClusterX.clear();
	buffClusterY.clear();
	buffClusterE.clear();
	buffClusterArea.clear();
	buffClusterFlag.clear();
	buffClusterTrigs.clear();
	buffClusterIsSplit.clear();
	buffClusterAnodeTiming.clear();
	eventDisplayFlag = kFALSE;
	if (enableEventMonitor)
	{
		if (nEventsDisplayed < nEvents2Display) eventDisplayFlag = kTRUE;
		if (nEvents2Display == -1) eventDisplayFlag = kTRUE;
		if (ievent == eventNumber2Display) eventDisplayFlag = kTRUE;
	}

	for (Int_t im = 0; im < nAMs; im++)
	{
		cathodeE[im] = 0;
		if (eventDisplayFlag)
		{
			imageE[im]->Reset();
			imageN[im]->Reset();
			imageT[im]->Reset();
			imageC[im]->Reset();
		}
		if (!analyseAM(ievent, im)) return kFALSE;
	}
	Bool_t plotEvent = kFALSE;
	
	if (saveSelectedComptonRootTree && savingPointLocation == 1 && buffClusterX[0].size() >= 2 && buffClusterX[1].size() >= 2)
	{
		cntEv2BSaved++;
		selEvt[evtTreeIdx] = cntEv2BSaved;
		for (Int_t im = 0; im < nAMs; im++)
		{
			for (Int_t kk = firstAMTreeEventIdx[im]; kk <= lastAMTreeEventIdx[im]; kk++) selEvtAM[im][kk] = cntEv2BSaved;
		}
	}
	
	// For now, make pairs just from two first clusters in cluster array sorted by energy
	// First Compton interaction - second cluster (second in energy), second interaction - first cluster (maximum energy), assuming forward Compton scattering
	for (Int_t im = 0; im < nAMs; im++)
	{
		if (buffClusterX[im].size() >= 2)
		{
			Float_t xc2 = buffClusterX[im][0]-1;
			Float_t yc2 = buffClusterY[im][0]-1;
			Float_t xc1 = buffClusterX[im][1]-1;
			Float_t yc1 = buffClusterY[im][1]-1;
			Float_t dist = TMath::Hypot(xc1-xc2,yc1-yc2);
			Float_t phiAng = getPhiAngleDeg(xc2-xc1,yc2-yc1);
			phiAngle[im]->Fill(phiAng);
		}
	}
	Float_t totalE0 = 0;
	Bool_t goodClusterSizes = kTRUE;
	Bool_t skipEvent = kFALSE;
	for (Int_t p = 0; p < buffClusterX[0].size(); p++)
	{
		totalE0 += buffClusterE[0][p];
		if (buffClusterArea[0][p] < minClusterSize4PairAnalysis || buffClusterArea[0][p] > maxClusterSize4PairAnalysis) goodClusterSizes = kFALSE;
		if (usePixelListToDisable && disabledPixels[Int_t(pixelPattern[0]->GetBinContent(buffClusterX[0][p]+1,buffClusterY[0][p]+1)+0.5)]) skipEvent = kTRUE;
	}
	Float_t totalE1 = 0;
	for (Int_t p = 0; p < buffClusterX[1].size(); p++)
	{
		totalE1 += buffClusterE[1][p];
		if (buffClusterArea[1][p] < minClusterSize4PairAnalysis || buffClusterArea[1][p] > maxClusterSize4PairAnalysis) goodClusterSizes = kFALSE;
		if (usePixelListToDisable && disabledPixels[Int_t(pixelPattern[1]->GetBinContent(buffClusterX[1][p]+1,buffClusterY[1][p]+1)+0.5)]) skipEvent = kTRUE;
	}
	if (minClusterSize4PairAnalysis < 0 || maxClusterSize4PairAnalysis < 0) goodClusterSizes = kTRUE;
	
	Bool_t goodNClusters = kTRUE;
	if(buffClusterX[0].size() < minNClusters4PairAnalysis || buffClusterX[0].size() > maxNClusters4PairAnalysis) goodNClusters = kFALSE;
	if(buffClusterX[1].size() < minNClusters4PairAnalysis || buffClusterX[1].size() > maxNClusters4PairAnalysis) goodNClusters = kFALSE;
	if (minNClusters4PairAnalysis < 0 || maxNClusters4PairAnalysis < 0) goodNClusters = kTRUE;
	
	if (buffClusterX[0].size() == 1 || buffClusterX[1].size() == 1)
	{
		if (buffClusterX[0].size() == 1 && buffClusterX[1].size() == 1) singleClusterSpecCorr->Fill(totalE0,totalE1);
		summedClusterSpecCorr->Fill(totalE0,totalE1);
	}
	if (buffClusterX[0].size() >= 2 && buffClusterX[1].size() >= 2 && goodClusterSizes && goodNClusters && !skipEvent)
	{
		if (saveSelectedComptonRootTree && savingPointLocation == 2)
		{
			cntEv2BSaved++;
			selEvt[evtTreeIdx] = cntEv2BSaved;
			for (Int_t im = 0; im < nAMs; im++)
			{
				for (Int_t kk = firstAMTreeEventIdx[im]; kk <= lastAMTreeEventIdx[im]; kk++) selEvtAM[im][kk] = cntEv2BSaved;
			}
		}
		if (cathodeE[0] > 1 && cathodeE[1] > 1)
		{
			cathodeSpecPair[0]->Fill(cathodeE[0]);
			cathodeSpecPair[1]->Fill(cathodeE[1]);
		}
		firstComptonClusterSize[0]->Fill(buffClusterArea[0][firstClusterIdx[0]]);
		secondComptonClusterSize[0]->Fill(buffClusterArea[0][secondClusterIdx[0]]);
		firstComptonClusterSize[1]->Fill(buffClusterArea[1][firstClusterIdx[1]]);
		secondComptonClusterSize[1]->Fill(buffClusterArea[1][secondClusterIdx[1]]);
		
		Float_t initialE1 = getScatteredEnergyFromAngle(Na22Energy, angleOfSecondHead);
		Float_t sum2clE0 = buffClusterE[0][firstClusterIdx[0]]+buffClusterE[0][secondClusterIdx[0]];
		Float_t sum2clE1 = buffClusterE[1][firstClusterIdx[1]]+buffClusterE[1][secondClusterIdx[1]];
		
		sortClusters(0);
		Float_t xcentre = 6, ycentre = 17;
		Float_t phiAng0 = -9999;
		Float_t xc2 = buffClusterX[0][secondClusterIdx[0]]-1;
		Float_t yc2 = buffClusterY[0][secondClusterIdx[0]]-1;
		Float_t xc1 = buffClusterX[0][firstClusterIdx[0]]-1;
		Float_t yc1 = buffClusterY[0][firstClusterIdx[0]]-1;
		Float_t dist0 = TMath::Hypot(xc1-xc2,yc1-yc2);
		Bool_t phi0_valid = kFALSE;
		if (dist0 < maxDistanceBetweenClusters4ComptonPair && dist0 > 0.01)
		{
			phiAng0 = getPhiAngleDeg(xc2-xc1,yc2-yc1);
			firstComptonCluster[0]->Fill(xc1+1,yc1+1);
			secondComptonCluster[0]->Fill(xc2+1,yc2+1);
			dPhiPlaneAM[0]->Fill(xcentre,ycentre-11);
			dPhiPlaneAM[0]->Fill(xc2+xcentre-xc1,yc2+ycentre-yc1-11);
			phi0_valid = kTRUE;
		}
		
		sortClusters(1);
		Float_t phiAng1 = -9999;
		Float_t phiAng1mirrored = -9999;
		Float_t phiAng1mirrored_shifted = -9999;
		xc2 = buffClusterX[1][secondClusterIdx[1]]-1;
		yc2 = buffClusterY[1][secondClusterIdx[1]]-1;
		xc1 = buffClusterX[1][firstClusterIdx[1]]-1;
		yc1 = buffClusterY[1][firstClusterIdx[1]]-1;
		Float_t dist1 = TMath::Hypot(xc1-xc2,yc1-yc2);
		Bool_t phi1_valid = kFALSE;
		if (dist1 < maxDistanceBetweenClusters4ComptonPair && dist1 > 0.01)
		{
			phiAng1 = getPhiAngleDeg(xc2-xc1,yc2-yc1);
			phiAng1mirrored = getPhiAngleDeg(xc1-xc2,yc2-yc1);
			phiAng1mirrored_shifted = phiAng1mirrored + relativePhiAngle;
			if (phiAng1mirrored_shifted > 360) phiAng1mirrored_shifted -= 360;
			firstComptonCluster[1]->Fill(xc1+1,yc1+1);
			secondComptonCluster[1]->Fill(xc2+1,yc2+1);
			dPhiPlaneAM[1]->Fill(xcentre,ycentre-11);
			dPhiPlaneAM[1]->Fill(xc2+xcentre-xc1,yc2+ycentre-yc1-11);
			phi1_valid = kTRUE;
		}
		if (phi0_valid && phi1_valid)
		{
			if (cathodeE[0] > 1 && cathodeE[1] > 1)
			{
				cathodeSpecSelPairs[0]->Fill(cathodeE[0]);
				cathodeSpecSelPairs[1]->Fill(cathodeE[1]);
			}
			Float_t theta0 = getThetaFromEnergy(Na22Energy, buffClusterE[0][firstClusterIdx[0]]);
			Float_t theta1 = getThetaFromEnergy(initialE1, buffClusterE[1][firstClusterIdx[1]]);
			Float_t theta0a = getThetaFromEnergy(Na22Energy, buffClusterE[0][secondClusterIdx[0]]);
			Float_t theta1a = getThetaFromEnergy(initialE1, buffClusterE[1][secondClusterIdx[1]]);
			
			phiAngleCorr->Fill(phiAng0,phiAng1mirrored_shifted);
			Float_t dPh = phiAng0-phiAng1mirrored_shifted;
			if (dPh > 180) dPh = 360 - dPh;
			if (dPh < -180) dPh = -360 - dPh;
			if (makeTimingCalibrationStuff && buffClusterX[0].size() == 2 && buffClusterX[1].size() == 2)
			{
				dAnodeTiming_vs_theta[0]->Fill(theta0,buffClusterAnodeTiming[0][firstClusterIdx[0]]-buffClusterAnodeTiming[0][secondClusterIdx[0]]);
				dAnodeTiming_vs_theta[1]->Fill(theta1,buffClusterAnodeTiming[1][firstClusterIdx[1]]-buffClusterAnodeTiming[1][secondClusterIdx[1]]);
				dAnodeTiming_vs_theta[0]->Fill(theta0a,buffClusterAnodeTiming[0][secondClusterIdx[0]]-buffClusterAnodeTiming[0][firstClusterIdx[0]]);
				dAnodeTiming_vs_theta[1]->Fill(theta1a,buffClusterAnodeTiming[1][secondClusterIdx[1]]-buffClusterAnodeTiming[1][firstClusterIdx[1]]);
				dAnodeTiming_vs_E2[0]->Fill(buffClusterE[0][secondClusterIdx[0]],buffClusterAnodeTiming[0][firstClusterIdx[0]]-buffClusterAnodeTiming[0][secondClusterIdx[0]]);
				dAnodeTiming_vs_E2[1]->Fill(buffClusterE[1][secondClusterIdx[1]],buffClusterAnodeTiming[1][firstClusterIdx[1]]-buffClusterAnodeTiming[1][secondClusterIdx[1]]);
				dAnodeTiming_vs_E1[0]->Fill(buffClusterE[0][firstClusterIdx[0]],buffClusterAnodeTiming[0][secondClusterIdx[0]]-buffClusterAnodeTiming[0][firstClusterIdx[0]]);
				dAnodeTiming_vs_E1[1]->Fill(buffClusterE[1][firstClusterIdx[1]],buffClusterAnodeTiming[1][secondClusterIdx[1]]-buffClusterAnodeTiming[1][firstClusterIdx[1]]);
				dAnodeTiming_vs_E12[0]->Fill(buffClusterE[0][firstClusterIdx[0]],buffClusterAnodeTiming[0][secondClusterIdx[0]]-buffClusterAnodeTiming[0][firstClusterIdx[0]]);
				dAnodeTiming_vs_E12[1]->Fill(buffClusterE[1][firstClusterIdx[1]],buffClusterAnodeTiming[1][secondClusterIdx[1]]-buffClusterAnodeTiming[1][firstClusterIdx[1]]);
				dAnodeTiming_vs_E12[0]->Fill(buffClusterE[0][secondClusterIdx[0]],buffClusterAnodeTiming[0][firstClusterIdx[0]]-buffClusterAnodeTiming[0][secondClusterIdx[0]]);
				dAnodeTiming_vs_E12[1]->Fill(buffClusterE[1][secondClusterIdx[1]],buffClusterAnodeTiming[1][firstClusterIdx[1]]-buffClusterAnodeTiming[1][secondClusterIdx[1]]);
				dAnodeTimingAbs_vs_E12[0]->Fill(buffClusterE[0][firstClusterIdx[0]],TMath::Abs(buffClusterAnodeTiming[0][secondClusterIdx[0]]-buffClusterAnodeTiming[0][firstClusterIdx[0]]));
				dAnodeTimingAbs_vs_E12[1]->Fill(buffClusterE[1][firstClusterIdx[1]],TMath::Abs(buffClusterAnodeTiming[1][secondClusterIdx[1]]-buffClusterAnodeTiming[1][firstClusterIdx[1]]));
				dAnodeTimingAbs_vs_E12[0]->Fill(buffClusterE[0][secondClusterIdx[0]],TMath::Abs(buffClusterAnodeTiming[0][firstClusterIdx[0]]-buffClusterAnodeTiming[0][secondClusterIdx[0]]));
				dAnodeTimingAbs_vs_E12[1]->Fill(buffClusterE[1][secondClusterIdx[1]],TMath::Abs(buffClusterAnodeTiming[1][firstClusterIdx[1]]-buffClusterAnodeTiming[1][secondClusterIdx[1]]));
				Float_t ddist1 = TMath::Sq(buffClusterX[0][secondClusterIdx[0]] - buffClusterX[0][firstClusterIdx[0]]);
				ddist1 += TMath::Sq(buffClusterY[0][secondClusterIdx[0]] - buffClusterY[0][firstClusterIdx[0]]);
				ddist1 += TMath::Sq((buffClusterAnodeTiming[0][secondClusterIdx[0]] - buffClusterAnodeTiming[0][firstClusterIdx[0]])/factor2ConvertAnodeTime2Distance);
				ddist1 = TMath::Sqrt(ddist1);
				dist12_vs_E1[0]->Fill(TMath::Power(buffClusterE[0][firstClusterIdx[0]],2.2),ddist1);
				dist12_vs_E2[0]->Fill(TMath::Power(buffClusterE[0][secondClusterIdx[0]],2.2),ddist1);
				dist12_vs_E12[0]->Fill(TMath::Power(buffClusterE[0][firstClusterIdx[0]],2.2),ddist1);
				dist12_vs_E12[0]->Fill(TMath::Power(buffClusterE[0][secondClusterIdx[0]],2.2),ddist1);
				Float_t ddist2 = TMath::Sq(buffClusterX[1][secondClusterIdx[1]] - buffClusterX[1][firstClusterIdx[1]]);
				ddist2 += TMath::Sq(buffClusterY[1][secondClusterIdx[1]] - buffClusterY[1][firstClusterIdx[1]]);
				ddist2 += TMath::Sq((buffClusterAnodeTiming[1][secondClusterIdx[1]] - buffClusterAnodeTiming[1][firstClusterIdx[1]])/factor2ConvertAnodeTime2Distance);
				ddist2 = TMath::Sqrt(ddist2);
				dist12_vs_E1[1]->Fill(TMath::Power(buffClusterE[1][firstClusterIdx[1]],2.2),ddist2);
				dist12_vs_E2[1]->Fill(TMath::Power(buffClusterE[1][secondClusterIdx[1]],2.2),ddist2);
				dist12_vs_E12[1]->Fill(TMath::Power(buffClusterE[1][firstClusterIdx[1]],2.2),ddist2);
				dist12_vs_E12[1]->Fill(TMath::Power(buffClusterE[1][secondClusterIdx[1]],2.2),ddist2);
				dist12_over_E12sq[0]->Fill(ddist1*1000/TMath::Power(buffClusterE[0][0],2.2),ddist1*1000/TMath::Power(buffClusterE[0][1],2.2));
				dist12_over_E12sq[1]->Fill(ddist2*1000/TMath::Power(buffClusterE[1][0],2.2),ddist2*1000/TMath::Power(buffClusterE[1][1],2.2));
			}
			if (theta0 > ThetaWindowFor4PhiAnalyis_min && theta0 < ThetaWindowFor4PhiAnalyis_max && theta1 > ThetaWindowFor4PhiAnalyis_min && theta1 < ThetaWindowFor4PhiAnalyis_max && cathodeE[0] < 540 && cathodeE[1] < 540) 
			{
				phiAngleSelEventsCorr->Fill(phiAng0,phiAng1mirrored_shifted);
				phiAngleSelEvents[0]->Fill(phiAng0);
				phiAngleSelEvents[1]->Fill(phiAng1mirrored_shifted);
				nClustersInSelEvent[0]->Fill(buffClusterX[0].size());
				nClustersInSelEvent[1]->Fill(buffClusterX[1].size());
				nClustersInSelEventCorr->Fill(buffClusterX[0].size(),buffClusterX[1].size());
				comptonHitsSpecClustersSizeCorr[0]->Fill(buffClusterArea[0][firstClusterIdx[0]],buffClusterArea[0][secondClusterIdx[0]]);
				comptonHitsSpecClustersSizeCorr[1]->Fill(buffClusterArea[1][firstClusterIdx[1]],buffClusterArea[1][secondClusterIdx[1]]);
				firstComptonClusterSizeSelEvents[0]->Fill(buffClusterArea[0][firstClusterIdx[0]]);
				secondComptonClusterSizeSelEvents[0]->Fill(buffClusterArea[0][secondClusterIdx[0]]);
				firstComptonClusterSizeSelEvents[1]->Fill(buffClusterArea[1][firstClusterIdx[1]]);
				secondComptonClusterSizeSelEvents[1]->Fill(buffClusterArea[1][secondClusterIdx[1]]);
				dPhiAngle->Fill(dPh);
				dPhiAngleNorm->Fill(dPh);
				dPhiAngle1->Fill(TMath::Abs(dPh));
				plotEvent = kTRUE;
				if (totalE0 >= EwindowFor4PhiAnalyis_min && totalE0 <= EwindowFor4PhiAnalyis_max && totalE1 >= EwindowFor4PhiAnalyis_min && totalE1 <= EwindowFor4PhiAnalyis_max)
				{
					dPhiAngle_Ewin->Fill(dPh);
					dPhiAngle1_Ewin->Fill(TMath::Abs(dPh));
					dPhiAngleNorm_Ewin->Fill(dPh);
				}
				if (makeTimingCalibrationStuff && buffClusterX[0].size() == 2 && buffClusterX[1].size() == 2)
				{
					anodeTimingTwoClustersAMCorr[0]->Fill(buffClusterAnodeTiming[0][firstClusterIdx[0]],buffClusterAnodeTiming[0][secondClusterIdx[0]]);
					anodeTimingTwoClustersAMCorr[1]->Fill(buffClusterAnodeTiming[1][firstClusterIdx[1]],buffClusterAnodeTiming[1][secondClusterIdx[1]]);
				}
			}
			dPhiAngle_2ClusterE->Fill(sum2clE0,dPh);
			dPhiAngle_2ClusterE->Fill(sum2clE1,dPh);
			thetaFromFirstClusterE[0]->Fill(theta0);
			thetaFromFirstClusterE[1]->Fill(theta1);
			thetaFromFirstClusterE_w[0]->Fill(theta0);
			thetaFromFirstClusterE_w[1]->Fill(theta1);
			thetaFromFirstClusterE_w[0]->Fill(theta0a);
			thetaFromFirstClusterE_w[1]->Fill(theta1a);
			if (buffClusterArea[0][firstClusterIdx[0]] == 1 && buffClusterArea[0][secondClusterIdx[0]] == 1)
			{
				thetaFromFirstClusterE_1pixOnly[0]->Fill(theta0);
				thetaFromFirstClusterE_1pixOnly_w[0]->Fill(theta0);
				thetaFromFirstClusterE_1pixOnly_w[0]->Fill(theta0a);
			}
			if (buffClusterArea[1][firstClusterIdx[1]] == 1 && buffClusterArea[1][secondClusterIdx[1]] == 1)
			{
				thetaFromFirstClusterE_1pixOnly[1]->Fill(theta1);
				thetaFromFirstClusterE_1pixOnly_w[1]->Fill(theta1);
				thetaFromFirstClusterE_1pixOnly_w[1]->Fill(theta1);
			}
			if (buffClusterArea[0][firstClusterIdx[0]] == 2 && buffClusterArea[0][secondClusterIdx[0]] == 2) thetaFromFirstClusterE_2pixOnly[0]->Fill(theta0);
			if (buffClusterArea[1][firstClusterIdx[1]] == 2 && buffClusterArea[1][secondClusterIdx[1]] == 2) thetaFromFirstClusterE_2pixOnly[1]->Fill(theta1);
			if (buffClusterArea[0][firstClusterIdx[0]] == 1 && buffClusterArea[0][secondClusterIdx[0]] == 2) thetaFromFirstClusterE_1_2pixOnly[0]->Fill(theta0);
			if (buffClusterArea[1][firstClusterIdx[1]] == 1 && buffClusterArea[1][secondClusterIdx[1]] == 2) thetaFromFirstClusterE_1_2pixOnly[1]->Fill(theta1);
			if (buffClusterArea[0][firstClusterIdx[0]] == 2 && buffClusterArea[0][secondClusterIdx[0]] == 1) thetaFromFirstClusterE_1_2pixOnly[0]->Fill(theta0);
			if (buffClusterArea[1][firstClusterIdx[1]] == 2 && buffClusterArea[1][secondClusterIdx[1]] == 1) thetaFromFirstClusterE_1_2pixOnly[1]->Fill(theta1);
			thetaFromFirstClusterE_dPhi[0]->Fill(theta0,dPh);
			thetaFromFirstClusterE_dPhi[1]->Fill(theta1,dPh);
			theta0_theta1_FromFirstClusterE->Fill(theta0,theta1);
			dPhiAngle_totalE->Fill(totalE0,dPh);
			dPhiAngle_totalE->Fill(totalE1,dPh);
		}
		comptonSummedSpec2ClustersCorr->Fill(sum2clE0,sum2clE1);
		comptonSummedSpecEventCorr->Fill(totalE0,totalE1);
		summedClusterSpecCorr->Fill(sum2clE0,sum2clE1);
		comptonSummedSpecEventClusters[0]->Fill(totalE0);
		comptonSummedSpecEventClusters[1]->Fill(totalE1);
	}

	for (Int_t im = 0; im < nAMs; im++)
	{
		if (buffClusterX[im].size() >= 2)
		{
			for (Int_t p = 0; p < 2; p++) comptonSpecClusters[im]->Fill(buffClusterE[im][p]);
			comptonSummedSpec2Clusters[im]->Fill(buffClusterE[im][0] + buffClusterE[im][1]);
			summedClusterEventsSpec[im]->Fill(buffClusterE[im][0] + buffClusterE[im][1]);
			comptonSpec2ClustersCorr[im]->Fill(buffClusterE[im][0],buffClusterE[im][1]);
			comptonHitsSpecClustersSizeCorr[im]->Fill(buffClusterArea[im][0],buffClusterArea[im][1]);
			comptonHitsClusters_nTrigsCorr[im]->Fill(buffClusterTrigs[im][0],buffClusterTrigs[im][1]);
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
			if (usePixelListToDisable && disabledPixels[Int_t(pixelPattern[im]->GetBinContent(buffClusterX[im][p]+1,buffClusterY[im][p]+1)+0.5)]) continue;
			allClustersCOGImage[im]->Fill(buffClusterX[im][p],buffClusterY[im][p],buffClusterE[im][p]);
			allClustersCOGFreqImage[im]->Fill(buffClusterX[im][p],buffClusterY[im][p]);
			allClustersFineCOGImage[im]->Fill(buffClusterX[im][p],buffClusterY[im][p],buffClusterE[im][p]);
			allClustersFineCOGFreqImage[im]->Fill(buffClusterX[im][p],buffClusterY[im][p]);
			if (plotPixelClusterSpectra) clusterSpec[im][Int_t(pixelPattern[im]->GetBinContent(buffClusterX[im][p]+1,buffClusterY[im][p]+1)+0.5)]->Fill(buffClusterE[im][p]);
		}
		if (buffClusterX[im].size() > 0) nClustersInEvent[im]->Fill(buffClusterX[im].size());
	}
	if (buffClusterX[0].size() > 0 && buffClusterX[1].size() > 0) nClustersInEventCorr->Fill(buffClusterX[0].size(),buffClusterX[1].size());

	if (plotEvent)
	{
		for (Int_t im = 0; im < nAMs; im++)
		{
			if (eventDisplayFlag)
			{
				c3->cd(im+1);
				if (buffClusterE[im].size() == 1) imageE[im]->SetTitle(Form("Event %d, 1 cluster",ievent));
				if (buffClusterE[im].size() > 1) imageE[im]->SetTitle(Form("Event %d, %ld clusters",ievent,buffClusterE[im].size()));
				imageE[im]->Draw("colz");
				imageT[im]->Draw("same");
				imageC[im]->Draw("same");
				imageN[im]->Draw("sametext");
				line2->DrawLine(imageE[im]->GetXaxis()->GetXmin(),imageE[im]->GetYaxis()->GetXmax()/2,imageE[im]->GetXaxis()->GetXmax(),imageE[im]->GetYaxis()->GetXmax()/2);
				line2->DrawLine(imageE[im]->GetXaxis()->GetXmax()/2,imageE[im]->GetYaxis()->GetXmin(),imageE[im]->GetXaxis()->GetXmax()/2,imageE[im]->GetYaxis()->GetXmax());
			}		
		}
		if (eventDisplayFlag)
		{
			c3->SaveAs(Form("%s/event%d.gif",evMonDir.Data(),ievent));
			nEventsDisplayed++;
		}
	}
	return kTRUE;
}

Bool_t analyseAM(const Int_t evt, const Int_t AM2do)
{
	Bool_t status = kFALSE;
	//Bool_t status = kTRUE;
	Int_t centrePixelX = -1;
	Int_t centrePixelY = -1;
	buffClusterXloc.clear();
	buffClusterYloc.clear();
	buffClusterEloc.clear();
	buffClusterArealoc.clear();
	buffClusterFlagloc.clear();
	buffClusterTrigsloc.clear();
	buffClusterIsSplitloc.clear();
	buffClusterAnodeTimingloc.clear();
	buffClusterXloc.resize(0,0);
	buffClusterYloc.resize(0,0);
	buffClusterEloc.resize(0,0);
	buffClusterArealoc.resize(0,0);
	buffClusterFlagloc.resize(0,0);
	buffClusterTrigsloc.resize(0,0);
	buffClusterIsSplitloc.resize(0,0);
	buffClusterAnodeTimingloc.resize(0,0);
	matrix_trigs.clear();
	matrix_trigs.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_flags.clear();
	matrix_flags.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_E.clear();
	matrix_E.resize(nPixXY, vector<Float_t>(nPixXY));
	matrix_Eneg.clear();
	matrix_Eneg.resize(nPixXY, vector<Float_t>(nPixXY));
	matrix_time.clear();
	matrix_time.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_timeCalib.clear();
	matrix_timeCalib.resize(nPixXY, vector<Float_t>(nPixXY));
	firstAMTreeEventIdx[AM2do] = 0;
	lastAMTreeEventIdx[AM2do] = 0;

	//cout << "=============Start=============" << endl;
	for (Int_t ii = 0; ii < Nev[AM2do] && event[AM2do] != evt; ii++)
	{
		cntEntryTree[AM2do]++;
		//cout << Form("cntEntryTree[%d] = %d",AM2do,cntEntryTree[AM2do]);
		chain1p[AM2do]->GetEntry(cntEntryTree[AM2do]);
		//cout << ", got it" << endl;
		if (event[AM2do] == evt)
		{
			firstAMTreeEventIdx[AM2do] = cntEntryTree[AM2do];
			break;
		}
	}
	while (event[AM2do] == evt && cntEntryTree[AM2do] <= Nev[AM2do])
	{
		if (pixel[AM2do] >= firstPixel && pixel[AM2do] <= lastPixel)
		{
			Int_t xx1 = -1, yy1 = -1;
			if (getPixel2Dcoordinates(pixel[AM2do],AM2do,xx1,yy1))
			{
				allPixelsImage[AM2do]->Fill(xx1-1,yy1-1,E[AM2do]);
				allPixelsFreqImage[AM2do]->Fill(xx1-1,yy1-1);
				if (eventDisplayFlag) imageN[AM2do]->SetBinContent(xx1,yy1,Int_t(E_neg[AM2do]+0.5));
				matrix_Eneg[xx1-1][yy1-1] = E_neg[AM2do];
				if (triggerFlag[AM2do] == 1)
				{
					matrix_time[xx1-1][yy1-1] = timeDetect[AM2do];
					matrix_timeCalib[xx1-1][yy1-1] = (timeDetect[AM2do]-tCalibAnodes[pixel[AM2do]][1]+rand3->Uniform(1.))/tCalibAnodes[pixel[AM2do]][0];				
				}
				//cout << Form("%d, %d - timing %d, trigFlag = %d",xx1,yy1,timeDetect[AM2do],triggerFlag[AM2do])<< endl;
				if (E[AM2do] > minPixelEnergyThr4ClusterReconstruction)
				{
					if ((useOnlyTriggeredPixelsInCluster && triggerFlag[AM2do] == 1) || !useOnlyTriggeredPixelsInCluster)
					{
						status = kTRUE;
						matrix_E[xx1-1][yy1-1] = E[AM2do];
						matrix_trigs[xx1-1][yy1-1] = triggerFlag[AM2do];
						matrix_flags[xx1-1][yy1-1] = newNeighbourFlag;
						if (eventDisplayFlag)
						{
							imageE[AM2do]->SetBinContent(xx1,yy1,E[AM2do]);
							imageN[AM2do]->SetBinContent(xx1,yy1,0);
							if (triggerFlag[AM2do] == 1) imageT[AM2do]->Fill(xx1-0.5,yy1-0.5);
						}
					}
				}
			}
		}
		if (pixel[AM2do] == getCath1Channel(AM2do,GM[AM2do],0))
		{
			cathodeSpecAllEvents[AM2do]->Fill(E[AM2do]);
			cathodeE[AM2do] = E[AM2do];
		}
		cntEntryTree[AM2do]++;
		chain1p[AM2do]->GetEntry(cntEntryTree[AM2do]);
	}
	cntEntryTree[AM2do]--;
	chain1p[AM2do]->GetEntry(cntEntryTree[AM2do]);
	lastAMTreeEventIdx[AM2do] = cntEntryTree[AM2do];
	if (!status) return status;
	
	Int_t clusterCounter = 0;
	Bool_t skipFrame = kFALSE;
	Float_t primaryClusterE = 0;
	Float_t E1 = -1, E2 = -1;
	Float_t xc1 = -1, yc1 = -1;
	Float_t xc2 = -1, yc2 = -1;
	Float_t aTime;
	for (Int_t ii = 0; ii < nPixXY&&(!skipFrame); ii++)
	{
		for (Int_t jj = 0; jj < nPixXY&&(!skipFrame); jj++)
		{
			if (matrix_flags[ii][jj] != newNeighbourFlag) continue;
			buffClusterXarr.clear();
			buffClusterYarr.clear();
			buffClusterTimearr.clear();
			buffClusterTrigsarr.clear();
			buffClusterPixelarr.clear();
			buffClusterAnodeTimearr.clear();
			buffClusterXarr.resize(0,0);
			buffClusterYarr.resize(0,0);
			buffClusterTimearr.resize(0,0);
			buffClusterTrigsarr.resize(0,0);
			buffClusterPixelarr.resize(0,0);
			buffClusterAnodeTimearr.resize(0,0);
			Float_t CoGx = 0;
			Float_t CoGy = 0;
			Int_t area = 0;
			Float_t totAmp = 0;
			Float_t totAmpEcal = 0;
			Float_t CoGxEcal = 0;
			Float_t CoGyEcal = 0;
			Bool_t isPrimaryCluster = kFALSE;
			Int_t nTrigsInCluster = 0;
			Float_t meanT = 0;
			
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
			buffClusterXarr.push_back(ii);
			buffClusterYarr.push_back(jj);
			buffClusterTimearr.push_back(matrix_time[ii][jj]);
			aTime = (matrix_time[ii][jj]-tCalibAnodes[Int_t(pixelPattern[AM2do]->GetBinContent(ii+1,jj+1)+0.5)][1]+rand3->Uniform(1.))/tCalibAnodes[Int_t(pixelPattern[AM2do]->GetBinContent(ii+1,jj+1)+0.5)][0];
			buffClusterAnodeTimearr.push_back(aTime);
			buffClusterTrigsarr.push_back(matrix_trigs[ii][jj]);
			buffClusterPixelarr.push_back(Int_t(pixelPattern[AM2do]->GetBinContent(ii+1,jj+1)+0.5));
			if (ii == centrePixelX-1 && jj == centrePixelY-1) isPrimaryCluster = kTRUE;
			if (matrix_trigs[ii][jj] == 1)
			{
				meanT += aTime;
				nTrigsInCluster++;
			}
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
						buffClusterXarr.push_back(buffXnew[ic]);
						buffClusterYarr.push_back(buffYnew[ic]);
						buffClusterTimearr.push_back(matrix_time[buffXnew[ic]][buffYnew[ic]]);
						aTime = (matrix_time[buffXnew[ic]][buffYnew[ic]]-tCalibAnodes[Int_t(pixelPattern[AM2do]->GetBinContent(buffXnew[ic]+1,buffYnew[ic]+1)+0.5)][1]+rand3->Uniform(1.))/tCalibAnodes[Int_t(pixelPattern[AM2do]->GetBinContent(buffXnew[ic]+1,buffYnew[ic]+1)+0.5)][0];
						buffClusterAnodeTimearr.push_back(aTime);
						buffClusterTrigsarr.push_back(matrix_trigs[buffXnew[ic]][buffYnew[ic]]);
						buffClusterPixelarr.push_back(Int_t(pixelPattern[AM2do]->GetBinContent(buffXnew[ic]+1,buffYnew[ic]+1)+0.5));
						if (buffXnew[ic] == centrePixelX-1 && buffYnew[ic] == centrePixelY-1) isPrimaryCluster = kTRUE;
						if (matrix_trigs[buffXnew[ic]][buffYnew[ic]] == 1)
						{
							meanT += aTime;
							nTrigsInCluster++;
						}
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
			Float_t xcen = CoGx/totAmp;
			Float_t ycen = CoGy/totAmp;
			if (nTrigsInCluster > 0) meanT /= nTrigsInCluster;
			Float_t namp;
			if (updateSinglePixelClusterCOGEneg && area == 1 && pixelStatusGlobal[AM2do][Int_t(xcen)][Int_t(ycen)] != -1)
			{
				namp = updateSinglePixelClusterCOGwithEneg(Int_t(xcen),Int_t(ycen),CoGx,CoGy);
				xcen = CoGx/namp;
				ycen = CoGy/namp;
			}
			Int_t stat = 0;
			if (updateTwoPixelClusterCOGEneg && area == 2)
			{
				updateTwoPixelClusterCOGwithEneg(AM2do,CoGx/totAmp,CoGy/totAmp,buffClusterXarr[0],buffClusterYarr[0],buffClusterXarr[1],buffClusterYarr[1],xcen,ycen);
				stat = 1;
				if (buffClusterXarr[0] == buffClusterXarr[1] || buffClusterYarr[0] == buffClusterYarr[1]) stat = 0;
			}
			buffClusterIsSplitloc.push_back(stat);
			if (eventDisplayFlag) imageC[AM2do]->Fill(xcen,ycen);
			if (makeRawTimingStuff || makeTimingCalibrationStuff)
			{
				if (nTrigsInCluster >= 2)
				{
					Int_t t1[2] = {0};
					Float_t t1a[2] = {0};
					Float_t t1b[2] = {0};
					for (Int_t it = 0; it < area; it++)
					{
						for (Int_t it2 = it+1; it2 < area; it2++)
						{
							t1b[0] = buffClusterAnodeTimearr[it];
							t1b[1] = buffClusterAnodeTimearr[it2];
							if (makeRawTimingStuff) rawTimingDiffAllClusters[AM2do]->Fill(TMath::Abs(buffClusterTimearr[it]-buffClusterTimearr[it2]));
							if (makeTimingCalibrationStuff) anodeTimingDiffAllClusters[AM2do]->Fill(TMath::Abs(t1b[0]-t1b[1]));
							if (nTrigsInCluster > 2)
							{
								if (makeRawTimingStuff) rawTimingDiff3PixClusters[AM2do]->Fill(TMath::Abs(buffClusterTimearr[it]-buffClusterTimearr[it2]));
								if (makeTimingCalibrationStuff)
								{
									anodeTimingDiff3PixClusters[AM2do]->Fill(TMath::Abs(t1b[0]-t1b[1]));
									anodeTiming3PixClustersCorr[AM2do]->Fill(t1b[0],t1b[1]);
								}
							}
							if (nTrigsInCluster == 2 && buffClusterTrigsarr[it] == 1 && buffClusterTrigsarr[it2] == 1) 
							{
								t1[0] = buffClusterTimearr[it];
								t1[1] = buffClusterTimearr[it2];
								t1a[0] = t1b[0];
								t1a[1] = t1b[1];
								if (makeRawTimingStuff) rawTimingDiff2PixClusters[AM2do]->Fill(TMath::Abs(buffClusterTimearr[it]-buffClusterTimearr[it2]));
								if (makeTimingCalibrationStuff)
								{
									anodeTimingDiff2PixClusters[AM2do]->Fill(TMath::Abs(t1b[0]-t1b[1]));
									anodeTiming2PixClustersCorr[AM2do]->Fill(t1b[0],t1b[1]);
								}
								if (buffClusterXarr[it] != buffClusterXarr[it2] && buffClusterYarr[it] != buffClusterYarr[it2])
								{
									if (!doNotUseCornerPixelsInPixelClusterCOGEneg)
									{
										if (makeRawTimingStuff) rawTimingDiff2PixDiagClusters[AM2do]->Fill(TMath::Abs(buffClusterTimearr[it]-buffClusterTimearr[it2]));
										if (makeTimingCalibrationStuff)	
										{
											anodeTimingDiff2PixDiagClusters[AM2do]->Fill(TMath::Abs(t1a[0]-t1a[1]));
											anodeTiming2PixDiagClustersCorr[AM2do]->Fill(t1a[0],t1a[1]);
										}
									}
								}
							}
						}
					}
					if (t1[0] > 0 && t1[1] > 0)
					{
						if (makeRawTimingStuff) rawTimingTwoPixelClustersCorr[AM2do]->Fill(t1[0],t1[1]);
						if (makeTimingCalibrationStuff) anodeTimingTwoPixelClustersCorr[AM2do]->Fill(t1a[0],t1a[1]);
					}
				}
			}
			buffClusterXloc.push_back(xcen);
			buffClusterYloc.push_back(ycen);
			buffClusterEloc.push_back(totAmp);
			buffClusterArealoc.push_back(area);
			if (isPrimaryCluster) buffClusterFlagloc.push_back(1);
			else buffClusterFlagloc.push_back(0);
			buffClusterTrigsloc.push_back(nTrigsInCluster);
			buffClusterAnodeTimingloc.push_back(meanT);
		}
	}
	if (clusterCounter > 0)
	{
		Int_t niter = 0;
		while (niter < buffClusterXloc.size())
		{
			Float_t maxE = -9999;
			Int_t idx = -1, maxArea, maxFlag, maxTrigs, maxIsSplit;
			Float_t maxX, maxY, maxT;
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
					maxIsSplit = buffClusterIsSplitloc[ic];
					maxT = buffClusterAnodeTimingloc[ic];
				}
			}
			buffClusterEloc.erase(buffClusterEloc.begin()+idx);
			buffClusterXloc.erase(buffClusterXloc.begin()+idx);
			buffClusterYloc.erase(buffClusterYloc.begin()+idx);
			buffClusterArealoc.erase(buffClusterArealoc.begin()+idx);
			buffClusterFlagloc.erase(buffClusterFlagloc.begin()+idx);
			buffClusterTrigsloc.erase(buffClusterTrigsloc.begin()+idx);
			buffClusterIsSplitloc.erase(buffClusterIsSplitloc.begin()+idx);
			buffClusterAnodeTimingloc.erase(buffClusterAnodeTimingloc.begin()+idx);
			buffClusterEloc.insert(buffClusterEloc.begin()+niter,maxE);
			buffClusterXloc.insert(buffClusterXloc.begin()+niter,maxX);
			buffClusterYloc.insert(buffClusterYloc.begin()+niter,maxY);
			buffClusterArealoc.insert(buffClusterArealoc.begin()+niter,maxArea);
			buffClusterFlagloc.insert(buffClusterFlagloc.begin()+niter,maxFlag);
			buffClusterTrigsloc.insert(buffClusterTrigsloc.begin()+niter,maxTrigs);
			buffClusterIsSplitloc.insert(buffClusterIsSplitloc.begin()+niter,maxIsSplit);
			buffClusterAnodeTimingloc.insert(buffClusterAnodeTimingloc.begin()+niter,maxT);
			niter++;
		}
		nTrigsInEvent_maxEcluster[AM2do]->Fill(buffClusterEloc[0],nTrigPixels[AM2do]);
		nTrigsInEvent[AM2do]->Fill(nTrigPixels[AM2do]);
		buffClusterX.push_back(buffClusterXloc);
		buffClusterY.push_back(buffClusterYloc);
		buffClusterE.push_back(buffClusterEloc);
		buffClusterArea.push_back(buffClusterArealoc);
		buffClusterFlag.push_back(buffClusterFlagloc);
		buffClusterTrigs.push_back(buffClusterTrigsloc);
		buffClusterIsSplit.push_back(buffClusterIsSplitloc);
		buffClusterAnodeTiming.push_back(buffClusterAnodeTimingloc);
	}
	//cout << "=============End=============" << endl;
	return status;
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
		
		if (sline.Contains("usePixelListToDisable ="))
		{
			sline.ReplaceAll("usePixelListToDisable =","");
			buf = sline.Atoi();
			usePixelListToDisable = kTRUE;
			if (buf == 0) usePixelListToDisable = kFALSE;
			if (printOutSetupFileParameters) cout << "usePixelListToDisable = " << usePixelListToDisable << endl;
		}		
	
		if (sline.Contains("pixelListFileName ="))
		{
			sline.ReplaceAll("pixelListFileName =","");
			pixelListFileName = sline;
			if (printOutSetupFileParameters) cout << "pixelListFileName = " << pixelListFileName << endl;
		}
	
		if (sline.Contains("enableEventMonitor ="))
		{
			sline.ReplaceAll("enableEventMonitor =","");
			buf = sline.Atoi();
			enableEventMonitor = kTRUE;
			if (buf == 0) enableEventMonitor = kFALSE;
			if (printOutSetupFileParameters) cout << "enableEventMonitor = " << enableEventMonitor << endl;
		}
		
		if (sline.Contains("nEvents2Display ="))
		{
			sline.ReplaceAll("nEvents2Display =","");
			nEvents2Display = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nEvents2Display = " << nEvents2Display << endl;
		}

		if (sline.Contains("eventNumber2Display ="))
		{
			sline.ReplaceAll("eventNumber2Display =","");
			eventNumber2Display = sline.Atoi();
			if (printOutSetupFileParameters) cout << "eventNumber2Display = " << eventNumber2Display << endl;
		}
		
		if (sline.Contains("makeTimingCalibrationStuff ="))
		{
			sline.ReplaceAll("makeTimingCalibrationStuff =","");
			buf = sline.Atoi();
			makeTimingCalibrationStuff = kTRUE;
			if (buf == 0) makeTimingCalibrationStuff = kFALSE;
			if (printOutSetupFileParameters) cout << "makeTimingCalibrationStuff = " << makeTimingCalibrationStuff << endl;
		}
		
		if (sline.Contains("makeRawTimingStuff ="))
		{
			sline.ReplaceAll("makeRawTimingStuff =","");
			buf = sline.Atoi();
			makeRawTimingStuff = kTRUE;
			if (buf == 0) makeRawTimingStuff = kFALSE;
			if (printOutSetupFileParameters) cout << "makeRawTimingStuff = " << makeRawTimingStuff << endl;
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
		
		if (sline.Contains("nEvents2Do ="))
		{
			sline.ReplaceAll("nEvents2Do =","");
			nEvents2Do = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nEvents2Do = " << nEvents2Do << endl;
		}
		
		if (sline.Contains("useOnlyTriggeredPixelsInCluster ="))
		{
			sline.ReplaceAll("useOnlyTriggeredPixelsInCluster =","");
			buf = sline.Atoi();
			useOnlyTriggeredPixelsInCluster = kTRUE;
			if (buf == 0) useOnlyTriggeredPixelsInCluster = kFALSE;
			if (printOutSetupFileParameters) cout << "useOnlyTriggeredPixelsInCluster = " << useOnlyTriggeredPixelsInCluster << endl;
		}
		
		if (sline.Contains("doNotUseCornerTouchingPixels4Clustering ="))
		{
			sline.ReplaceAll("doNotUseCornerTouchingPixels4Clustering =","");
			buf = sline.Atoi();
			doNotUseCornerTouchingPixels4Clustering = kTRUE;
			if (buf == 0) doNotUseCornerTouchingPixels4Clustering = kFALSE;
			if (printOutSetupFileParameters) cout << "doNotUseCornerTouchingPixels4Clustering = " << doNotUseCornerTouchingPixels4Clustering << endl;
		}
		
		if (sline.Contains("makeCosineFit ="))
		{
			sline.ReplaceAll("makeCosineFit =","");
			buf = sline.Atoi();
			makeCosineFit = kTRUE;
			if (buf == 0) makeCosineFit = kFALSE;
			if (printOutSetupFileParameters) cout << "makeCosineFit = " << makeCosineFit << endl;
		}
		
		if (sline.Contains("minClusterSize4PairAnalysis ="))
		{
			sline.ReplaceAll("minClusterSize4PairAnalysis =","");
			minClusterSize4PairAnalysis = sline.Atoi();
			if (printOutSetupFileParameters) cout << "minClusterSize4PairAnalysis = " << minClusterSize4PairAnalysis << endl;
		}
		
		if (sline.Contains("maxClusterSize4PairAnalysis ="))
		{
			sline.ReplaceAll("maxClusterSize4PairAnalysis =","");
			maxClusterSize4PairAnalysis = sline.Atoi();
			if (printOutSetupFileParameters) cout << "maxClusterSize4PairAnalysis = " << maxClusterSize4PairAnalysis << endl;
		}
		
		if (sline.Contains("updateSinglePixelClusterCOGEneg ="))
		{
			sline.ReplaceAll("updateSinglePixelClusterCOGEneg =","");
			buf = sline.Atoi();
			updateSinglePixelClusterCOGEneg = kTRUE;
			if (buf == 0) updateSinglePixelClusterCOGEneg = kFALSE;
			if (printOutSetupFileParameters) cout << "updateSinglePixelClusterCOGEneg = " << updateSinglePixelClusterCOGEneg << endl;
		}
		
		if (sline.Contains("updateTwoPixelClusterCOGEneg ="))
		{
			sline.ReplaceAll("updateTwoPixelClusterCOGEneg =","");
			buf = sline.Atoi();
			updateTwoPixelClusterCOGEneg = kTRUE;
			if (buf == 0) updateTwoPixelClusterCOGEneg = kFALSE;
			if (printOutSetupFileParameters) cout << "updateTwoPixelClusterCOGEneg = " << updateTwoPixelClusterCOGEneg << endl;
		}
		
		if (sline.Contains("updateTwoPixelClusterCOG_1Donly ="))
		{
			sline.ReplaceAll("updateTwoPixelClusterCOG_1Donly =","");
			buf = sline.Atoi();
			updateTwoPixelClusterCOG_1Donly = kTRUE;
			if (buf == 0) updateTwoPixelClusterCOG_1Donly = kFALSE;
			if (printOutSetupFileParameters) cout << "updateTwoPixelClusterCOG_1Donly = " << updateTwoPixelClusterCOG_1Donly << endl;
		}
		
		if (sline.Contains("minRawTiming4H ="))
		{
			sline.ReplaceAll("minRawTiming4H =","");
			minRawTiming4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "minRawTiming4H = " << minRawTiming4H << endl;
		}

		if (sline.Contains("maxRawTiming4H ="))
		{
			sline.ReplaceAll("maxRawTiming4H =","");
			maxRawTiming4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "maxRawTiming4H = " << maxRawTiming4H << endl;
		}

		if (sline.Contains("nBinsRawTiming4H ="))
		{
			sline.ReplaceAll("nBinsRawTiming4H =","");
			nBinsRawTiming4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nBinsRawTiming4H = " << nBinsRawTiming4H << endl;
		}
		
		if (sline.Contains("minDeltaRawTiming4H ="))
		{
			sline.ReplaceAll("minDeltaRawTiming4H =","");
			minDeltaRawTiming4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "minDeltaRawTiming4H = " << minDeltaRawTiming4H << endl;
		}

		if (sline.Contains("maxDeltaRawTiming4H ="))
		{
			sline.ReplaceAll("maxDeltaRawTiming4H =","");
			maxDeltaRawTiming4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "maxDeltaRawTiming4H = " << maxDeltaRawTiming4H << endl;
		}

		if (sline.Contains("nBinsDeltaRawTiming4H ="))
		{
			sline.ReplaceAll("nBinsDeltaRawTiming4H =","");
			nBinsDeltaRawTiming4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nBinsDeltaRawTiming4H = " << nBinsDeltaRawTiming4H << endl;
		}
		
		if (sline.Contains("minAnodeTiming4H ="))
		{
			sline.ReplaceAll("minAnodeTiming4H =","");
			minAnodeTiming4H = sline.Atof();
			if (printOutSetupFileParameters) cout << "minAnodeTiming4H = " << minAnodeTiming4H << endl;
		}

		if (sline.Contains("maxAnodeTiming4H ="))
		{
			sline.ReplaceAll("maxAnodeTiming4H =","");
			maxAnodeTiming4H = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxAnodeTiming4H = " << maxAnodeTiming4H << endl;
		}

		if (sline.Contains("nBinsAnodeTiming4H ="))
		{
			sline.ReplaceAll("nBinsAnodeTiming4H =","");
			nBinsAnodeTiming4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nBinsAnodeTiming4H = " << nBinsAnodeTiming4H << endl;
		}
		
		if (sline.Contains("minDeltaAnodeTiming4H ="))
		{
			sline.ReplaceAll("minDeltaAnodeTiming4H =","");
			minDeltaAnodeTiming4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "minDeltaAnodeTiming4H = " << minDeltaAnodeTiming4H << endl;
		}

		if (sline.Contains("maxDeltaAnodeTiming4H ="))
		{
			sline.ReplaceAll("maxDeltaAnodeTiming4H =","");
			maxDeltaAnodeTiming4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "maxDeltaAnodeTiming4H = " << maxDeltaAnodeTiming4H << endl;
		}

		if (sline.Contains("nBinsDeltaAnodeTiming4H ="))
		{
			sline.ReplaceAll("nBinsDeltaAnodeTiming4H =","");
			nBinsDeltaAnodeTiming4H = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nBinsDeltaAnodeTiming4H = " << nBinsDeltaAnodeTiming4H << endl;
		}

		if (sline.Contains("maxDtBetweenPixel4Clustering ="))
		{
			sline.ReplaceAll("maxDtBetweenPixel4Clustering =","");
			maxDtBetweenPixel4Clustering = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxDtBetweenPixel4Clustering = " << maxDtBetweenPixel4Clustering << endl;
		}
		
		if (sline.Contains("applyMaxDtCut4Clustering ="))
		{
			sline.ReplaceAll("applyMaxDtCut4Clustering =","");
			buf = sline.Atoi();
			applyMaxDtCut4Clustering = kTRUE;
			if (buf == 0) applyMaxDtCut4Clustering = kFALSE;
			if (printOutSetupFileParameters) cout << "applyMaxDtCut4Clustering = " << applyMaxDtCut4Clustering << endl;
		}
		
		if (sline.Contains("doNotUseCornerPixelsInPixelClusterCOGEneg ="))
		{
			sline.ReplaceAll("doNotUseCornerPixelsInPixelClusterCOGEneg =","");
			buf = sline.Atoi();
			doNotUseCornerPixelsInPixelClusterCOGEneg = kTRUE;
			if (buf == 0) doNotUseCornerPixelsInPixelClusterCOGEneg = kFALSE;
			if (printOutSetupFileParameters) cout << "doNotUseCornerPixelsInPixelClusterCOGEneg = " << doNotUseCornerPixelsInPixelClusterCOGEneg << endl;
		}

		if (sline.Contains("minNClusters4PairAnalysis ="))
		{
			sline.ReplaceAll("minNClusters4PairAnalysis =","");
			minNClusters4PairAnalysis = sline.Atoi();
			if (printOutSetupFileParameters) cout << "minNClusters4PairAnalysis = " << minNClusters4PairAnalysis << endl;
		}

		if (sline.Contains("maxNClusters4PairAnalysis ="))
		{
			sline.ReplaceAll("maxNClusters4PairAnalysis =","");
			maxNClusters4PairAnalysis = sline.Atoi();
			if (printOutSetupFileParameters) cout << "maxNClusters4PairAnalysis = " << maxNClusters4PairAnalysis << endl;
		}		
		
		if (sline.Contains("plotPixelClusterSpectra ="))
		{
			sline.ReplaceAll("plotPixelClusterSpectra =","");
			buf = sline.Atoi();
			plotPixelClusterSpectra = kTRUE;
			if (buf == 0) plotPixelClusterSpectra = kFALSE;
			if (printOutSetupFileParameters) cout << "plotPixelClusterSpectra = " << plotPixelClusterSpectra << endl;
		}
		
		if (sline.Contains("saveSelectedComptonRootTree ="))
		{
			sline.ReplaceAll("saveSelectedComptonRootTree =","");
			buf = sline.Atoi();
			saveSelectedComptonRootTree = kTRUE;
			if (buf == 0) saveSelectedComptonRootTree = kFALSE;
			if (printOutSetupFileParameters) cout << "saveSelectedComptonRootTree = " << saveSelectedComptonRootTree << endl;
		}
		
		if (sline.Contains("savingPointLocation ="))
		{
			sline.ReplaceAll("savingPointLocation =","");
			savingPointLocation = sline.Atoi();
			if (printOutSetupFileParameters) cout << "savingPointLocation = " << savingPointLocation << endl;
		}
		
		if (sline.Contains("factor2ConvertAnodeTime2Distance ="))
		{
			sline.ReplaceAll("factor2ConvertAnodeTime2Distance =","");
			factor2ConvertAnodeTime2Distance = sline.Atof();
			if (printOutSetupFileParameters) cout << "factor2ConvertAnodeTime2Distance = " << factor2ConvertAnodeTime2Distance << endl;
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
			if (doNotUseCornerTouchingPixels4Clustering && (k-ix == m-iy || k-ix == iy-m)) continue;
			if (applyMaxDtCut4Clustering && TMath::Abs(matrix_timeCalib[ix][iy] - matrix_timeCalib[k][m]) > maxDtBetweenPixel4Clustering) continue;
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

Float_t updateSinglePixelClusterCOGwithEneg(const Int_t ix, const Int_t iy, Float_t &x1, Float_t &y1)
{
	x1 = 0;
	y1 = 0;
	Float_t tamp = 0;
	for (Int_t k = ix-1; k <= ix+1; k++)
	{
		if (k < 0) continue;
		if (k > 10) continue;
		if (ix == 0 && k == 1) continue;
		if (ix == 10 && k == 9) continue;
		for (Int_t m = iy-1; m <= iy+1; m++)
		{
			if (m < 11) continue;
			if (m > 21) continue;
			if (iy == 11 && m == 12) continue;
			if (iy == 21 && m == 20) continue;
			if (k == ix && m == iy) continue;
			if (matrix_Eneg[k][m] < 0) continue;
			tamp += matrix_Eneg[k][m];
			x1 += (k+0.5)*matrix_Eneg[k][m];
			y1 += (m+0.5)*matrix_Eneg[k][m];
		}
	}
	return tamp;
}

void updateTwoPixelClusterCOGwithEneg(const Int_t am, const Float_t ixc, const Float_t iyc, const Int_t ix1, const Int_t iy1, const Int_t ix2, const Int_t iy2, Float_t &x1, Float_t &y1)
{
	x1 = ixc;
	y1 = iyc;
	if (pixelStatusGlobal[am][ix1][iy1] < 1 && pixelStatusGlobal[am][ix2][iy2] < 1) return;
	Float_t totNE1 = 0;
	Float_t totNE2 = 0;
	if (ix1 == ix2) // vertical cluster
	{
		if (iy1 > iy2)
		{
			if (!updateTwoPixelClusterCOG_1Donly && iy1 <= 20) totNE1 += matrix_Eneg[ix1-1][iy1+1];
			totNE1 += matrix_Eneg[ix1-1][iy1];
			totNE1 += matrix_Eneg[ix1-1][iy2];
			if (!updateTwoPixelClusterCOG_1Donly && iy2 >= 11) totNE1 += matrix_Eneg[ix1-1][iy2-1];
			if (!updateTwoPixelClusterCOG_1Donly && iy1 <= 20) totNE2 += matrix_Eneg[ix1+1][iy1+1];
			totNE2 += matrix_Eneg[ix1+1][iy1];
			totNE2 += matrix_Eneg[ix1+1][iy2];
			if (!updateTwoPixelClusterCOG_1Donly && iy2 >= 11) totNE2 += matrix_Eneg[ix1+1][iy2-1];
		}
		if (iy1 < iy2)
		{
			if (!updateTwoPixelClusterCOG_1Donly && iy2 <= 20) totNE1 += matrix_Eneg[ix1-1][iy2+1];
			totNE1 += matrix_Eneg[ix1-1][iy2];
			totNE1 += matrix_Eneg[ix1-1][iy1];
			if (!updateTwoPixelClusterCOG_1Donly && iy1 >= 11) totNE1 += matrix_Eneg[ix1-1][iy1-1];
			if (!updateTwoPixelClusterCOG_1Donly && iy2 <= 20) totNE2 += matrix_Eneg[ix1+1][iy2+1];
			totNE2 += matrix_Eneg[ix1+1][iy2];
			totNE2 += matrix_Eneg[ix1+1][iy1];
			if (!updateTwoPixelClusterCOG_1Donly && iy1 >= 11) totNE2 += matrix_Eneg[ix1+1][iy1-1];
		}
		x1 = (totNE1*(ix1-0.5)+totNE2*(ix1+1.5))/(totNE1+totNE2);
	}
	if (iy1 == iy2) // horizontal cluster
	{
		if (ix1 > ix2)
		{
			if (!updateTwoPixelClusterCOG_1Donly && ix2 >= 1) totNE1 += matrix_Eneg[ix2-1][iy1-1];
			totNE1 += matrix_Eneg[ix2][iy1-1];
			totNE1 += matrix_Eneg[ix1][iy1-1];
			if (!updateTwoPixelClusterCOG_1Donly && ix1 <= 9) totNE1 += matrix_Eneg[ix1+1][iy1-1];
			if (!updateTwoPixelClusterCOG_1Donly && ix2 >= 1) totNE2 += matrix_Eneg[ix2-1][iy1+1];
			totNE2 += matrix_Eneg[ix2][iy1+1];
			totNE2 += matrix_Eneg[ix1][iy1+1];
			if (!updateTwoPixelClusterCOG_1Donly && ix1 <= 9) totNE2 += matrix_Eneg[ix1+1][iy1+1];
		}
		if (ix1 < ix2)
		{
			if (!updateTwoPixelClusterCOG_1Donly && ix1 >= 1) totNE1 += matrix_Eneg[ix1-1][iy1-1];
			totNE1 += matrix_Eneg[ix1][iy1-1];
			totNE1 += matrix_Eneg[ix2][iy1-1];
			if (!updateTwoPixelClusterCOG_1Donly && ix2 <= 9) totNE1 += matrix_Eneg[ix2+1][iy1-1];
			if (!updateTwoPixelClusterCOG_1Donly && ix1 >= 1) totNE2 += matrix_Eneg[ix1-1][iy1+1];
			totNE2 += matrix_Eneg[ix1][iy1+1];
			totNE2 += matrix_Eneg[ix2][iy1+1];
			if (!updateTwoPixelClusterCOG_1Donly && ix2 <= 9) totNE2 += matrix_Eneg[ix2+1][iy1+1];
		}
		y1 = (totNE1*(iy1-0.5)+totNE2*(iy1+1.5))/(totNE1+totNE2);
	}
}

Double_t fitfun(Double_t *x, Double_t *par)
{
	return par[0]*TMath::Cos(2*x[0]*TMath::DegToRad()-TMath::Pi())+par[1];
}

void sortClusters(const Int_t am)
{
	firstClusterIdx[am] = 1;
	secondClusterIdx[am] = 0;
	if (buffClusterTrigs[am][firstClusterIdx[am]] >= 3)
	{
		firstClusterIdx[am] = 0;
		secondClusterIdx[am] = 1;
	}
	if (buffClusterIsSplit[am][firstClusterIdx[am]] == 1 && buffClusterIsSplit[am][secondClusterIdx[am]] == 0)
	{
		firstClusterIdx[am] = 0;
		secondClusterIdx[am] = 1;
	}
	if (buffClusterX[am].size() > 2)
	{
		firstClusterIdx[am] = 0;
		Float_t minDist = -9999;
		Int_t idx = -1;
		for (Int_t p = 1; p < buffClusterX[am].size(); p++)
		{
			Float_t dist = TMath::Hypot(buffClusterX[am][0]-buffClusterX[am][p],buffClusterY[am][0]-buffClusterY[am][p]);
			if (dist > minDist)
			{
				minDist = dist;
				idx = p;
			}
		}
		secondClusterIdx[am] = idx;
	}
}