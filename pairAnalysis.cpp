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

const TString ver = "test";
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
Float_t pixelPitch = 0.8;
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
Float_t XRF_Cd = 23.17;
Float_t XRF_Te = 27.47;
Float_t angleOfSecondHead = 90;
Float_t relativePhiAngle = 90;
Float_t Ewindow4PhiAnalyis_min = 530;
Float_t Ewindow4PhiAnalyis_max = 530;
Float_t Theta1WindowFor4PhiAnalyis_min = 60;
Float_t Theta1WindowFor4PhiAnalyis_max = 90;
Float_t Theta2WindowFor4PhiAnalyis_min = 60;
Float_t Theta2WindowFor4PhiAnalyis_max = 90;
Bool_t makeAsymmetricThetaWindow = kFALSE;
Bool_t makeCosineFit = kFALSE;
Float_t *factor2ConvertAnodeTime2Distance;
Float_t minZ = -5;
Float_t maxZ = 10;
Int_t nBinsZ = 200;

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
Bool_t doNotUseCornerTouchingTriggeredPixels4Clustering = kTRUE;
Bool_t doNotUseCornerTouchingPixels4Clustering = kTRUE;
Int_t maxClusterSizeToSkipFrame = 25;
Bool_t useOnlyTriggeredPixelsInCluster = kFALSE;
Int_t minClusterSize4PairAnalysis = -1;
Int_t maxClusterSize4PairAnalysis = -1;
Float_t maxDtBetweenPixel4Clustering = 200;
Bool_t applyMaxDtCut4Clustering = kFALSE;
Bool_t doNotUseCornerPixelsInPixelClusterCOGEneg = kFALSE;
Int_t minNClusters4PairAnalysis = 1;
Int_t maxNClusters4PairAnalysis = 2;

Bool_t enableUpdateEnergyClusters = kFALSE;

// Single pixel clusters
Float_t maxDtBetweenSinglePixClusters4updatingEnergy_low = 200;
Float_t maxDtBetweenSinglePixClusters4updatingEnergy_up = 200;
// Not required for updating the energy
Bool_t updateSinglePixelClusterCOG = kFALSE;

// Two pixel clusters
Float_t maxDtBetweenTwoPixClusters4updatingEnergy_low = 200;
Float_t maxDtBetweenTwoPixClusters4updatingEnergy_up = 200;
// Not required for updating the energy
Bool_t updateTwoPixelClusterCOG = kFALSE;
Bool_t updateTwoPixelClusterCOGshort_1D = kFALSE;
Bool_t updateTwoPixelClusterCOGlong_1D = kFALSE;

Bool_t eventHasToHaveCathodeSignal = 1;
Float_t maxCathodeEnergy = 540;

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
TH2I **pixelPattern;
TH2I **pixelTypePattern;
TH2I **usedCentrePixelPattern;

std::vector<std::vector<Int_t>> matrix_trigs;
std::vector<std::vector<Int_t>> matrix_flags;
std::vector<std::vector<Int_t>> matrix_time;
std::vector<std::vector<Float_t>> matrix_timeCalib;
std::vector<std::vector<Float_t>> matrix_E;
std::vector<std::vector<Int_t>> matrix_trigs_4sum;
std::vector<std::vector<Float_t>> matrix_Eall; // pixel energy + assigned neighbour
std::vector<std::vector<Float_t>> matrix_Eneg;
std::vector<Int_t> buffX;
std::vector<Int_t> buffY;
std::vector<Int_t> buffXnew;
std::vector<Int_t> buffYnew;

std::vector<std::vector<Float_t>> buffClusterX;
std::vector<std::vector<Float_t>> buffClusterY;
std::vector<std::vector<Float_t>> buffClusterE;
std::vector<std::vector<Float_t>> buffClusterEorig;
std::vector<std::vector<Float_t>> buffClusterArea;
std::vector<std::vector<Int_t>> buffClusterFlag;
std::vector<std::vector<Int_t>> buffClusterTrigs;
std::vector<std::vector<Int_t>> buffClusterIsSplit;
std::vector<std::vector<Int_t>> buffClusterOrientation;
std::vector<std::vector<Float_t>> buffClusterAnodeTiming;
std::vector<Float_t> buffClusterXloc;
std::vector<Float_t> buffClusterYloc;
std::vector<Float_t> buffClusterEloc;
std::vector<Float_t> buffClusterEorigloc;
std::vector<Float_t> buffClusterArealoc;
std::vector<Int_t> buffClusterFlagloc;
std::vector<Int_t> buffClusterTrigsloc;
std::vector<Int_t> buffClusterXarr;
std::vector<Int_t> buffClusterYarr;
std::vector<Int_t> buffClusterTimearr;
std::vector<Int_t> buffClusterAnodeTimearr;
std::vector<Int_t> buffClusterTrigsarr;
std::vector<Int_t> buffClusterPixelarr;
std::vector<Int_t> buffClusterIsSplitloc;
std::vector<Float_t> buffClusterAnodeTimingloc;
std::vector<Float_t> buff_E_neighb_4sum;
std::vector<Int_t> buffClusterOrientationloc;
std::vector<Int_t> buff_X_neighb_4sum;
std::vector<Int_t> buff_Y_neighb_4sum;

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
TLine *line1;
Bool_t *disabledPixels;
Int_t nDisabledPixels = 0;

TH1F **cathodeSpecAllEvents;
TH1F **cathodeSpecPair;
TH1F **cathodeSpecSelPairs;
TH1F **allClustersSpec;
TH1F **nTrigsInEvent;
TH1F **comptonSummedSpec2Clusters;
TH1F **comptonSummedSpec2Clusters_1PixClusters;
TH1F **comptonSummedSpec2Clusters_2PixClusters;
TH1F **comptonSummedSpec2Clusters_1_2PixClusters;
TH1F **comptonSummedSpec2Clusters_Eunsummed;
TH1F **comptonSummedSpec2Clusters_1PixClusters_Eunsummed;
TH1F **comptonSummedSpec2Clusters_2PixClusters_Eunsummed;
TH1F **comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed;
TH1F **comptonSpecClusters;
TH1F **comptonSpecClusters_Eunsummed;
TH2F **comptonSpec2ClustersCorr;
TH2F **clusterSpec2ClusterEventsCorr;
TH2F *comptonSummedSpec2ClustersCorr;
TH1F *comptonSummedSpec2Clusters2Heads;
TH2F *comptonSummedSpec2ClustersSelEventsCorr;
TH1F *comptonSummedSpec2ClustersSelEvents2Heads;
TH1F **phiAngle;
TH1F **phiAngleSelEvents;
TH1F *dPhiAngle;
TH1F *dPhiAngle_Ewin;
TH1F *dPhiAngle_w;
TH1F *dPhiAngle_Ewin_w;
TH1F *dPhiAngleNorm;
TH1F *dPhiAngleNorm_Ewin;
TH1F *dPhiAngleNorm_w;
TH1F *dPhiAngleNorm_Ewin_w;
TH1F *dPhiAngle1;
TH1F *dPhiAngle1_Ewin;
TH1F *dPhiAngle1norm;
TH1F *dPhiAngle1norm_Ewin;
TH2F *dPhiAngle_2ClusterE;
TH1F *dPhiXYtAngle;
TH1F *dPhiXYtAngle_Ewin;
TH1F *dPhiXYtAngleNorm;
TH1F *dPhiXYtAngleNorm_Ewin;
TH1F *dPhiXYtAngle1;
TH1F *dPhiXYtAngle1_Ewin;
TH1F *dPhiXYtAngle1Norm;
TH1F *dPhiXYtAngle1Norm_Ewin;
TH2F *phiAngleCorr;
TH2F *phiAngleSelEventsCorr;
TH2F **dPhiPlaneAM;
TH2F **firstComptonCluster;
TH2F **secondComptonCluster;
TH2F **clusterE_vsSize;
TH1F **thetaFromFirstClusterE;
TH2F **thetaFromFirstClusterE_dPhi;
TH2F **thetaFromFirstClusterE_dPhi_w;
TH2F **thetaFromXYtClusterE_dPhi;
TH2F **thetaFromXYtClusterE_dPhi_w;
TH1F **thetaFromXYtE;
TH1F **thetaFromXYtE_w;
TH2F **theta_firstCluster_XYtE_Corr;
TH2F **theta_firstCluster_XYtE_Corr_w;
TH1F **clusterSize;
TH2F *theta0_theta1_FromFirstClusterE;
TH2F **clusterSize_nClusters;
TH1F **numberOfClusters;
TH2F **clusterMaxE_vsSize;
TH2F **clusterSize2ClusterEventsCorr;
TH2F **imageE;
TH2F **imageT;
TH2F **imageC;
TH2I **imageN;
TH1F **nClustersInEvent;
TH1F **nClustersInSelEvent;
TH1F **firstComptonClusterSize;
TH1F **secondComptonClusterSize;
TH1F **firstComptonClusterSizeSelEvents;
TH1F **secondComptonClusterSizeSelEvents;
TH1F **thetaFromFirstClusterE_w;
TH1F **thetaFromFirstClusterE_1pixOnly;
TH1F **thetaFromFirstClusterE_1pixOnly_w;
TH1F **thetaFromFirstClusterE_2pixOnly;
TH1F **thetaFromFirstClusterE_1_2pixOnly;
TH1F **anodeTimingDiff2PixDiagClusters;
TH1F **anodeTimingDiff2PixClusters;
TH1F **anodeTimingDiffAllClusters;
TH2F **anodeTiming2PixDiagClustersCorr;
TH2F **anodeTiming2PixClustersCorr;
TH2F **allClustersCOGImage;
TH2F **allPixelsImage;
TH2F **allClustersFineCOGImage;
TH2F **allPixelsFreqImage;
TH2F **allClustersCOGFreqImage;
TH2F **allClustersFineCOGFreqImage;
TH1F ***clusterSpec;
TH2F **dAnodeTiming_vs_theta;
TH2F **dAnodeTiming_vs_E1;
TH2F **dAnodeTiming_vs_E2;
TH2F **dAnodeTiming_vs_E12;
TH2F **dAnodeTimingAbs_vs_E12;
TH2F **anodeTiming1PixelClustersNeigbCorr;
TH2F **anodeTiming2PixelClustersNeigbCorr;
TH2F **dAnodeTiming_vs_nE_1PixelClustersNeigb;
TH2F **dAnodeTiming_vs_sumE_1PixelClustersNeigb;
TH2F **dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters;
TH1F **pixelSpecNeighbour1PixClusters;
TH2F **energySpecCorr1PixClusters;
TH2F **energySpecCorr2PixClusters;
TH1F **clusterSpec_1ClusterEvents_summed;
TH1F **clusterSpec_1ClusterEvents;
TH1F **firstComptonClusterE;
TH1F **secondComptonClusterE;
TH2F **anodeTimingNeigb_nE_1PixelClusters;
TH2F **anodeTimingNeigb_nE_2PixelClusters;
TH1F **clusterZfromAnodeTiming;
TH1F **comptonEventsDeltaZpix;
TH1F **comptonEventsDeltaZdist;
TH1F **clusterSummingStats1PixClusters;
TH1F **clusterSummingStats2PixClusters;

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
void sortClusters(const Int_t);
Double_t fitMAC(Double_t *, Double_t *);
void updateClustersWithNeighboursEnergy(const Int_t, const Int_t);
Float_t getThetaAngleFromHitsXYt(const Int_t, const Float_t, const Float_t, const Float_t, const Float_t, const Float_t, const Float_t);

const Int_t nMAC = 11;
const Int_t nMACbd = 5;
const Double_t e4MAC[nMAC] = {40,50,60,80,100,150,200,300,400,500,600}; 
const Double_t MAC[nMAC] = {19.3,10.67,6.542,3.019,1.671,0.6071,0.3246,0.1628,0.1147,0.09291,0.08042}; 
TH1F * prob_h1;
TH1F * prob2_h1;

TFile *hfile1;//hfile1 = new TFile(Form("%s%s/histos.root",pathRoot.Data(),subdir.Data()),"recreate");
TFile *hfile2;//hfile1 = new TFile(Form("%s%s/histos.root",pathRoot.Data(),subdir.Data()),"recreate");
TTree* tree = new TTree("tree","Output events tree");

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
	
	if ((updateTwoPixelClusterCOG && updateTwoPixelClusterCOGshort_1D) || (updateTwoPixelClusterCOG && updateTwoPixelClusterCOGlong_1D) || (updateTwoPixelClusterCOGshort_1D && updateTwoPixelClusterCOGlong_1D))
	{
		cerr << "ERROR: Wrong settings for two pixel cluster COG update, only one of the options could be enabled:" << endl;
		cerr << "updateTwoPixelClusterCOG = " << updateTwoPixelClusterCOG << endl; 
		cerr << "updateTwoPixelClusterCOGshort_1D = " << updateTwoPixelClusterCOGshort_1D << endl; 
		cerr << "updateTwoPixelClusterCOGlong_1D = " << updateTwoPixelClusterCOGlong_1D << endl;
		cerr << "Exiting." << endl;
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
	
	TString subdir = "Pair_analysis_" + ver;
	subdir = Form("%s_pixType%d",subdir.Data(),typeOfPixelsToUse);
	if (useOnlyTriggeredPixelsInCluster) subdir += "_onlyTrigPixels";
	if (minClusterSize4PairAnalysis > 0 && maxClusterSize4PairAnalysis > 0)
	{
		if (minClusterSize4PairAnalysis == maxClusterSize4PairAnalysis) subdir += Form("_%d-pixCls",minClusterSize4PairAnalysis);
		if (minClusterSize4PairAnalysis < maxClusterSize4PairAnalysis) subdir += Form("_%d-%d-pixCls",minClusterSize4PairAnalysis,maxClusterSize4PairAnalysis);
	}
	if (minNClusters4PairAnalysis > 0 && maxNClusters4PairAnalysis > 0)
	{
		if (minNClusters4PairAnalysis == maxNClusters4PairAnalysis) subdir += Form("_%d-nCls",minNClusters4PairAnalysis);
		if (minNClusters4PairAnalysis < maxNClusters4PairAnalysis) subdir += Form("_%d-%d-nCls",minNClusters4PairAnalysis,maxNClusters4PairAnalysis);
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
	
	TFile* f_ref = new TFile("probabilites_ref.root", "read");
	prob_h1 = (TH1F *) f_ref->Get("probThetaRec_E1_vs_E1");
	prob2_h1 = (TH1F *) f_ref->Get("prob2ThetaRec_E1_vs_E1");
	prob_h1->SetName("prob_h1");
	prob_h1->SetDirectory(0);
	prob2_h1->SetName("prob2_h1");
	prob2_h1->SetDirectory(0);
	f_ref->Close();
	
	hfile1 = new TFile(Form("%s%s/histos.root",pathRoot.Data(),subdir.Data()),"recreate");
	hfile2 = new TFile(Form("%s%s/studentOutput.root",pathRoot.Data(),subdir.Data()),"recreate");
//	tree->Branch("dPhi",&dPh,"dPhi/F");
//	tree->Branch("theta0",&theta0,"theta0/F");
//	tree->Branch("theta1",&theta1,"theta1/F");
//	tree->Branch("energy0",&sum2clE0,"energy0/F");
//	tree->Branch("energy1",&sum2clE1,"energy1/F");
//	tree->Branch("thetaCut",&thetaWindowCut,"thetaCut/O");
//	tree->Branch("energyCut",&EnergyWindowCut,"energyCut/O");
	
	image = new TH2F("image","Event image",nPixXY,0,nPixXY,nPixXY,0,nPixXY);
	image->GetXaxis()->SetTitle("Pixel number");
	image->GetXaxis()->SetTitleOffset(1.2);
	image->GetYaxis()->SetTitle("Pixel number");
	image->GetYaxis()->SetTitleOffset(1.2);
	
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
	numberOfClusters = new TH1F*[nAMs];
	clusterSize2ClusterEventsCorr = new TH2F*[nAMs];
	clusterSpec_1ClusterEvents = new TH1F*[nAMs];
	cathodeSpecAllEvents = new TH1F*[nAMs];
	cathodeSpecPair = new TH1F*[nAMs];
	cathodeSpecSelPairs = new TH1F*[nAMs];
	allClustersSpec = new TH1F*[nAMs];
	pixelSpecNeighbour1PixClusters = new TH1F*[nAMs];
	clusterSpec_1ClusterEvents_summed = new TH1F*[nAMs];
	comptonSpecClusters_Eunsummed = new TH1F*[nAMs];
	comptonSummedSpec2Clusters_Eunsummed = new TH1F*[nAMs];
	comptonSummedSpec2Clusters_1PixClusters_Eunsummed = new TH1F*[nAMs];
	comptonSummedSpec2Clusters_2PixClusters_Eunsummed = new TH1F*[nAMs];
	comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed = new TH1F*[nAMs];
	if (enableUpdateEnergyClusters)
	{		
		comptonSpecClusters = new TH1F*[nAMs];
		comptonSummedSpec2Clusters = new TH1F*[nAMs];
		comptonSummedSpec2Clusters_1PixClusters = new TH1F*[nAMs];
		comptonSummedSpec2Clusters_2PixClusters = new TH1F*[nAMs];
		comptonSummedSpec2Clusters_1_2PixClusters = new TH1F*[nAMs];
		energySpecCorr1PixClusters = new TH2F*[nAMs];
		energySpecCorr2PixClusters = new TH2F*[nAMs];
		clusterSummingStats1PixClusters = new TH1F*[nAMs];
		clusterSummingStats2PixClusters = new TH1F*[nAMs];
	}
	comptonSpec2ClustersCorr = new TH2F*[nAMs];
	clusterSpec2ClusterEventsCorr = new TH2F*[nAMs];
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
	thetaFromFirstClusterE_dPhi_w = new TH2F*[nAMs];
	thetaFromXYtClusterE_dPhi = new TH2F*[nAMs];
	thetaFromXYtClusterE_dPhi_w = new TH2F*[nAMs];
	thetaFromXYtE = new TH1F*[nAMs];
	thetaFromXYtE_w = new TH1F*[nAMs];
	theta_firstCluster_XYtE_Corr = new TH2F*[nAMs];
	theta_firstCluster_XYtE_Corr_w = new TH2F*[nAMs];
	clusterSize = new TH1F*[nAMs];
	clusterSize_nClusters = new TH2F*[nAMs];
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
	firstComptonClusterE = new TH1F*[nAMs];
	secondComptonClusterE = new TH1F*[nAMs];
	comptonEventsDeltaZpix = new TH1F*[nAMs];
	if (plotPixelClusterSpectra) clusterSpec = new TH1F**[nAMs];
	comptonEventsDeltaZdist = new TH1F*[nAMs];
	if (makeTimingCalibrationStuff)
	{
		clusterZfromAnodeTiming = new TH1F*[nAMs];
		anodeTimingDiff2PixClusters = new TH1F*[nAMs];
		anodeTimingDiffAllClusters = new TH1F*[nAMs];
		if (!doNotUseCornerPixelsInPixelClusterCOGEneg)
		{
			anodeTiming2PixDiagClustersCorr = new TH2F*[nAMs];
			anodeTimingDiff2PixDiagClusters = new TH1F*[nAMs];
		}
		anodeTiming2PixClustersCorr = new TH2F*[nAMs];
		dAnodeTiming_vs_theta = new TH2F*[nAMs];
		dAnodeTiming_vs_E1 = new TH2F*[nAMs];
		dAnodeTiming_vs_E2 = new TH2F*[nAMs];
		dAnodeTiming_vs_E12 = new TH2F*[nAMs];
		dAnodeTimingAbs_vs_E12 = new TH2F*[nAMs];
		if (enableUpdateEnergyClusters)
		{		
			dAnodeTiming_vs_nE_1PixelClustersNeigb = new TH2F*[nAMs];
			dAnodeTiming_vs_sumE_1PixelClustersNeigb = new TH2F*[nAMs];
			dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters = new TH2F*[nAMs];
			anodeTiming1PixelClustersNeigbCorr = new TH2F*[nAMs];
			anodeTiming2PixelClustersNeigbCorr = new TH2F*[nAMs];
			anodeTimingNeigb_nE_1PixelClusters = new TH2F*[nAMs];
			anodeTimingNeigb_nE_2PixelClusters = new TH2F*[nAMs];
		}
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

		thetaFromFirstClusterE_dPhi_w[im] = new TH2F(Form("thetaFromFirstClusterE_dPhi_w_AM%d",im),Form("%s, #theta calculated from first interaction energy, weighted, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180,nBins_dPhi,-180,180);
		thetaFromFirstClusterE_dPhi_w[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromFirstClusterE_dPhi_w[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));
		thetaFromFirstClusterE_dPhi_w[im]->GetYaxis()->SetTitleOffset(1.4);
		thetaFromFirstClusterE_dPhi_w[im]->GetYaxis()->SetTitle(Form("#Delta#varphi AM%d, degrees",im));

		thetaFromXYtClusterE_dPhi[im] = new TH2F(Form("thetaFromXYtClusterE_dPhi_AM%d",im),Form("%s, #theta calculated from first interaction energy, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180,nBins_dPhi,-180,180);
		thetaFromXYtClusterE_dPhi[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromXYtClusterE_dPhi[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));
		thetaFromXYtClusterE_dPhi[im]->GetYaxis()->SetTitleOffset(1.4);
		thetaFromXYtClusterE_dPhi[im]->GetYaxis()->SetTitle(Form("#Delta#varphi AM%d, degrees",im));

		thetaFromXYtClusterE_dPhi_w[im] = new TH2F(Form("thetaFromXYtClusterE_dPhi_w_AM%d",im),Form("%s, #theta calculated from first interaction energy, weighted, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180,nBins_dPhi,-180,180);
		thetaFromXYtClusterE_dPhi_w[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromXYtClusterE_dPhi_w[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));
		thetaFromXYtClusterE_dPhi_w[im]->GetYaxis()->SetTitleOffset(1.4);
		thetaFromXYtClusterE_dPhi_w[im]->GetYaxis()->SetTitle(Form("#Delta#varphi AM%d, degrees",im));
		
		thetaFromXYtE[im] = new TH1F(Form("thetaFromXYtE_AM%d",im),Form("%s, #theta calculated from X-Y-t, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180);
		thetaFromXYtE[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromXYtE[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));

		thetaFromXYtE_w[im] = new TH1F(Form("thetaFromXYtE_w_AM%d",im),Form("%s, #theta calculated from X-Y-t, weighted, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180);
		thetaFromXYtE_w[im]->GetXaxis()->SetTitleOffset(1.2);
		thetaFromXYtE_w[im]->GetXaxis()->SetTitle(Form("#theta in AM%d, degrees",im));
	
		theta_firstCluster_XYtE_Corr[im] = new TH2F(Form("theta_firstCluster_XYtE_Corr_AM%d",im),Form("%s, #theta calculated from X-Y-t, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180,nBinsTheta4H,0,180);
		theta_firstCluster_XYtE_Corr[im]->GetXaxis()->SetTitleOffset(1.2);
		theta_firstCluster_XYtE_Corr[im]->GetXaxis()->SetTitle(Form("#theta using first cluster energy, degrees, AM%d",im));
		theta_firstCluster_XYtE_Corr[im]->GetYaxis()->SetTitleOffset(1.4);
		theta_firstCluster_XYtE_Corr[im]->GetYaxis()->SetTitle(Form("#theta using X-Y-t, degrees, AM%d",im));
	
		theta_firstCluster_XYtE_Corr_w[im] = new TH2F(Form("theta_firstCluster_XYtE_Corr_w_AM%d",im),Form("%s, #theta calculated from X-Y-t, weighted, AM%d",spectrumName.Data(),im),nBinsTheta4H,0,180,nBinsTheta4H,0,180);
		theta_firstCluster_XYtE_Corr_w[im]->GetXaxis()->SetTitleOffset(1.2);
		theta_firstCluster_XYtE_Corr_w[im]->GetXaxis()->SetTitle(Form("#theta using first cluster energy, degrees, AM%d",im));
		theta_firstCluster_XYtE_Corr_w[im]->GetYaxis()->SetTitleOffset(1.4);
		theta_firstCluster_XYtE_Corr_w[im]->GetYaxis()->SetTitle(Form("#theta using X-Y-t, degrees, AM%d",im));
	
		nTrigsInEvent[im] = new TH1F(Form("nTrigsInEvent_AM%d",im),Form("%s, number of triggers per event, AM%d",spectrumName.Data(),im),maxNTrigs4H-minNTrigs4H,minNTrigs4H,maxNTrigs4H);
		nTrigsInEvent[im]->GetXaxis()->SetTitleOffset(1.2);
		nTrigsInEvent[im]->GetXaxis()->SetTitle(Form("Number of triggers in AM%d",im));

		nClustersInEvent[im] = new TH1F(Form("nClustersInEvent_AM%d",im),Form("%s, number of clusters per event, AM%d",spectrumName.Data(),im),maxNTrigs4H-minNTrigs4H,minNTrigs4H,maxNTrigs4H);
		nClustersInEvent[im]->GetXaxis()->SetTitleOffset(1.2);
		nClustersInEvent[im]->GetXaxis()->SetTitle(Form("Number of clusters in AM%d",im));

		nClustersInSelEvent[im] = new TH1F(Form("nClustersInSelEvent_AM%d",im),Form("%s, number of clusters/event in selected events, AM%d",spectrumName.Data(),im),maxNTrigs4H-minNTrigs4H,minNTrigs4H,maxNTrigs4H);
		nClustersInSelEvent[im]->GetXaxis()->SetTitleOffset(1.2);
		nClustersInSelEvent[im]->GetXaxis()->SetTitle(Form("Number of clusters in AM%d",im));
		
		clusterSize_nClusters[im] = new TH2F(Form("clusterSize_nClusters_AM%d",im),Form("%s, number of triggers in cluster vs cluster size, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H,maxNTrigs4H_2D-minNTrigs4H_2D,minNTrigs4H_2D,maxNTrigs4H_2D);
		clusterSize_nClusters[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterSize_nClusters[im]->GetXaxis()->SetTitle(Form("Cluster size, pixels, AM%d",im));	
		clusterSize_nClusters[im]->GetYaxis()->SetTitleOffset(1.4);
		clusterSize_nClusters[im]->GetYaxis()->SetTitle("Number of clusters/event");
			
		comptonEventsDeltaZpix[im] = new TH1F(Form("comptonEventsDeltaZpix_AM%d",im),Form("%s, #Deltaz between Compton event clusters, AM%d",spectrumName.Data(),im),12,0,12);
		comptonEventsDeltaZpix[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonEventsDeltaZpix[im]->GetXaxis()->SetTitle(Form("Distance, pixels, AM%d",im));
		
		if (makeTimingCalibrationStuff)
		{
			if (!doNotUseCornerPixelsInPixelClusterCOGEneg) 
			{
				anodeTimingDiff2PixDiagClusters[im] = new TH1F(Form("anodeTimingDiff2PixDiagClusters_AM%d",im),Form("%s, Difference in triggered pixel anode timing, two triggered diagonal pixel clusters only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H);
				anodeTimingDiff2PixDiagClusters[im]->GetXaxis()->SetTitleOffset(1.2);
				anodeTimingDiff2PixDiagClusters[im]->GetXaxis()->SetTitle(Form("Anode timing difference between pixels, ns, AM%d",im));	
			
				anodeTiming2PixDiagClustersCorr[im] = new TH2F(Form("anodeTiming2PixDiagClustersCorr_AM%d",im),Form("%s, Correlation between triggered pixel anode timing, two triggered diagonal pixel clusters only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H);
				anodeTiming2PixDiagClustersCorr[im]->GetXaxis()->SetTitleOffset(1.2);
				anodeTiming2PixDiagClustersCorr[im]->GetXaxis()->SetTitle(Form("Anode timing, first pixel, ns, AM%d",im));	
				anodeTiming2PixDiagClustersCorr[im]->GetYaxis()->SetTitleOffset(1.4);
				anodeTiming2PixDiagClustersCorr[im]->GetYaxis()->SetTitle(Form("Anode timing, second pixel, ns, AM%d",im));	
			}
			
			anodeTimingDiff2PixClusters[im] = new TH1F(Form("anodeTimingDiff2PixClusters_AM%d",im),Form("%s, Difference in triggered pixel anode timing, 2-pix clusters only, before cut on #Deltat, AM%d",spectrumName.Data(),im),nBinsDeltaAnodeTiming4H,minDeltaAnodeTiming4H,maxDeltaAnodeTiming4H);
			anodeTimingDiff2PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			anodeTimingDiff2PixClusters[im]->GetXaxis()->SetTitle(Form("Anode timing difference between pixels, ns, AM%d",im));	
			
			anodeTimingDiffAllClusters[im] = new TH1F(Form("anodeTimingDiffAllClusters_AM%d",im),Form("%s, Difference in triggered pixel anode timing, all clusters, AM%d",spectrumName.Data(),im),nBinsDeltaAnodeTiming4H,minDeltaAnodeTiming4H,maxDeltaAnodeTiming4H);
			anodeTimingDiffAllClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			anodeTimingDiffAllClusters[im]->GetXaxis()->SetTitle(Form("Anode timing difference between pixels, ns, AM%d",im));	
			
			if (enableUpdateEnergyClusters)
			{
				dAnodeTiming_vs_nE_1PixelClustersNeigb[im] = new TH2F(Form("dAnodeTiming_vs_nE_1PixelClustersNeigb_AM%d",im),Form("%s, #Deltat(triggered pixel - neighbour) vs its energy, 1-pix clusters only, AM%d",spectrumName.Data(),im),nBinsE,minE,50,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H+200);
				dAnodeTiming_vs_nE_1PixelClustersNeigb[im]->GetXaxis()->SetTitleOffset(1.2);
				dAnodeTiming_vs_nE_1PixelClustersNeigb[im]->GetXaxis()->SetTitle(Form("Neighbour pixel energy, keV, AM%d",im));	
				dAnodeTiming_vs_nE_1PixelClustersNeigb[im]->GetYaxis()->SetTitleOffset(1.5);
				dAnodeTiming_vs_nE_1PixelClustersNeigb[im]->GetYaxis()->SetTitle(Form("#Deltat(triggered pixel - neighbour), ns, AM%d",im));	
				
				dAnodeTiming_vs_sumE_1PixelClustersNeigb[im] = new TH2F(Form("dAnodeTiming_vs_sumE_1PixelClustersNeigb_AM%d",im),Form("%s, #Deltat(triggered pixel - neighbour) vs  neighbour pixel + triggered pixel energy, 1-pix clusters, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H+200);
				dAnodeTiming_vs_sumE_1PixelClustersNeigb[im]->GetXaxis()->SetTitleOffset(1.2);
				dAnodeTiming_vs_sumE_1PixelClustersNeigb[im]->GetXaxis()->SetTitle(Form("Neighbour pixel + triggered pixel energy, keV, AM%d",im));	
				dAnodeTiming_vs_sumE_1PixelClustersNeigb[im]->GetYaxis()->SetTitleOffset(1.5);
				dAnodeTiming_vs_sumE_1PixelClustersNeigb[im]->GetYaxis()->SetTitle(Form("#Deltat(triggered pixel - neighbour), ns, AM%d",im));
				
				dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters[im] = new TH2F(Form("dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters_AM%d",im),Form("%s, #Deltat(triggered pixel - neighbour) vs neighbour anode timing, 1-pix clusters only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H+200,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H+200);
				dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters[im]->GetXaxis()->SetTitleOffset(1.2);
				dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters[im]->GetXaxis()->SetTitle(Form("Anode timing, neighbour pixel, ns, AM%d",im));	
				dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters[im]->GetYaxis()->SetTitleOffset(1.5);
				dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters[im]->GetYaxis()->SetTitle(Form("#Deltat(triggered pixel - neighbour), ns, AM%d",im));	

				anodeTiming1PixelClustersNeigbCorr[im] = new TH2F(Form("anodeTiming1PixelClustersNeigbCorr_AM%d",im),Form("%s, Correlation between triggered pixel anode timing and its neighbours, 1-pix clusters only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H+200,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H+200);
				anodeTiming1PixelClustersNeigbCorr[im]->GetXaxis()->SetTitleOffset(1.2);
				anodeTiming1PixelClustersNeigbCorr[im]->GetXaxis()->SetTitle(Form("Anode timing, main pixel, ns, AM%d",im));	
				anodeTiming1PixelClustersNeigbCorr[im]->GetYaxis()->SetTitleOffset(1.5);
				anodeTiming1PixelClustersNeigbCorr[im]->GetYaxis()->SetTitle(Form("Anode timing, neighbour pixel, ns, AM%d",im));	

				anodeTiming2PixelClustersNeigbCorr[im] = new TH2F(Form("anodeTiming2PixelClustersNeigbCorr_AM%d",im),Form("%s, Correlation between triggered pixel anode timing and its neighbours, 2-pix clusters only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H+200,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H+200);
				anodeTiming2PixelClustersNeigbCorr[im]->GetXaxis()->SetTitleOffset(1.2);
				anodeTiming2PixelClustersNeigbCorr[im]->GetXaxis()->SetTitle(Form("Anode timing, main pixel, ns, AM%d",im));	
				anodeTiming2PixelClustersNeigbCorr[im]->GetYaxis()->SetTitleOffset(1.5);
				anodeTiming2PixelClustersNeigbCorr[im]->GetYaxis()->SetTitle(Form("Anode timing, neighbour pixel, ns, AM%d",im));
				
				anodeTimingNeigb_nE_1PixelClusters[im] = new TH2F(Form("anodeTimingNeigb_nE_1PixelClusters_AM%d",im),Form("%s, Correlation between neighbour energy and anode timing, 1-pix clusters only, AM%d",spectrumName.Data(),im),nBinsE,minE,50,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H+200);
				anodeTimingNeigb_nE_1PixelClusters[im]->GetXaxis()->SetTitleOffset(1.2);
				anodeTimingNeigb_nE_1PixelClusters[im]->GetXaxis()->SetTitle(Form("Neighbour energy, AM%d, keV",im));	
				anodeTimingNeigb_nE_1PixelClusters[im]->GetYaxis()->SetTitleOffset(1.5);
				anodeTimingNeigb_nE_1PixelClusters[im]->GetYaxis()->SetTitle(Form("Anode timing, neighbour pixel, ns, AM%d",im));	
				
				anodeTimingNeigb_nE_2PixelClusters[im] = new TH2F(Form("anodeTimingNeigb_nE_2PixelClusters_AM%d",im),Form("%s, Correlation between neighbour energy and anode timing, 2-pix clusters only, AM%d",spectrumName.Data(),im),nBinsE,minE,50,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H+200);
				anodeTimingNeigb_nE_2PixelClusters[im]->GetXaxis()->SetTitleOffset(1.2);
				anodeTimingNeigb_nE_2PixelClusters[im]->GetXaxis()->SetTitle(Form("Neighbour energy, AM%d, keV",im));	
				anodeTimingNeigb_nE_2PixelClusters[im]->GetYaxis()->SetTitleOffset(1.5);
				anodeTimingNeigb_nE_2PixelClusters[im]->GetYaxis()->SetTitle(Form("Anode timing, neighbour pixel, ns, AM%d",im));					
			}
			
			anodeTiming2PixClustersCorr[im] = new TH2F(Form("anodeTiming2PixClustersCorr_AM%d",im),Form("%s, Correlation between triggered pixel anode timing, two triggered pixel clusters only, AM%d",spectrumName.Data(),im),nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H,nBinsAnodeTiming4H,minAnodeTiming4H,maxAnodeTiming4H);
			anodeTiming2PixClustersCorr[im]->GetXaxis()->SetTitleOffset(1.2);
			anodeTiming2PixClustersCorr[im]->GetXaxis()->SetTitle(Form("Anode timing, first pixel, ns, AM%d",im));	
			anodeTiming2PixClustersCorr[im]->GetYaxis()->SetTitleOffset(1.4);
			anodeTiming2PixClustersCorr[im]->GetYaxis()->SetTitle(Form("Anode timing, second pixel, ns, AM%d",im));	

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
			dAnodeTiming_vs_E2[im]->GetYaxis()->SetTitle(Form("#Delta(Anode timing), ns, AM%d",im));
			
			dAnodeTiming_vs_E12[im] = new TH2F(Form("dAnodeTiming_vs_E12_AM%d",im),Form("%s, #Delta(Anode timing) = t_{1} - t_{2} vs both cluster energies, two cluster events only, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,100,-500,500);
			dAnodeTiming_vs_E12[im]->GetXaxis()->SetTitleOffset(1.2);
			dAnodeTiming_vs_E12[im]->GetXaxis()->SetTitle(Form("Cluster energy, AM%d, keV",im));
			dAnodeTiming_vs_E12[im]->GetYaxis()->SetTitleOffset(1.5);
			dAnodeTiming_vs_E12[im]->GetYaxis()->SetTitle(Form("#Delta(Anode timing), ns, AM%d",im));
			
			dAnodeTimingAbs_vs_E12[im] = new TH2F(Form("dAnodeTimingAbs_vs_E12_AM%d",im),Form("%s, |#Delta(Anode timing)| vs both cluster energies, two cluster events only, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,100,0,500);
			dAnodeTimingAbs_vs_E12[im]->GetXaxis()->SetTitleOffset(1.2);
			dAnodeTimingAbs_vs_E12[im]->GetXaxis()->SetTitle(Form("Cluster energy, AM%d, keV",im));
			dAnodeTimingAbs_vs_E12[im]->GetYaxis()->SetTitleOffset(1.5);
			dAnodeTimingAbs_vs_E12[im]->GetYaxis()->SetTitle(Form("|#Delta(Anode timing)|, ns, AM%d",im));

			clusterZfromAnodeTiming[im] = new TH1F(Form("clusterZfromAnodeTiming_AM%d",im),Form("%s, Cluster Z-coordinate calculated from anode timing, 1-2-pix clusters, AM%d",spectrumName.Data(),im),nBinsZ,minZ,maxZ);
			clusterZfromAnodeTiming[im]->GetXaxis()->SetTitleOffset(1.2);
			clusterZfromAnodeTiming[im]->GetXaxis()->SetTitle(Form("Cluster Z, mm, AM%d",im));			
		}
		comptonEventsDeltaZdist[im] = new TH1F(Form("comptonEventsDeltaZdist_AM%d",im),Form("%s, #Deltaz between Compton event clusters using anode timing, AM%d",spectrumName.Data(),im),40,0,10);
		comptonEventsDeltaZdist[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonEventsDeltaZdist[im]->GetXaxis()->SetTitle(Form("Distance, mm, AM%d",im));
		
		clusterSize[im] = new TH1F(Form("clusterSize_AM%d",im),Form("%s, cluster size, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		clusterSize[im]->GetXaxis()->SetTitleOffset(1.4);
		clusterSize[im]->GetXaxis()->SetTitle(Form("Cluster size, pixels, AM%d",im));		
		
		numberOfClusters[im] = new TH1F(Form("numberOfClusters_AM%d",im),Form("%s, number of clusters per event, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		numberOfClusters[im]->GetXaxis()->SetTitleOffset(1.4);
		numberOfClusters[im]->GetXaxis()->SetTitle(Form("Number of clusters/event, AM%d",im));		
		
		clusterE_vsSize[im] = new TH2F(Form("clusterE_vsSize_AM%d",im),Form("%s, cluster size vs cluster energy, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		clusterE_vsSize[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterE_vsSize[im]->GetXaxis()->SetTitle(Form("Cluster energy, keV, AM%d",im));		
		clusterE_vsSize[im]->GetYaxis()->SetTitleOffset(1.4);
		clusterE_vsSize[im]->GetYaxis()->SetTitle(Form("Cluster size, pixels, AM%d",im));			
		
		clusterMaxE_vsSize[im] = new TH2F(Form("clusterMaxE_vsSize_AM%d",im),Form("%s, max energy cluster size vs its size, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		clusterMaxE_vsSize[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterMaxE_vsSize[im]->GetXaxis()->SetTitle(Form("Cluster energy, keV, AM%d",im));
		clusterMaxE_vsSize[im]->GetYaxis()->SetTitleOffset(1.4);
		clusterMaxE_vsSize[im]->GetYaxis()->SetTitle(Form("Cluster size, pixels, AM%d",im));		
	
		clusterSize2ClusterEventsCorr[im] = new TH2F(Form("clusterSize2ClusterEventsCorr_AM%d",im),Form("%s, first vs second cluster size, 2 cluster events, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H,maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		clusterSize2ClusterEventsCorr[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterSize2ClusterEventsCorr[im]->GetXaxis()->SetTitle(Form("First cluster size, pixels, AM%d",im));		
		clusterSize2ClusterEventsCorr[im]->GetYaxis()->SetTitleOffset(1.4);
		clusterSize2ClusterEventsCorr[im]->GetYaxis()->SetTitle(Form("Second cluster size, pixels, AM%d",im));		
		
		allClustersSpec[im] = new TH1F(Form("allClustersSpec_AM%d",im),Form("%s, all cluster spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		allClustersSpec[im]->GetXaxis()->SetTitleOffset(1.2);
		allClustersSpec[im]->GetXaxis()->SetTitle(Form("Cluster energy, keV, AM%d",im));	

		pixelSpecNeighbour1PixClusters[im] = new TH1F(Form("pixelSpecNeighbour1PixClusters_AM%d",im),Form("%s, energy spectrum of neighbour with maximum energy, single pixel cluster, AM%d",spectrumName.Data(),im),120,0,60);
		pixelSpecNeighbour1PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
		pixelSpecNeighbour1PixClusters[im]->GetXaxis()->SetTitle(Form("Cluster energy, keV, AM%d",im));	
		
		clusterSpec_1ClusterEvents_summed[im] = new TH1F(Form("clusterSpec_1ClusterEvents_summed_AM%d",im),Form("%s, Summed energy spectrum (max neighbour + cluster), 1 cluster events, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		clusterSpec_1ClusterEvents_summed[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterSpec_1ClusterEvents_summed[im]->GetXaxis()->SetTitle(Form("Summed energy, keV, AM%d",im));		
		
		firstComptonClusterE[im] = new TH1F(Form("firstComptonClusterE_AM%d",im),Form("%s, First Compton cluster spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		firstComptonClusterE[im]->GetXaxis()->SetTitleOffset(1.2);
		firstComptonClusterE[im]->GetXaxis()->SetTitle(Form("Energy, keV, AM%d",im));		
		
		secondComptonClusterE[im] = new TH1F(Form("secondComptonClusterE_AM%d",im),Form("%s, Second Compton cluster spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		secondComptonClusterE[im]->GetXaxis()->SetTitleOffset(1.2);
		secondComptonClusterE[im]->GetXaxis()->SetTitle(Form("Energy, keV, AM%d",im));		
		
		clusterSpec_1ClusterEvents[im] = new TH1F(Form("clusterSpec_1ClusterEvents_AM%d",im),Form("%s, unsummed energy spectrum, 1 cluster events, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		clusterSpec_1ClusterEvents[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterSpec_1ClusterEvents[im]->GetXaxis()->SetTitle(Form("Unsummed energy, keV, AM%d",im));		
		
		cathodeSpecAllEvents[im] = new TH1F(Form("cathodeSpecAllEvents_AM%d",im),Form("%s, cathode spectrum, all evetns, AM%d",spectrumName.Data(),im),200,0,600);
		cathodeSpecAllEvents[im]->GetXaxis()->SetTitleOffset(1.2);
		cathodeSpecAllEvents[im]->GetXaxis()->SetTitle(Form("Cathode energy, keV, AM%d",im));		

		cathodeSpecPair[im] = new TH1F(Form("cathodeSpecPair_AM%d",im),Form("%s, cathode spectrum, Compton pair events, AM%d",spectrumName.Data(),im),200,0,600);
		cathodeSpecPair[im]->GetXaxis()->SetTitleOffset(1.2);
		cathodeSpecPair[im]->GetXaxis()->SetTitle(Form("Cathode energy, keV, AM%d",im));		

		cathodeSpecSelPairs[im] = new TH1F(Form("cathodeSpecSelPairs_AM%d",im),Form("%s, cathode spectrum, selected Compton pair events, AM%d",spectrumName.Data(),im),200,0,600);
		cathodeSpecSelPairs[im]->GetXaxis()->SetTitleOffset(1.2);
		cathodeSpecSelPairs[im]->GetXaxis()->SetTitle(Form("Cathode energy, keV, AM%d",im));		

		comptonSpec2ClustersCorr[im] = new TH2F(Form("comptonSpec2ClustersCorr_AM%d",im),Form("%s, unsummed Compton interactions, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,nBinsE,minE,maxE);
		comptonSpec2ClustersCorr[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonSpec2ClustersCorr[im]->GetXaxis()->SetTitle(Form("First cluster energy, keV, AM%d",im));		
		comptonSpec2ClustersCorr[im]->GetYaxis()->SetTitleOffset(1.4);
		comptonSpec2ClustersCorr[im]->GetYaxis()->SetTitle(Form("Second cluster energy, keV, AM%d",im));		

		clusterSpec2ClusterEventsCorr[im] = new TH2F(Form("clusterSpec2ClusterEventsCorr_AM%d",im),Form("%s, cluster energies, two cluster events, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,nBinsE,minE,maxE);
		clusterSpec2ClusterEventsCorr[im]->GetXaxis()->SetTitleOffset(1.2);
		clusterSpec2ClusterEventsCorr[im]->GetXaxis()->SetTitle(Form("First cluster energy, keV, AM%d",im));		
		clusterSpec2ClusterEventsCorr[im]->GetYaxis()->SetTitleOffset(1.4);
		clusterSpec2ClusterEventsCorr[im]->GetYaxis()->SetTitle(Form("Second cluster energy, keV, AM%d",im));		

		firstComptonClusterSize[im] = new TH1F(Form("firstComptonClusterSize_AM%d",im),Form("%s, cluster size of first Compton cluster, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		firstComptonClusterSize[im]->GetXaxis()->SetTitleOffset(1.4);
		firstComptonClusterSize[im]->GetXaxis()->SetTitle(Form("Cluster size, pixels, AM%d",im));		

		secondComptonClusterSize[im] = new TH1F(Form("secondComptonClusterSize_AM%d",im),Form("%s, cluster size of second Compton cluster, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		secondComptonClusterSize[im]->GetXaxis()->SetTitleOffset(1.4);
		secondComptonClusterSize[im]->GetXaxis()->SetTitle(Form("Cluster size, pixels, AM%d",im));	

		firstComptonClusterSizeSelEvents[im] = new TH1F(Form("firstComptonClusterSizeSelEvents_AM%d",im),Form("%s, cluster size of first Compton cluster, #Delta#theta selected events, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		firstComptonClusterSizeSelEvents[im]->GetXaxis()->SetTitleOffset(1.4);
		firstComptonClusterSizeSelEvents[im]->GetXaxis()->SetTitle(Form("Cluster size, pixels, AM%d",im));		

		secondComptonClusterSizeSelEvents[im] = new TH1F(Form("secondComptonClusterSizeSelEvents_AM%d",im),Form("%s, cluster size of second Compton cluster, #Delta#theta selected events, AM%d",spectrumName.Data(),im),maxClusterSize4H-minClusterSize4H,minClusterSize4H,maxClusterSize4H);
		secondComptonClusterSizeSelEvents[im]->GetXaxis()->SetTitleOffset(1.4);
		secondComptonClusterSizeSelEvents[im]->GetXaxis()->SetTitle(Form("Cluster size, pixels, AM%d",im));		
		
		if (enableUpdateEnergyClusters)
		{
			comptonSpecClusters[im] = new TH1F(Form("comptonSpecClusters_AM%d",im),Form("%s, unsummed Compton interaction cluster spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
			comptonSpecClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			comptonSpecClusters[im]->GetXaxis()->SetTitle(Form("Cluster energy, keV, AM%d",im));		

			comptonSummedSpec2Clusters[im] = new TH1F(Form("comptonSummedSpec2Clusters_AM%d",im),Form("%s, two interaction summed Compton spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
			comptonSummedSpec2Clusters[im]->GetXaxis()->SetTitleOffset(1.2);
			comptonSummedSpec2Clusters[im]->GetXaxis()->SetTitle(Form("Summed event energy, keV, AM%d",im));		
			
			comptonSummedSpec2Clusters_1PixClusters[im] = new TH1F(Form("comptonSummedSpec2Clusters_1PixClusters_AM%d",im),Form("%s, two interaction summed Compton spectrum, 1-pix clusters, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
			comptonSummedSpec2Clusters_1PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			comptonSummedSpec2Clusters_1PixClusters[im]->GetXaxis()->SetTitle(Form("Summed event energy, keV, AM%d",im));		
			
			comptonSummedSpec2Clusters_2PixClusters[im] = new TH1F(Form("comptonSummedSpec2Clusters_2PixClusters_AM%d",im),Form("%s, two interaction summed Compton spectrum, 2-pix clusters, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
			comptonSummedSpec2Clusters_2PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			comptonSummedSpec2Clusters_2PixClusters[im]->GetXaxis()->SetTitle(Form("Summed event energy, keV, AM%d",im));		
			
			comptonSummedSpec2Clusters_1_2PixClusters[im] = new TH1F(Form("comptonSummedSpec2Clusters_1_2PixClusters_AM%d",im),Form("%s, two interaction summed Compton spectrum, 1-2-pix clusters, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
			comptonSummedSpec2Clusters_1_2PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			comptonSummedSpec2Clusters_1_2PixClusters[im]->GetXaxis()->SetTitle(Form("Summed event energy, keV, AM%d",im));		
		
			energySpecCorr1PixClusters[im] = new TH2F(Form("energySpecCorr1PixClusters_AM%d",im),Form("%s, single pixel cluster energy vs neighbour pixel energy, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,120,0,60);
			energySpecCorr1PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			energySpecCorr1PixClusters[im]->GetXaxis()->SetTitle(Form("Cluster energy, keV, AM%d",im));		
			energySpecCorr1PixClusters[im]->GetYaxis()->SetTitleOffset(1.4);
			energySpecCorr1PixClusters[im]->GetYaxis()->SetTitle(Form("Neighbour pixel energy, keV, AM%d",im));		
		
			energySpecCorr2PixClusters[im] = new TH2F(Form("energySpecCorr2PixClusters_AM%d",im),Form("%s, two pixel cluster energy vs total neighbours pixel energy, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE,120,0,60);
			energySpecCorr2PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			energySpecCorr2PixClusters[im]->GetXaxis()->SetTitle(Form("Cluster energy, keV, AM%d",im));		
			energySpecCorr2PixClusters[im]->GetYaxis()->SetTitleOffset(1.4);
			energySpecCorr2PixClusters[im]->GetYaxis()->SetTitle(Form("Total neighbour pixel energy, keV, AM%d",im));

			clusterSummingStats1PixClusters[im] = new TH1F(Form("clusterSummingStats1PixClusters_AM%d",im),Form("%s, cluster summing statistics, 1-pix clusters, AM%d",spectrumName.Data(),im),5,0,5);
			clusterSummingStats1PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			clusterSummingStats1PixClusters[im]->GetXaxis()->SetTitle(Form("Number of neighbours added, AM%d",im));		

			clusterSummingStats2PixClusters[im] = new TH1F(Form("clusterSummingStats2PixClusters_AM%d",im),Form("%s, cluster summing statistics, 2-pix clusters, AM%d",spectrumName.Data(),im),5,0,5);
			clusterSummingStats2PixClusters[im]->GetXaxis()->SetTitleOffset(1.2);
			clusterSummingStats2PixClusters[im]->GetXaxis()->SetTitle(Form("Number of neighbours added, AM%d",im));		
		}
		
		comptonSummedSpec2Clusters_Eunsummed[im] = new TH1F(Form("comptonSummedSpec2Clusters_Eunsummed_AM%d",im),Form("%s, two interaction summed Compton spectrum, unsummed clusters, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		comptonSummedSpec2Clusters_Eunsummed[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonSummedSpec2Clusters_Eunsummed[im]->GetXaxis()->SetTitle(Form("Summed event energy, keV, AM%d",im));		
		
		comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im] = new TH1F(Form("comptonSummedSpec2Clusters_1PixClusters_Eunsummed_AM%d",im),Form("%s, two interaction summed Compton spectrum, 1-pix unsummed clusters, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im]->GetXaxis()->SetTitle(Form("Summed event energy, keV, AM%d",im));		
		
		comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im] = new TH1F(Form("comptonSummedSpec2Clusters_2PixClusters_Eunsummed_AM%d",im),Form("%s, two interaction summed Compton spectrum, 2-pix unsummed clusters, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im]->GetXaxis()->SetTitle(Form("Summed event energy, keV, AM%d",im));			
		
		comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im] = new TH1F(Form("comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed_AM%d",im),Form("%s, two interaction summed Compton spectrum, 1-2-pix unsummed clusters, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im]->GetXaxis()->SetTitle(Form("Summed event energy, keV, AM%d",im));			

		comptonSpecClusters_Eunsummed[im] = new TH1F(Form("comptonSpecClusters_Eunsummed_AM%d",im),Form("%s, unsummed Compton interaction cluster spectrum, AM%d",spectrumName.Data(),im),nBinsE,minE,maxE);
		comptonSpecClusters_Eunsummed[im]->GetXaxis()->SetTitleOffset(1.2);
		comptonSpecClusters_Eunsummed[im]->GetXaxis()->SetTitle(Form("Unsummed cluster energy, keV, AM%d",im));		
		
		phiAngle[im] = new TH1F(Form("phiAngle_AM%d",im),Form("%s, cluster azimuthal angle, AM%d",spectrumName.Data(),im),360,0,360);
		phiAngle[im]->GetXaxis()->SetTitleOffset(1.2);
		phiAngle[im]->GetXaxis()->SetTitle(Form("#varphi AM%d, degrees",im));

		phiAngleSelEvents[im] = new TH1F(Form("phiAngleSelEvents_AM%d",im),Form("%s, cluster azimuthal angle, selected events, #Delta#theta selected events, AM%d",spectrumName.Data(),im),360,0,360);
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
	
	comptonSummedSpec2Clusters2Heads = new TH1F("comptonSummedSpec2Clusters2Heads",Form("%s, two interaction summed Compton spectrum, sum of two heads",spectrumName.Data()),nBinsE*2,minE,maxE*2);
	comptonSummedSpec2Clusters2Heads->GetXaxis()->SetTitleOffset(1.2);
	comptonSummedSpec2Clusters2Heads->GetXaxis()->SetTitle("Summed total event energy, keV");		
		
	comptonSummedSpec2ClustersCorr = new TH2F("comptonSummedSpec2ClustersCorr",Form("%s, two interaction summed Compton spectra in AMs",spectrumName.Data()),nBinsE,minE,maxE,nBinsE,minE,maxE);
	comptonSummedSpec2ClustersCorr->GetXaxis()->SetTitleOffset(1.2);
	comptonSummedSpec2ClustersCorr->GetXaxis()->SetTitle("AM0 summed cluster energy, keV");
	comptonSummedSpec2ClustersCorr->GetYaxis()->SetTitleOffset(1.4);
	comptonSummedSpec2ClustersCorr->GetYaxis()->SetTitle("AM1 summed cluster energy, keV");	
	
	comptonSummedSpec2ClustersSelEvents2Heads = new TH1F("comptonSummedSpec2ClustersSelEvents2Heads",Form("%s, two interaction summed Compton spectrum, sum of two heads, #Delta#theta selected events",spectrumName.Data()),nBinsE*2,minE,maxE*2);
	comptonSummedSpec2ClustersSelEvents2Heads->GetXaxis()->SetTitleOffset(1.2);
	comptonSummedSpec2ClustersSelEvents2Heads->GetXaxis()->SetTitle("Summed total event energy, keV");		
		
	comptonSummedSpec2ClustersSelEventsCorr = new TH2F("comptonSummedSpec2ClustersSelEventsCorr",Form("%s, two interaction summed Compton spectra in AMs, #Delta#theta selected events",spectrumName.Data()),nBinsE,minE,maxE,nBinsE,minE,maxE);
	comptonSummedSpec2ClustersSelEventsCorr->GetXaxis()->SetTitleOffset(1.2);
	comptonSummedSpec2ClustersSelEventsCorr->GetXaxis()->SetTitle("AM0 summed cluster energy, keV");
	comptonSummedSpec2ClustersSelEventsCorr->GetYaxis()->SetTitleOffset(1.4);
	comptonSummedSpec2ClustersSelEventsCorr->GetYaxis()->SetTitle("AM1 summed cluster energy, keV");	

	dPhiXYtAngle = new TH1F("dPhiXYtAngle",Form("%s, #Delta#varphi using X-Y-t #theta window",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiXYtAngle->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngle->GetXaxis()->SetTitle("#Delta#varphi, degrees");		

	dPhiXYtAngle_Ewin = new TH1F("dPhiXYtAngle_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, using X-Y-t #theta window",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi,-180,180);
	dPhiXYtAngle_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngle_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		

	dPhiXYtAngleNorm = new TH1F("dPhiXYtAngleNorm",Form("%s, #Delta#varphi, normalised at #Delta#varphi = 0, using X-Y-t #theta window",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiXYtAngleNorm->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngleNorm->GetXaxis()->SetTitle("#Delta#varphi, degrees");		

	dPhiXYtAngleNorm_Ewin = new TH1F("dPhiXYtAngleNorm_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0, using X-Y-t #theta window",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi,-180,180);
	dPhiXYtAngleNorm_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngleNorm_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		

	dPhiXYtAngle1 = new TH1F("dPhiXYtAngle1",Form("%s, #Delta#varphi using X-Y-t #theta window",spectrumName.Data()),nBins_dPhi/2,0,180);
	dPhiXYtAngle1->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngle1->GetXaxis()->SetTitle("#Delta#varphi, degrees");		

	dPhiXYtAngle1_Ewin = new TH1F("dPhiXYtAngle1_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, using X-Y-t #theta window",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi/2,0,180);
	dPhiXYtAngle1_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngle1_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		

	dPhiXYtAngle1Norm = new TH1F("dPhiXYtAngle1Norm",Form("%s, #Delta#varphi, normalised at #Delta#varphi = 0, using X-Y-t #theta window",spectrumName.Data()),nBins_dPhi/2,0,180);
	dPhiXYtAngle1Norm->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngle1Norm->GetXaxis()->SetTitle("#Delta#varphi, degrees");		

	dPhiXYtAngle1Norm_Ewin = new TH1F("dPhiXYtAngle1Norm_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0, using X-Y-t #theta window",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi/2,0,180);
	dPhiXYtAngle1Norm_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngle1Norm_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle = new TH1F("dPhiAngle",Form("%s, #Delta#varphi",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiAngle->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle_w = new TH1F("dPhiAngle_w",Form("%s, #Delta#varphi, weighted",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiAngle_w->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_w->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngleNorm = new TH1F("dPhiAngleNorm",Form("%s, #Delta#varphi, normalised at #Delta#varphi = 0",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiAngleNorm->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngleNorm->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngleNorm_w = new TH1F("dPhiAngleNorm_w",Form("%s, #Delta#varphi, normalised at #Delta#varphi = 0, weighted",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiAngleNorm_w->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngleNorm_w->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle1 = new TH1F("dPhiAngle1",Form("%s, |#Delta#varphi|",spectrumName.Data()),nBins_dPhi/2,0,180);
	dPhiAngle1->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle1->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle1norm = new TH1F("dPhiAngle1norm",Form("%s, |#Delta#varphi|, normalised at #Delta#varphi = 0",spectrumName.Data()),nBins_dPhi/2,0,180);
	dPhiAngle1norm->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle1norm->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle_Ewin = new TH1F("dPhiAngle_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi,-180,180);
	dPhiAngle_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle_Ewin_w = new TH1F("dPhiAngle_Ewin_w",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, weighted",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi,-180,180);
	dPhiAngle_Ewin_w->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_Ewin_w->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngleNorm_Ewin = new TH1F("dPhiAngleNorm_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi,-180,180);
	dPhiAngleNorm_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngleNorm_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngleNorm_Ewin_w = new TH1F("dPhiAngleNorm_Ewin_w",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0, weighted",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi,-180,180);
	dPhiAngleNorm_Ewin_w->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngleNorm_Ewin_w->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle1_Ewin = new TH1F("dPhiAngle1_Ewin",Form("%s, |#Delta#varphi| within  %.0f keV < E < %.0f keV",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi/2,0,180);
	dPhiAngle1_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle1_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle1norm_Ewin = new TH1F("dPhiAngle1norm_Ewin",Form("%s, |#Delta#varphi| within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi/2,0,180);
	dPhiAngle1norm_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle1norm_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle_2ClusterE = new TH2F("dPhiAngle_2ClusterE",Form("%s, #Delta#varphi vs summed energy of two first Compton clusters",spectrumName.Data()),nBinsE,minE,maxE,nBins_dPhi,-180,180);
	dPhiAngle_2ClusterE->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_2ClusterE->GetXaxis()->SetTitle("Total summed energy, keV");
	dPhiAngle_2ClusterE->GetYaxis()->SetTitleOffset(1.4);
	dPhiAngle_2ClusterE->GetYaxis()->SetTitle("#Delta#varphi, degrees");		
	
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
	matrix_Eall.clear();
	matrix_Eall.resize(nPixXY, vector<Float_t>(nPixXY));
	matrix_Eneg.clear();
	matrix_Eneg.resize(nPixXY, vector<Float_t>(nPixXY));
	matrix_time.clear();
	matrix_time.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_timeCalib.clear();
	matrix_timeCalib.resize(nPixXY, vector<Float_t>(nPixXY));
	matrix_trigs_4sum.clear();
	matrix_trigs_4sum.resize(nPixXY, vector<Int_t>(nPixXY));
	buff_E_neighb_4sum.clear();
	buff_E_neighb_4sum.resize(0,0);
	buff_X_neighb_4sum.clear();
	buff_X_neighb_4sum.resize(0,0);
	buff_Y_neighb_4sum.clear();
	buff_Y_neighb_4sum.resize(0,0);
		
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
	line1 = new TLine();
	line1->SetLineColor(1);
	line1->SetLineWidth(1);
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
	TH1F **numberOfClusters_norm = new TH1F*[nAMs];
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
		line1->DrawLine(comptonSpec2ClustersCorr[im]->GetXaxis()->GetXmin(),Na22Energy-pixelPattern[im]->GetXaxis()->GetXmin(),Na22Energy-pixelPattern[im]->GetYaxis()->GetXmin(),comptonSpec2ClustersCorr[im]->GetYaxis()->GetXmin());
		line1->DrawLine(comptonSpec2ClustersCorr[im]->GetXaxis()->GetXmin(),Na22Energy-XRF_Cd-pixelPattern[im]->GetXaxis()->GetXmin(),Na22Energy-XRF_Cd-pixelPattern[im]->GetYaxis()->GetXmin(),comptonSpec2ClustersCorr[im]->GetYaxis()->GetXmin());
		line1->DrawLine(comptonSpec2ClustersCorr[im]->GetXaxis()->GetXmin(),Na22Energy-XRF_Te-pixelPattern[im]->GetXaxis()->GetXmin(),Na22Energy-XRF_Te-pixelPattern[im]->GetYaxis()->GetXmin(),comptonSpec2ClustersCorr[im]->GetYaxis()->GetXmin());
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(comptonSpec2ClustersCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(comptonSpec2ClustersCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
					
		clusterSpec2ClusterEventsCorr[im]->Draw("colz");
		line1->DrawLine(clusterSpec2ClusterEventsCorr[im]->GetXaxis()->GetXmin(),Na22Energy-pixelPattern[im]->GetXaxis()->GetXmin(),Na22Energy-pixelPattern[im]->GetYaxis()->GetXmin(),clusterSpec2ClusterEventsCorr[im]->GetYaxis()->GetXmin());
		line1->DrawLine(clusterSpec2ClusterEventsCorr[im]->GetXaxis()->GetXmin(),Na22Energy-XRF_Cd-pixelPattern[im]->GetXaxis()->GetXmin(),Na22Energy-XRF_Cd-pixelPattern[im]->GetYaxis()->GetXmin(),clusterSpec2ClusterEventsCorr[im]->GetYaxis()->GetXmin());
		line1->DrawLine(clusterSpec2ClusterEventsCorr[im]->GetXaxis()->GetXmin(),Na22Energy-XRF_Te-pixelPattern[im]->GetXaxis()->GetXmin(),Na22Energy-XRF_Te-pixelPattern[im]->GetYaxis()->GetXmin(),clusterSpec2ClusterEventsCorr[im]->GetYaxis()->GetXmin());
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(clusterSpec2ClusterEventsCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(clusterSpec2ClusterEventsCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
	
		dPhiPlaneAM[im]->Draw("colz");
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dPhiPlaneAM[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dPhiPlaneAM[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
	
		clusterSize2ClusterEventsCorr[im]->Draw("colz");
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(clusterSize2ClusterEventsCorr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(clusterSize2ClusterEventsCorr[im]->GetName()).Data()));
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
		
		theta_firstCluster_XYtE_Corr[im]->Draw("colz");
		line1->DrawLine(0,0,180,180);
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(theta_firstCluster_XYtE_Corr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(theta_firstCluster_XYtE_Corr[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);

		theta_firstCluster_XYtE_Corr_w[im]->Draw("colz");
		line1->DrawLine(0,0,180,180);
		c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(theta_firstCluster_XYtE_Corr_w[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(1);
		c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(theta_firstCluster_XYtE_Corr_w[im]->GetName()).Data()));
		c2->GetPad(0)->SetLogz(0);
		
		c1->cd(0);
		allClustersSpec[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(allClustersSpec[im]->GetName()).Data()));

		pixelSpecNeighbour1PixClusters[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(pixelSpecNeighbour1PixClusters[im]->GetName()).Data()));

		clusterSpec_1ClusterEvents_summed[im]->Draw();
		Y_text_top = 0.82;
		getHistoPeak(clusterSpec_1ClusterEvents_summed[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(clusterSpec_1ClusterEvents_summed[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/clusterSpec_1ClusterEvents_summed[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(clusterSpec_1ClusterEvents_summed[im]->GetName()).Data()));
	
		firstComptonClusterE[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(firstComptonClusterE[im]->GetName()).Data()));
		
		secondComptonClusterE[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(secondComptonClusterE[im]->GetName()).Data()));

		firstComptonClusterE[im]->Draw();
		secondComptonClusterE[im]->SetLineColor(2);
		secondComptonClusterE[im]->Draw("same");
		c1->SaveAs(Form("%s/%s_overlay.gif",outputPathAM[im].Data(),TString(firstComptonClusterE[im]->GetName()).Data()));
		
		clusterSpec_1ClusterEvents[im]->Draw();
		Y_text_top = 0.82;
		getHistoPeak(clusterSpec_1ClusterEvents[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(clusterSpec_1ClusterEvents[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/clusterSpec_1ClusterEvents[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(clusterSpec_1ClusterEvents[im]->GetName()).Data()));
	
		cathodeSpecAllEvents[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(cathodeSpecAllEvents[im]->GetName()).Data()));
		
		cathodeSpecPair[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(cathodeSpecPair[im]->GetName()).Data()));
		
		cathodeSpecSelPairs[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(cathodeSpecSelPairs[im]->GetName()).Data()));
				
		comptonSpecClusters_Eunsummed[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSpecClusters_Eunsummed[im]->GetName()).Data()));
		
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
		
		clusterSize_nClusters[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(clusterSize_nClusters[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(clusterSize_nClusters[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		if (makeTimingCalibrationStuff)
		{			
			c2->cd(0);			
			if (enableUpdateEnergyClusters)
			{
				anodeTiming1PixelClustersNeigbCorr[im]->Draw("colz");
				if (enableUpdateEnergyClusters)
				{
					line2->DrawLine(anodeTiming1PixelClustersNeigbCorr[im]->GetXaxis()->GetXmin(),
						anodeTiming1PixelClustersNeigbCorr[im]->GetXaxis()->GetXmin()+maxDtBetweenSinglePixClusters4updatingEnergy_up,
							anodeTiming1PixelClustersNeigbCorr[im]->GetXaxis()->GetXmax()-maxDtBetweenSinglePixClusters4updatingEnergy_up,
								anodeTiming1PixelClustersNeigbCorr[im]->GetYaxis()->GetXmax());
					line2->DrawLine(anodeTiming1PixelClustersNeigbCorr[im]->GetXaxis()->GetXmin()+maxDtBetweenSinglePixClusters4updatingEnergy_low,
						anodeTiming1PixelClustersNeigbCorr[im]->GetYaxis()->GetXmin(),
							anodeTiming1PixelClustersNeigbCorr[im]->GetXaxis()->GetXmax(),
								anodeTiming1PixelClustersNeigbCorr[im]->GetXaxis()->GetXmax()-maxDtBetweenSinglePixClusters4updatingEnergy_low);
				}
				c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(anodeTiming1PixelClustersNeigbCorr[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(1);
				c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(anodeTiming1PixelClustersNeigbCorr[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(0);

				anodeTiming2PixelClustersNeigbCorr[im]->Draw("colz");
				if (enableUpdateEnergyClusters)
				{
					line2->DrawLine(anodeTiming2PixelClustersNeigbCorr[im]->GetXaxis()->GetXmin(),
						anodeTiming2PixelClustersNeigbCorr[im]->GetXaxis()->GetXmin()+maxDtBetweenTwoPixClusters4updatingEnergy_up,
							anodeTiming2PixelClustersNeigbCorr[im]->GetXaxis()->GetXmax()-maxDtBetweenTwoPixClusters4updatingEnergy_up,
								anodeTiming2PixelClustersNeigbCorr[im]->GetYaxis()->GetXmax());
					line2->DrawLine(anodeTiming2PixelClustersNeigbCorr[im]->GetXaxis()->GetXmin()+maxDtBetweenTwoPixClusters4updatingEnergy_low,
						anodeTiming2PixelClustersNeigbCorr[im]->GetYaxis()->GetXmin(),
							anodeTiming2PixelClustersNeigbCorr[im]->GetXaxis()->GetXmax(),
								anodeTiming2PixelClustersNeigbCorr[im]->GetXaxis()->GetXmax()-maxDtBetweenTwoPixClusters4updatingEnergy_low);
				}
				c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(anodeTiming2PixelClustersNeigbCorr[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(1);
				c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(anodeTiming2PixelClustersNeigbCorr[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(0);
				
				dAnodeTiming_vs_nE_1PixelClustersNeigb[im]->Draw("colz");
				c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_nE_1PixelClustersNeigb[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(1);
				c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_nE_1PixelClustersNeigb[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(0);
				
				dAnodeTiming_vs_sumE_1PixelClustersNeigb[im]->Draw("colz");
				c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_sumE_1PixelClustersNeigb[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(1);
				c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_sumE_1PixelClustersNeigb[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(0);
				
				dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters[im]->Draw("colz");
				c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(1);
				c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(0);
				
				anodeTimingNeigb_nE_1PixelClusters[im]->Draw("colz");
				c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(anodeTimingNeigb_nE_1PixelClusters[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(1);
				c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(anodeTimingNeigb_nE_1PixelClusters[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(0);
				
				anodeTimingNeigb_nE_2PixelClusters[im]->Draw("colz");
				c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(anodeTimingNeigb_nE_2PixelClusters[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(1);
				c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(anodeTimingNeigb_nE_2PixelClusters[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(0);
			}
			
			if (!doNotUseCornerPixelsInPixelClusterCOGEneg) 
			{
				anodeTiming2PixDiagClustersCorr[im]->Draw("colz");
				c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(anodeTiming2PixDiagClustersCorr[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(1);
				c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(anodeTiming2PixDiagClustersCorr[im]->GetName()).Data()));
				c2->GetPad(0)->SetLogz(0);
				
				c1->cd(0);
				anodeTimingDiff2PixDiagClusters[im]->Draw();
				c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(anodeTimingDiff2PixDiagClusters[im]->GetName()).Data()));
			}
			
			c2->cd(0);
			anodeTiming2PixClustersCorr[im]->Draw("colz");
			if (applyMaxDtCut4Clustering)
			{
				line2->DrawLine(anodeTiming2PixClustersCorr[im]->GetXaxis()->GetXmin(),
					anodeTiming2PixClustersCorr[im]->GetYaxis()->GetXmin()+maxDtBetweenPixel4Clustering,
						anodeTiming2PixClustersCorr[im]->GetXaxis()->GetXmax()-maxDtBetweenPixel4Clustering,
							anodeTiming2PixClustersCorr[im]->GetYaxis()->GetXmax());
				line2->DrawLine(anodeTiming2PixClustersCorr[im]->GetXaxis()->GetXmin()+maxDtBetweenPixel4Clustering,
					anodeTiming2PixClustersCorr[im]->GetYaxis()->GetXmin(),
						anodeTiming2PixClustersCorr[im]->GetXaxis()->GetXmax(),
							anodeTiming2PixClustersCorr[im]->GetXaxis()->GetXmax()-maxDtBetweenPixel4Clustering);
			}
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(anodeTiming2PixClustersCorr[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(anodeTiming2PixClustersCorr[im]->GetName()).Data()));
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
			
			dAnodeTiming_vs_E1[im]->Draw("colz");
			c2->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_E1[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(1);
			c2->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(dAnodeTiming_vs_E1[im]->GetName()).Data()));
			c2->GetPad(0)->SetLogz(0);
						
			c1->cd(0);
			anodeTimingDiff2PixClusters[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(anodeTimingDiff2PixClusters[im]->GetName()).Data()));
			
			anodeTimingDiffAllClusters[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(anodeTimingDiffAllClusters[im]->GetName()).Data()));			
					
			clusterZfromAnodeTiming[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(clusterZfromAnodeTiming[im]->GetName()).Data()));
			c1->GetPad(0)->SetLogy(1);
			c1->SaveAs(Form("%s/%s_logy.gif",outputPathAM[im].Data(),TString(clusterZfromAnodeTiming[im]->GetName()).Data()));
			c1->GetPad(0)->SetLogy(0);
		}
		
		comptonEventsDeltaZdist[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonEventsDeltaZdist[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogy(1);
		c1->SaveAs(Form("%s/%s_logy.gif",outputPathAM[im].Data(),TString(comptonEventsDeltaZdist[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogy(0);

		comptonEventsDeltaZpix[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonEventsDeltaZpix[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogy(1);
		c1->SaveAs(Form("%s/%s_logy.gif",outputPathAM[im].Data(),TString(comptonEventsDeltaZpix[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogy(0);

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
		
		if (enableUpdateEnergyClusters)
		{
			clusterSummingStats1PixClusters[im]->SetNormFactor(1);
			clusterSummingStats1PixClusters[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(clusterSummingStats1PixClusters[im]->GetName()).Data()));
			c1->GetPad(0)->SetLogy(1);
			c1->SaveAs(Form("%s/%s_logy.png",outputPathAM[im].Data(),TString(clusterSummingStats1PixClusters[im]->GetName()).Data()));
			c1->GetPad(0)->SetLogy(0);		

			clusterSummingStats2PixClusters[im]->SetNormFactor(1);
			clusterSummingStats2PixClusters[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(clusterSummingStats2PixClusters[im]->GetName()).Data()));
			c1->GetPad(0)->SetLogy(1);
			c1->SaveAs(Form("%s/%s_logy.png",outputPathAM[im].Data(),TString(clusterSummingStats2PixClusters[im]->GetName()).Data()));
			c1->GetPad(0)->SetLogy(0);		

			comptonSpecClusters[im]->Draw();
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSpecClusters[im]->GetName()).Data()));

			comptonSummedSpec2Clusters[im]->Draw();
			line1->DrawLine(Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters[im]->GetMinimum(),Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters[im]->GetMaximum()/3);
			line1->DrawLine(Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters[im]->GetMinimum(),Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters[im]->GetMaximum()/3);
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
		
			comptonSummedSpec2Clusters_1PixClusters[im]->Draw();
			line1->DrawLine(Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_1PixClusters[im]->GetMinimum(),Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_1PixClusters[im]->GetMaximum()/3);
			line1->DrawLine(Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_1PixClusters[im]->GetMinimum(),Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_1PixClusters[im]->GetMaximum()/3);
			Y_text_top = 0.82;
			getHistoPeak(comptonSummedSpec2Clusters_1PixClusters[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
			Peak_FWHM_raw = getFWHM(comptonSummedSpec2Clusters_1PixClusters[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
			txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
			Y_text_top -= Y_text_step;
			txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
			Y_text_top -= Y_text_step;
			txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
			Y_text_top -= Y_text_step;
			txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpec2Clusters_1PixClusters[im]->GetBinWidth(1)));
			line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSummedSpec2Clusters_1PixClusters[im]->GetName()).Data()));
		
			comptonSummedSpec2Clusters_2PixClusters[im]->Draw();
			line1->DrawLine(Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_2PixClusters[im]->GetMinimum(),Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_2PixClusters[im]->GetMaximum()/3);
			line1->DrawLine(Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_2PixClusters[im]->GetMinimum(),Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_2PixClusters[im]->GetMaximum()/3);
			Y_text_top = 0.82;
			getHistoPeak(comptonSummedSpec2Clusters_2PixClusters[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
			Peak_FWHM_raw = getFWHM(comptonSummedSpec2Clusters_2PixClusters[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
			txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
			Y_text_top -= Y_text_step;
			txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
			Y_text_top -= Y_text_step;
			txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
			Y_text_top -= Y_text_step;
			txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpec2Clusters_2PixClusters[im]->GetBinWidth(1)));
			line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSummedSpec2Clusters_2PixClusters[im]->GetName()).Data()));
		
			comptonSummedSpec2Clusters_1_2PixClusters[im]->Draw();
			line1->DrawLine(Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_1_2PixClusters[im]->GetMinimum(),Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_1_2PixClusters[im]->GetMaximum()/3);
			line1->DrawLine(Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_1_2PixClusters[im]->GetMinimum(),Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_1_2PixClusters[im]->GetMaximum()/3);
			Y_text_top = 0.82;
			getHistoPeak(comptonSummedSpec2Clusters_1_2PixClusters[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
			Peak_FWHM_raw = getFWHM(comptonSummedSpec2Clusters_1_2PixClusters[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
			txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
			Y_text_top -= Y_text_step;
			txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
			Y_text_top -= Y_text_step;
			txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
			Y_text_top -= Y_text_step;
			txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpec2Clusters_1_2PixClusters[im]->GetBinWidth(1)));
			line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
			c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSummedSpec2Clusters_1_2PixClusters[im]->GetName()).Data()));
			
			energySpecCorr1PixClusters[im]->Draw("colz");
			c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(energySpecCorr1PixClusters[im]->GetName()).Data()));
			c1->GetPad(0)->SetLogz(1);
			c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(energySpecCorr1PixClusters[im]->GetName()).Data()));
			c1->GetPad(0)->SetLogz(0);		
			
			energySpecCorr2PixClusters[im]->Draw("colz");
			c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(energySpecCorr2PixClusters[im]->GetName()).Data()));
			c1->GetPad(0)->SetLogz(1);
			c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(energySpecCorr2PixClusters[im]->GetName()).Data()));
			c1->GetPad(0)->SetLogz(0);		
		}
		
		c1->cd(0);
		comptonSummedSpec2Clusters_Eunsummed[im]->Draw();
		line1->DrawLine(Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_Eunsummed[im]->GetMinimum(),Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_Eunsummed[im]->GetMaximum()/3);
		line1->DrawLine(Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_Eunsummed[im]->GetMinimum(),Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_Eunsummed[im]->GetMaximum()/3);
		Y_text_top = 0.82;
		getHistoPeak(comptonSummedSpec2Clusters_Eunsummed[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(comptonSummedSpec2Clusters_Eunsummed[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpec2Clusters_Eunsummed[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSummedSpec2Clusters_Eunsummed[im]->GetName()).Data()));
		
		comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im]->Draw();
		line1->DrawLine(Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im]->GetMinimum(),Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im]->GetMaximum()/3);
		line1->DrawLine(Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im]->GetMinimum(),Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im]->GetMaximum()/3);
		Y_text_top = 0.82;
		getHistoPeak(comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im]->GetName()).Data()));
		
		comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im]->Draw();
		line1->DrawLine(Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im]->GetMinimum(),Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im]->GetMaximum()/3);
		line1->DrawLine(Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im]->GetMinimum(),Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im]->GetMaximum()/3);
		Y_text_top = 0.82;
		getHistoPeak(comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im]->GetName()).Data()));
		
		comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im]->Draw();
		line1->DrawLine(Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im]->GetMinimum(),Ewindow4PhiAnalyis_min,comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im]->GetMaximum()/3);
		line1->DrawLine(Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im]->GetMinimum(),Ewindow4PhiAnalyis_max,comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im]->GetMaximum()/3);
		Y_text_top = 0.82;
		getHistoPeak(comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im], minPhotopeakE4PeakSearch, maxPhotopeakE4PeakSearch, locMaxPeak, locMaxPeakHeight);
		Peak_FWHM_raw = getFWHM(comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im], locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
		txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
		Y_text_top -= Y_text_step;
		txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im]->GetBinWidth(1)));
		line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im]->GetName()).Data()));
			
		phiAngle[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(phiAngle[im]->GetName()).Data()));
		
		phiAngleSelEvents[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(phiAngleSelEvents[im]->GetName()).Data()));
		
		thetaFromFirstClusterE_dPhi[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_dPhi[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_dPhi[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		thetaFromFirstClusterE_dPhi_w[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_dPhi_w[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(thetaFromFirstClusterE_dPhi_w[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		thetaFromXYtClusterE_dPhi[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(thetaFromXYtClusterE_dPhi[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(thetaFromXYtClusterE_dPhi[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
		
		thetaFromXYtClusterE_dPhi_w[im]->Draw("colz");
		c1->SaveAs(Form("%s/%s.png",outputPathAM[im].Data(),TString(thetaFromXYtClusterE_dPhi_w[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(1);
		c1->SaveAs(Form("%s/%s_logz.png",outputPathAM[im].Data(),TString(thetaFromXYtClusterE_dPhi_w[im]->GetName()).Data()));
		c1->GetPad(0)->SetLogz(0);
			
		thetaFromXYtE[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(thetaFromXYtE[im]->GetName()).Data()));
		
		thetaFromXYtE_w[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(thetaFromXYtE_w[im]->GetName()).Data()));
		
		thetaFromFirstClusterE[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(thetaFromFirstClusterE[im]->GetName()).Data()));
		
		thetaFromFirstClusterE_w[im]->Draw("hist");
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
		
		numberOfClusters[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(numberOfClusters[im]->GetName()).Data()));
		numberOfClusters_norm[im] = (TH1F *) numberOfClusters[im]->Clone();
		numberOfClusters_norm[im]->SetName(Form("numberOfClusters_norm_AM%d",im));
		numberOfClusters_norm[im]->SetNormFactor(1);
		numberOfClusters_norm[im]->GetYaxis()->SetTitle("Normalised counts");
		numberOfClusters_norm[im]->GetYaxis()->SetTitleOffset(1.3);
		numberOfClusters_norm[im]->Draw();
		c1->SaveAs(Form("%s/%s.gif",outputPathAM[im].Data(),TString(numberOfClusters_norm[im]->GetName()).Data()));
		
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
	
	dPhiAngle_w->SetMinimum(0);
	dPhiAngle_w->Draw("hist");
	if (makeCosineFit)
	{
		TF1 *ffunc = new TF1("fit",fitfun,dPhiAngle_w->GetXaxis()->GetXmin(),dPhiAngle_w->GetXaxis()->GetXmax(),2);
		dPhiAngle_w->Fit(ffunc,"RQ0");
		ffunc->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle_w->GetName()).Data()));
	
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
	
	dPhiAngleNorm_w->SetMinimum(0);
	b1 = dPhiAngleNorm_w->GetXaxis()->FindBin(-0.01);
	b2 = dPhiAngleNorm_w->GetXaxis()->FindBin(0.01);
	nf = (dPhiAngleNorm_w->GetBinContent(b1) + dPhiAngleNorm_w->GetBinContent(b2))/2;
	for (Int_t ib = 1; ib <= dPhiAngleNorm_w->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiAngleNorm_w->GetBinContent(ib)/nf;
		dPhiAngleNorm_w->SetBinContent(ib,bc);
	}
	dPhiAngleNorm_w->Draw("hist");
	if (makeCosineFit)
	{
		TF1 *ffunc = new TF1("fit",fitfun,dPhiAngleNorm_w->GetXaxis()->GetXmin(),dPhiAngleNorm_w->GetXaxis()->GetXmax(),2);
		dPhiAngleNorm_w->Fit(ffunc,"RQ0");
		ffunc->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngleNorm_w->GetName()).Data()));
	
	dPhiAngle1norm->SetMinimum(0);
	Float_t nf1 = (dPhiAngle1norm->GetBinContent(1) + dPhiAngle1norm->GetBinContent(dPhiAngle1norm->GetNbinsX()))/2;
	for (Int_t ib = 1; ib <= dPhiAngle1norm->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiAngle1norm->GetBinContent(ib)/nf1;
		dPhiAngle1norm->SetBinContent(ib,bc);
	}
	dPhiAngle1norm->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc1 = new TF1("fit",fitfun,dPhiAngle1norm->GetXaxis()->GetXmin(),dPhiAngle1norm->GetXaxis()->GetXmax(),2);
		dPhiAngle1norm->Fit(ffunc1,"RQ0");
		ffunc1->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle1norm->GetName()).Data()));

	dPhiAngle_Ewin->SetMinimum(0);
	dPhiAngle_Ewin->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc_Ewin = new TF1("fit_Ewin",fitfun,dPhiAngle_Ewin->GetXaxis()->GetXmin(),dPhiAngle_Ewin->GetXaxis()->GetXmax(),2);
		dPhiAngle_Ewin->Fit(ffunc_Ewin,"RQ0");
		ffunc_Ewin->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle_Ewin->GetName()).Data()));

	dPhiAngle_Ewin_w->SetMinimum(0);
	dPhiAngle_Ewin_w->Draw("hist");
	if (makeCosineFit)
	{
		TF1 *ffunc_Ewin = new TF1("fit_Ewin",fitfun,dPhiAngle_Ewin_w->GetXaxis()->GetXmin(),dPhiAngle_Ewin_w->GetXaxis()->GetXmax(),2);
		dPhiAngle_Ewin_w->Fit(ffunc_Ewin,"RQ0");
		ffunc_Ewin->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle_Ewin_w->GetName()).Data()));

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

	dPhiAngleNorm_Ewin_w->SetMinimum(0);
	b1w = dPhiAngleNorm_Ewin_w->GetXaxis()->FindBin(-0.01);
	b2w = dPhiAngleNorm_Ewin_w->GetXaxis()->FindBin(0.01);
	nf_w = (dPhiAngleNorm_Ewin_w->GetBinContent(b1w) + dPhiAngleNorm_Ewin_w->GetBinContent(b2w))/2;
	for (Int_t ib = 1; ib <= dPhiAngleNorm_Ewin_w->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiAngleNorm_Ewin_w->GetBinContent(ib)/nf_w;
		dPhiAngleNorm_Ewin_w->SetBinContent(ib,bc);
	}
	dPhiAngleNorm_Ewin_w->Draw("hist");
	if (makeCosineFit)
	{
		TF1 *ffunc_Ewin = new TF1("fit_Ewin",fitfun,dPhiAngleNorm_Ewin_w->GetXaxis()->GetXmin(),dPhiAngleNorm_Ewin_w->GetXaxis()->GetXmax(),2);
		dPhiAngleNorm_Ewin_w->Fit(ffunc_Ewin,"RQ0");
		ffunc_Ewin->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngleNorm_Ewin_w->GetName()).Data()));
	
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

	dPhiAngle1norm_Ewin->SetMinimum(0);
	Float_t nf1_w = (dPhiAngle1norm_Ewin->GetBinContent(1) + dPhiAngle1norm_Ewin->GetBinContent(dPhiAngle1norm_Ewin->GetNbinsX()))/2;
	for (Int_t ib = 1; ib <= dPhiAngle1norm_Ewin->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiAngle1norm_Ewin->GetBinContent(ib)/nf1_w;
		dPhiAngle1norm_Ewin->SetBinContent(ib,bc);
	}
	dPhiAngle1norm_Ewin->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc2_Ewin = new TF1("fit_Ewin",fitfun,dPhiAngle1norm_Ewin->GetXaxis()->GetXmin(),dPhiAngle1norm_Ewin->GetXaxis()->GetXmax(),2);
		dPhiAngle1norm_Ewin->Fit(ffunc2_Ewin,"RQ0");
		ffunc2_Ewin->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiAngle1norm_Ewin->GetName()).Data()));
	
	dPhiXYtAngle->SetMinimum(0);
	dPhiXYtAngle->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc = new TF1("fit",fitfun,dPhiXYtAngle->GetXaxis()->GetXmin(),dPhiXYtAngle->GetXaxis()->GetXmax(),2);
		dPhiXYtAngle->Fit(ffunc,"RQ0");
		ffunc->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiXYtAngle->GetName()).Data()));	
	
	dPhiXYtAngle_Ewin->SetMinimum(0);
	dPhiXYtAngle_Ewin->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc_Ewin = new TF1("fit_Ewin",fitfun,dPhiXYtAngle_Ewin->GetXaxis()->GetXmin(),dPhiXYtAngle_Ewin->GetXaxis()->GetXmax(),2);
		dPhiXYtAngle_Ewin->Fit(ffunc_Ewin,"RQ0");
		ffunc_Ewin->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiXYtAngle_Ewin->GetName()).Data()));
	
	dPhiXYtAngleNorm->SetMinimum(0);
	Int_t b1x = dPhiXYtAngleNorm->GetXaxis()->FindBin(-0.01);
	Int_t b2x = dPhiXYtAngleNorm->GetXaxis()->FindBin(0.01);
	Float_t nfx = (dPhiXYtAngleNorm->GetBinContent(b1x) + dPhiXYtAngleNorm->GetBinContent(b2x))/2;
	for (Int_t ib = 1; ib <= dPhiXYtAngleNorm->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiXYtAngleNorm->GetBinContent(ib)/nfx;
		dPhiXYtAngleNorm->SetBinContent(ib,bc);
	}
	dPhiXYtAngleNorm->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc = new TF1("fit",fitfun,dPhiXYtAngleNorm->GetXaxis()->GetXmin(),dPhiXYtAngleNorm->GetXaxis()->GetXmax(),2);
		dPhiXYtAngleNorm->Fit(ffunc,"RQ0");
		ffunc->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiXYtAngleNorm->GetName()).Data()));
	
	dPhiXYtAngleNorm_Ewin->SetMinimum(0);
	Int_t b1wx = dPhiXYtAngleNorm_Ewin->GetXaxis()->FindBin(-0.01);
	Int_t b2wx = dPhiXYtAngleNorm_Ewin->GetXaxis()->FindBin(0.01);
	Float_t nf_wx = (dPhiXYtAngleNorm_Ewin->GetBinContent(b1wx) + dPhiXYtAngleNorm_Ewin->GetBinContent(b2wx))/2;
	for (Int_t ib = 1; ib <= dPhiXYtAngleNorm_Ewin->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiXYtAngleNorm_Ewin->GetBinContent(ib)/nf_wx;
		dPhiXYtAngleNorm_Ewin->SetBinContent(ib,bc);
	}
	dPhiXYtAngleNorm_Ewin->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc_Ewin = new TF1("fit_Ewin",fitfun,dPhiXYtAngleNorm_Ewin->GetXaxis()->GetXmin(),dPhiXYtAngleNorm_Ewin->GetXaxis()->GetXmax(),2);
		dPhiXYtAngleNorm_Ewin->Fit(ffunc_Ewin,"RQ0");
		ffunc_Ewin->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiXYtAngleNorm_Ewin->GetName()).Data()));
	
	dPhiXYtAngle1->SetMinimum(0);
	dPhiXYtAngle1->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc = new TF1("fit",fitfun,dPhiXYtAngle1->GetXaxis()->GetXmin(),dPhiXYtAngle1->GetXaxis()->GetXmax(),2);
		dPhiXYtAngle1->Fit(ffunc,"RQ0");
		ffunc->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiXYtAngle1->GetName()).Data()));
	
	dPhiXYtAngle1_Ewin->SetMinimum(0);
	dPhiXYtAngle1_Ewin->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc = new TF1("fit",fitfun,dPhiXYtAngle1_Ewin->GetXaxis()->GetXmin(),dPhiXYtAngle1_Ewin->GetXaxis()->GetXmax(),2);
		dPhiXYtAngle1_Ewin->Fit(ffunc,"RQ0");
		ffunc->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiXYtAngle1_Ewin->GetName()).Data()));
	
	dPhiXYtAngle1Norm->SetMinimum(0);
	Float_t nf1x = (dPhiXYtAngle1Norm->GetBinContent(1) + dPhiXYtAngle1Norm->GetBinContent(dPhiXYtAngle1Norm->GetNbinsX()))/2;
	for (Int_t ib = 1; ib <= dPhiXYtAngle1Norm->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiXYtAngle1Norm->GetBinContent(ib)/nf1x;
		dPhiXYtAngle1Norm->SetBinContent(ib,bc);
	}
	dPhiXYtAngle1Norm->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc = new TF1("fit",fitfun,dPhiXYtAngle1Norm->GetXaxis()->GetXmin(),dPhiXYtAngle1Norm->GetXaxis()->GetXmax(),2);
		dPhiXYtAngle1Norm->Fit(ffunc,"RQ0");
		ffunc->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiXYtAngle1Norm->GetName()).Data()));
	
	dPhiXYtAngle1Norm_Ewin->SetMinimum(0);
	nf1x = (dPhiXYtAngle1Norm_Ewin->GetBinContent(1) + dPhiXYtAngle1Norm_Ewin->GetBinContent(dPhiXYtAngle1Norm_Ewin->GetNbinsX()))/2;
	for (Int_t ib = 1; ib <= dPhiXYtAngle1Norm_Ewin->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiXYtAngle1Norm_Ewin->GetBinContent(ib)/nf1x;
		dPhiXYtAngle1Norm_Ewin->SetBinContent(ib,bc);
	}
	dPhiXYtAngle1Norm_Ewin->Draw();
	if (makeCosineFit)
	{
		TF1 *ffunc = new TF1("fit",fitfun,dPhiXYtAngle1Norm_Ewin->GetXaxis()->GetXmin(),dPhiXYtAngle1Norm_Ewin->GetXaxis()->GetXmax(),2);
		dPhiXYtAngle1Norm_Ewin->Fit(ffunc,"RQ0");
		ffunc->Draw("same");
	}
	c1->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(dPhiXYtAngle1Norm_Ewin->GetName()).Data()));
	
	comptonSummedSpec2Clusters2Heads->Draw();
	line1->DrawLine(Ewindow4PhiAnalyis_min*2,comptonSummedSpec2Clusters2Heads->GetMinimum(),Ewindow4PhiAnalyis_min*2,comptonSummedSpec2Clusters2Heads->GetMaximum()/3);
	line1->DrawLine(Ewindow4PhiAnalyis_max*2,comptonSummedSpec2Clusters2Heads->GetMinimum(),Ewindow4PhiAnalyis_max*2,comptonSummedSpec2Clusters2Heads->GetMaximum()/3);
	Y_text_top = 0.82;
	getHistoPeak(comptonSummedSpec2Clusters2Heads, minPhotopeakE4PeakSearch*2, maxPhotopeakE4PeakSearch*2, locMaxPeak, locMaxPeakHeight);
	Peak_FWHM_raw = getFWHM(comptonSummedSpec2Clusters2Heads, locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch*2, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
	txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
	Y_text_top -= Y_text_step;
	txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
	Y_text_top -= Y_text_step;
	txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
	Y_text_top -= Y_text_step;
	txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpec2Clusters2Heads->GetBinWidth(1)));
	line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
	c1->SaveAs(Form("%s/%s.gif",outputPathPairs.Data(),TString(comptonSummedSpec2Clusters2Heads->GetName()).Data()));
	
	comptonSummedSpec2ClustersSelEvents2Heads->Draw();
	line1->DrawLine(Ewindow4PhiAnalyis_min*2,comptonSummedSpec2ClustersSelEvents2Heads->GetMinimum(),Ewindow4PhiAnalyis_min*2,comptonSummedSpec2ClustersSelEvents2Heads->GetMaximum()/3);
	line1->DrawLine(Ewindow4PhiAnalyis_max*2,comptonSummedSpec2ClustersSelEvents2Heads->GetMinimum(),Ewindow4PhiAnalyis_max*2,comptonSummedSpec2ClustersSelEvents2Heads->GetMaximum()/3);
	Y_text_top = 0.82;
	getHistoPeak(comptonSummedSpec2ClustersSelEvents2Heads, minPhotopeakE4PeakSearch*2, maxPhotopeakE4PeakSearch*2, locMaxPeak, locMaxPeakHeight);
	Peak_FWHM_raw = getFWHM(comptonSummedSpec2ClustersSelEvents2Heads, locMaxPeak, locMaxPeakHeight, 200, 200, minPhotopeakE4PeakSearch*2, 4000, Peak_FWHM_raw_leftEdge, Peak_FWHM_raw_rightEdge);
	txt->DrawLatex(X_text,Y_text_top,Form("peak = %.1f keV",locMaxPeak));
	Y_text_top -= Y_text_step;
	txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f keV",Peak_FWHM_raw));
	Y_text_top -= Y_text_step;
	txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.2f %%",Peak_FWHM_raw/locMaxPeak*100));
	Y_text_top -= Y_text_step;
	txt->DrawLatex(X_text,Y_text_top,Form("FWHM = %.1f bins",Peak_FWHM_raw/comptonSummedSpec2ClustersSelEvents2Heads->GetBinWidth(1)));
	line2->DrawLine(Peak_FWHM_raw_leftEdge,locMaxPeakHeight/2,Peak_FWHM_raw_rightEdge,locMaxPeakHeight/2);
	c1->SaveAs(Form("%s/%s.gif",outputPathPairs.Data(),TString(comptonSummedSpec2ClustersSelEvents2Heads->GetName()).Data()));

	c2->cd(0);
	comptonSummedSpec2ClustersCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(comptonSummedSpec2ClustersCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(comptonSummedSpec2ClustersCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);

	comptonSummedSpec2ClustersSelEventsCorr->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(comptonSummedSpec2ClustersSelEventsCorr->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(comptonSummedSpec2ClustersSelEventsCorr->GetName()).Data()));
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

	theta0_theta1_FromFirstClusterE->Draw("colz");
	c2->SaveAs(Form("%s/%s.png",outputPathPairs.Data(),TString(theta0_theta1_FromFirstClusterE->GetName()).Data()));
	c2->GetPad(0)->SetLogz(1);
	c2->SaveAs(Form("%s/%s_logz.png",outputPathPairs.Data(),TString(theta0_theta1_FromFirstClusterE->GetName()).Data()));
	c2->GetPad(0)->SetLogz(0);
	
	delete c0;
	delete c1;
	delete c2;
	//logfile.close();
	
	if (totalTimeElapced < 60) cout << Form("Total running time is %.2f sec",totalTimeElapced) << endl;
	if (totalTimeElapced >= 60 && totalTimeElapced < 3600) cout << Form("Total running time is %dm:%.0f sec",Int_t(totalTimeElapced/60),totalTimeElapced - Int_t(totalTimeElapced/60)*60) << endl;
	if (totalTimeElapced >= 3600) cout << Form("Total running time is %dh:%dm:%.0fs",Int_t(totalTimeElapced/3600),Int_t((totalTimeElapced - Int_t(totalTimeElapced/3600)*3600)/60),
		totalTimeElapced - Int_t(totalTimeElapced/3600)*3600 - Int_t((totalTimeElapced - Int_t(totalTimeElapced/3600)*3600)/60)*60) << endl;
	
	hfile1->cd();
	for (Int_t im = 0; im < nAMs; im++)
	{
		cathodeSpecAllEvents[im]->Write();
		cathodeSpecPair[im]->Write();
		cathodeSpecSelPairs[im]->Write();
		allClustersSpec[im]->Write();
		pixelSpecNeighbour1PixClusters[im]->Write();
		nTrigsInEvent[im]->Write();
		comptonSummedSpec2Clusters_Eunsummed[im]->Write();
		comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im]->Write();
		comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im]->Write();
		comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im]->Write();
		if (enableUpdateEnergyClusters)
		{
			comptonSpecClusters[im]->Write();
			clusterSummingStats1PixClusters[im]->Write();
			clusterSummingStats2PixClusters[im]->Write();
			comptonSummedSpec2Clusters[im]->Write();
			comptonSummedSpec2Clusters_1PixClusters[im]->Write();
			comptonSummedSpec2Clusters_2PixClusters[im]->Write();
			comptonSummedSpec2Clusters_1_2PixClusters[im]->Write();
			energySpecCorr1PixClusters[im]->Write();
			energySpecCorr2PixClusters[im]->Write();
		}
		comptonSpecClusters_Eunsummed[im]->Write();
		comptonSpec2ClustersCorr[im]->Write();
		clusterSpec2ClusterEventsCorr[im]->Write();
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
		thetaFromXYtE[im]->Write();
		thetaFromXYtE_w[im]->Write();
		thetaFromFirstClusterE[im]->Write();
		thetaFromFirstClusterE_dPhi[im]->Write();
		thetaFromFirstClusterE_dPhi_w[im]->Write();
		thetaFromXYtClusterE_dPhi[im]->Write();
		thetaFromXYtClusterE_dPhi_w[im]->Write();
		nTrigsInEvent_norm[im]->Write();
		clusterSize[im]->Write();
		clusterSize2ClusterEventsCorr[im]->Write();
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
		clusterSpec_1ClusterEvents[im]->Write();
		firstComptonClusterE[im]->Write();
		secondComptonClusterE[im]->Write();
		clusterSize_norm[im]->Write();
		numberOfClusters[im]->Write();
		numberOfClusters_norm[im]->Write();
		clusterSize_nClusters[im]->Write();
		theta_firstCluster_XYtE_Corr[im]->Write();
		theta_firstCluster_XYtE_Corr_w[im]->Write();
		comptonEventsDeltaZpix[im]->Write();
		comptonEventsDeltaZdist[im]->Write();
		if (makeTimingCalibrationStuff)
		{
			if (!doNotUseCornerPixelsInPixelClusterCOGEneg)
			{
				anodeTimingDiff2PixDiagClusters[im]->Write();
				anodeTiming2PixDiagClustersCorr[im]->Write();
			}
			clusterZfromAnodeTiming[im]->Write();
			anodeTimingDiff2PixClusters[im]->Write();
			anodeTiming2PixClustersCorr[im]->Write();
			anodeTimingDiffAllClusters[im]->Write();
			dAnodeTiming_vs_theta[im]->Write();
			dAnodeTiming_vs_E1[im]->Write();
			dAnodeTiming_vs_E2[im]->Write();
			dAnodeTiming_vs_E12[im]->Write();
			dAnodeTimingAbs_vs_E12[im]->Write();
			if (enableUpdateEnergyClusters)
			{
				dAnodeTiming_vs_sumE_1PixelClustersNeigb[im]->Write();
				dAnodeTiming_vs_nE_1PixelClustersNeigb[im]->Write();
				dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters[im]->Write();
				anodeTiming1PixelClustersNeigbCorr[im]->Write();
				anodeTiming2PixelClustersNeigbCorr[im]->Write();
				anodeTimingNeigb_nE_1PixelClusters[im]->Write();
				anodeTimingNeigb_nE_2PixelClusters[im]->Write();
			}
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
	dPhiXYtAngle->Write();
	dPhiXYtAngle_Ewin->Write();
	dPhiXYtAngleNorm->Write();
	dPhiXYtAngleNorm_Ewin->Write();
	dPhiXYtAngle1->Write();
	dPhiXYtAngle1_Ewin->Write();
	dPhiXYtAngle1Norm->Write();
	dPhiXYtAngle1Norm_Ewin->Write();
	dPhiAngle->Write();
	dPhiAngleNorm->Write();
	dPhiAngle_Ewin->Write();
	dPhiAngleNorm_Ewin->Write();
	dPhiAngle1->Write();
	dPhiAngle1_Ewin->Write();
	dPhiAngle1norm->Write();
	dPhiAngle1norm_Ewin->Write();
	dPhiAngle_2ClusterE->Write();
	phiAngleCorr->Write();
	phiAngleSelEventsCorr->Write();
	comptonSummedSpec2ClustersCorr->Write();
	comptonSummedSpec2Clusters2Heads->Write();
	comptonSummedSpec2ClustersSelEventsCorr->Write();
	comptonSummedSpec2ClustersSelEvents2Heads->Write();
	theta0_theta1_FromFirstClusterE->Write();

	hfile1->Close();
	hfile2->cd();
	tree->Write();
	hfile2->Close();
}

Bool_t analyseNextEvent(const Int_t ievent)
{
	buffClusterX.clear();
	buffClusterY.clear();
	buffClusterE.clear();
	buffClusterEorig.clear();
	buffClusterArea.clear();
	buffClusterFlag.clear();
	buffClusterTrigs.clear();
	buffClusterIsSplit.clear();
	buffClusterAnodeTiming.clear();
	buffClusterOrientation.clear();
	eventDisplayFlag = kFALSE;
	if (enableEventMonitor)
	{
		if (nEventsDisplayed < nEvents2Display) eventDisplayFlag = kTRUE;
		if (nEvents2Display == -1) eventDisplayFlag = kTRUE;
		if (ievent == eventNumber2Display) eventDisplayFlag = kTRUE;
	}

	for (Int_t im = 0; im < nAMs; im++)
	{
		cathodeE[im] = -9999;
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
	
	Bool_t skipEvent = kFALSE;
	for (Int_t p = 0; p < buffClusterX[0].size(); p++)
	{
		if (usePixelListToDisable && disabledPixels[Int_t(pixelPattern[0]->GetBinContent(buffClusterX[0][p]+1,buffClusterY[0][p]+1)+0.5)]) skipEvent = kTRUE;
	}
	for (Int_t p = 0; p < buffClusterX[1].size(); p++)
	{
		if (usePixelListToDisable && disabledPixels[Int_t(pixelPattern[1]->GetBinContent(buffClusterX[1][p]+1,buffClusterY[1][p]+1)+0.5)]) skipEvent = kTRUE;
	}
	if (buffClusterX[0].size() >= 2 && buffClusterX[1].size() >= 2 && !skipEvent)
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
			phiAngle[0]->Fill(phiAng0);
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
			phiAngle[1]->Fill(phiAng1mirrored_shifted);
			firstComptonCluster[1]->Fill(xc1+1,yc1+1);
			secondComptonCluster[1]->Fill(xc2+1,yc2+1);
			dPhiPlaneAM[1]->Fill(xcentre,ycentre-11);
			dPhiPlaneAM[1]->Fill(xc2+xcentre-xc1,yc2+ycentre-yc1-11);
			phi1_valid = kTRUE;
		}
		
		for (Int_t im = 0; im < nAMs; im++)
		{
			if (buffClusterX[im].size() != 2) continue;
			Float_t ddist = TMath::Hypot(buffClusterX[im][firstClusterIdx[im]]-buffClusterX[im][secondClusterIdx[im]],buffClusterY[im][firstClusterIdx[im]]-buffClusterY[im][secondClusterIdx[im]]);
			comptonEventsDeltaZpix[im]->Fill(ddist);
			comptonEventsDeltaZdist[im]->Fill(ddist*pixelPitch);
		}
		
		if (phi0_valid && phi1_valid)
		{
			Float_t www1 = prob_h1->GetBinContent(prob_h1->GetXaxis()->FindBin(buffClusterE[0][firstClusterIdx[0]]));
			Float_t www2 = prob2_h1->GetBinContent(prob2_h1->GetXaxis()->FindBin(buffClusterE[0][secondClusterIdx[0]]));
			Float_t www1a = prob_h1->GetBinContent(prob_h1->GetXaxis()->FindBin(buffClusterE[1][firstClusterIdx[1]]));
			Float_t www2a = prob2_h1->GetBinContent(prob2_h1->GetXaxis()->FindBin(buffClusterE[1][secondClusterIdx[1]]));			
			if (cathodeE[0] > 1 && cathodeE[1] > 1)
			{
				cathodeSpecSelPairs[0]->Fill(cathodeE[0]);
				cathodeSpecSelPairs[1]->Fill(cathodeE[1]);
			}
			Float_t theta0 = getThetaFromEnergy(Na22Energy, buffClusterE[0][firstClusterIdx[0]]);
			Float_t theta1 = getThetaFromEnergy(initialE1, buffClusterE[1][firstClusterIdx[1]]);
			Float_t theta0a = getThetaFromEnergy(Na22Energy, buffClusterE[0][secondClusterIdx[0]]);
			Float_t theta1a = getThetaFromEnergy(initialE1, buffClusterE[1][secondClusterIdx[1]]);
			
			Float_t theta0a_xyt = getThetaAngleFromHitsXYt(0,buffClusterX[0][firstClusterIdx[0]],buffClusterX[0][secondClusterIdx[0]],buffClusterY[0][firstClusterIdx[0]],buffClusterY[0][secondClusterIdx[0]],buffClusterAnodeTiming[0][firstClusterIdx[0]],buffClusterAnodeTiming[0][secondClusterIdx[0]]);
			Float_t theta1a_xyt = getThetaAngleFromHitsXYt(1,buffClusterX[1][firstClusterIdx[1]],buffClusterX[1][secondClusterIdx[1]],buffClusterY[1][firstClusterIdx[1]],buffClusterY[1][secondClusterIdx[1]],buffClusterAnodeTiming[1][firstClusterIdx[1]],buffClusterAnodeTiming[1][secondClusterIdx[1]]);
			Float_t theta0_xyt = 180 - theta0a_xyt;
			Float_t theta1_xyt = 180 - theta1a_xyt;
			
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
			}
			Bool_t skipEvent_cathode = kFALSE;
			if (eventHasToHaveCathodeSignal)
			{
				if (cathodeE[0] < -9000) skipEvent_cathode = kTRUE;
				if (cathodeE[1] < -9000) skipEvent_cathode = kTRUE;
			}
			if (cathodeE[0] > maxCathodeEnergy || cathodeE[1] > maxCathodeEnergy) skipEvent_cathode = kTRUE;
		
			// true is good
			Bool_t thetaWindowCut = kTRUE;
			if (makeAsymmetricThetaWindow)
			{
				if (!((theta0 > Theta1WindowFor4PhiAnalyis_min && theta0 < Theta1WindowFor4PhiAnalyis_max && theta1 > Theta2WindowFor4PhiAnalyis_min && theta1 < Theta2WindowFor4PhiAnalyis_max)
						|| (theta0 > Theta2WindowFor4PhiAnalyis_min && theta0 < Theta2WindowFor4PhiAnalyis_max && theta1 > Theta1WindowFor4PhiAnalyis_min && theta1 < Theta1WindowFor4PhiAnalyis_max)))
									thetaWindowCut = kFALSE;
			}
			else
			{
				if (theta0 <= Theta1WindowFor4PhiAnalyis_min || theta0 >= Theta1WindowFor4PhiAnalyis_max) thetaWindowCut = kFALSE;
				if (theta1 <= Theta2WindowFor4PhiAnalyis_min || theta1 >= Theta2WindowFor4PhiAnalyis_max) thetaWindowCut = kFALSE;
			}
			
			Bool_t thetaXYtWindowCut = kTRUE;
			if (makeAsymmetricThetaWindow)
			{
				if (!((theta0_xyt > Theta1WindowFor4PhiAnalyis_min && theta0_xyt < Theta1WindowFor4PhiAnalyis_max && theta1_xyt > Theta2WindowFor4PhiAnalyis_min && theta1_xyt < Theta2WindowFor4PhiAnalyis_max)
						|| (theta0_xyt > Theta2WindowFor4PhiAnalyis_min && theta0_xyt < Theta2WindowFor4PhiAnalyis_max && theta1_xyt > Theta1WindowFor4PhiAnalyis_min && theta1_xyt < Theta1WindowFor4PhiAnalyis_max)))
									thetaXYtWindowCut = kFALSE;
			}
			else
			{
				if (theta0_xyt <= Theta1WindowFor4PhiAnalyis_min || theta0_xyt >= Theta1WindowFor4PhiAnalyis_max) thetaXYtWindowCut = kFALSE;
				if (theta1_xyt <= Theta2WindowFor4PhiAnalyis_min || theta1_xyt >= Theta2WindowFor4PhiAnalyis_max) thetaXYtWindowCut = kFALSE;
			}
		
			Bool_t EnergyWindowCut = kTRUE;
			if (sum2clE0 <= Ewindow4PhiAnalyis_min || sum2clE0 >= Ewindow4PhiAnalyis_max) EnergyWindowCut = kFALSE;
			if (sum2clE1 <= Ewindow4PhiAnalyis_min || sum2clE1 >= Ewindow4PhiAnalyis_max) EnergyWindowCut = kFALSE;
		
			if (!skipEvent_cathode)
			{
				if (thetaWindowCut) // These are good events, that fill dPhi plot etc.
				{
					comptonSummedSpec2ClustersSelEventsCorr->Fill(sum2clE0,sum2clE1);
					comptonSummedSpec2ClustersSelEvents2Heads->Fill(sum2clE0+sum2clE1);
					phiAngleSelEventsCorr->Fill(phiAng0,phiAng1mirrored_shifted);
					phiAngleSelEvents[0]->Fill(phiAng0);
					phiAngleSelEvents[1]->Fill(phiAng1mirrored_shifted);
					nClustersInSelEvent[0]->Fill(buffClusterX[0].size());
					nClustersInSelEvent[1]->Fill(buffClusterX[1].size());
					firstComptonClusterSizeSelEvents[0]->Fill(buffClusterArea[0][firstClusterIdx[0]]);
					secondComptonClusterSizeSelEvents[0]->Fill(buffClusterArea[0][secondClusterIdx[0]]);
					firstComptonClusterSizeSelEvents[1]->Fill(buffClusterArea[1][firstClusterIdx[1]]);
					secondComptonClusterSizeSelEvents[1]->Fill(buffClusterArea[1][secondClusterIdx[1]]);
					dPhiAngle->Fill(dPh);
					dPhiAngleNorm->Fill(dPh);
					dPhiAngle1->Fill(TMath::Abs(dPh));
					dPhiAngle1norm->Fill(TMath::Abs(dPh));
					plotEvent = kTRUE;
					if (EnergyWindowCut)
					{
						dPhiAngle_Ewin->Fill(dPh);
						dPhiAngle1_Ewin->Fill(TMath::Abs(dPh));
						dPhiAngle1norm_Ewin->Fill(TMath::Abs(dPh));
						dPhiAngleNorm_Ewin->Fill(dPh);
					}
				}
				if (thetaWindowCut) 
				{
					dPhiXYtAngle->Fill(dPh);
					dPhiXYtAngleNorm->Fill(dPh);
					dPhiXYtAngle1->Fill(dPh);
					dPhiXYtAngle1Norm->Fill(dPh);
					if (EnergyWindowCut)
					{
						dPhiXYtAngle_Ewin->Fill(dPh);
						dPhiXYtAngleNorm_Ewin->Fill(dPh);
						dPhiXYtAngle1_Ewin->Fill(dPh);
						dPhiXYtAngle1Norm_Ewin->Fill(dPh);
					}
				}				
				if (theta0 > Theta1WindowFor4PhiAnalyis_min && theta0 < Theta1WindowFor4PhiAnalyis_max && theta1 > Theta1WindowFor4PhiAnalyis_min && theta1 < Theta1WindowFor4PhiAnalyis_max)	
				{
					dPhiAngle_w->Fill(dPh,www1*www1a);
					dPhiAngleNorm_w->Fill(dPh,www1*www1a);
					if (EnergyWindowCut)
					{
						dPhiAngle_Ewin_w->Fill(dPh,www1*www1a);
						dPhiAngleNorm_Ewin_w->Fill(dPh,www1*www1a);
					}
				}
				if (theta0 > Theta1WindowFor4PhiAnalyis_min && theta0 < Theta1WindowFor4PhiAnalyis_max && theta1a > Theta1WindowFor4PhiAnalyis_min && theta1a < Theta1WindowFor4PhiAnalyis_max)	
				{
					dPhiAngle_w->Fill(dPh,www1*www2a);
					dPhiAngleNorm_w->Fill(dPh,www1*www2a);
					if (EnergyWindowCut)
					{
						dPhiAngle_Ewin_w->Fill(dPh,www1*www2a);
						dPhiAngleNorm_Ewin_w->Fill(dPh,www1*www2a);
					}
				}
				if (theta0a > Theta1WindowFor4PhiAnalyis_min && theta0a < Theta1WindowFor4PhiAnalyis_max && theta1 > Theta1WindowFor4PhiAnalyis_min && theta1 < Theta1WindowFor4PhiAnalyis_max)	
				{
					dPhiAngle_w->Fill(dPh,www2*www1a);
					dPhiAngleNorm_w->Fill(dPh,www2*www1a);
					if (EnergyWindowCut)
					{
						dPhiAngle_Ewin_w->Fill(dPh,www2*www1a);
						dPhiAngleNorm_Ewin_w->Fill(dPh,www2*www1a);
					}
				}
				if (theta0a > Theta1WindowFor4PhiAnalyis_min && theta0a < Theta1WindowFor4PhiAnalyis_max && theta1a > Theta1WindowFor4PhiAnalyis_min && theta1a < Theta1WindowFor4PhiAnalyis_max)	
				{
					dPhiAngle_w->Fill(dPh,www1*www2a);
					dPhiAngleNorm_w->Fill(dPh,www1*www2a);
					if (EnergyWindowCut)
					{
						dPhiAngle_Ewin_w->Fill(dPh,www1*www2a);
						dPhiAngleNorm_Ewin_w->Fill(dPh,www1*www2a);
					}
				}
			
				dPhiAngle_2ClusterE->Fill(sum2clE0,dPh);
				dPhiAngle_2ClusterE->Fill(sum2clE1,dPh);
				if (EnergyWindowCut)
				{
					thetaFromFirstClusterE[0]->Fill(theta0);
					thetaFromFirstClusterE[1]->Fill(theta1);
					thetaFromFirstClusterE_w[0]->Fill(theta0,www1);
					thetaFromFirstClusterE_w[0]->Fill(theta0a,www2);
					thetaFromFirstClusterE_w[1]->Fill(theta1,www1a);
					thetaFromFirstClusterE_w[1]->Fill(theta1a,www2a);
					thetaFromXYtE[0]->Fill(theta0_xyt);
					thetaFromXYtE[1]->Fill(theta1_xyt);
					thetaFromXYtE_w[0]->Fill(theta0_xyt,www1);
					thetaFromXYtE_w[1]->Fill(theta1_xyt,www1a);
					thetaFromXYtE_w[0]->Fill(theta0a_xyt,www2);
					thetaFromXYtE_w[1]->Fill(theta1a_xyt,www2a);
					
					theta_firstCluster_XYtE_Corr[0]->Fill(theta0,theta0_xyt);
					theta_firstCluster_XYtE_Corr[1]->Fill(theta1,theta1_xyt);
					
					theta_firstCluster_XYtE_Corr_w[0]->Fill(theta0,theta0_xyt,www1);
					theta_firstCluster_XYtE_Corr_w[0]->Fill(theta0a,theta0a_xyt,www2);
					theta_firstCluster_XYtE_Corr_w[1]->Fill(theta1,theta1_xyt,www1a);
					theta_firstCluster_XYtE_Corr_w[1]->Fill(theta1a,theta1a_xyt,www2a);
					
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
					
					thetaFromXYtClusterE_dPhi[0]->Fill(theta0_xyt,dPh);
					thetaFromXYtClusterE_dPhi[1]->Fill(theta1_xyt,dPh);
					
					if (theta0 > 60) thetaFromFirstClusterE_dPhi_w[0]->Fill(theta0,dPh,www1);
					if (theta0a > 60) thetaFromFirstClusterE_dPhi_w[0]->Fill(theta0a,dPh,www2);
					if (theta1 > 60) thetaFromFirstClusterE_dPhi_w[1]->Fill(theta1,dPh,www1a);
					if (theta1a > 60) thetaFromFirstClusterE_dPhi_w[1]->Fill(theta1a,dPh,www2a);
					
					if (theta0_xyt > 60) thetaFromXYtClusterE_dPhi_w[0]->Fill(theta0_xyt,dPh,www1);
					if (theta0a_xyt > 60) thetaFromXYtClusterE_dPhi_w[0]->Fill(theta0a_xyt,dPh,www2);
					if (theta1_xyt > 60) thetaFromXYtClusterE_dPhi_w[1]->Fill(theta1_xyt,dPh,www1a);
					if (theta1a_xyt > 60) thetaFromXYtClusterE_dPhi_w[1]->Fill(theta1a_xyt,dPh,www2a);					
					
					theta0_theta1_FromFirstClusterE->Fill(theta0,theta1);
				} // end of: if (EnergyWindowCut)
				tree->Fill();
			} // end of: if (!skipEvent_cathode)
		} // end of: if (phi0_valid && phi1_valid)
		comptonSummedSpec2ClustersCorr->Fill(sum2clE0,sum2clE1);
		comptonSummedSpec2Clusters2Heads->Fill(sum2clE0+sum2clE1);
	}

	Float_t tEE1 = buffClusterE[0][firstClusterIdx[0]]+buffClusterE[0][secondClusterIdx[0]];
	Float_t tEE2 = buffClusterE[1][firstClusterIdx[1]]+buffClusterE[1][secondClusterIdx[1]];
	Bool_t EnergyWindowCut2 = kTRUE;
	if (tEE1 <= Ewindow4PhiAnalyis_min || tEE1 >= Ewindow4PhiAnalyis_max) EnergyWindowCut2 = kFALSE;
	if (tEE2 <= Ewindow4PhiAnalyis_min || tEE2 >= Ewindow4PhiAnalyis_max) EnergyWindowCut2 = kFALSE;

	if (buffClusterX[0].size() == 2 && buffClusterX[1].size() == 2 && EnergyWindowCut2)
	{
		firstComptonClusterE[0]->Fill(buffClusterE[0][firstClusterIdx[0]]);
		secondComptonClusterE[0]->Fill(buffClusterE[0][secondClusterIdx[0]]);
		firstComptonClusterE[1]->Fill(buffClusterE[1][firstClusterIdx[1]]);
		secondComptonClusterE[1]->Fill(buffClusterE[1][secondClusterIdx[1]]);
	}

	for (Int_t im = 0; im < nAMs; im++)
	{
		if (buffClusterX[im].size() == 2)
		{
			for (Int_t p = 0; p < 2; p++)
			{
				if (enableUpdateEnergyClusters) comptonSpecClusters[im]->Fill(buffClusterE[im][p]);
				comptonSpecClusters_Eunsummed[im]->Fill(buffClusterEorig[im][p]);
			}
			if (enableUpdateEnergyClusters) comptonSummedSpec2Clusters[im]->Fill(buffClusterE[im][0] + buffClusterE[im][1]);
			comptonSummedSpec2Clusters_Eunsummed[im]->Fill(buffClusterEorig[im][0] + buffClusterEorig[im][1]);
			if (buffClusterArea[im][0] == 1 && buffClusterArea[im][1] == 1)
			{
				if (enableUpdateEnergyClusters) comptonSummedSpec2Clusters_1PixClusters[im]->Fill(buffClusterE[im][0] + buffClusterE[im][1]);
				comptonSummedSpec2Clusters_1PixClusters_Eunsummed[im]->Fill(buffClusterEorig[im][0] + buffClusterEorig[im][1]);
			}
			if (buffClusterArea[im][0] == 2 && buffClusterArea[im][1] == 2)
			{
				if (enableUpdateEnergyClusters) comptonSummedSpec2Clusters_2PixClusters[im]->Fill(buffClusterE[im][0] + buffClusterE[im][1]);
				comptonSummedSpec2Clusters_2PixClusters_Eunsummed[im]->Fill(buffClusterEorig[im][0] + buffClusterEorig[im][1]);
			}
			if ((buffClusterArea[im][0] == 1 && buffClusterArea[im][1] == 2) || (buffClusterArea[im][0] == 2 && buffClusterArea[im][1] == 1))
			{
				if (enableUpdateEnergyClusters) comptonSummedSpec2Clusters_1_2PixClusters[im]->Fill(buffClusterE[im][0] + buffClusterE[im][1]);
				comptonSummedSpec2Clusters_1_2PixClusters_Eunsummed[im]->Fill(buffClusterEorig[im][0] + buffClusterEorig[im][1]);
			}
			comptonSpec2ClustersCorr[im]->Fill(buffClusterE[im][0],buffClusterE[im][1]);
		}
		for (Int_t p = 0; p < buffClusterX[im].size(); p++)
		{
			clusterE_vsSize[im]->Fill(buffClusterE[im][p],buffClusterArea[im][p]);
			allClustersSpec[im]->Fill(buffClusterE[im][p]);
			if (usePixelListToDisable && disabledPixels[Int_t(pixelPattern[im]->GetBinContent(buffClusterX[im][p]+1,buffClusterY[im][p]+1)+0.5)]) continue;
			allClustersCOGImage[im]->Fill(buffClusterX[im][p],buffClusterY[im][p],buffClusterE[im][p]);
			allClustersCOGFreqImage[im]->Fill(buffClusterX[im][p],buffClusterY[im][p]);
			allClustersFineCOGImage[im]->Fill(buffClusterX[im][p],buffClusterY[im][p],buffClusterE[im][p]);
			allClustersFineCOGFreqImage[im]->Fill(buffClusterX[im][p],buffClusterY[im][p]);
			if (plotPixelClusterSpectra) clusterSpec[im][Int_t(pixelPattern[im]->GetBinContent(buffClusterX[im][p]+1,buffClusterY[im][p]+1)+0.5)]->Fill(buffClusterE[im][p]);
		}
		if (buffClusterX[im].size() > 0) nClustersInEvent[im]->Fill(buffClusterX[im].size());
	}

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
	buffClusterEorigloc.clear();
	buffClusterArealoc.clear();
	buffClusterFlagloc.clear();
	buffClusterTrigsloc.clear();
	buffClusterIsSplitloc.clear();
	buffClusterOrientationloc.clear();
	buffClusterAnodeTimingloc.clear();
	buffClusterXloc.resize(0,0);
	buffClusterYloc.resize(0,0);
	buffClusterEloc.resize(0,0);
	buffClusterEorigloc.resize(0,0);
	buffClusterArealoc.resize(0,0);
	buffClusterFlagloc.resize(0,0);
	buffClusterTrigsloc.resize(0,0);
	buffClusterIsSplitloc.resize(0,0);
	buffClusterOrientationloc.resize(0,0);
	buffClusterAnodeTimingloc.resize(0,0);
	matrix_trigs.clear();
	matrix_trigs.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_flags.clear();
	matrix_flags.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_E.clear();
	matrix_E.resize(nPixXY, vector<Float_t>(nPixXY));
	matrix_Eall.clear();
	matrix_Eall.resize(nPixXY, vector<Float_t>(nPixXY));
	matrix_Eneg.clear();
	matrix_Eneg.resize(nPixXY, vector<Float_t>(nPixXY));
	matrix_time.clear();
	matrix_time.resize(nPixXY, vector<Int_t>(nPixXY));
	matrix_timeCalib.clear();
	matrix_timeCalib.resize(nPixXY, vector<Float_t>(nPixXY));
	matrix_trigs_4sum.clear();
	matrix_trigs_4sum.resize(nPixXY, vector<Int_t>(nPixXY));
	buff_E_neighb_4sum.clear();
	buff_E_neighb_4sum.resize(0,0);
	buff_X_neighb_4sum.clear();
	buff_X_neighb_4sum.resize(0,0);
	buff_Y_neighb_4sum.clear();
	buff_Y_neighb_4sum.resize(0,0);

	firstAMTreeEventIdx[AM2do] = 0;
	lastAMTreeEventIdx[AM2do] = 0;

	//cout << Form("========Start==AM%d===event %d=====",AM2do,evt) << endl;
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
					//matrix_timeCalib[xx1-1][yy1-1] = (timeDetect[AM2do]-tCalibAnodes[pixel[AM2do]][1]+rand3->Uniform(1.))/tCalibAnodes[pixel[AM2do]][0];				
				}
				matrix_timeCalib[xx1-1][yy1-1] = (timeDetect[AM2do]-tCalibAnodes[pixel[AM2do]][1]+rand3->Uniform(1.))/tCalibAnodes[pixel[AM2do]][0];
				if (E[AM2do] > minPixelEnergyThr4ClusterReconstruction)
				{
					if ((useOnlyTriggeredPixelsInCluster && triggerFlag[AM2do] == 1) || !useOnlyTriggeredPixelsInCluster)
					{
						status = kTRUE;
						matrix_E[xx1-1][yy1-1] = E[AM2do];
						matrix_trigs[xx1-1][yy1-1] = triggerFlag[AM2do];
						//matrix_trigs_4sum[xx1-1][yy1-1] = triggerFlag[AM2do];
						matrix_flags[xx1-1][yy1-1] = newNeighbourFlag;
						if (eventDisplayFlag)
						{
							imageE[AM2do]->SetBinContent(xx1,yy1,E[AM2do]);
							imageN[AM2do]->SetBinContent(xx1,yy1,0);
							if (triggerFlag[AM2do] == 1) imageT[AM2do]->Fill(xx1-0.5,yy1-0.5);
						}
					}
					if (triggerFlag[AM2do] != 1)
					{
						buff_E_neighb_4sum.push_back(E[AM2do]);
						buff_X_neighb_4sum.push_back(xx1-1);
						buff_Y_neighb_4sum.push_back(yy1-1);
					}
				}
				matrix_Eall[xx1-1][yy1-1] = E[AM2do];
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
			
			for (Int_t ia = 0; ia < buffClusterXarr.size(); ia++)
			{
				matrix_trigs_4sum[buffClusterXarr[ia]][buffClusterYarr[ia]] = clusterCounter;
			}
			
			Int_t orient = 0; // 0 - horizontal (X) orientation, 1 - vertical (Y) orientation
			if (area == 2)
				if (buffClusterYarr[0] == buffClusterYarr[1]) orient = 1;
			
			if (makeTimingCalibrationStuff)
			{
				if (nTrigsInCluster <= 2) clusterZfromAnodeTiming[AM2do]->Fill(meanT/factor2ConvertAnodeTime2Distance[AM2do]);
				Int_t stat = 0;
				if (nTrigsInCluster == 2 && applyMaxDtCut4Clustering)
					if (TMath::Abs(buffClusterAnodeTimearr[0] - buffClusterAnodeTimearr[1]) > maxDtBetweenPixel4Clustering) stat = 1;
				buffClusterIsSplitloc.push_back(stat);
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
							if (makeTimingCalibrationStuff) anodeTimingDiffAllClusters[AM2do]->Fill(TMath::Abs(t1b[0]-t1b[1]));
							if (nTrigsInCluster == 2 && buffClusterTrigsarr[it] == 1 && buffClusterTrigsarr[it2] == 1) 
							{
								t1[0] = buffClusterTimearr[it];
								t1[1] = buffClusterTimearr[it2];
								t1a[0] = t1b[0];
								t1a[1] = t1b[1];
								if (makeTimingCalibrationStuff)
								{
									anodeTimingDiff2PixClusters[AM2do]->Fill(TMath::Abs(t1b[0]-t1b[1]));
									anodeTiming2PixClustersCorr[AM2do]->Fill(t1b[0],t1b[1]);
								}
								if (buffClusterXarr[it] != buffClusterXarr[it2] && buffClusterYarr[it] != buffClusterYarr[it2])
								{
									if (!doNotUseCornerPixelsInPixelClusterCOGEneg)
									{
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
				}
			}
			buffClusterXloc.push_back(xcen);
			buffClusterYloc.push_back(ycen);
			buffClusterEloc.push_back(totAmp);
			buffClusterEorigloc.push_back(totAmp);
			buffClusterArealoc.push_back(area);
			if (isPrimaryCluster) buffClusterFlagloc.push_back(1);
			else buffClusterFlagloc.push_back(0);
			buffClusterTrigsloc.push_back(nTrigsInCluster);
			buffClusterAnodeTimingloc.push_back(meanT);
			buffClusterOrientationloc.push_back(orient);
		}
	}
	if (clusterCounter > 0)
	{
		numberOfClusters[AM2do]->Fill(clusterCounter);
		Bool_t clusterSizesGood = kTRUE;
		Bool_t hasSplitClisters = kFALSE;
		Float_t maxE = -9999;
		Int_t maxEidx = -1;
		for (Int_t ic = 0; ic < clusterCounter; ic++)
		{
			if (buffClusterEloc[ic] > maxE)
			{
				maxE = buffClusterEloc[ic];
				maxEidx = ic;
			}
			clusterE_vsSize[AM2do]->Fill(buffClusterEloc[ic],buffClusterArealoc[ic]);
			clusterSize[AM2do]->Fill(buffClusterArealoc[ic]);
			clusterSize_nClusters[AM2do]->Fill(buffClusterArealoc[ic],clusterCounter);
			if (buffClusterArealoc[ic] < minClusterSize4PairAnalysis || buffClusterArealoc[ic] > maxClusterSize4PairAnalysis) clusterSizesGood = kFALSE;
			if (applyMaxDtCut4Clustering && buffClusterArealoc[ic] == 2)
				if (buffClusterIsSplitloc[ic] == 1) hasSplitClisters = kTRUE;
		}
		clusterMaxE_vsSize[AM2do]->Fill(maxE,buffClusterArealoc[maxEidx]);
		if (minClusterSize4PairAnalysis == -1 || maxClusterSize4PairAnalysis == -1) clusterSizesGood = kTRUE;
		Bool_t goodNClusters = kTRUE;
		if(clusterCounter < minNClusters4PairAnalysis || clusterCounter > maxNClusters4PairAnalysis) goodNClusters = kFALSE;
		if (minNClusters4PairAnalysis < 0 || maxNClusters4PairAnalysis < 0) goodNClusters = kTRUE;
		
		if (clusterCounter == 2)
		{
			clusterSize2ClusterEventsCorr[AM2do]->Fill(buffClusterArealoc[0],buffClusterArealoc[1]);
			clusterSpec2ClusterEventsCorr[AM2do]->Fill(buffClusterEloc[0],buffClusterEloc[1]);
		}
		if (clusterCounter == 1) clusterSpec_1ClusterEvents[AM2do]->Fill(buffClusterEloc[0]);
		if ((clusterCounter == 1 || (clusterCounter == 2 && !hasSplitClisters)) && clusterSizesGood)
		{
			updateClustersWithNeighboursEnergy(AM2do,evt);
		}
		if (clusterCounter == 1) clusterSpec_1ClusterEvents_summed[AM2do]->Fill(buffClusterEloc[0]);
		if (!clusterSizesGood || !goodNClusters || hasSplitClisters) return kFALSE;
		Int_t niter = 0;
		while (niter < buffClusterXloc.size())
		{
			Float_t maxE = -9999;
			Float_t maxEorig = -9999;
			Int_t idx = -1, maxArea, maxFlag, maxTrigs, maxIsSplit, maxO;
			Float_t maxX, maxY, maxT;
			for (Int_t ic = niter; ic < buffClusterXloc.size(); ic++)
			{
				if (eventDisplayFlag && niter == 0) imageC[AM2do]->Fill(buffClusterXloc[ic],buffClusterYloc[ic]);
				if (buffClusterEloc[ic] > maxE)
				{
					maxE = buffClusterEloc[ic];
					maxEorig = buffClusterEorigloc[ic];
					idx = ic;
					maxX = buffClusterXloc[ic];
					maxY = buffClusterYloc[ic];
					maxArea = buffClusterArealoc[ic];
					maxFlag = buffClusterFlagloc[ic];
					maxTrigs = buffClusterTrigsloc[ic];
					maxIsSplit = buffClusterIsSplitloc[ic];
					maxT = buffClusterAnodeTimingloc[ic];
					maxO = buffClusterOrientationloc[ic];
				}
			}
			buffClusterEloc.erase(buffClusterEloc.begin()+idx);
			buffClusterEorigloc.erase(buffClusterEorigloc.begin()+idx);
			buffClusterXloc.erase(buffClusterXloc.begin()+idx);
			buffClusterYloc.erase(buffClusterYloc.begin()+idx);
			buffClusterArealoc.erase(buffClusterArealoc.begin()+idx);
			buffClusterFlagloc.erase(buffClusterFlagloc.begin()+idx);
			buffClusterTrigsloc.erase(buffClusterTrigsloc.begin()+idx);
			buffClusterIsSplitloc.erase(buffClusterIsSplitloc.begin()+idx);
			buffClusterAnodeTimingloc.erase(buffClusterAnodeTimingloc.begin()+idx);
			buffClusterOrientationloc.erase(buffClusterOrientationloc.begin()+idx);
			buffClusterEloc.insert(buffClusterEloc.begin()+niter,maxE);
			buffClusterEorigloc.insert(buffClusterEorigloc.begin()+niter,maxEorig);
			buffClusterXloc.insert(buffClusterXloc.begin()+niter,maxX);
			buffClusterYloc.insert(buffClusterYloc.begin()+niter,maxY);
			buffClusterArealoc.insert(buffClusterArealoc.begin()+niter,maxArea);
			buffClusterFlagloc.insert(buffClusterFlagloc.begin()+niter,maxFlag);
			buffClusterTrigsloc.insert(buffClusterTrigsloc.begin()+niter,maxTrigs);
			buffClusterIsSplitloc.insert(buffClusterIsSplitloc.begin()+niter,maxIsSplit);
			buffClusterAnodeTimingloc.insert(buffClusterAnodeTimingloc.begin()+niter,maxT);
			buffClusterOrientationloc.insert(buffClusterOrientationloc.begin()+niter,maxO);
			niter++;
		}
		nTrigsInEvent[AM2do]->Fill(nTrigPixels[AM2do]);
		buffClusterX.push_back(buffClusterXloc);
		buffClusterY.push_back(buffClusterYloc);
		buffClusterE.push_back(buffClusterEloc);
		buffClusterEorig.push_back(buffClusterEorigloc);
		buffClusterArea.push_back(buffClusterArealoc);
		buffClusterFlag.push_back(buffClusterFlagloc);
		buffClusterTrigs.push_back(buffClusterTrigsloc);
		buffClusterIsSplit.push_back(buffClusterIsSplitloc);
		buffClusterAnodeTiming.push_back(buffClusterAnodeTimingloc);
		buffClusterOrientation.push_back(buffClusterOrientationloc);
	}
	else status = kFALSE;
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
		if (sline.Contains("pixelPitch ="))
		{
			sline.ReplaceAll("pixelPitch =","");
			pixelPitch = sline.Atof();
			if (printOutSetupFileParameters) cout << "pixelPitch = " << pixelPitch << endl;
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
		
		if (sline.Contains("XRF_Cd ="))
		{
			sline.ReplaceAll("XRF_Cd =","");
			XRF_Cd = sline.Atof();
			if (printOutSetupFileParameters) cout << "XRF_Cd = " << XRF_Cd << endl;
		}
		
		if (sline.Contains("XRF_Te ="))
		{
			sline.ReplaceAll("XRF_Te =","");
			XRF_Te = sline.Atof();
			if (printOutSetupFileParameters) cout << "XRF_Te = " << XRF_Te << endl;
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
		
		if (sline.Contains("Ewindow4PhiAnalyis_min ="))
		{
			sline.ReplaceAll("Ewindow4PhiAnalyis_min =","");
			Ewindow4PhiAnalyis_min = sline.Atof();
			if (printOutSetupFileParameters) cout << "Ewindow4PhiAnalyis_min = " << Ewindow4PhiAnalyis_min << endl;
		}
		
		if (sline.Contains("Ewindow4PhiAnalyis_max ="))
		{
			sline.ReplaceAll("Ewindow4PhiAnalyis_max =","");
			Ewindow4PhiAnalyis_max = sline.Atof();
			if (printOutSetupFileParameters) cout << "Ewindow4PhiAnalyis_max = " << Ewindow4PhiAnalyis_max << endl;
		}
		
		if (sline.Contains("Theta1WindowFor4PhiAnalyis_min ="))
		{
			sline.ReplaceAll("Theta1WindowFor4PhiAnalyis_min =","");
			Theta1WindowFor4PhiAnalyis_min = sline.Atof();
			if (printOutSetupFileParameters) cout << "Theta1WindowFor4PhiAnalyis_min = " << Theta1WindowFor4PhiAnalyis_min << endl;
		}
		
		if (sline.Contains("Theta1WindowFor4PhiAnalyis_max ="))
		{
			sline.ReplaceAll("Theta1WindowFor4PhiAnalyis_max =","");
			Theta1WindowFor4PhiAnalyis_max = sline.Atof();
			if (printOutSetupFileParameters) cout << "Theta1WindowFor4PhiAnalyis_max = " << Theta1WindowFor4PhiAnalyis_max << endl;
		}
		
		if (sline.Contains("Theta2WindowFor4PhiAnalyis_min ="))
		{
			sline.ReplaceAll("Theta2WindowFor4PhiAnalyis_min =","");
			Theta2WindowFor4PhiAnalyis_min = sline.Atof();
			if (printOutSetupFileParameters) cout << "Theta2WindowFor4PhiAnalyis_min = " << Theta2WindowFor4PhiAnalyis_min << endl;
		}
		
		if (sline.Contains("Theta2WindowFor4PhiAnalyis_max ="))
		{
			sline.ReplaceAll("Theta2WindowFor4PhiAnalyis_max =","");
			Theta2WindowFor4PhiAnalyis_max = sline.Atof();
			if (printOutSetupFileParameters) cout << "Theta2WindowFor4PhiAnalyis_max = " << Theta2WindowFor4PhiAnalyis_max << endl;
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
		
		if (sline.Contains("doNotUseCornerTouchingTriggeredPixels4Clustering ="))
		{
			sline.ReplaceAll("doNotUseCornerTouchingTriggeredPixels4Clustering =","");
			buf = sline.Atoi();
			doNotUseCornerTouchingTriggeredPixels4Clustering = kTRUE;
			if (buf == 0) doNotUseCornerTouchingTriggeredPixels4Clustering = kFALSE;
			if (printOutSetupFileParameters) cout << "doNotUseCornerTouchingTriggeredPixels4Clustering = " << doNotUseCornerTouchingTriggeredPixels4Clustering << endl;
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
			factor2ConvertAnodeTime2Distance = new Float_t[nAMs];
			sline.ReplaceAll("factor2ConvertAnodeTime2Distance =","");
			TObjArray *tt = sline.Tokenize(" ");
			for (Int_t ii = 0; ii < nAMs; ii++)
			{
				factor2ConvertAnodeTime2Distance[ii] = (((TObjString *) (tt->At(ii)))->String()).Atof();
				if (printOutSetupFileParameters) cout << "factor2ConvertAnodeTime2Distance[" << ii << "] = " << factor2ConvertAnodeTime2Distance[ii] << endl;
			}
		}
		
		if (sline.Contains("minZ ="))
		{
			sline.ReplaceAll("minZ =","");
			minZ = sline.Atof();
			if (printOutSetupFileParameters) cout << "minZ = " << minZ << endl;
		}
		
		if (sline.Contains("maxZ ="))
		{
			sline.ReplaceAll("maxZ =","");
			maxZ = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxZ = " << maxZ << endl;
		}
		
		if (sline.Contains("nBinsZ ="))
		{
			sline.ReplaceAll("nBinsZ =","");
			nBinsZ = sline.Atoi();
			if (printOutSetupFileParameters) cout << "nBinsZ = " << nBinsZ << endl;
		}
		
		if (sline.Contains("maxCathodeEnergy ="))
		{
			sline.ReplaceAll("maxCathodeEnergy =","");
			maxCathodeEnergy = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxCathodeEnergy = " << maxCathodeEnergy << endl;
		}
		
		if (sline.Contains("eventHasToHaveCathodeSignal ="))
		{
			sline.ReplaceAll("eventHasToHaveCathodeSignal =","");
			buf = sline.Atoi();
			eventHasToHaveCathodeSignal = kTRUE;
			if (buf == 0) eventHasToHaveCathodeSignal = kFALSE;
			if (printOutSetupFileParameters) cout << "eventHasToHaveCathodeSignal = " << eventHasToHaveCathodeSignal << endl;
		}
		
		if (sline.Contains("makeAsymmetricThetaWindow ="))
		{
			sline.ReplaceAll("makeAsymmetricThetaWindow =","");
			buf = sline.Atoi();
			makeAsymmetricThetaWindow = kTRUE;
			if (buf == 0) makeAsymmetricThetaWindow = kFALSE;
			if (printOutSetupFileParameters) cout << "makeAsymmetricThetaWindow = " << makeAsymmetricThetaWindow << endl;
		}
		
		if (sline.Contains("enableUpdateEnergyClusters ="))
		{
			sline.ReplaceAll("enableUpdateEnergyClusters =","");
			buf = sline.Atoi();
			enableUpdateEnergyClusters = kTRUE;
			if (buf == 0) enableUpdateEnergyClusters = kFALSE;
			if (printOutSetupFileParameters) cout << "enableUpdateEnergyClusters = " << enableUpdateEnergyClusters << endl;
		}
		
		if (sline.Contains("maxDtBetweenSinglePixClusters4updatingEnergy_low ="))
		{
			sline.ReplaceAll("maxDtBetweenSinglePixClusters4updatingEnergy_low =","");
			maxDtBetweenSinglePixClusters4updatingEnergy_low = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxDtBetweenSinglePixClusters4updatingEnergy_low = " << maxDtBetweenSinglePixClusters4updatingEnergy_low << endl;
		}
		
		if (sline.Contains("maxDtBetweenSinglePixClusters4updatingEnergy_up ="))
		{
			sline.ReplaceAll("maxDtBetweenSinglePixClusters4updatingEnergy_up =","");
			maxDtBetweenSinglePixClusters4updatingEnergy_up = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxDtBetweenSinglePixClusters4updatingEnergy_up = " << maxDtBetweenSinglePixClusters4updatingEnergy_up << endl;
		}
		
		if (sline.Contains("maxDtBetweenTwoPixClusters4updatingEnergy_low ="))
		{
			sline.ReplaceAll("maxDtBetweenTwoPixClusters4updatingEnergy_low =","");
			maxDtBetweenTwoPixClusters4updatingEnergy_low = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxDtBetweenTwoPixClusters4updatingEnergy_low = " << maxDtBetweenTwoPixClusters4updatingEnergy_low << endl;
		}
		
		if (sline.Contains("maxDtBetweenTwoPixClusters4updatingEnergy_up ="))
		{
			sline.ReplaceAll("maxDtBetweenTwoPixClusters4updatingEnergy_up =","");
			maxDtBetweenTwoPixClusters4updatingEnergy_up = sline.Atof();
			if (printOutSetupFileParameters) cout << "maxDtBetweenTwoPixClusters4updatingEnergy_up = " << maxDtBetweenTwoPixClusters4updatingEnergy_up << endl;
		}
		
		if (sline.Contains("updateSinglePixelClusterCOG ="))
		{
			sline.ReplaceAll("updateSinglePixelClusterCOG =","");
			buf = sline.Atoi();
			updateSinglePixelClusterCOG = kTRUE;
			if (buf == 0) updateSinglePixelClusterCOG = kFALSE;
			if (printOutSetupFileParameters) cout << "updateSinglePixelClusterCOG = " << updateSinglePixelClusterCOG << endl;
		}
		
		if (sline.Contains("updateTwoPixelClusterCOG ="))
		{
			sline.ReplaceAll("updateTwoPixelClusterCOG =","");
			buf = sline.Atoi();
			updateTwoPixelClusterCOG = kTRUE;
			if (buf == 0) updateTwoPixelClusterCOG = kFALSE;
			if (printOutSetupFileParameters) cout << "updateTwoPixelClusterCOG = " << updateTwoPixelClusterCOG << endl;
		}
		
		if (sline.Contains("updateTwoPixelClusterCOGshort_1D ="))
		{
			sline.ReplaceAll("updateTwoPixelClusterCOGshort_1D =","");
			buf = sline.Atoi();
			updateTwoPixelClusterCOGshort_1D = kTRUE;
			if (buf == 0) updateTwoPixelClusterCOGshort_1D = kFALSE;
			if (printOutSetupFileParameters) cout << "updateTwoPixelClusterCOGshort_1D = " << updateTwoPixelClusterCOGshort_1D << endl;
		}
		
		if (sline.Contains("updateTwoPixelClusterCOGlong_1D ="))
		{
			sline.ReplaceAll("updateTwoPixelClusterCOGlong_1D =","");
			buf = sline.Atoi();
			updateTwoPixelClusterCOGlong_1D = kTRUE;
			if (buf == 0) updateTwoPixelClusterCOGlong_1D = kFALSE;
			if (printOutSetupFileParameters) cout << "updateTwoPixelClusterCOGlong_1D = " << updateTwoPixelClusterCOGlong_1D << endl;
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
			//if (matrix_trigs[k][m] == 1 && doNotUseCornerTouchingTriggeredPixels4Clustering && (k-ix == m-iy || k-ix == iy-m)) continue;
			//if (matrix_trigs[k][m] == 0 && doNotUseCornerTouchingPixels4Clustering && (k-ix == m-iy || k-ix == iy-m)) continue;
			//if (matrix_trigs[k][m] == 1 && doNotUseCornerTouchingTriggeredPixels4Clustering && (k-ix == m-iy || k-ix == iy-m || ix-k == iy-m || ix-k == m-iy)) continue;
			//if (matrix_trigs[k][m] == 0 && doNotUseCornerTouchingPixels4Clustering && (k-ix == m-iy || k-ix == iy-m || ix-k == iy-m || ix-k == m-iy)) continue;
			if (matrix_trigs[k][m] == 1 && doNotUseCornerTouchingTriggeredPixels4Clustering && k == ix-win && m == iy-win) continue;
			if (matrix_trigs[k][m] == 1 && doNotUseCornerTouchingTriggeredPixels4Clustering && k == ix+win && m == iy+win) continue;
			if (matrix_trigs[k][m] == 1 && doNotUseCornerTouchingTriggeredPixels4Clustering && k == ix-win && m == iy+win) continue;
			if (matrix_trigs[k][m] == 1 && doNotUseCornerTouchingTriggeredPixels4Clustering && k == ix+win && m == iy-win) continue;
			if (matrix_trigs[k][m] == 0 && doNotUseCornerTouchingPixels4Clustering && k == ix-win && m == iy-win) continue;
			if (matrix_trigs[k][m] == 0 && doNotUseCornerTouchingPixels4Clustering && k == ix+win && m == iy+win) continue;
			if (matrix_trigs[k][m] == 0 && doNotUseCornerTouchingPixels4Clustering && k == ix-win && m == iy+win) continue;
			if (matrix_trigs[k][m] == 0 && doNotUseCornerTouchingPixels4Clustering && k == ix+win && m == iy-win) continue;
			//if (applyMaxDtCut4Clustering && TMath::Abs(matrix_timeCalib[ix][iy] - matrix_timeCalib[k][m]) > maxDtBetweenPixel4Clustering) continue;
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

Double_t fitfun(Double_t *x, Double_t *par)
{
	return par[0]*TMath::Cos(2*x[0]*TMath::DegToRad()-TMath::Pi())+par[1];
}

void sortClusters(const Int_t am)
{
	firstClusterIdx[am] = 0;
	secondClusterIdx[am] = 1;
	/*
	Float_t xx = 6;
	Float_t yy = 17;
	if (am == 1)
	{
		xx = 4;
		yy = 16;
	}
	Float_t dist1 = TMath::Hypot(buffClusterX[am][0]-xx,buffClusterY[am][0]-yy);
	Float_t dist2 = TMath::Hypot(buffClusterX[am][1]-xx,buffClusterY[am][1]-yy);
	if (dist1 > dist2)
	{
		firstClusterIdx[am] = 1;
		secondClusterIdx[am] = 0;
	}
	*/
}

void updateClustersWithNeighboursEnergy(const Int_t am, const Int_t eve2)
{
	Int_t vsize = buff_E_neighb_4sum.size();
	Int_t *idx = new Int_t[vsize];
	Float_t *arr = new Float_t[vsize];
	for (Int_t ia = 0; ia < vsize; ia++) arr[ia] = buff_E_neighb_4sum[ia];
	TMath::Sort(vsize,arr,idx,1);
	Float_t x2new[2] = {0}, y2new[2] = {0}, totE2pix[2] = {0};
	Int_t idx_2pix[2] = {-1,-1};
	Bool_t COGupdated = kFALSE;
	Int_t nbAdded[2] = {0,0};
	
	for (Int_t ia = 0; ia < vsize; ia++) 
	{
		Int_t ix = buff_X_neighb_4sum[idx[ia]];
		Int_t iy = buff_Y_neighb_4sum[idx[ia]];
		Float_t minDTime = 9999;
		Int_t minDTime_x = -1;
		Int_t minDTime_y = -1;
		Float_t maxEn = -9999;
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
				if (matrix_trigs_4sum[k][m] < 1) continue;
				Float_t dist = TMath::Abs(matrix_timeCalib[ix][iy]-matrix_timeCalib[k][m]);
				if (dist < minDTime)
				{
					minDTime = dist;
					minDTime_x = k;
					minDTime_y = m;
					maxEn = arr[idx[ia]];
				}
			}
		}
		if (maxEn > minPixelEnergyThr4ClusterReconstruction)
		{
			Int_t clusterN = matrix_trigs_4sum[minDTime_x][minDTime_y]-1;
			Float_t dt = matrix_timeCalib[ix][iy]-matrix_timeCalib[minDTime_x][minDTime_y];
			Bool_t isGood = kFALSE;
			if (buffClusterArealoc[clusterN] == 1 && enableUpdateEnergyClusters)
			{	
				if (dt < 0 && -dt < maxDtBetweenSinglePixClusters4updatingEnergy_up) isGood = kTRUE;
				if (dt >= 0 && dt < maxDtBetweenSinglePixClusters4updatingEnergy_low) isGood = kTRUE;
				pixelSpecNeighbour1PixClusters[am]->Fill(maxEn);
				energySpecCorr1PixClusters[am]->Fill(buffClusterEloc[clusterN],maxEn);
				anodeTiming1PixelClustersNeigbCorr[am]->Fill(matrix_timeCalib[minDTime_x][minDTime_y],matrix_timeCalib[ix][iy]);				
				dAnodeTiming_vs_nE_1PixelClustersNeigb[am]->Fill(maxEn,matrix_timeCalib[minDTime_x][minDTime_y]-matrix_timeCalib[ix][iy]);
				dAnodeTiming_vs_sumE_1PixelClustersNeigb[am]->Fill(buffClusterEloc[clusterN]+maxEn,matrix_timeCalib[minDTime_x][minDTime_y]-matrix_timeCalib[ix][iy]);
				dAnodeTiming_vs_NeigbAnodeTiming_1PixelClusters[am]->Fill(matrix_timeCalib[ix][iy],matrix_timeCalib[minDTime_x][minDTime_y]-matrix_timeCalib[ix][iy]);
				anodeTimingNeigb_nE_1PixelClusters[am]->Fill(maxEn,matrix_timeCalib[ix][iy]);
				if (updateSinglePixelClusterCOG)
				{
					Float_t x1new = (minDTime_x+0.5)*buffClusterEloc[clusterN];
					Float_t y1new = (minDTime_y+0.5)*buffClusterEloc[clusterN];
					x1new += (ix+0.5)*maxEn;
					y1new += (iy+0.5)*maxEn;
					buffClusterXloc[clusterN] = x1new/(buffClusterEloc[clusterN]+maxEn);
					buffClusterYloc[clusterN] = y1new/(buffClusterEloc[clusterN]+maxEn);
				}
				if (isGood) nbAdded[clusterN]++;
			}
			if (buffClusterArealoc[clusterN] == 2 && enableUpdateEnergyClusters)
			{
				if (dt < 0 && -dt < maxDtBetweenTwoPixClusters4updatingEnergy_up) isGood = kTRUE;
				if (dt >= 0 && dt < maxDtBetweenTwoPixClusters4updatingEnergy_low) isGood = kTRUE;
				anodeTiming2PixelClustersNeigbCorr[am]->Fill(matrix_timeCalib[minDTime_x][minDTime_y],matrix_timeCalib[ix][iy]);
				anodeTimingNeigb_nE_1PixelClusters[am]->Fill(maxEn,matrix_timeCalib[ix][iy]);
				if (updateTwoPixelClusterCOG)
				{
					x2new[clusterN] += (ix+0.5)*maxEn;
					y2new[clusterN] += (iy+0.5)*maxEn;
					COGupdated = kTRUE;
					totE2pix[clusterN] += maxEn;
				}
				if (updateTwoPixelClusterCOGshort_1D)
				{
					if (buffClusterOrientationloc[clusterN] == 0)
					{
						y2new[clusterN] += (iy+0.5)*maxEn;
						COGupdated = kTRUE;
					}
					if (buffClusterOrientationloc[clusterN] == 1)
					{
						x2new[clusterN] += (ix+0.5)*maxEn;
						COGupdated = kTRUE;
					}
				}
				if (updateTwoPixelClusterCOGlong_1D)
				{
					if (buffClusterOrientationloc[clusterN] == 0)
					{
						x2new[clusterN] += (ix+0.5)*maxEn;
						COGupdated = kTRUE;
					}
					if (buffClusterOrientationloc[clusterN] == 1)
					{
						y2new[clusterN] += (iy+0.5)*maxEn;
						COGupdated = kTRUE;
					}
				}		
				if (isGood) nbAdded[clusterN]++;
			}
			if (isGood)
			{
				buffClusterEloc[clusterN] += maxEn;
				matrix_trigs_4sum[minDTime_x][minDTime_y] = 0;
				if (buffClusterArealoc.size() == 1) clusterSpec_1ClusterEvents_summed[am]->Fill(buffClusterEloc[0]);
				if (eventDisplayFlag) imageE[am]->SetBinContent(ix+1,iy+1,maxEn);
			}
		}
	}
	if (COGupdated)
	{
		for (Int_t ia = 0; ia < buffClusterEloc.size(); ia++)
		{
			if (buffClusterArealoc[ia] != 2) continue;
			totE2pix[ia] += buffClusterEorigloc[ia];
			x2new[ia] += buffClusterXloc[ia]*buffClusterEorigloc[ia];
			y2new[ia] += buffClusterYloc[ia]*buffClusterEorigloc[ia];
		}
		for (Int_t ia = 0; ia < buffClusterEloc.size(); ia++)
		{
			if (buffClusterArealoc[ia] != 2 || TMath::Abs(totE2pix[ia]) < 1) continue;
			buffClusterXloc[ia] = x2new[ia]/totE2pix[ia];
			buffClusterYloc[ia] = y2new[ia]/totE2pix[ia];
		}
	}
	for (Int_t ia = 0; ia < buffClusterEloc.size(); ia++)
	{
		if (buffClusterArealoc[ia] == 1) clusterSummingStats1PixClusters[am]->Fill(nbAdded[ia]);
		if (buffClusterArealoc[ia] == 2) clusterSummingStats2PixClusters[am]->Fill(nbAdded[ia]);
	}
	
	delete arr, idx;
}

Double_t fitMAC(Double_t *x, Double_t *par)
{
	return par[0]*TMath::Power(x[0]-par[1],par[2]);
}

Float_t getThetaAngleFromHitsXYt(const Int_t am, const Float_t x1, const Float_t x2, const Float_t y1, const Float_t y2, const Float_t t1, const Float_t t2)
{
	return TMath::ATan(TMath::Hypot(x1-x2,y1-y2)*pixelPitch/(TMath::Abs(t1-t2)/factor2ConvertAnodeTime2Distance[am]))*TMath::RadToDeg();
}
