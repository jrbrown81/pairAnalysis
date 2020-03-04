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

const TString ver = "v13";
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

const Float_t eRestMass = 511;
Float_t Na22Energy = 511;
Float_t XRF_Cd = 23.17;
Float_t XRF_Te = 27.47;
Float_t angleOfSecondHead = 90;
Float_t relativePhiAngle = 90;
Float_t Ewindow4PhiAnalysis_min[2] = {530,530};
Float_t Ewindow4PhiAnalysis_max[2] = {530,530};

Float_t AMcoords[2][3] = {{0,0,-80},{0,0,80}};

Float_t Theta1WindowFor4PhiAnalysis_min = 60;
Float_t Theta1WindowFor4PhiAnalysis_max = 90;
Float_t Theta2WindowFor4PhiAnalysis_min = 60;
Float_t Theta2WindowFor4PhiAnalysis_max = 90;

Float_t bestTheta1WindowFor4PhiAnalysis_min = 60;
Float_t bestTheta1WindowFor4PhiAnalysis_max = 90;
Float_t bestTheta2WindowFor4PhiAnalysis_min = 60;
Float_t bestTheta2WindowFor4PhiAnalysis_max = 90;

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
TH1F *dPhiAngle_bestTheta;
TH1F *dPhiAngle_bestTheta_Ewin;
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
TH2F *theta0_theta1_FromFirstClusterE_Ewin;
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

void setupModuleHistos(Int_t im)
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

void setupGlobalHistos()
{
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

	dPhiXYtAngle_Ewin = new TH1F("dPhiXYtAngle_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, using X-Y-t #theta window",spectrumName.Data(),Ewindow4PhiAnalysis_min[0],Ewindow4PhiAnalysis_max[0]),nBins_dPhi,-180,180);
	dPhiXYtAngle_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngle_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiXYtAngleNorm = new TH1F("dPhiXYtAngleNorm",Form("%s, #Delta#varphi, normalised at #Delta#varphi = 0, using X-Y-t #theta window",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiXYtAngleNorm->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngleNorm->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiXYtAngleNorm_Ewin = new TH1F("dPhiXYtAngleNorm_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0, using X-Y-t #theta window",spectrumName.Data(),Ewindow4PhiAnalysis_min[0],Ewindow4PhiAnalysis_max[0]),nBins_dPhi,-180,180);
	dPhiXYtAngleNorm_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngleNorm_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiXYtAngle1 = new TH1F("dPhiXYtAngle1",Form("%s, #Delta#varphi using X-Y-t #theta window",spectrumName.Data()),nBins_dPhi/2,0,180);
	dPhiXYtAngle1->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngle1->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiXYtAngle1_Ewin = new TH1F("dPhiXYtAngle1_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, using X-Y-t #theta window",spectrumName.Data(),Ewindow4PhiAnalysis_min[0],Ewindow4PhiAnalysis_max[0]),nBins_dPhi/2,0,180);
	dPhiXYtAngle1_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngle1_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiXYtAngle1Norm = new TH1F("dPhiXYtAngle1Norm",Form("%s, #Delta#varphi, normalised at #Delta#varphi = 0, using X-Y-t #theta window",spectrumName.Data()),nBins_dPhi/2,0,180);
	dPhiXYtAngle1Norm->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngle1Norm->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiXYtAngle1Norm_Ewin = new TH1F("dPhiXYtAngle1Norm_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0, using X-Y-t #theta window",spectrumName.Data(),Ewindow4PhiAnalysis_min[0],Ewindow4PhiAnalysis_max[0]),nBins_dPhi/2,0,180);
	dPhiXYtAngle1Norm_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiXYtAngle1Norm_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiAngle = new TH1F("dPhiAngle",Form("%s, #Delta#varphi",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiAngle->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiAngle_bestTheta = new TH1F("dPhiAngle_bestTheta",Form("%s, #Delta#varphi",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiAngle_bestTheta->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_bestTheta->GetXaxis()->SetTitle("#Delta#varphi, degrees");

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

	dPhiAngle_Ewin = new TH1F("dPhiAngle_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV",spectrumName.Data(),Ewindow4PhiAnalysis_min[0],Ewindow4PhiAnalysis_max[0]),nBins_dPhi,-180,180);
	dPhiAngle_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");
	
	dPhiAngle_bestTheta_Ewin = new TH1F("dPhiAngle_bestTheta_Ewin",Form("%s, #Delta#varphi",spectrumName.Data()),nBins_dPhi,-180,180);
	dPhiAngle_bestTheta_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_bestTheta_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiAngle_Ewin_w = new TH1F("dPhiAngle_Ewin_w",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, weighted",spectrumName.Data(),Ewindow4PhiAnalysis_min[0],Ewindow4PhiAnalysis_max[0]),nBins_dPhi,-180,180);
	dPhiAngle_Ewin_w->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_Ewin_w->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiAngleNorm_Ewin = new TH1F("dPhiAngleNorm_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0",spectrumName.Data(),Ewindow4PhiAnalysis_min[0],Ewindow4PhiAnalysis_max[0]),nBins_dPhi,-180,180);
	dPhiAngleNorm_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngleNorm_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiAngleNorm_Ewin_w = new TH1F("dPhiAngleNorm_Ewin_w",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0, weighted",spectrumName.Data(),Ewindow4PhiAnalysis_min[0],Ewindow4PhiAnalysis_max[0]),nBins_dPhi,-180,180);
	dPhiAngleNorm_Ewin_w->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngleNorm_Ewin_w->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiAngle1_Ewin = new TH1F("dPhiAngle1_Ewin",Form("%s, |#Delta#varphi| within  %.0f keV < E < %.0f keV",spectrumName.Data(),Ewindow4PhiAnalysis_min[0],Ewindow4PhiAnalysis_max[0]),nBins_dPhi/2,0,180);
	dPhiAngle1_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle1_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");

	dPhiAngle1norm_Ewin = new TH1F("dPhiAngle1norm_Ewin",Form("%s, |#Delta#varphi| within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0",spectrumName.Data(),Ewindow4PhiAnalysis_min[0],Ewindow4PhiAnalysis_max[0]),nBins_dPhi/2,0,180);
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

	theta0_theta1_FromFirstClusterE_Ewin = new TH2F("theta0_theta1_FromFirstClusterE_Ewin",Form("%s, #theta_{0} - #theta_{1} calculated from first interaction energy, with E window",spectrumName.Data()),nBinsTheta4H,0,180,nBinsTheta4H,0,180);
	theta0_theta1_FromFirstClusterE_Ewin->GetXaxis()->SetTitleOffset(1.2);
	theta0_theta1_FromFirstClusterE_Ewin->GetXaxis()->SetTitle("#theta in AM0, degrees");
	theta0_theta1_FromFirstClusterE_Ewin->GetYaxis()->SetTitleOffset(1.4);
	theta0_theta1_FromFirstClusterE_Ewin->GetYaxis()->SetTitle("#theta in AM1, degrees");
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
		
		if (sline.Contains("E1window4PhiAnalysis_min ="))
		{
			sline.ReplaceAll("E1window4PhiAnalysis_min =","");
			Ewindow4PhiAnalysis_min[0] = sline.Atof();
			if (printOutSetupFileParameters) cout << "E1window4PhiAnalysis_min = " << Ewindow4PhiAnalysis_min[0] << endl;
		}
		if (sline.Contains("E1window4PhiAnalysis_max ="))
		{
			sline.ReplaceAll("E1window4PhiAnalysis_max =","");
			Ewindow4PhiAnalysis_max[0] = sline.Atof();
			if (printOutSetupFileParameters) cout << "E1window4PhiAnalysis_max = " << Ewindow4PhiAnalysis_max[0] << endl;
		}
		if (sline.Contains("E2window4PhiAnalysis_min ="))
		{
			sline.ReplaceAll("E2window4PhiAnalysis_min =","");
			Ewindow4PhiAnalysis_min[1] = sline.Atof();
			if (printOutSetupFileParameters) cout << "E2window4PhiAnalysis_min = " << Ewindow4PhiAnalysis_min[1] << endl;
		}
		if (sline.Contains("E2window4PhiAnalysis_max ="))
		{
			sline.ReplaceAll("E2window4PhiAnalysis_max =","");
			Ewindow4PhiAnalysis_max[1] = sline.Atof();
			if (printOutSetupFileParameters) cout << "E2window4PhiAnalysis_max = " << Ewindow4PhiAnalysis_max[1] << endl;
		}
		
		if (sline.Contains("Theta1WindowFor4PhiAnalysis_min ="))
		{
			sline.ReplaceAll("Theta1WindowFor4PhiAnalysis_min =","");
			Theta1WindowFor4PhiAnalysis_min = sline.Atof();
			if (printOutSetupFileParameters) cout << "Theta1WindowFor4PhiAnalysis_min = " << Theta1WindowFor4PhiAnalysis_min << endl;
		}
		if (sline.Contains("Theta1WindowFor4PhiAnalysis_max ="))
		{
			sline.ReplaceAll("Theta1WindowFor4PhiAnalysis_max =","");
			Theta1WindowFor4PhiAnalysis_max = sline.Atof();
			if (printOutSetupFileParameters) cout << "Theta1WindowFor4PhiAnalysis_max = " << Theta1WindowFor4PhiAnalysis_max << endl;
		}
		if (sline.Contains("Theta2WindowFor4PhiAnalysis_min ="))
		{
			sline.ReplaceAll("Theta2WindowFor4PhiAnalysis_min =","");
			Theta2WindowFor4PhiAnalysis_min = sline.Atof();
			if (printOutSetupFileParameters) cout << "Theta2WindowFor4PhiAnalysis_min = " << Theta2WindowFor4PhiAnalysis_min << endl;
		}
		if (sline.Contains("Theta2WindowFor4PhiAnalysis_max ="))
		{
			sline.ReplaceAll("Theta2WindowFor4PhiAnalysis_max =","");
			Theta2WindowFor4PhiAnalysis_max = sline.Atof();
			if (printOutSetupFileParameters) cout << "Theta2WindowFor4PhiAnalysis_max = " << Theta2WindowFor4PhiAnalysis_max << endl;
		}
		
		if (sline.Contains("theta1Best4PhiAnalysis_min ="))
		{
			sline.ReplaceAll("theta1Best4PhiAnalysis_min =","");
			bestTheta1WindowFor4PhiAnalysis_min = sline.Atof();
			if (printOutSetupFileParameters) cout << "theta1Best4PhiAnalysis_min = " << bestTheta1WindowFor4PhiAnalysis_min << endl;
		}
		if (sline.Contains("theta1Best4PhiAnalysis_max ="))
		{
			sline.ReplaceAll("theta1Best4PhiAnalysis_max =","");
			bestTheta1WindowFor4PhiAnalysis_max = sline.Atof();
			if (printOutSetupFileParameters) cout << "theta1Best4PhiAnalysis_max = " << bestTheta1WindowFor4PhiAnalysis_max << endl;
		}
		if (sline.Contains("theta2Best4PhiAnalysis_min ="))
		{
			sline.ReplaceAll("theta2Best4PhiAnalysis_min =","");
			bestTheta2WindowFor4PhiAnalysis_min = sline.Atof();
			if (printOutSetupFileParameters) cout << "theta2Best4PhiAnalysis_min = " << bestTheta2WindowFor4PhiAnalysis_min << endl;
		}
		if (sline.Contains("theta2Best4PhiAnalysis_max ="))
		{
			sline.ReplaceAll("theta2Best4PhiAnalysis_max =","");
			bestTheta2WindowFor4PhiAnalysis_max = sline.Atof();
			if (printOutSetupFileParameters) cout << "theta2Best4PhiAnalysis_max = " << bestTheta2WindowFor4PhiAnalysis_max << endl;
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
