//#define MATLAB_MEX_FILE 
#include "SeqSim.h"		// DSP header
#include "fftw3.h"		// The Fastest Fourier Transform in the West
#include <algorithm>    // std::min_element, std::max_element
#include <iostream>
#include <fstream>

#include <cstring>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

namespace SEQSIM
{

	#ifdef MATLAB_MEX_FILE
		#define PRINT mexPrintf
	#else
		#define PRINT std::printf
		#define ERROR throw std::invalid_argument
		void mexEvalString(const char *str) {}
		void mexErrMsgTxt(const char *err_msg) {
			printf ("%s\n", err_msg );
			exit(1);
		}
	#endif
		
	// sprintf
	#pragma warning(disable:4996)

	//*******************************************************************
	// ENUMS/DEFINITIONS/STRUCTURES/TEMPLATES                           
	//******************************************************************* 
	static const double sfGRT		= 10.0e-6;			// Gradient raster time in s
	static const double sfPI		= 3.14159265359;	// Mathematical constant
	static const double sfGAMMA_1H	= 42.575575;		// Gyromagnetic ratio of 1H in MHz/T

	static double sfEmpiricalFactor	= 1.0e-2;			// The amplitude definition of the eddy currents is not clear


	#define DECAY_TIME_FACTOR 5 // Simulate the exponential decay of the eddy currents for a time period of DECAY_TIME_FACTOR*Longest_Time_Constant

	template <class T> T& max( T& a,  T& b) {
		return (a<b) ? b : a;     // or: return comp(a,b)?b:a; for version (2)
	}

	template <class T>  T& min( T& a,  T& b) {
		return !(b<a) ? a : b;     // or: return !comp(b,a)?a:b; for version (2)
	}

	void mexErrMsgOriginAndTxt(const char* pModule, const char* Meassage) {

        string full_text = string(pModule) + string(Meassage);

		mexErrMsgTxt(full_text.c_str());
	}

	void mexErrMsgOriginAndTxt(const char* pModule, const char* Meassage, const char* Dump) {

        string full_text = string(pModule) + string(Meassage) + string(Dump);	

		mexErrMsgTxt(full_text.c_str());
	}


	// Splits a full path into component file parts
	FileParts fileparts(string &fullpath)
	{

		size_t idxSlash = fullpath.rfind("/");
		if (idxSlash == string::npos) {
			idxSlash = fullpath.rfind("\\");
		}
		size_t idxDot = fullpath.rfind(".");

		FileParts fp;
		if (idxSlash != string::npos && idxDot != string::npos) {
			fp.path = fullpath.substr(0, idxSlash + 1);
			fp.name = fullpath.substr(idxSlash + 1, idxDot - idxSlash - 1);
			fp.ext = fullpath.substr(idxDot);
		}
		else if (idxSlash == string::npos && idxDot == string::npos) {
			fp.name = fullpath;
		}
		else if ( idxSlash == string::npos) {
			fp.name = fullpath.substr(0, idxDot);
			fp.ext = fullpath.substr(idxDot);
		}
		else { // only idxDot == string::npos
			fp.path = fullpath.substr(0, idxSlash + 1);
			fp.name = fullpath.substr(idxSlash + 1);
		}
		return fp;
	}


	//*******************************************************************
	// CONSTRUCTOR                                                       
	//******************************************************************* 
	DSP::DSP() :
		m_lGradientShapeLength(0),
		m_afMULTI_PURPOSEX(NULL),
		m_afMULTI_PURPOSEY(NULL),
		m_afMULTI_PURPOSEZ(NULL),
		m_afMULTI_PURPOSE_INTERPX(NULL),
		m_afMULTI_PURPOSE_INTERPY(NULL),
		m_afMULTI_PURPOSE_INTERPZ(NULL),
		m_uiVerbose(DISPLAY_NONE),
		m_fLargestTau(0),
		m_lExponentialLength(0),
		m_lConvolutionlLength(0),
		m_afExponentialX(NULL),
		m_afExponentialY(NULL),
		m_afExponentialZ(NULL),
		m_lRXSampleLength(0),
		m_lRXEvents(0),
		m_lTXEvents(0),
		m_lTrigEvents(0),
		m_adRXTimes(NULL),
		m_adTrigTimes(NULL),
		m_lCurrentRXSampleLength(0),
		m_adTXCenterTimes(NULL),
		m_bDebugMode(false),
		m_aiMatrixA(NULL),
		m_aiMatrixB(NULL),
		m_uiRXEventLength(NULL),
		m_lCurrentGCSampleLength(0),
		m_bEccCompensationAvailable(false),
		m_lCurrentTXNumber(0),
		m_lCurrentRXNumber(0),
		m_lDataType(DataType::GRADIENT),
		m_lOutputMode(OutputMode::FULL)
	{  
		m_sDSVFileNamePrefix = string("");
		
		// Get start time
		m_tstart = time(0);
	
		// Create matrices
		m_aiMatrixA = new double*[3];
		m_aiMatrixB = new double*[3];

		for (int i = 0; i < 3; ++i) {
			m_aiMatrixA[i] = new double[3];
			m_aiMatrixB[i] = new double[3];
		}

		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				m_aiMatrixA[i][j] = 0.0;
			}
		}
	}

	//*******************************************************************
	// DESTRUCTOR                                                        
	//******************************************************************* 
	DSP::~DSP(){
    
		// Delete matrices
		for (int i = 0; i < 3; ++i) {
			if (m_aiMatrixA[i] != NULL) {
				delete[] m_aiMatrixA[i];
			}
			if (m_aiMatrixB[i] != NULL) {
				delete[] m_aiMatrixB[i];
			}
		}

		if (m_aiMatrixA != NULL) {
			delete[] m_aiMatrixA;
			m_aiMatrixA = NULL;
		}

		if (m_aiMatrixB != NULL) {
			delete[] m_aiMatrixB;
			m_aiMatrixB = NULL;
		}
		
				
		// Delete exponentials
		if (m_afExponentialX != NULL) {
			delete[] m_afExponentialX;
			m_afExponentialX = NULL;
		}

		if (m_afExponentialY != NULL) {
			delete[] m_afExponentialY;
			m_afExponentialY = NULL;
		}

		if (m_afExponentialZ != NULL) {
			delete[] m_afExponentialZ;
			m_afExponentialZ = NULL;
		}


		#ifdef MATLAB_MEX_FILE 

			if (m_lOutputMode == OutputMode::FULL) {
				// The full data is returned to Matlab
				// Delete the interpolated data

				if (m_afMULTI_PURPOSE_INTERPX != NULL) {
					delete[] m_afMULTI_PURPOSE_INTERPX;
					m_afMULTI_PURPOSE_INTERPX = NULL;
				}
				if (m_afMULTI_PURPOSE_INTERPY != NULL) {
					delete[] m_afMULTI_PURPOSE_INTERPY;
					m_afMULTI_PURPOSE_INTERPY = NULL;
				}
				if (m_afMULTI_PURPOSE_INTERPZ != NULL) {
					delete[] m_afMULTI_PURPOSE_INTERPZ;
					m_afMULTI_PURPOSE_INTERPZ = NULL;
				}

			}
			else {

				// The interpolated data is returned to Matlab
				// Delete the full shape

				if (m_afMULTI_PURPOSEX != NULL) {
					delete[] m_afMULTI_PURPOSEX;
					m_afMULTI_PURPOSEX = NULL;
				}
				if (m_afMULTI_PURPOSEY != NULL) {
					delete[] m_afMULTI_PURPOSEY;
					m_afMULTI_PURPOSEY = NULL;
				}
				if (m_afMULTI_PURPOSEZ != NULL) {
					delete[] m_afMULTI_PURPOSEZ;
					m_afMULTI_PURPOSEZ = NULL;
				}
			}

		#else	

			// If the arrays are not returned to Matlab
			// throw everything into the trash bin

			// Delete the full shape
			if (m_afMULTI_PURPOSEX != NULL) {
				delete[] m_afMULTI_PURPOSEX;
				m_afMULTI_PURPOSEX = NULL;
			}
			if (m_afMULTI_PURPOSEY != NULL) {
				delete[] m_afMULTI_PURPOSEY;
				m_afMULTI_PURPOSEY = NULL;
			}
			if (m_afMULTI_PURPOSEZ != NULL) {
				delete[] m_afMULTI_PURPOSEZ;
				m_afMULTI_PURPOSEZ = NULL;
			}

			// Delete the interpolated data
			if (m_afMULTI_PURPOSE_INTERPX != NULL) {
				delete[] m_afMULTI_PURPOSE_INTERPX;
				m_afMULTI_PURPOSE_INTERPX = NULL;
			}
			if (m_afMULTI_PURPOSE_INTERPY != NULL) {
				delete[] m_afMULTI_PURPOSE_INTERPY;
				m_afMULTI_PURPOSE_INTERPY = NULL;
			}
			if (m_afMULTI_PURPOSE_INTERPZ != NULL) {
				delete[] m_afMULTI_PURPOSE_INTERPZ;
				m_afMULTI_PURPOSE_INTERPZ = NULL;
			}

			// Delete RX start indices
			if (m_uiRXEventLength != NULL) {
				delete[] m_uiRXEventLength;
				m_uiRXEventLength = NULL;
			}

			// Delete RX times
			if (m_adRXTimes != NULL) {
				delete[] m_adRXTimes;
				m_adRXTimes = NULL;
			}

			// Delete TX times
			if (m_adTXCenterTimes != NULL) {
				delete[] m_adTXCenterTimes;
				m_adTXCenterTimes = NULL;
			}

			// Delete Trigger times
			if (m_adTrigTimes != NULL) {
				delete[] m_adTrigTimes;
				m_adTrigTimes = NULL;
			}
			
		#endif	

		// Get end time
		m_tend = time(0);   

		char timestr[4096];
		getDurationString(m_tend - m_tstart, timestr);

		if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) { PRINT("Computation finished in %s! \n", timestr); mexEvalString("drawnow"); }
	}


	//*******************************************************************
	// SET FUNCTIONS                                                     
	//******************************************************************* 
	void DSP::setFileName(const char * pFileName){

		string fullpath(pFileName);
		m_fpXMLFile = fileparts(fullpath);
		 
		cout << "Setting XML file: " << fullpath << endl;

	}
	 
	// The class assumes that the input files are names GX.dsv etc
	// A prefix for the files can be defines to also read files called
	// e.g. temp_GX.dsv, i.e. setDSVFileNamePrefix("temp_")
	void DSP::setDSVFileNamePrefix(const char * pDSVFileNamePrefix){
		m_sDSVFileNamePrefix = string(pDSVFileNamePrefix);
		cout << "Using the following prefix for the DSV files: " << m_sDSVFileNamePrefix << endl;
	}
	 

	 void DSP::setDSVFolderPath(const char * pFolderPath) {
		 m_sDSVFolderPath = string(pFolderPath);
		 cout << "Setting DSV folder: " << m_sDSVFolderPath << endl;
	 }


	 // See enum Verbose in types.h 
	 void DSP::setVerboseMode(int verbose) {
		 m_uiVerbose = verbose; 
	 }
 

	 void  DSP::setDataType(DataType type) {
		 m_lDataType = type;
	 }


	 void  DSP::setOutputMode(OutputMode mode) {
		 m_lOutputMode = mode;
	 }
	 

	 void DSP::setB0CorrCoeff(ECC_Coeff &CoeffX, ECC_Coeff &CoeffY, ECC_Coeff &CoeffZ) {
	 
		 m_ECC_Coeff_X = CoeffX;
		 m_ECC_Coeff_Y = CoeffY;
		 m_ECC_Coeff_Z = CoeffZ;

	 }

	void DSP::setB0CorrCoeff(	vector<double> &vfCoeffXamp, vector<double> &vfCoeffXtau,
								vector<double> &vfCoeffYamp, vector<double> &vfCoeffYtau,
								vector<double> &vfCoeffZamp, vector<double> &vfCoeffZtau     ) {

		 ECC_Coeff CoeffX, CoeffY, CoeffZ;

		 CoeffX.vfAmp = vfCoeffXamp;
		 CoeffX.vfTau = vfCoeffXtau;

		 CoeffY.vfAmp = vfCoeffYamp;
		 CoeffY.vfTau = vfCoeffYtau;

		 CoeffZ.vfAmp = vfCoeffZamp;
		 CoeffZ.vfTau = vfCoeffZtau;

		 setB0CorrCoeff(CoeffX, CoeffY, CoeffZ);
	}


	//*******************************************************************
	// GETDURATIONSTRING()                                               
	//******************************************************************* 

	// Convert duration given in seconds into a time string: X hours, Y minutes and Z seconds
	void DSP::getDurationString(long duration, char *str) {

		int seconds, hours, minutes;
		minutes = duration / 60;
		seconds = duration % 60;
		hours = minutes / 60;
		minutes = minutes % 60;

		char hs[2], ms[2], ss[2];
		str[0] = '\0';
		if (hours > 1)   {hs[0] = 's'; hs[1] = '\0';}	else hs[0] = '\0';
		if (minutes > 1) {ms[0] = 's'; ms[1] = '\0';}	else ms[0] = '\0';
		if (seconds > 1) {ss[0] = 's'; ss[1] = '\0';}	else ss[0] = '\0';

		if (hours > 0)
			sprintf(str, "%d hour%s, %d minute%s and %d second%s", hours, hs, minutes, ms, seconds, ss);
		else if (minutes > 0)
			sprintf(str, "%d minute%s and %d second%s", minutes, ms, seconds, ss);
		else {
			if (seconds < 1) seconds = 1;
				sprintf(str, "%d second%s", seconds, ss);
			}
		}
 
	void DSP::writeToFile() {

		PRINT("Writing data to file...\n");

		m_fpXMLFile.path;

		ofstream outputFileX, outputFileY, outputFileZ;
		string fileNameX(m_fpXMLFile.path);
		string fileNameY(m_fpXMLFile.path);
		string fileNameZ(m_fpXMLFile.path);


		switch (m_lDataType) {
		case DataType::GRADIENT: 
			fileNameX.append("GX.txt");
			fileNameY.append("GY.txt");
			fileNameZ.append("GZ.txt");
			break;      
		case DataType::KSPACE:
			fileNameX.append("KX.txt");
			fileNameY.append("KY.txt");
			fileNameZ.append("KZ.txt");
			break;
		case DataType::SLEWRATE:
			fileNameX.append("SX.txt");
			fileNameY.append("SY.txt");
			fileNameZ.append("SZ.txt");
			break;
		case DataType::EDDYCURRENT:
			fileNameX.append("ECX.txt");
			fileNameY.append("ECY.txt");
			fileNameZ.append("ECZ.txt");
			break;
		case DataType::EDDYPHASE:
			fileNameX.append("EP.txt");
			break;
		}


		// Write data
		if (m_lDataType == DataType::EDDYPHASE) {

			// Open file
			outputFileX.open(fileNameX);
			
			if (outputFileX.is_open())
			{
				if (m_lOutputMode == OutputMode::FULL)
					for (long t = 0; t < m_lConvolutionlLength; t++) {
						outputFileX << m_afMULTI_PURPOSEX[t] << "\n";
					}
				else {
					for (long t = 0; t < m_lRXSampleLength; t++) {
						outputFileX << m_afMULTI_PURPOSE_INTERPX[t] << "\n";
					}
				}

				// Close file
				outputFileX.close();

			}
			else {
				throw "Unable to open file";
			}
		}
		else {

			// Open files
			outputFileX.open(fileNameX);
			outputFileY.open(fileNameY);
			outputFileZ.open(fileNameZ);

			if (outputFileX.is_open() && outputFileY.is_open() && outputFileZ.is_open())
			{
				if (m_lOutputMode == OutputMode::FULL)
					for (long t = 0; t < m_lConvolutionlLength; t++) {
						outputFileX << m_afMULTI_PURPOSEX[t] << "\n";
						outputFileY << m_afMULTI_PURPOSEY[t] << "\n";
						outputFileZ << m_afMULTI_PURPOSEZ[t] << "\n";
					}
				else {
					for (long t = 0; t < m_lRXSampleLength; t++) {
						outputFileX << m_afMULTI_PURPOSE_INTERPX[t] << "\n";
						outputFileY << m_afMULTI_PURPOSE_INTERPY[t] << "\n";
						outputFileZ << m_afMULTI_PURPOSE_INTERPZ[t] << "\n";
					}
				}

				// Close files
				outputFileX.close();
				outputFileY.close();
				outputFileZ.close();

			}
			else {
				throw "Unable to open file";
			}
		}
		
		//---------------------------------------------------------------------
		// RX Times & RX Samples                                                         
		//---------------------------------------------------------------------
		ofstream outputFile;
		outputFile.open(m_fpXMLFile.path + "RXTimes.txt");
		if (outputFile.is_open()) {
			for (long t = 0; t < m_lRXSampleLength; t++) {
				outputFile << m_adRXTimes[t] << "\n";
			}
			outputFile.close();
		}

		outputFile.open(m_fpXMLFile.path + "RXSamp.txt");
		if (outputFile.is_open()) {
			for (long t = 0; t < m_lRXEvents; t++) {
				outputFile << m_uiRXEventLength[t] << "\n";
			}
			outputFile.close();
		}

		//---------------------------------------------------------------------
		// TX Times                                                         
		//---------------------------------------------------------------------
		outputFile.open(m_fpXMLFile.path + "TXTimes.txt");
		if (outputFile.is_open()) {
			for (long t = 0; t < m_lTXEvents; t++) {
				outputFile << m_adTXCenterTimes[t] << "\n";
			}
			outputFile.close();
		}

		//---------------------------------------------------------------------
		// Trigger Times                                                         
		//---------------------------------------------------------------------
		outputFile.open(m_fpXMLFile.path + "m_lTrigEvents.txt");
		if (outputFile.is_open()) {
			for (long t = 0; t < m_lTXEvents; t++) {
				outputFile << m_adTrigTimes[t] << "\n";
			}
			outputFile.close();
		}

	}

	//*******************************************************************
	// CALMEMORYREQUIREMENT()                                            
	//******************************************************************* 

	// Loops through the text file and counts the number of gradient and receiver events
	// Thereafter the text file is rewound till the first Halt gradient instruction
	void DSP::calcMemoryRequirement(){

		static const char *ptModule = { "DSP::calcMemoryRequirement(): " };
	
		// Gradient shape length
		m_lGradientShapeLength = 0;

		// Number of readout events
		m_lRXEvents = 0;

		// Number of transmit events
		m_lTXEvents = 0;

		// Number of trigger events
		m_lTrigEvents = 0;

		// Number of ADC samples (must be the same for all imaging scans)
		long ADCSamples = -1;
		

		//---------------------------------------------------------------------
		// Get size of gradient shape
		//---------------------------------------------------------------------
		string sfileName;
		double dHoriDelta = 0.0;
		double dVertFactor = 0.0;
		long lLengthGY = 0;
		long lLengthGZ = 0;


		// GX
		sfileName = m_sDSVFolderPath + "/" + m_sDSVFileNamePrefix + "GRX.dsv";
		{
			ifstream in(sfileName.c_str());

			if (!in) {
				cout << "Cannot open input file: " << sfileName << endl;
				mexErrMsgOriginAndTxt(ptModule, "Check provided folder. DSV file needs to be named: GRX.dsv, GRY.dsv, GRZ.dsv, and INF.dsv (+ an optional prefix). \n");
			}
			this->getDSVParams(in, m_lGradientShapeLength, dHoriDelta, dVertFactor);
		}

		// GY
		sfileName = m_sDSVFolderPath + "/" + m_sDSVFileNamePrefix + "GRY.dsv";
		{
			ifstream in(sfileName.c_str());

			if (!in) {
				cout << "Cannot open input file: " << sfileName << endl;
				mexErrMsgOriginAndTxt(ptModule, "Check provided folder. DSV file needs to be named: GRX.dsv, GRY.dsv, GRZ.dsv, and INF.dsv (+ an optional prefix). \n");
			}
			this->getDSVParams(in, lLengthGY, dHoriDelta, dVertFactor);
		}

		// GZ
		sfileName = m_sDSVFolderPath + "/" + m_sDSVFileNamePrefix + "GRZ.dsv";
		{
			ifstream in(sfileName.c_str());

			if (!in) {
				cout << "Cannot open input file: " << sfileName << endl;
				mexErrMsgOriginAndTxt(ptModule, "Check provided folder. DSV file needs to be named: GRX.dsv, GRY.dsv, GRZ.dsv, and INF.dsv (+ an optional prefix). \n");
			}
			this->getDSVParams(in, lLengthGZ, dHoriDelta, dVertFactor);
		}

		if (m_lGradientShapeLength != lLengthGY)
			mexErrMsgOriginAndTxt(ptModule, "The length of GX and GY are not consistent.\n");

		if (m_lGradientShapeLength != lLengthGZ)
			mexErrMsgOriginAndTxt(ptModule, "The length of GX and GZ are not consistent.\n");


		//---------------------------------------------------------------------
		// Get file Parameters
		//---------------------------------------------------------------------
		sfileName = m_sDSVFolderPath + "/" + m_sDSVFileNamePrefix + "INF.dsv";
		{
			ifstream in(sfileName.c_str());

			if (!in) {
				cout << "Cannot open input file: " << sfileName << endl;
				mexErrMsgOriginAndTxt(ptModule, "Check provided folder. DSV file needs to be named: GRX.dsv, GRY.dsv, GRZ.dsv, and INF.dsv (+ an optional prefix). \n");
			}
			this->getInfoParams(in, m_lRXSampleLength, m_lRXEvents, m_lTXEvents);
		}


		//---------------------------------------------------------------------
		// Display some information                              
		//--------------------------------------------------------------------- 
		if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {

			long total = int(m_lGradientShapeLength * sfGRT);

			char str[4096];
			getDurationString(total, str);

			PRINT("Found DSP gradient instructions for %s!\n", str);
			PRINT("Found %ld RX samples in %ld RX events.\n", m_lRXSampleLength, m_lRXEvents);
			PRINT("Found %ld TX events.\n", m_lTXEvents);
			PRINT("Found %ld Trigger events (EXTR).\n", m_lTrigEvents);
			mexEvalString("drawnow");
		}

	}
 


	//*******************************************************************
	// COMPUTEEXPONENTIALS()                                             
	//******************************************************************* 

	// Calculate the exponential decay curves
	 void DSP::computeExponentials(){

		 if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {
			 PRINT("Computing exponentials...  \n");
			 mexEvalString("drawnow");
		 }

		 for (long t = 0; t < m_lExponentialLength; t++) {
			// X
			m_afExponentialX[t] = 0.0;
			for (int j = 0; j < m_ECC_Coeff_X.vfTau.size(); j++) {
				if (m_ECC_Coeff_X.vfTau.at(j)!=0)
					m_afExponentialX[t] += sfEmpiricalFactor*m_ECC_Coeff_X.vfAmp.at(j)*exp(-sfGRT*t / m_ECC_Coeff_X.vfTau.at(j));
			}

			// Y
			m_afExponentialY[t] = 0.0;
			for (int j = 0; j < m_ECC_Coeff_Y.vfTau.size(); j++) {
				if (m_ECC_Coeff_Y.vfTau.at(j) != 0)
					m_afExponentialY[t] += sfEmpiricalFactor*m_ECC_Coeff_Y.vfAmp.at(j)*exp(-sfGRT*t / m_ECC_Coeff_Y.vfTau.at(j));
			 }

			// Z
			m_afExponentialZ[t] = 0.0;
			for (int j = 0; j < m_ECC_Coeff_Z.vfTau.size(); j++) {
				 if (m_ECC_Coeff_Z.vfTau.at(j) != 0)
					 m_afExponentialZ[t] += sfEmpiricalFactor*m_ECC_Coeff_Z.vfAmp.at(j)*exp(-sfGRT*t / m_ECC_Coeff_Z.vfTau.at(j));
			}
		 }
		 
		// Set the rest to zero
		for (long t = m_lExponentialLength; t < m_lConvolutionlLength; t++) {
			m_afExponentialX[t] = 0.0;
			m_afExponentialY[t] = 0.0;
			m_afExponentialZ[t] = 0.0;
		}
	 }

	//*******************************************************************
	// DETERMINELONGESTTIMECONSTANT()                                    
	//******************************************************************* 

	// Determines the longest decay time on all axes
	void DSP::determineLongestTimeConstant() {

		if(	m_ECC_Coeff_X.vfTau.size() > 0 && m_ECC_Coeff_X.vfTau.size() > 0 &&
			m_ECC_Coeff_Y.vfTau.size() > 0 && m_ECC_Coeff_Y.vfTau.size() > 0 &&
			m_ECC_Coeff_Z.vfTau.size() > 0 && m_ECC_Coeff_Z.vfTau.size() > 0) {

			// Longest decay constant per axis
			double largest_tau_X = *max_element(m_ECC_Coeff_X.vfTau.begin(), m_ECC_Coeff_X.vfTau.end());
			double largest_tau_Y = *max_element(m_ECC_Coeff_Y.vfTau.begin(), m_ECC_Coeff_Y.vfTau.end());
			double largest_tau_Z = *max_element(m_ECC_Coeff_Z.vfTau.begin(), m_ECC_Coeff_Z.vfTau.end());

			// Longest overall decay constant
			m_fLargestTau = max(largest_tau_X, max(largest_tau_Y, largest_tau_Z));

			if (m_fLargestTau == 0) {
				PRINT("Longest decay constant is zero.\n");
			}
			if (m_fLargestTau > 5) {
				PRINT("Longest decay constant is longer than 5 seconds. This is unlikely.\n");
			}
			if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {
				PRINT("Longest decay constant %f seconds.\n", m_fLargestTau);
			}
			m_lExponentialLength = (long)(ceil(m_fLargestTau) * 1e5)*DECAY_TIME_FACTOR; // Gradient raster time is 10 us
			if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {
				PRINT("Exponential shape length = %ld.\n", m_lExponentialLength);
			}
		 }
		 else {
			// No coefficients have been provided
			PRINT("No ECC coefficients were provided.\n");
			m_fLargestTau = 0;
			m_lExponentialLength = 1;
		 }

	 }


	 //*******************************************************************
	 // COMPUTEECC()                                                      
	 //******************************************************************* 

	 // Computes the ECC via convolution for the specifed axis
	 // If afS is a Null pointer, the first input array will used as input and output array 
	 // If afS is a not a Null pointer, afS will used as input and afECC as output array 
	 void DSP::computeECC(double *adSlewRate, double *afExponential) {

		static const char *ptModule = { "DSP::computeECC(): " };

		// Allocate memory for FFT result
		fftw_complex *outGrad = fftw_alloc_complex(m_lConvolutionlLength);
		fftw_complex *outExp  = fftw_alloc_complex(m_lConvolutionlLength);

		// Create plan for FFT and execute
		fftw_plan plan;

		// FFT slew rate
		plan = fftw_plan_dft_r2c_1d(m_lConvolutionlLength, adSlewRate, outGrad, FFTW_ESTIMATE);
		fftw_execute(plan);
	 
		// FFT exponential
		plan = fftw_plan_dft_r2c_1d(m_lConvolutionlLength, afExponential, outExp, FFTW_ESTIMATE);
		fftw_execute(plan);

	 	 
		// Multiply both FFTs
		for (int i = 0; i < m_lConvolutionlLength; i++) {

			double real = outGrad[i][0] * outExp[i][0] - outGrad[i][1] * outExp[i][1];
			double imag = outGrad[i][1] * outExp[i][0] + outGrad[i][0] * outExp[i][1];

			outGrad[i][0] = real;
			outGrad[i][1] = imag;

		}
	 
		 // Transform back
		plan = fftw_plan_dft_c2r_1d(m_lConvolutionlLength, outGrad, adSlewRate, FFTW_ESTIMATE);
		fftw_execute(plan);
	
		// Correct amplitude
		double fCorrectionFactor = 1.0 / m_lConvolutionlLength;

		for (int i = 0; i < m_lConvolutionlLength; i++) {
			adSlewRate[i] *= fCorrectionFactor;
		}
	
		fftw_free(outGrad);
		fftw_free(outExp);
		fftw_destroy_plan(plan);

	 }

	 // Computes the ECC on each axis
	 void DSP::computeECC(){

		 static const char *ptModule = { "DSP::computeECC(): " };
		 	
		 // Compute ECC on all axes
		 if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {PRINT("Computing eddy current compensation for X axis...  \n"); mexEvalString("drawnow"); }
		 computeECC(m_afMULTI_PURPOSEX, m_afExponentialX);
	 
		 if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) { PRINT("Computing eddy current compensation for Y axis...  \n"); mexEvalString("drawnow"); }
		 computeECC(m_afMULTI_PURPOSEY, m_afExponentialY);

		 if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) { PRINT("Computing eddy current compensation for Z axis...  \n"); mexEvalString("drawnow"); }
		 computeECC(m_afMULTI_PURPOSEZ, m_afExponentialZ);
	 	 
	 }
 

	 //*******************************************************************
	 // CALCUALTEDERIVATIVE()                                             
	 //******************************************************************* 

	 // Calculates the derivatives of the provided data
	 void DSP::calculateDerivative(double *adData) {

			// Calculate derivative and destroy input
			for (long t = 0; t < m_lConvolutionlLength - 1; t++) {
				adData[t] = (adData[t + 1] - adData[t]) / sfGRT*1.0e-3; // Convert to T/m/s
			}
			adData[m_lConvolutionlLength - 1] = 0;
	 
	 }

	 // Calculates the derivatives for the gradient on each axis
	 void DSP::calculateDerivative() {
	 
			 if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) { PRINT("Computing derivative of GX...  \n"); mexEvalString("drawnow"); }
			 calculateDerivative(m_afMULTI_PURPOSEX);

			 if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) { PRINT("Computing derivative of GY...  \n"); mexEvalString("drawnow"); }
			 calculateDerivative(m_afMULTI_PURPOSEY);

			 if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) { PRINT("Computing derivative of GZ...  \n"); mexEvalString("drawnow"); }
			 calculateDerivative(m_afMULTI_PURPOSEZ);
		 
	 }

	 // Calculates the derivatives for the gradient on each axis
	 void DSP::calculateIntegral() {

		 if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) { PRINT("Computing integral of GX...  \n"); mexEvalString("drawnow"); }
		 calculateIntegral(m_afMULTI_PURPOSEX, m_afMULTI_PURPOSEX);

		 if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) { PRINT("Computing integral of GY...  \n"); mexEvalString("drawnow"); }
		 calculateIntegral(m_afMULTI_PURPOSEY, m_afMULTI_PURPOSEY);

		 if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) { PRINT("Computing integral of GZ...  \n"); mexEvalString("drawnow"); }
		 calculateIntegral(m_afMULTI_PURPOSEZ, m_afMULTI_PURPOSEZ);
	 
	 }


	 void DSP::calculateIntegral(double *adData, double *adIntegral, bool bNullAtTXCenter) {
		 

		 static const char *ptModule = { "DSP::calculateIntegral(): " };

		// Calculate integral and keep input
		adIntegral[0] = 0.0;
		
		// Find index of first TX center
		long cTXPulse = 0;
		long  x0 = 0;
		
		if (m_lTXEvents ==0  && bNullAtTXCenter){
			cout << "Warning: No TX events were detected (not implemented for version prior to VE). Hence the integral cannot be nulled at the TX center." << endl;
		}
						
		if (m_lTXEvents > 0){
			x0 = (long)ceil((m_adTXCenterTimes[cTXPulse] - sfGRT) / sfGRT);
		}		

		// Loop over array
		for (long t = 1; t < m_lConvolutionlLength; t++) {

			adIntegral[t] = adIntegral[t - 1] + adData[t] * sfGRT;
			
			// Null integral at center of RF pulse
			if (bNullAtTXCenter) {
				if (t == x0) {
					adIntegral[t] = 0;

					// Get next RF pulse center
					cTXPulse++;
					if (cTXPulse < m_lTXEvents) {
						x0 = (long)ceil((m_adTXCenterTimes[cTXPulse] - sfGRT) / sfGRT);
					}
				}
			}
		}

		return;
	 }


	 //*******************************************************************
	 // INTERPOLATEPHASE()                                                
	 //******************************************************************* 

	 // Interpolate the data at the positions of the readout samples
	 void DSP::interpolateData(double *afData, double *afDataInterp) {


		static const char *ptModule = { "DSP::interpolateData(): " };


		 for (long t = 0; t < m_lRXSampleLength; t++) {

			 // Requested sampling point
			 // Convert to raster units
			 double x = (m_adRXTimes[t]-sfGRT)/sfGRT; // First time sample is at 1 GRT, the first sample index is 0
		 		
			 // Find closest sampling raster time
			 long x0 = (long)floor( x );
			 long x1 = (long)ceil ( x );
		 
			 // Check indices
			 if(x0<0)
				 mexErrMsgOriginAndTxt(ptModule,"Invalid array index.\n");

			 if (x1 >= m_lConvolutionlLength)
				 x1 = m_lConvolutionlLength - 1;

			 // Get values at gradient raster times
			 double y0,y1;
		 		
			 y0 = afData[x0];
			 y1 = afData[x1];

		 		 
			 // Interpolate
			 if (x0 == x1) {
				// Sampling point is already on raster time
				 afDataInterp[t] = y0;
			 }
			 else {
				 // Perform linear interpolation
				 afDataInterp[t] = y0 + (x - (double)x0)*(y1 - y0) / ((double)(x1 - x0));
			 }

		 }


	 }

	 //*******************************************************************
	 // ALLOCATEMEMORY()                                                  
	 //******************************************************************* 


	// Allocates memory. In debug mode all arrays are handed over to MATLAB at the end of 
	// the calculation
	#ifdef MATLAB_MEX_FILE 
		void DSP::allocateMemory(mxArray *plhs[]) {

	 
		//---------------------------------------------------------------------
		// RX Samples                                                      
		//---------------------------------------------------------------------
		plhs[0] = mxCreateNumericMatrix(m_lRXEvents, 1, mxINT32_CLASS, mxREAL);
		m_uiRXEventLength = (uint32_t*)mxGetPr(plhs[0]);

		//---------------------------------------------------------------------
		// RX Times                                                         
		//---------------------------------------------------------------------
		plhs[1] = mxCreateNumericMatrix(m_lRXSampleLength, 1, mxDOUBLE_CLASS, mxREAL);
		m_adRXTimes = (double*)mxGetPr(plhs[1]);

		//---------------------------------------------------------------------
		// TX Times                                                         
		//---------------------------------------------------------------------
		plhs[2] = mxCreateNumericMatrix(m_lTXEvents, 1, mxDOUBLE_CLASS, mxREAL);
		m_adTXCenterTimes = (double*)mxGetPr(plhs[2]);

		//---------------------------------------------------------------------
		// Trigger Times                                                   
		//---------------------------------------------------------------------
		plhs[3] = mxCreateNumericMatrix(m_lTrigEvents, 1, mxDOUBLE_CLASS, mxREAL);
		m_adTrigTimes = (double*)mxGetPr(plhs[3]);

	 	 
		// The length of the linear convolution of two vectors of length,
		// M and L is M + L - 1, so we will extend our two vectors to that length before
		// computing the circular convolution using the DFT.
		m_lConvolutionlLength = m_lGradientShapeLength + m_lExponentialLength - 1;

		if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {
			PRINT("Gradient shape length = %d.\n", m_lGradientShapeLength);
			PRINT("Convolution length = %d.\n", m_lConvolutionlLength);
		}

		double dMemRequ = (double)m_lRXEvents * sizeof(uint32_t) + (double)(m_lRXSampleLength + m_lTXEvents + m_lTrigEvents) * sizeof(double);

		if (m_lOutputMode == OutputMode::FULL) {
			// Data + exponentials
			dMemRequ += (double)6.0 * m_lConvolutionlLength * sizeof(double);
		}
		else {
			// Data + exponentials + interpolated data
			dMemRequ += (double)6.0 * m_lConvolutionlLength * sizeof(double);
			dMemRequ += (double)3.0 * m_lRXSampleLength * sizeof(double);
		}


		if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {
			PRINT("Allocating %.2f MB \n\n", 1.0 * dMemRequ / 1e6);
			mexEvalString("drawnow");
		}
	
		//---------------------------------------------------------------------
		// Gradients                                                        
		//---------------------------------------------------------------------
		if (m_lOutputMode == OutputMode::FULL) {

			// Return the full data

			plhs[4] = mxCreateNumericMatrix(1, m_lConvolutionlLength, mxDOUBLE_CLASS, mxREAL);
			plhs[5] = mxCreateNumericMatrix(1, m_lConvolutionlLength, mxDOUBLE_CLASS, mxREAL);
			plhs[6] = mxCreateNumericMatrix(1, m_lConvolutionlLength, mxDOUBLE_CLASS, mxREAL);

			// Assign
			m_afMULTI_PURPOSEX = (double*)mxGetPr(plhs[4]);
			m_afMULTI_PURPOSEY = (double*)mxGetPr(plhs[5]);
			m_afMULTI_PURPOSEZ = (double*)mxGetPr(plhs[6]);

		}
		else {

			// Return the interpolated data

			plhs[4] = mxCreateNumericMatrix(1, m_lRXSampleLength, mxDOUBLE_CLASS, mxREAL);
			plhs[5] = mxCreateNumericMatrix(1, m_lRXSampleLength, mxDOUBLE_CLASS, mxREAL);
			plhs[6] = mxCreateNumericMatrix(1, m_lRXSampleLength, mxDOUBLE_CLASS, mxREAL);

			// Assign
			m_afMULTI_PURPOSE_INTERPX = (double*)mxGetPr(plhs[4]);
			m_afMULTI_PURPOSE_INTERPY = (double*)mxGetPr(plhs[5]);
			m_afMULTI_PURPOSE_INTERPZ = (double*)mxGetPr(plhs[6]);

		 
			// Allocate memory for the full data
			if (m_afMULTI_PURPOSEX != NULL) { delete[] m_afMULTI_PURPOSEX; m_afMULTI_PURPOSEX = NULL; }
			if (m_afMULTI_PURPOSEY != NULL) { delete[] m_afMULTI_PURPOSEY; m_afMULTI_PURPOSEY = NULL; }
			if (m_afMULTI_PURPOSEZ != NULL) { delete[] m_afMULTI_PURPOSEZ; m_afMULTI_PURPOSEZ = NULL; }

			m_afMULTI_PURPOSEX = new double[m_lConvolutionlLength];
			m_afMULTI_PURPOSEY = new double[m_lConvolutionlLength];
			m_afMULTI_PURPOSEZ = new double[m_lConvolutionlLength];
		 
			memset(m_afMULTI_PURPOSEX, 0, m_lConvolutionlLength * sizeof(double));
			memset(m_afMULTI_PURPOSEY, 0, m_lConvolutionlLength * sizeof(double));
			memset(m_afMULTI_PURPOSEZ, 0, m_lConvolutionlLength * sizeof(double));													 
		}

		//---------------------------------------------------------------------
		// Exponentials                                                      
		//---------------------------------------------------------------------
		if (m_afExponentialX != NULL) { delete[] m_afExponentialX; m_afExponentialX = NULL; }
		if (m_afExponentialY != NULL) { delete[] m_afExponentialY; m_afExponentialY = NULL; }
		if (m_afExponentialZ != NULL) { delete[] m_afExponentialZ; m_afExponentialZ = NULL; }

		m_afExponentialX = new double[m_lConvolutionlLength];
		m_afExponentialY = new double[m_lConvolutionlLength];
		m_afExponentialZ = new double[m_lConvolutionlLength];

		memset(m_afExponentialX, 0, m_lConvolutionlLength * sizeof(double));
		memset(m_afExponentialY, 0, m_lConvolutionlLength * sizeof(double));
		memset(m_afExponentialZ, 0, m_lConvolutionlLength * sizeof(double));															  
	 
		}

	#else

		void DSP::allocateMemory() {

			// The length of the linear convolution of two vectors of length,
			// M and L is M + L - 1, so we will extend our two vectors to that length before
			// computing the circular convolution using the DFT.
			m_lConvolutionlLength = m_lGradientShapeLength + m_lExponentialLength - 1;

			if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {
				PRINT("Gradient shape length = %ld.\n", m_lGradientShapeLength);
				PRINT("Convolution length = %ld.\n", m_lConvolutionlLength);
			}

			double dMemRequ = (double)m_lRXEvents * sizeof(uint32_t) + (double)(m_lRXSampleLength + m_lTXEvents + m_lTrigEvents) * sizeof(double);

			if (m_lOutputMode == OutputMode::FULL) {
				// Data + exponentials
				dMemRequ += (double)6.0*m_lConvolutionlLength * sizeof(double);
			}
			else {
				// Data + exponentials + interpolated data
				dMemRequ += (double)6.0*m_lConvolutionlLength * sizeof(double);
				dMemRequ += (double)3.0*m_lRXSampleLength * sizeof(double);
			}


			if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {
				PRINT("Allocating %.2f MB \n\n", 1.0 * dMemRequ / 1.0e6);
				mexEvalString("drawnow");
			}


			//---------------------------------------------------------------------
			// Gradients                                                        
			//---------------------------------------------------------------------
			if (m_afMULTI_PURPOSEX != NULL) { delete[] m_afMULTI_PURPOSEX; m_afMULTI_PURPOSEX = NULL; }
			if (m_afMULTI_PURPOSEY != NULL) { delete[] m_afMULTI_PURPOSEY; m_afMULTI_PURPOSEY = NULL; }
			if (m_afMULTI_PURPOSEZ != NULL) { delete[] m_afMULTI_PURPOSEZ; m_afMULTI_PURPOSEZ = NULL; }

			m_afMULTI_PURPOSEX = new double[m_lConvolutionlLength];
			m_afMULTI_PURPOSEY = new double[m_lConvolutionlLength];
			m_afMULTI_PURPOSEZ = new double[m_lConvolutionlLength];

			memset(m_afMULTI_PURPOSEX, 0, m_lConvolutionlLength * sizeof(double));
			memset(m_afMULTI_PURPOSEY, 0, m_lConvolutionlLength * sizeof(double));
			memset(m_afMULTI_PURPOSEZ, 0, m_lConvolutionlLength * sizeof(double));															 

			if (m_afMULTI_PURPOSE_INTERPX != NULL) { delete[] m_afMULTI_PURPOSE_INTERPX; m_afMULTI_PURPOSE_INTERPX = NULL; }
			if (m_afMULTI_PURPOSE_INTERPY != NULL) { delete[] m_afMULTI_PURPOSE_INTERPY; m_afMULTI_PURPOSE_INTERPY = NULL; }
			if (m_afMULTI_PURPOSE_INTERPZ != NULL) { delete[] m_afMULTI_PURPOSE_INTERPZ; m_afMULTI_PURPOSE_INTERPZ = NULL; }

			m_afMULTI_PURPOSE_INTERPX = new double[m_lRXSampleLength];
			m_afMULTI_PURPOSE_INTERPY = new double[m_lRXSampleLength];
			m_afMULTI_PURPOSE_INTERPZ = new double[m_lRXSampleLength];

			memset(m_afMULTI_PURPOSE_INTERPX, 0, m_lRXSampleLength * sizeof(double));
			memset(m_afMULTI_PURPOSE_INTERPY, 0, m_lRXSampleLength * sizeof(double));
			memset(m_afMULTI_PURPOSE_INTERPZ, 0, m_lRXSampleLength * sizeof(double));																
			//---------------------------------------------------------------------
			// Exponentials                                                      
			//---------------------------------------------------------------------
			if (m_afExponentialX != NULL) { delete[] m_afExponentialX; m_afExponentialX = NULL; }
			if (m_afExponentialY != NULL) { delete[] m_afExponentialY; m_afExponentialY = NULL; }
			if (m_afExponentialZ != NULL) { delete[] m_afExponentialZ; m_afExponentialZ = NULL; }

			m_afExponentialX = new double[m_lConvolutionlLength];
			m_afExponentialY = new double[m_lConvolutionlLength];
			m_afExponentialZ = new double[m_lConvolutionlLength];

			memset(m_afExponentialX, 0, m_lConvolutionlLength * sizeof(double));
			memset(m_afExponentialY, 0, m_lConvolutionlLength * sizeof(double));
			memset(m_afExponentialZ, 0, m_lConvolutionlLength * sizeof(double));														   
			//---------------------------------------------------------------------
			// RX Times                                                         
			//---------------------------------------------------------------------
			if (m_adRXTimes != NULL) { delete[] m_adRXTimes; m_adRXTimes = NULL; }
			m_adRXTimes = new double[m_lRXSampleLength];

			if (m_uiRXEventLength != NULL) { delete[] m_uiRXEventLength; m_uiRXEventLength = NULL; }
			m_uiRXEventLength = new uint32_t[m_lRXEvents];

			//---------------------------------------------------------------------
			// TX Times                                                         
			//---------------------------------------------------------------------
			if (m_adTXCenterTimes != NULL) { delete[] m_adTXCenterTimes; m_adTXCenterTimes = NULL; }
			m_adTXCenterTimes = new double[m_lTXEvents];

			//---------------------------------------------------------------------
			// Trigger Times                                                         
			//---------------------------------------------------------------------
			if (m_adTrigTimes != NULL) { delete[] m_adTrigTimes; m_adTrigTimes = NULL; }
			m_adTrigTimes = new double[m_lTrigEvents];

				
		}

	#endif


	//*******************************************************************
	// OPENFILE()                                                        
	//******************************************************************* 

	// Open the text file
	void DSP::openFile() {

		static const char *ptModule = { "DSP::openFile(): " };
		string fullpath = m_fpXMLFile.path + m_fpXMLFile.name + m_fpXMLFile.ext;

		if (fullpath.size()==0){
			mexErrMsgOriginAndTxt(ptModule,"File name not set.\n");
		}

		if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {
			PRINT("Loading from file: %s\n", fullpath.c_str());
		}	 
	}

	//*******************************************************************
	// READSHAPES()                                                      
	//******************************************************************* 

	void DSP::readGCShapes() {

		static const char *ptModule = { "DSP::readGCShapes(): " };


		string sfileName;
		double dHoriDelta = 0.0;
		double dVertFactor = 0.0;
		long lLengthGY = 0;
		long lLengthGZ = 0;

		// GX
		sfileName = m_sDSVFolderPath + "/" + m_sDSVFileNamePrefix + "GRX.dsv";
		{
			ifstream in(sfileName.c_str());

			if (!in) {
			cout << "Cannot open input file: " << sfileName << endl;
			mexErrMsgOriginAndTxt(ptModule, "Check provided folder. DSV file needs to be named: GRX.dsv, GRY.dsv, GRZ.dsv, and INF.dsv (+ an optional prefix). \n");
			}

			// Get the gradient shape
			this->getDSV(in, m_afMULTI_PURPOSEX);
			in.close();
		}

		// GY
		sfileName = m_sDSVFolderPath + "/" + m_sDSVFileNamePrefix + "GRY.dsv";
		{
			ifstream in(sfileName.c_str());

			if (!in) {
			cout << "Cannot open input file: " << sfileName << endl;
			mexErrMsgOriginAndTxt(ptModule, "Check provided folder. DSV file needs to be named: GRX.dsv, GRY.dsv, GRZ.dsv, and INF.dsv (+ an optional prefix). \n");
			}

			// Get the gradient shape
			this->getDSV(in, m_afMULTI_PURPOSEY);
			in.close();
		}

		// GZ
		sfileName = m_sDSVFolderPath + "/" + m_sDSVFileNamePrefix + "GRZ.dsv";
		{
			ifstream in(sfileName.c_str());

			if (!in) {
			cout << "Cannot open input file: " << sfileName << endl;
			mexErrMsgOriginAndTxt(ptModule, "Check provided folder. DSV file needs to be named: GRX.dsv, GRY.dsv, GRZ.dsv, and INF.dsv (+ an optional prefix). \n");
			}

			// Get the gradient shape
			this->getDSV(in, m_afMULTI_PURPOSEZ);
			in.close();
		}


		//---------------------------------------------------------------------
		// Get file Parameters
		//---------------------------------------------------------------------
		sfileName = m_sDSVFolderPath + "/" + m_sDSVFileNamePrefix + "INF.dsv";
		{
			ifstream in(sfileName.c_str());

			if (!in) {
			cout << "Cannot open input file: " << sfileName << endl;
			mexErrMsgOriginAndTxt(ptModule, "Check provided folder. DSV file needs to be named: GRX.dsv, GRY.dsv, GRZ.dsv, and INF.dsv (+ an optional prefix). \n");
			}

			this->getInfo(in);
		}

	}


	 void DSP::applyPhaseModulation(complex<float> *data, unsigned short numberOfSamples, unsigned int ulScanCounter) {

		// Scan counter should start at zero

		static const char *ptModule = { "DSP::applyPhaseModulation(): " };
			     
		 if (m_bEccCompensationAvailable) {
         
			 // Checks
			 if (ulScanCounter == m_lRXEvents) {
				 PRINT(ptModule, "ACQEND will not be processed.\n");
				 // ACQEND is not processed by DSP.cpp because it appears in the XML file after the GC-<Halt> instruction.
				 return;
			 }
			 else if (ulScanCounter > m_lRXEvents - 1) {
				 mexErrMsgOriginAndTxt(ptModule, "Scan counter exceeds number of stored ADC events.\n");
			 }
			 else if(ulScanCounter < 0) {
				 mexErrMsgOriginAndTxt(ptModule, "Scan counter cannot be smaller than zero.\n");
			 }
         
			 // Expected number of samples for current readout
			 long lExpectedNumberOfSamples = 0;
			 if (ulScanCounter == 0) {
				 lExpectedNumberOfSamples = m_uiRXEventLength[ulScanCounter];
			 }
			 else {
				 lExpectedNumberOfSamples = m_uiRXEventLength[ulScanCounter] - m_uiRXEventLength[ulScanCounter-1];
			 }
			 
			 if (numberOfSamples != lExpectedNumberOfSamples) {
				 static char tLine[1000];
				 tLine[0] = '\0';
				 printf(tLine, "Number of RX samples (%d) for current scan counter (%d) is not equal to expected number (%d).\n", numberOfSamples, ulScanCounter + 1, lExpectedNumberOfSamples);  //Scan counter (Siemens raw data) starts orginally at 1. Therefor we add here 1 again.
				 mexErrMsgOriginAndTxt(ptModule, tLine);
			 }
			  
			 long lStartIndex = 0;
			 if (ulScanCounter > 0)
				 lStartIndex = m_uiRXEventLength[ulScanCounter-1];

			 // Apply phase modulation	
			 for (long i = 0; i < numberOfSamples; i++) {
				   data[i] =  data[i] * polar((float)1.0, (float)(-1.0*m_afMULTI_PURPOSE_INTERPX[lStartIndex + i]));
				   //data[i] = complex<float>(m_afMULTI_PURPOSE_INTERPX[lStartIndex + i],0); // For  debug only
			 }
        
		 }
		 
	 }



	 double DSP::getADCStartTime(unsigned int ulScanCounter) {
		 static const char *ptModule = { "DSP::getADCStartTime(): " };

		 throw "Not implemented.";
		 
	 }

	 void DSP::getInfoParams(ifstream &in, long &lRXSampleLength, long &lRXEvents, long &lTXEvents) {

		 string line;

		 while (std::getline(in, line))
		 {	
	 
			// Extract TX events for VE
			if (line.find("| Rel.Time |") != std::string::npos) {
				
				// Read block info
				while (std::getline(in, line))
				{	
					
					if (line.find("#INFO-END") != std::string::npos) {
						break;
					}
										
					// Get positions of all '|' characters
					vector<size_t> positions; 

					size_t pos = line.find("|", 0);
					while(pos != string::npos)
					{
						positions.push_back(pos);
						pos = line.find("|",pos+1);
					}
					
					if (positions.size()>3){
						string txString = line.substr(positions.at(1)+1,positions.at(2)-positions.at(1));
						if (txString.find(":") != std::string::npos) {
							// Count TX events
							lTXEvents++;
						}	
					}
				}
								
				continue;
			}
			// VD
			else if (line.find("Rel.Time") != std::string::npos) {
				// Read block info
				while (std::getline(in, line))
				{

					if (line.find("#INFO-END") != std::string::npos) {
						break;
					}

					int s = line.size()>23 ? 23 : line.size();
					

					string st = line.substr(0, s - 1);

					// Get position of '/' characters - Separator between TXFlip and Duration
					vector<size_t> positions;

					size_t pos = st.find("/", 0);
					if (pos != string::npos) {
						// Count TX events
						lTXEvents++;
					}					

				}

				continue;
			}

			 
	 
	 			 
			 if (line.find("adc MeasHeader") != std::string::npos) {
				 lRXEvents++;
				 continue;
			 }
			 
			 
			 
			 if (line.find("adc MeasHeader") != std::string::npos) {
				 lRXEvents++;
				 continue;
			 }

			 if (line.find("ushSamplesInScan") != std::string::npos) {

				 char * pch;
				 pch = strtok((char *)line.c_str(), ":;");
				 pch = strtok(NULL, "=;");

				 lRXSampleLength += stod(pch);
			 }
		 }

		 // The last RX event can be disregarded. It has zero samples
		 lRXEvents--;

		 in.clear();
		 in.seekg(0);

	 };

	 

	 void DSP::getDSVParams(ifstream &in, long &lSamples, double &dHoriDelta, double &dVertFactor) {
		 
		 string line;
		 while (std::getline(in, line))
		 {

			 // Stop reading when [VALUES] array starts
			 std::size_t found = line.find("[VALUES]");
			 if (found != std::string::npos) {
				 break;
			 }
			 
			 char * pch;
			 pch = strtok((char *)line.c_str(), "=");

			 while (pch != NULL)
			 {
				 if (strcmp(pch, "SAMPLES") == 0) {

					 // Get value
					 pch = strtok(NULL, "=");

					 lSamples = stoi(pch);

				 }
				 else if (strcmp(pch, "HORIDELTA") == 0) {

					 // Get value
					 pch = strtok(NULL, "=");

					 dHoriDelta = stod(pch);

				 }
				 else if (strcmp(pch, "VERTFACTOR") == 0) {

					 // Get value
					 pch = strtok(NULL, "=");

					 dVertFactor = 1.0 / stod(pch);

				 }

				 pch = strtok(NULL, "=");
			 }

		 }
	 }


	 void DSP::getInfo(ifstream &in) {
		 static const char *ptModule = { "DSP::getInfo(): " };

		 //---------------------------------------------------------------------
		 // Read the values
		 //---------------------------------------------------------------------
		 vector<double> vDSV;
		 string line;
		 
		 long lStartTimeEventBlock;
		 
		 long lStartTime;
		 double dADCDuration;
		 long lADCSamples;

		 m_lCurrentRXSampleLength = 0;
		 m_lCurrentTXNumber = 0;

		 long cRXEvent = 0;
		 while (std::getline(in, line))
		 {
	
			// Event block
			if (line.find("EventBlock") != std::string::npos) {
				 char * pch;
				 pch = strtok((char *)line.c_str(), " ");
				 pch = strtok(NULL, " ");

				 // Start time
				 pch = strtok(NULL, " ");
				 lStartTimeEventBlock = stod(pch);
				 			 
				 
				// Read block info
				bool isVE = true;
				while (std::getline(in, line))
				{

					if (line.find("#INFO-END") != std::string::npos) {
						break;
					}

					if (line.find("| Rel.Time |") != std::string::npos) {
						isVE = true;
						continue;
					}
					else if (line.find("Rel.Time") != std::string::npos) {
						isVE = false;
						continue;
					}

					if (isVE) {
						// Get positions of all '|' characters
						vector<size_t> positions;

						size_t pos = line.find("|", 0);
						while (pos != string::npos)
						{
							positions.push_back(pos);
							pos = line.find("|", pos + 1);
						}

						if (positions.size() > 3) {


							string txString = line.substr(positions.at(1) + 1, positions.at(2) - positions.at(1) - 1);
							if (txString.find(":") != std::string::npos) {

								string flipDurString = txString.substr(txString.find(":") + 1, string::npos);


								char * pch;
								pch = strtok((char *)flipDurString.c_str(), "/");

								// Duration
								pch = strtok(NULL, "/");
								long lDur = stod(pch);


								// Relative time
								string relTimeString = line.substr(positions.at(0) + 1, positions.at(1) - positions.at(0) - 1);
								long lrelTime = stod(relTimeString);

								// Set TX center (assuming it is a symmetric pulse)						
								m_adTXCenterTimes[m_lCurrentTXNumber] = 1.0e-6*double(lStartTimeEventBlock + lrelTime + lDur / 2);
								m_lCurrentTXNumber++;

							}
						}
					} 
					else {// is VD

						int s = line.size() > 23 ? 23 : line.size();
						
						string flipDurString = line.substr(8, 12);

						// Get position of '/' characters - Separator between TXFlip and Duration
						vector<size_t> positions;

						size_t pos = flipDurString.find("/", 0);
						if (pos != string::npos) {
							char * pch;
							pch = strtok((char *)flipDurString.c_str(), "/");

							// Duration
							pch = strtok(NULL, "/");
							long lDur = stod(pch);

							// Relative time
							string relTimeString = line.substr(0, 8);
							long lrelTime = stod(relTimeString);

							// Set TX center (assuming it is a symmetric pulse)						
							m_adTXCenterTimes[m_lCurrentTXNumber] = 1.0e-6*double(lStartTimeEventBlock + lrelTime + lDur / 2);
							m_lCurrentTXNumber++;
							
						}
					}
				}
								
				continue;

			}
	
			// Measurement header
			if (line.find("MeasHeader") != std::string::npos) {
				 
				 char * pch;
				 pch = strtok((char *)line.c_str(), " ");
				 pch = strtok(NULL, " ");

				 // Start time
				 pch = strtok(NULL, " ");
				 lStartTime = stod(pch);

				 // ADC duration
				 pch = strtok(NULL, " ");
				 dADCDuration = stod(pch);				 
			
			}

			if (line.find("ushSamplesInScan") != std::string::npos) {

				char * pch;
				pch = strtok((char *)line.c_str(), ":;");
				pch = strtok(NULL, "=;");

				lADCSamples = stoi(pch);
				 
				if (lADCSamples > 0 && cRXEvent < m_lRXEvents) {

					double dDwellTime = dADCDuration / lADCSamples;

					for (long t = 0; t < lADCSamples; t++) {
						if (m_lCurrentRXSampleLength + t < m_lRXSampleLength) {
							m_adRXTimes[m_lCurrentRXSampleLength + t] = ((double)lStartTime + dDwellTime * (double)t)*1e-6; // Final unit is seconds
						}
					}

					// Update shape length
					m_lCurrentRXSampleLength += lADCSamples;

					// Set number of samples per ADC
					m_uiRXEventLength[cRXEvent] = m_lCurrentRXSampleLength;


					// Update RX event counter
					cRXEvent++;
				}

			}
		 }

	 }

	 bool is_number(const std::string& s)
	 {
		 char* end = 0;
		 double val = strtod(s.c_str(), &end);
		 return end != s.c_str() && val != HUGE_VAL;
	 }

	 bool DSP::getDSV(ifstream &in, double *y) {

		 static const char *ptModule = { "DSP::getDSV(): " };

		 //---------------------------------------------------------------------
		 // Get file name
		 //---------------------------------------------------------------------
		 if (!in) {
			 mexErrMsgOriginAndTxt(ptModule, "Cannot open input file stream. \n");
		 }

		 //---------------------------------------------------------------------
		 // Get file Parameters
		 //---------------------------------------------------------------------
		 long lSamples = 0;
		 double dHoriDelta = 0.0;
		 double dVertFactor = 0.0;

		 this->getDSVParams(in, lSamples, dHoriDelta, dVertFactor);

		 //---------------------------------------------------------------------
		 // Read the values
		 //---------------------------------------------------------------------
		 vector<double> vDSV;
		 string line;
		 while (std::getline(in, line))
		 {
			 if (is_number(line)) {
				 vDSV.push_back(stod(line));
			 }
			 else {
				 break;
			 }
		 }

		 //---------------------------------------------------------------------
		 // Compute shape
		 //---------------------------------------------------------------------
		 // Set length of DSV file
		 long lNdsv = vDSV.size();

		 // Init gradient shape
		 for (int n = 0; n < lSamples; n++) {
			 y[n] = 0;
		 }

		 // Init counters
		 long i = 0;
		 long n = 0;

		 // First sample 
		 y[0] = vDSV.at(0);
		 double prev = vDSV.at(0);

		 while (i < lNdsv - 1 && n < lSamples - 2)
		 {
			 i++;
			 y[n + 1] = y[n] + vDSV.at(i);

			 if (vDSV.at(i) != prev) {
				 n = n + 1;
				 prev = vDSV.at(i);
			 }
			 else {
				 long nrep = vDSV.at(i + 1);
				 if (nrep == 0) {
					 i = i + 2;
					 n = n + 1;
				 }
				 else {
					 n = n + 1;


					 for (int a = 0; a < nrep; a++) {
						 y[n + 1 + a] = y[n] + (a + 1)*prev;
					 }

					 n = n + nrep;
					 i = i + 2;
				 }

				 if (i < lNdsv) {
					 y[n + 1] = y[n] + vDSV.at(i);
					 n = n + 1;
					 prev = vDSV.at(i);
				 }
			 }

		 }

		 //---------------------------------------------------------------------
		 // Scaling
		 //---------------------------------------------------------------------
		 for (int n = 0; n < lSamples; n++) {
			 y[n] = y[n]*dVertFactor;
		 }

		 return true;
	 }


	//*******************************************************************
	// RUN()                                                             
	//******************************************************************* 

	// Read text file and simulate ECC

	#ifdef MATLAB_MEX_FILE
	 void DSP::run(mxArray *plhs[])
	#else
	 void DSP::run() 
	#endif

	{
		//---------------------------------------------------------------------
		// Print Info                                                         
		//---------------------------------------------------------------------
		if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {
			PRINT("\n#####################################\n");
			PRINT("#      Gradient reader (VD/VE)      #\n");
			PRINT("#              and                  #\n");
			PRINT("#         ECC calculator            #\n");
			PRINT("#       Version 2.0 (2020)          #\n");
			PRINT("#####################################\n\n");
		}


		//---------------------------------------------------------------------
		// Open file                                                         
		//---------------------------------------------------------------------
		// DSV folder has not been set / check if xml file is present
		if (m_sDSVFolderPath.size() == 0) {
			this->openFile();
		}

		//---------------------------------------------------------------------
		// Calculate memory requirement (set m_lGradientShapeLength)                              
		//---------------------------------------------------------------------
		this->calcMemoryRequirement();

		//---------------------------------------------------------------------
		// Determine largest decay constant (set m_lExponentialLength)                                              
		//---------------------------------------------------------------------
		this->determineLongestTimeConstant();

		//---------------------------------------------------------------------
		// Allocate memory                                                   
		//---------------------------------------------------------------------
		#ifdef MATLAB_MEX_FILE
			this->allocateMemory(plhs);
		#else
			this->allocateMemory();
		#endif

		//---------------------------------------------------------------------
		// Read shapes                                                 
		//---------------------------------------------------------------------
		this->readGCShapes();

		//---------------------------------------------------------------------
		// Calculate integral of gradient                                           
		//---------------------------------------------------------------------
		if (m_lDataType == DataType::KSPACE) {

			// Calculate integral of gradient
			this->calculateIntegral();

			for (long t = 0; t < m_lConvolutionlLength; t++) {
				m_afMULTI_PURPOSEX[t] *= 2.0*sfPI*sfGAMMA_1H*1000.0;
				m_afMULTI_PURPOSEY[t] *= 2.0*sfPI*sfGAMMA_1H*1000.0;
				m_afMULTI_PURPOSEZ[t] *= 2.0*sfPI*sfGAMMA_1H*1000.0;
			}

		}

		//---------------------------------------------------------------------
		// Calculate derivative of gradient                                           
		//---------------------------------------------------------------------
		if (m_lDataType == DataType::SLEWRATE || m_lDataType == DataType::EDDYCURRENT || m_lDataType == DataType::EDDYPHASE) {
			this->calculateDerivative();
		}


		//---------------------------------------------------------------------
		// Calculate B0 eddy currents                                 
		//---------------------------------------------------------------------
		if (m_lDataType == DataType::EDDYCURRENT || m_lDataType == DataType::EDDYPHASE) {

			//---------------------------------------------------------------------
			// Compute exponentials                                                 
			//---------------------------------------------------------------------
			this->computeExponentials();


			//---------------------------------------------------------------------
			// Compute convolution via DFT ( conv(x,y) = ifft(fft(x).*fft(y)) )
			//---------------------------------------------------------------------
			this->computeECC();
		}


		//---------------------------------------------------------------------
		// Calculate phase caused by B0 eddy currents                                 
		//---------------------------------------------------------------------
		if ( m_lDataType == DataType::EDDYPHASE ) {

			// Compute sum
			for (long t = 0; t < m_lConvolutionlLength; t++) {
				m_afMULTI_PURPOSEX[t] += m_afMULTI_PURPOSEY[t];
				m_afMULTI_PURPOSEX[t] += m_afMULTI_PURPOSEZ[t];
			}

			for (long t = 0; t < 10; t++) {
				double a = m_afMULTI_PURPOSEX[t];
			}

			//---------------------------------------------------------------------
			// Calculate phase by trapezoidal integration
			//---------------------------------------------------------------------
			// Do null phase at center of TX pulse
			calculateIntegral(m_afMULTI_PURPOSEX, m_afMULTI_PURPOSEX, true);

			for (long t = 0; t < m_lConvolutionlLength; t++) {
				m_afMULTI_PURPOSEX[t] *= 2.0*sfPI*sfGAMMA_1H;
			}
						
		}

		for (long t = 0; t < 10; t++) {
			double a = m_afMULTI_PURPOSEX[t];
		}


		//---------------------------------------------------------------------
		// Interpolate to RX events                      
		//---------------------------------------------------------------------
		if (m_lOutputMode == OutputMode::INTERPOLATED_TO_RX) {
						
			this->interpolateData(m_afMULTI_PURPOSEX, m_afMULTI_PURPOSE_INTERPX);
					   
			if (m_lDataType != DataType::EDDYPHASE) {
				this->interpolateData(m_afMULTI_PURPOSEY, m_afMULTI_PURPOSE_INTERPY);
				this->interpolateData(m_afMULTI_PURPOSEZ, m_afMULTI_PURPOSE_INTERPZ);
			

				if (m_lTXEvents > 1) {

					// Set RX times relative to RF center
					long txPulse = -1;

					double x0 = 0.0;

					// Loop over array
					for (long t = 0; t < m_lRXSampleLength; t++) {

						// Null integral at center of RF pulse
						if (txPulse + 1 < m_lTXEvents) {

							if ( m_adRXTimes[t] >= m_adTXCenterTimes[txPulse + 1] ) {

								txPulse++;
								x0 = m_adTXCenterTimes[txPulse];
							}

						}

						m_adRXTimes[t] -= x0;  

					}
				}

			}
		}

		if (m_lOutputMode == OutputMode::INTERPOLATED_TO_RX && m_lDataType == DataType::EDDYPHASE) {
			m_bEccCompensationAvailable = true;
		}	
	}

}
