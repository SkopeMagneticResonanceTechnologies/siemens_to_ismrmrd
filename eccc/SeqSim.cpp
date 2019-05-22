//	----------------------------------------------------------------------------------------
//	  Copyright (C) Skope AG 2018  All Rights Reserved.
//	----------------------------------------------------------------------------------------
//
//	  Project:	DSP READER
//	     File:	SeqSim.cpp
//	  Version:	1.1
//	   Author:	MC^2
//	     Date:	08.11.2018
//	     Lang:	C++
//
//Description:	This class can generate gradient/k-space/slew rate/B0 eddy currents/B0 eddy current phase 
//				data based on an XML file containing gradient instructions. The XML must be generated in 
//				IDEA with the function call defined below. The class is intended to be compiled as a dll
//				to be integrated later in the SiemensToMRD converter. A wrapper for Matlab (SeqSim_mex.cpp)
//				and a console application (SeqSimRun.cpp) are also provided. The compilation instructions
//				can be found within the source file. Additional inputs to the class are the name of the 
//				XML file, the type of data to be returned, the outputmode and (if eddy currents should be
//				calculated) the ECC coefficients that can be extracted from the Siemens raw data header.
//
//	 Options:   setDataType defines the type to be returned by the class: 
//					- gradient
//					- k-space
//					- slew rate
//					- B0 eddy currents
//					- B0 eddy current phase 	
//
//				setOutputMode:
//					FULL : return the full simulated data shape on the gradient raster time (GRT)
//					INTERPOLATED_TO_RX : return the interpolated data shape at the RX sample locations
//
//		Units:	
//				gradient [mT/m/s]			
//				slew rate [T/m/s]			
//				k-space [1/m]				
//				eddy currents [uT]			
//				eddy current phase [rad]	
//				TX/RX/Trig times [s]
//
// Important notes:
//		-	The first time sample is at 10e-6 s and not at zero! This is important if you want to determine the k-space position at a specific 
//			RX time for example. In Matlab you would use: kspaceInterp = interp1([1:length(kspaceOnGRT)]*10e-6, kspaceOnGRT, RXTimes)
//
//		-	If the outputMode is 
//									* 0: the RX time will be calculated relative to the start of the sequence
//									* 1: the RX time will be calculated relative to the center of the last RF pulse
//
//	----------------------------------------------------------------------------------------
// To create an xml file in IDEA call:
//
//    	MriSeqSimRunSeq -p prot.pro DSP.xml
//		MriSeqSimRunSeq -p prot.pro -c PAC_1 DSP.xml  (to use custom coil environment'%MEASDAT%/MeasCoil_PAC_1.dat')
//
//	----------------------------------------------------------------------------------------


//#include "stdafx.h"

//#define MATLAB_MEX_FILE 
#include "SeqSim.h"		// DSP header
#include "fftw3.h"		// The Fastest Fourier Transform in the West
#include <algorithm>    // std::min_element, std::max_element
//#include <string.h> 
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
	// HELPER CLASSES                                                    
	//******************************************************************* 

	// Axis properties in event block


	// Readout properties in event block
	Readout::Readout() {
		bIsValidScan = false;
		fDwellTime	 = 0.0;
		lSamples	 = 0;
	}

	Axis::Axis() {
		Increment_PE.fOffset = 0.0;
		Increment_PE.fValue = 0.0;

		Increment_RO.fOffset = 0.0;
		Increment_RO.fValue = 0.0;

		Increment_SL.fOffset = 0.0;
		Increment_SL.fValue = 0.0;

		Shape_PE.ID = -1;

		Shape_RO.ID = -1;

		Shape_SL.ID = -1;
	}

	// TX properties in event block
	TX::TX() {
		TxType txType    = Undefined;
		double fAsymmetry = 0.5;
		double fDwellTime = 0.0;
		long lSamples    = 0;
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
			// throw everthing into the trash bin

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
			
			if (outputFileX.is_open() && outputFileY.is_open() && outputFileZ.is_open())
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

		static const char *ptModule = { "DSP::calcMemoryRequirement()" };
	
		//---------------------------------------------------------------------
		// Initialise
		//---------------------------------------------------------------------

		pugi::xml_node root;
	
		// Gradient shape length
		m_lGradientShapeLength	= 0;

		// Number of readout events
		m_lRXEvents				= 0;

		// Number of transmit events
		m_lTXEvents             = 0;

		// Number of trigger events
		m_lTrigEvents			= 0;

		// Number of ADC samples (must be the same for all imaging scans)
		long ADCSamples			= -1;

		//---------------------------------------------------------------------
		// Load document
		//---------------------------------------------------------------------
		string fullpath = m_fpXMLFile.path + m_fpXMLFile.name + m_fpXMLFile.ext;

		if (!m_doc.load_file(fullpath.c_str())) {
			mexErrMsgOriginAndTxt(ptModule, "Error loading xml file.\n");
		}

		//---------------------------------------------------------------------
		// Read event blocks
		//---------------------------------------------------------------------
		START:

		// Delete full path for now. If an additional xml file need to be loaded
		// then it will be set during the next for-loop
		fullpath.clear();
	
		// Load root
		root = m_doc.child("NUMARIS4_DSP_SIMULATION");
	
		// Loop through children of NUMARIS4_DSP_SIMULATION
		for (pugi::xml_node node = root.first_child(); node; node = node.next_sibling())
		{
		
			// Convert name to string
			std::string sNodeName(node.name());

			//---------------------------------------------------------------------
			// Event block
			//---------------------------------------------------------------------
			if (sNodeName.compare("EventBlock") == 0) {

				//---------------------------------------------------------------------
				// GC
				//---------------------------------------------------------------------
				pugi::xml_node nodeGC = node.child("GC");

				// Loop through children of GC
				for (pugi::xml_node child = nodeGC.first_child(); child; child = child.next_sibling())
				{

					// Convert name to string
					std::string sChildName(child.name());

					if (sChildName.compare("Control") == 0) {
					
						// Get attribute 'Ticks'
						if (pugi::xml_attribute attr = child.attribute("Ticks")) // attribute really exists
						{
							// Get number of samples
							m_lGradientShapeLength += attr.as_int();
						}

						// Check ig control has a child called 'Trigger'
						for (pugi::xml_node childControl = child.first_child(); childControl; childControl = childControl.next_sibling())
						{
							// Convert name to string
							std::string sChildName2(childControl.name());

							if (sChildName2.compare("Trigger") == 0) {

								// Get text 
								string sTriggerType = childControl.text().as_string();

								if (sTriggerType.compare("TX") == 0) {
									m_lTXEvents++;
								}
							}

							if (sChildName2.compare("Sync") == 0) {

								// Only consider EXT_TRIGGER
								if (pugi::xml_attribute attr = childControl.attribute("Ext")) 
								{
									// One could also differentiate between EXTR0 and EXTR1 if needed
									// by getting the attribute value
									m_lTrigEvents++;
								}

							}

						}

					}
					else if (sChildName.compare("Halt") == 0) {
						if ((m_uiVerbose & DISPLAY_ADVANCED) == DISPLAY_ADVANCED) {
							PRINT("%s: Halt instruction received.\n", ptModule);
						}
						goto FINISH;
					}
				}

				//---------------------------------------------------------------------
				// RX
				//---------------------------------------------------------------------
				pugi::xml_node nodeRX = node.child("RX");

				// Loop through children of RX
				for (pugi::xml_node child = nodeRX.first_child(); child; child = child.next_sibling())
				{

					// Convert name to string
					std::string sChildName(child.name());

					if (sChildName.compare("Readout") == 0) {

						// Get attribute 'NumberOfPoints'
						if (pugi::xml_attribute attr = child.attribute("NumberOfPoints")) // attribute really exists
						{

							pugi::xml_node header = child.child("Info").child("Header");
							pugi::xml_attribute evalMask = header.attribute("aulEvalInfoMask0");
						
							unsigned long ulEvalInfoMask = evalMask.as_uint();
						
							// Only consider imaging scans
							if (
							//	   (ulEvalInfoMask & MdhBitField::MDH_FLAG_ONLINE) == MdhBitField::MDH_FLAG_ONLINE &&
							//	 (ulEvalInfoMask & MdhBitField::MDH_FLAG_PHASCOR) != MdhBitField::MDH_FLAG_PHASCOR &&
							//	 (ulEvalInfoMask & MdhBitField::MDH_FLAG_NOISEADJSCAN) != MdhBitField::MDH_FLAG_NOISEADJSCAN &&
								 (ulEvalInfoMask & MdhBitField::MDH_FLAG_ACQEND) != MdhBitField::MDH_FLAG_ACQEND) {

								// Get number of samples
								long Samples = attr.as_int();
							
								// Set ADC samples if was not set yet
								if (ADCSamples == -1)
									ADCSamples = Samples;

								//if (ADCSamples != Samples)
								//	mexErrMsgOriginAndTxt(ptModule,"Number of ADC samples must be the same for all acqusitions. (Possible reason: Seperate Ref scans have an oversampling factor equal to 1.)\n");

								// Add samples to counter
								m_lRXSampleLength += Samples;
								m_lRXEvents++;

							} // Flag check
						} // Attribute Number of points
					
					} // Readout

				} // Loop over children
		
			} // Loop over Event blocks
			else if (sNodeName.compare("Continue") == 0) {

				// Get attribute 'File'
				if (pugi::xml_attribute attr = node.attribute("File")) // attribute really exists
				{
					// Get next file name
					pugi::char_t *cFilename = (char*)attr.value();

					fullpath = m_fpXMLFile.path + string(cFilename); // The extension is already included in cFilename

				}
				else {
					mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'File' in Continue.\n");
				}
			
			}
		} // root
	
		FINISH:
	
		//---------------------------------------------------------------------
		// Load next xml file is necessary
		//---------------------------------------------------------------------
		if (fullpath.size()>0) {

			if (!m_doc.load_file(fullpath.c_str())) {
				mexErrMsgOriginAndTxt(ptModule, "Error loading next xml file: \n", fullpath.c_str());
			}
			else {
				//PRINT("%s: Loading next XML file: %s.\n", ptModule, fullpath.c_str());
				// Jump back and parse the next file
				goto START;
			}
		}
	
		//---------------------------------------------------------------------
		// Check if the needed data is present                              
		//--------------------------------------------------------------------- 
		if (m_lGradientShapeLength == 0)
			mexErrMsgOriginAndTxt(ptModule,"Could not find any gradient DSP information!\n");

		if (m_lRXEvents == 0)
			mexErrMsgOriginAndTxt(ptModule,"Could not find any receiver DSP information!\n");

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
	// RUNINSTRUCTIONS()                                                 
	//******************************************************************* 

	// Runs all the gradient and receiver instructions specified in the text file
	 void DSP::runInstructions(){

		 static const char *ptModule = { "DSP::runInstructions()" };
	 

		 m_lCurrentRXNumber   = 0;
		 m_lCurrentTXNumber   = 0;
		 m_lCurrentTrigNumber = 0;

		 //---------------------------------------------------------------------
		 // Initialise
		 //---------------------------------------------------------------------
		 // Storage for filename
		 pugi::xml_node root;
	 
		 //---------------------------------------------------------------------
		 // Load document
		 //---------------------------------------------------------------------
		 string fullpath = m_fpXMLFile.path + m_fpXMLFile.name + m_fpXMLFile.ext;

		 if (!m_doc.load_file(fullpath.c_str())) {
			 mexErrMsgOriginAndTxt(ptModule, "Error loading xml file.");
		 }

		 //---------------------------------------------------------------------
		 // Read event blocks
		 //---------------------------------------------------------------------
		 START:

		 // Load root
		 root = m_doc.child("NUMARIS4_DSP_SIMULATION");

		 // Set fullpath to zero for now. If an additional XML file needs to be loaded
		 // then it will be set during the next for-loop
		 fullpath.clear();

		 // Loop through children of NUMARIS4_DSP_SIMULATION
		 for (pugi::xml_node node = root.first_child(); node; node = node.next_sibling())
		 {

			 // Convert name to string
			 std::string sNodeName(node.name());

			 //---------------------------------------------------------------------
			 // Event block
			 //---------------------------------------------------------------------
			 if (sNodeName.compare("EventBlock") == 0) {

			 
				 //---------------------------------------------------------------------
				 // Initialise
				 //---------------------------------------------------------------------
				 vector<Readout> vReadout;
				 vector<TX> vTx;
			 			 
				 //---------------------------------------------------------------------
				 // RX
				 //---------------------------------------------------------------------
				 pugi::xml_node nodeRX = node.child("RX");

				 // Loop through children of RX
				 for (pugi::xml_node child = nodeRX.first_child(); child; child = child.next_sibling())
				 {
			
					 // Convert name to string
					 std::string sChildName(child.name());

					 if (sChildName.compare("Readout") == 0) {

						 // Create local readout object
						 Readout readout;
					 
						 // Get attribute 'NumberOfPoints'
						 if (pugi::xml_attribute attr = child.attribute("NumberOfPoints")) {
							 readout.lSamples = attr.as_int();
						 } 
						 else {
							 mexErrMsgOriginAndTxt(ptModule,"Could not find attribute 'NumberOfPoints' in Readout.");
						 }

						 // Get attribute dwelltime
						 if (pugi::xml_attribute attr = child.attribute("Decimation")) {
							 readout.fDwellTime = attr.as_double()/10;
						 }
						 else {
							 mexErrMsgOriginAndTxt(ptModule,"Could not find attribute 'Decimation' in Readout.");
						 }
					 

						 // Retrieve evaluation mask
						 pugi::xml_node header = child.child("Info").child("Header");
						 pugi::xml_attribute evalMask = header.attribute("aulEvalInfoMask0");

						 unsigned long ulEvalInfoMask = evalMask.as_uint();


						 // Only consider imaging scans
						 if (
						 //	  (ulEvalInfoMask & MdhBitField::MDH_FLAG_ONLINE) == MdhBitField::MDH_FLAG_ONLINE &&
						 //	  (ulEvalInfoMask & MdhBitField::MDH_FLAG_PHASCOR) != MdhBitField::MDH_FLAG_PHASCOR &&
						 //	  (ulEvalInfoMask & MdhBitField::MDH_FLAG_NOISEADJSCAN) != MdhBitField::MDH_FLAG_NOISEADJSCAN &&
					 		  (ulEvalInfoMask & MdhBitField::MDH_FLAG_ACQEND) != MdhBitField::MDH_FLAG_ACQEND) {

							 // Set flag
							 readout.bIsValidScan = true;

						  } // Flag check

						  // Add readout to vector
						 vReadout.push_back(readout);

					 } // Readout
				 } // Loop over children
			 

				 //---------------------------------------------------------------------
				 // TX
				 //---------------------------------------------------------------------
				 pugi::xml_node nodeTX = node.child("TX");
		
				 // Loop through children of GC
				 for (pugi::xml_node child = nodeTX.first_child(); child; child = child.next_sibling())
				 {

					 // Convert name to string
					 std::string sChildName(child.name());

					 if (sChildName.compare("RfShape") == 0) {

						 // Create local readout object
						 TX tx;

						 // Get attribute 'Type'
						if (pugi::xml_attribute attr = child.attribute("Type")) // attribute really exists
						{
							string sExcitationType = attr.as_string();
							if (sExcitationType.compare("Excitation") == 0) {
								tx.txType = Excitation;
							}
							else {
								mexErrMsgOriginAndTxt(ptModule, "Excitation type unkown or not yet implemented.");
							}

						}

						// Get attribute 'Asymmetry'
						if (pugi::xml_attribute attr = child.attribute("Asymmetry")) // attribute really exists
						{
							tx.fAsymmetry = attr.as_double();
						}

						pugi::xml_node nodeTrigger = child.child("Trigger");
					
						// Get attribute 'Type'
						if (pugi::xml_attribute attr = nodeTrigger.attribute("Time")) // attribute really exists
						{
							tx.fDwellTime = attr.as_double()/10.0;
						}

						// Get attribute 'NumberOfSamples'
						if (pugi::xml_attribute attr = nodeTrigger.attribute("NumberOfSamples")) // attribute really exists
						{
							tx.lSamples = attr.as_int();
						}

						// Add readout to vector
						vTx.push_back(tx);

					 }
				
				 }

				 //---------------------------------------------------------------------
				 // GC
				 //---------------------------------------------------------------------
				 pugi::xml_node nodeGC = node.child("GC");

				 // Loop through children of GC
				 for (pugi::xml_node child = nodeGC.first_child(); child; child = child.next_sibling())
				 {

					 // Convert name to string
					 std::string sChildName(child.name());

					 if (sChildName.compare("Control") == 0) {
						 processControl(child, vReadout, vTx);
					 }
					 else if (sChildName.compare("Rotation") == 0) {
						 processRotation(child);
					 }
					 else if (sChildName.compare("Halt") == 0) {
						 goto FINISH;
					 }
					 else {
						 PRINT("Cannot handle: %s\n", child.name());
						 mexErrMsgOriginAndTxt(ptModule,"Unknown GC instruction.");
					 }
				 }



			 }
			 else if (sNodeName.compare("Continue") == 0) {

				 // Get attribute 'File'
				 if (pugi::xml_attribute attr = node.attribute("File")) // attribute really exists
				 {
					 // Get next file name
					 pugi::char_t *cFilename = (char*)attr.value();

					 fullpath = m_fpXMLFile.path + string(cFilename); // The extension is already included in cFilename

				 }
				 else {
					 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'File' in Continue.");
				 }
			 }
		 }

		FINISH:

		 //---------------------------------------------------------------------
		 // Load next xml file is necessary
		 //---------------------------------------------------------------------
		 if (fullpath.size()>0) {

			 if (!m_doc.load_file(fullpath.c_str())) {
				 mexErrMsgOriginAndTxt(ptModule, "Error loading next xml file: ", fullpath.c_str());
			 }
			 else {
				 // Jump back and parse the next file
				 goto START;
			 }
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

		 if (m_ECC_Coeff_X.vfTau.size() > 0 && m_ECC_Coeff_X.vfTau.size() > 0 &&
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
			 // No coefficents have been provided
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
		long  x0 = (long)ceil((m_adTXCenterTimes[cTXPulse] - sfGRT) / sfGRT);

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
				 mexErrMsgOriginAndTxt(ptModule,"DSP::interpolateData(): Invalid array index.\n");

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

		long lMemRequ = m_lRXEvents * sizeof(uint32_t) + (m_lRXSampleLength + m_lTXEvents + m_lTrigEvents) * sizeof(double);

		if (m_lOutputMode == OutputMode::FULL) {
			// Data + exponentials
			lMemRequ += 6 * m_lConvolutionlLength * sizeof(double);
		}
		else {
			// Data + exponentials + interpolated data
			lMemRequ += 6 * m_lConvolutionlLength * sizeof(double);
			lMemRequ += 3 * m_lRXSampleLength * sizeof(double);
		}


		if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {
			PRINT("Allocating %.2f MB \n\n", 1.0 * lMemRequ / 1e6);
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

			long lMemRequ = m_lRXEvents * sizeof(uint32_t) + (m_lRXSampleLength + m_lTXEvents + m_lTrigEvents) * sizeof(double);

			if (m_lOutputMode == OutputMode::FULL) {
				// Data + exponentials
				lMemRequ += 6 * m_lConvolutionlLength * sizeof(double);
			}
			else {
				// Data + exponentials + interpolated data
				lMemRequ += 6 * m_lConvolutionlLength * sizeof(double);
				lMemRequ += 3 * m_lRXSampleLength * sizeof(double);
			}


			if ((m_uiVerbose & DISPLAY_BASIC) == DISPLAY_BASIC) {
				PRINT("Allocating %.2f MB \n\n", 1.0 * lMemRequ / 1e6);
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


			if (m_afMULTI_PURPOSE_INTERPX != NULL) { delete[] m_afMULTI_PURPOSE_INTERPX; m_afMULTI_PURPOSE_INTERPX = NULL; }
			if (m_afMULTI_PURPOSE_INTERPY != NULL) { delete[] m_afMULTI_PURPOSE_INTERPY; m_afMULTI_PURPOSE_INTERPY = NULL; }
			if (m_afMULTI_PURPOSE_INTERPZ != NULL) { delete[] m_afMULTI_PURPOSE_INTERPZ; m_afMULTI_PURPOSE_INTERPZ = NULL; }

			m_afMULTI_PURPOSE_INTERPX = new double[m_lRXSampleLength];
			m_afMULTI_PURPOSE_INTERPY = new double[m_lRXSampleLength];
			m_afMULTI_PURPOSE_INTERPZ = new double[m_lRXSampleLength];

			//---------------------------------------------------------------------
			// Exponentials                                                      
			//---------------------------------------------------------------------
			if (m_afExponentialX != NULL) { delete[] m_afExponentialX; m_afExponentialX = NULL; }
			if (m_afExponentialY != NULL) { delete[] m_afExponentialY; m_afExponentialY = NULL; }
			if (m_afExponentialZ != NULL) { delete[] m_afExponentialZ; m_afExponentialZ = NULL; }

			m_afExponentialX = new double[m_lConvolutionlLength];
			m_afExponentialY = new double[m_lConvolutionlLength];
			m_afExponentialZ = new double[m_lConvolutionlLength];

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
	 
		 if (!m_doc.load_file(fullpath.c_str()) ){
			 mexErrMsgOriginAndTxt(ptModule,"Error loading xml file.\n");
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

		 //---------------------------------------------------------------------
		 // Read shapes
		 //---------------------------------------------------------------------
		 pugi::xml_node shapes = m_doc.child("NUMARIS4_DSP_SIMULATION").child("Shapes");

		 // Loop through shapes
		 for (pugi::xml_node nodeShape = shapes.first_child(); nodeShape; nodeShape = nodeShape.next_sibling())
		 {

			 // Convert name to string
			 std::string sShapeName(nodeShape.name());

			 //---------------------------------------------------------------------
			 // TX shapes
			 //---------------------------------------------------------------------
			 if (sShapeName.compare("TXShape") == 0) {
				 // Nothing to do
			 }
			 //---------------------------------------------------------------------
			 // GC Shape
			 //---------------------------------------------------------------------
			 else if (sShapeName.compare("GCShape") == 0) {

				 long lSamples;

				 // Get attribute 'Samples'
				 if (pugi::xml_attribute attr = nodeShape.attribute("Samples")) // attribute really exists
				 {
					 // Get number of samples
					 lSamples = attr.as_int();

					 // Create vector for current shape
					 vector<double> vRow;
					 vRow.reserve(lSamples);

					 // Collect all data
					 for (pugi::xml_node data = nodeShape.first_child(); data; data = data.next_sibling()) {
						 vRow.push_back(data.text().as_double());
					 }

					 // Append to vector to vector
					 m_vGCShapes.push_back(vRow);

				 }
				 else {
					 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'Samples' in GCShape.");
				 }

			 }
			 //---------------------------------------------------------------------
			 // Otherwise
			 //---------------------------------------------------------------------
			 else {
				 mexErrMsgOriginAndTxt(ptModule, "Error unknown shape name.");
			 }

		 } // loop Shapes children
	 }


	 //*******************************************************************
	 // PROCESSROTATION()                                                 
	 //******************************************************************* 

	 void DSP::processRotation(pugi::xml_node &node) {

		 static const char *ptModule = { "DSP::processRotation(): " };

		 //---------------------------------------------------------------------
		 // Init
		 //---------------------------------------------------------------------
		 LogicalAxis logicalAxis;
		 GradientAxis gradientAxis;

		 //---------------------------------------------------------------------
		 // Read attributes
		 //---------------------------------------------------------------------

		 // Get attribute 'Log'
		 if (pugi::xml_attribute attr = node.attribute("Log")) // attribute really exists
		 {
			 const pugi::char_t *cLog = attr.value();
			 if (cLog[0] == 'A') {
				 logicalAxis = LogicalAxis::LOG_A;
			 }
			 else if (cLog[0] == 'B') {
				 logicalAxis = LogicalAxis::LOG_B;
			 }
			 else {
				 mexErrMsgOriginAndTxt(ptModule, "Unknown 'Log' in 'Axis'.\n");
			 }
		 }
		 else {
			 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'Log' in 'Axis'.\n");
		 }

		 // Get attribute 'Gradient'
		 if (pugi::xml_attribute attr = node.attribute("Gradient")) // attribute really exists
		 {
			 // Convert name to string
			 std::string sAttrStr(attr.value());

			 if (sAttrStr.compare("PE") == 0) {
				 gradientAxis = GradientAxis::GRAD_PE;
			 }
			 else if (sAttrStr.compare("RO") == 0) {
				 gradientAxis = GradientAxis::GRAD_RO;
			 }
			 else if (sAttrStr.compare("SL") == 0) {
				 gradientAxis = GradientAxis::GRAD_SL;
			 }
			 else {
				 mexErrMsgOriginAndTxt(ptModule, "Unknown value for attribute 'Gradient' in 'Axis': ", attr.value());
			 }
		 }
		 else {
			 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'Gradient' in 'Axis'.");
		 }

		 //---------------------------------------------------------------------
		 // Get child = vector
		 //---------------------------------------------------------------------
		 pugi::xml_node nodeVector = node.child("Vector");

		 double afrotvector[3];

		 // Get attribute 'Value'
		 if (pugi::xml_attribute attr = nodeVector.attribute("Rot0")) {
			 afrotvector[0] = attr.as_double();
		 }
		 else {
			 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'Rot0' in 'Vector'.");
		 }

		 if (pugi::xml_attribute attr = nodeVector.attribute("Rot1")) {
			 afrotvector[1] = attr.as_double();
		 }
		 else {
			 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'Rot1' in 'Vector'.");
		 }

		 if (pugi::xml_attribute attr = nodeVector.attribute("Rot2")) {
			 afrotvector[2] = attr.as_double();
		 }
		 else {
			 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'Rot2' in 'Vector'.");
		 }

		 //---------------------------------------------------------------------
		 // Update rotation matrix
		 //---------------------------------------------------------------------
		 switch (logicalAxis) {
		 case LOG_A:
			 m_aiMatrixA[0][gradientAxis] = afrotvector[0];
			 m_aiMatrixA[1][gradientAxis] = afrotvector[1];
			 m_aiMatrixA[2][gradientAxis] = afrotvector[2];
			 break;
		 case LOG_B:
			 m_aiMatrixB[0][gradientAxis] = afrotvector[0];
			 m_aiMatrixB[1][gradientAxis] = afrotvector[1];
			 m_aiMatrixB[2][gradientAxis] = afrotvector[2];
			 break;
		 }

		 if ((m_uiVerbose & DISPLAY_ROTMAT) == DISPLAY_ROTMAT) {
			 PRINT("Matrix A\n");
			 PRINT("|%.2f %.2f %.2f|\n", m_aiMatrixA[0][0], m_aiMatrixA[0][1], m_aiMatrixA[0][2]);
			 PRINT("|%.2f %.2f %.2f|\n", m_aiMatrixA[1][0], m_aiMatrixA[1][1], m_aiMatrixA[1][2]);
			 PRINT("|%.2f %.2f %.2f|\n", m_aiMatrixA[2][0], m_aiMatrixA[2][1], m_aiMatrixA[2][2]);

			 PRINT("Matrix B\n");
			 PRINT("|%.2f %.2f %.2f|\n", m_aiMatrixB[0][0], m_aiMatrixB[0][1], m_aiMatrixB[0][2]);
			 PRINT("|%.2f %.2f %.2f|\n", m_aiMatrixB[1][0], m_aiMatrixB[1][1], m_aiMatrixB[1][2]);
			 PRINT("|%.2f %.2f %.2f|\n", m_aiMatrixB[2][0], m_aiMatrixB[2][1], m_aiMatrixB[2][2]);
		 }

	 }

	 //*******************************************************************
	 // PROCESSSHAPE()                                                    
	 //******************************************************************* 

	 void DSP::processShape(pugi::xml_node &node, Axis *axis) {

		 static const char *ptModule = { "DSP::processShape(): " };

		 //---------------------------------------------------------------------
		 // Init
		 //---------------------------------------------------------------------
		 LogicalAxis logicalAxis;
		 int iID = -1;

		 //---------------------------------------------------------------------
		 // Read attributes
		 //---------------------------------------------------------------------
		 // Get attribute 'Value'
		 if (pugi::xml_attribute attr = node.attribute("ID")) // attribute really exists
		 {
			 iID = attr.as_int();
		 }
		 else {
			 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'ID' in 'Shape'.");
		 }

		 // Get attribute 'Log'
		 if (pugi::xml_attribute attr = node.attribute("Log")) // attribute really exists
		 {
			 const pugi::char_t *cLog = attr.value();
			 if (cLog[0] == 'A') {
				 logicalAxis = LogicalAxis::LOG_A;
			 }
			 else if (cLog[0] == 'B') {
				 logicalAxis = LogicalAxis::LOG_B;
			 }
			 else {
				 mexErrMsgOriginAndTxt(ptModule, "Unknown 'Log' in 'Shape': ", attr.value());
			 }
		 }
		 else {
			 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'Log' in 'Shape'.");
		 }

		 //---------------------------------------------------------------------
		 // Read attributes and set properites of shape
		 //---------------------------------------------------------------------
		 // Get attribute 'Gradient'
		 if (pugi::xml_attribute attr = node.attribute("Gradient")) // attribute really exists
		 {

			 // Convert name to string
			 std::string sAttrStr(attr.value());

			 if (sAttrStr.compare("PE") == 0) {
				 axis[logicalAxis].Shape_PE.ID = iID;
			 }
			 else if (sAttrStr.compare("RO") == 0) {
				 axis[logicalAxis].Shape_RO.ID = iID;
			 }
			 else if (sAttrStr.compare("SL") == 0) {
				 axis[logicalAxis].Shape_SL.ID = iID;
			 }
			 else {
				 mexErrMsgOriginAndTxt(ptModule, "Unknown value for attribute 'Gradient' in 'Shape': ", attr.value());
			 }
		 }
		 else {
			 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'Gradient' in 'Shape'.");
		 }

	 }

	 //*******************************************************************
	 // PROCESSINCREMENT()                                                
	 //******************************************************************* 

	 void DSP::processIncrement(pugi::xml_node &node, Axis *axis) {

		 static const char *ptModule = { "DSP::processIncrement(): " };

		 //---------------------------------------------------------------------
		 // Init
		 //---------------------------------------------------------------------
		 LogicalAxis logicalAxis;
		 double fOffset = 0.0;
		 double fValue = 0.0;

		 //---------------------------------------------------------------------
		 // Read attributes
		 //---------------------------------------------------------------------
		 // Get attribute 'Value'
		 if (pugi::xml_attribute attr = node.attribute("Value")) // attribute really exists
		 {
			 fValue = attr.as_double();
		 }

		 // Get attribute 'Offset'
		 if (pugi::xml_attribute attr = node.attribute("Offset")) // attribute really exists
		 {
			 fOffset = attr.as_double();
		 }

		 // Get attribute 'Log'
		 if (pugi::xml_attribute attr = node.attribute("Log")) // attribute really exists
		 {
			 const pugi::char_t *cLog = attr.value();
			 if (cLog[0] == 'A') {
				 logicalAxis = LogicalAxis::LOG_A;
			 }
			 else if (cLog[0] == 'B') {
				 logicalAxis = LogicalAxis::LOG_B;
			 }
			 else {
				 mexErrMsgOriginAndTxt(ptModule, "Unknown 'Log' in 'Increment': ", attr.value());
			 }
		 }
		 else {
			 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'Log' in 'Increment'.");
		 }

		 //---------------------------------------------------------------------
		 // Read attributes and set properties of axis
		 //---------------------------------------------------------------------
		 // Get attribute 'Gradient'
		 if (pugi::xml_attribute attr = node.attribute("Gradient")) // attribute really exists
		 {

			 // Convert name to string
			 std::string sAttrStr(attr.value());

			 if (sAttrStr.compare("PE") == 0) {
				 axis[logicalAxis].Increment_PE.fOffset = fOffset;
				 axis[logicalAxis].Increment_PE.fValue  = fValue;
			 }
			 else if (sAttrStr.compare("RO") == 0) {
				 axis[logicalAxis].Increment_RO.fOffset = fOffset;
				 axis[logicalAxis].Increment_RO.fValue  = fValue;
			 }
			 else if (sAttrStr.compare("SL") == 0) {
				 axis[logicalAxis].Increment_SL.fOffset = fOffset;
				 axis[logicalAxis].Increment_SL.fValue  = fValue;
			 }
			 else {
				 mexErrMsgOriginAndTxt(ptModule, "Unknown value for attribute 'Gradient' in 'Increment'.");
			 }
		 }
		 else {
			 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'Gradient' in 'Increment'.");
		 }

	 }



	 //*******************************************************************
	 // PROCESSAXIS()                                                     
	 //******************************************************************* 

	 void DSP::processAxis(pugi::xml_node &node, Axis *axis) {

		 static const char *ptModule = { "DSP::processAxis(): " };

		 //---------------------------------------------------------------------
		 // Loop through children
		 //---------------------------------------------------------------------
		 for (pugi::xml_node child = node.first_child(); child; child = child.next_sibling())
		 {
			 // Convert name to string
			 std::string sChildName(child.name());

			 //---------------------------------------------------------------------
			 // Shape
			 //---------------------------------------------------------------------
			 if (sChildName.compare("Shape") == 0) {
				 processShape(child, axis);
			 }
			 //---------------------------------------------------------------------
			 // Increment
			 //---------------------------------------------------------------------
			 else if (sChildName.compare("Increment") == 0) {
				 processIncrement(child, axis);
			 }
			 //---------------------------------------------------------------------
			 // Otherwise
			 //---------------------------------------------------------------------
			 else {
				 mexErrMsgOriginAndTxt(ptModule, "Unknown Axis instruction: ", child.name());
			 }
		 }
	 }

	 //*******************************************************************
	 // APPLYAXIS()                                                       
	 //******************************************************************* 

	 // Combine shape and apply rotation
	 void DSP::applyAxis(Axis *axis, double *GX, double *GY, double *GZ, long lTicks) {

		 //---------------------------------------------------------------------
		 // Allocate space for gradient shapes on both logical axes
		 //---------------------------------------------------------------------
		 double *Phase_A = new double[lTicks];
		 double *Read_A = new double[lTicks];
		 double *Slice_A = new double[lTicks];

		 double *Phase_B = new double[lTicks];
		 double *Read_B = new double[lTicks];
		 double *Slice_B = new double[lTicks];

		 //---------------------------------------------------------------------
		 // Loop over samples
		 //---------------------------------------------------------------------
		 for (long i = 0; i < lTicks; i++) {

			 //---------------------------------------------------------------------
			 // Log axis A
			 //---------------------------------------------------------------------
			 if (axis[LOG_A].Shape_PE.ID >= 0) {
				 Phase_A[i] = (m_vGCShapes[axis[LOG_A].Shape_PE.ID]).at(i);
			 }
			 else {
				 Phase_A[i] = i*axis[LOG_A].Increment_PE.fValue + axis[LOG_A].Increment_PE.fOffset;
			 }

			 if (axis[LOG_A].Shape_RO.ID >= 0) {
				 Read_A[i] = (m_vGCShapes[axis[LOG_A].Shape_RO.ID]).at(i);
			 }
			 else {
				 Read_A[i] = i*axis[LOG_A].Increment_RO.fValue + axis[LOG_A].Increment_RO.fOffset;
			 }

			 if (axis[LOG_A].Shape_SL.ID >= 0) {
				 Slice_A[i] = (m_vGCShapes[axis[LOG_A].Shape_SL.ID]).at(i);
			 }
			 else {
				 Slice_A[i] = i*axis[LOG_A].Increment_SL.fValue + axis[LOG_A].Increment_SL.fOffset;
			 }

			 //---------------------------------------------------------------------
			 // Log axis B
			 //---------------------------------------------------------------------
			 if (axis[LOG_B].Shape_PE.ID >= 0) {
				 Phase_B[i] = (m_vGCShapes[axis[LOG_B].Shape_PE.ID]).at(i);
			 }
			 else {
				 Phase_B[i] = i*axis[LOG_B].Increment_PE.fValue + axis[LOG_B].Increment_PE.fOffset;
			 }

			 if (axis[LOG_B].Shape_RO.ID >= 0) {
				 Read_B[i] = (m_vGCShapes[axis[LOG_B].Shape_RO.ID]).at(i);
			 }
			 else {
				 Read_B[i] = i*axis[LOG_B].Increment_RO.fValue + axis[LOG_B].Increment_RO.fOffset;
			 }

			 if (axis[LOG_B].Shape_SL.ID >= 0) {
				 Slice_B[i] = (m_vGCShapes[axis[LOG_B].Shape_SL.ID]).at(i);
			 }
			 else {
				 Slice_B[i] = i*axis[LOG_B].Increment_SL.fValue + axis[LOG_B].Increment_SL.fOffset;
			 }

			 //---------------------------------------------------------------------
			 // Rotate to xyz and sum
			 //---------------------------------------------------------------------
			 GX[i] += m_aiMatrixA[0][0] * Phase_A[i] + m_aiMatrixA[0][1] * Read_A[i] + m_aiMatrixA[0][2] * Slice_A[i];
			 GY[i] += m_aiMatrixA[1][0] * Phase_A[i] + m_aiMatrixA[1][1] * Read_A[i] + m_aiMatrixA[1][2] * Slice_A[i];
			 GZ[i] += m_aiMatrixA[2][0] * Phase_A[i] + m_aiMatrixA[2][1] * Read_A[i] + m_aiMatrixA[2][2] * Slice_A[i];

			 GX[i] += m_aiMatrixB[0][0] * Phase_B[i] + m_aiMatrixB[0][1] * Read_B[i] + m_aiMatrixB[0][2] * Slice_B[i];
			 GY[i] += m_aiMatrixB[1][0] * Phase_B[i] + m_aiMatrixB[1][1] * Read_B[i] + m_aiMatrixB[1][2] * Slice_B[i];
			 GZ[i] += m_aiMatrixB[2][0] * Phase_B[i] + m_aiMatrixB[2][1] * Read_B[i] + m_aiMatrixB[2][2] * Slice_B[i];

		 }

		 //---------------------------------------------------------------------
		 // Free memory
		 //---------------------------------------------------------------------
		 delete[] Phase_A;
		 delete[] Read_A;
		 delete[] Slice_A;

		 delete[] Phase_B;
		 delete[] Read_B;
		 delete[] Slice_B;

	 }

	 //*******************************************************************
	 // PROCESSTRIGGER()                                                  
	 //******************************************************************* 

	 void DSP::processTrigger(pugi::xml_node &node, vector<Readout> &vReadout, vector<TX> &vTx) {

		 static const char *ptModule = { "DSP::processTrigger(): " };

		 //---------------------------------------------------------------------
		 // Get first child which contains the data (RX, TX, or FreqPhase)                                                        
		 //---------------------------------------------------------------------
		 pugi::xml_node data = node.first_child();

		 // Convert name to string
		 std::string sTrigger(data.value());

		 //---------------------------------------------------------------------
		 // RX trigger                                                       
		 //---------------------------------------------------------------------
		 if (sTrigger.compare("RX") == 0) {

			 if (vReadout.size() == 0) {
				 mexErrMsgOriginAndTxt(ptModule, "DSP::processTrigger: Size of vReadout is zero.");
			 }

			 // Get RX trigger delay
			 double dDelay = 0.0;
			 if (pugi::xml_attribute attr = node.attribute("Delay")) // attribute really exists
			 {
				 dDelay = attr.as_double();
			 }

			 //---------------------------------------------------------------------
			 // Is it an imaging scan?                                                   
			 //---------------------------------------------------------------------
			 if (vReadout.at(0).bIsValidScan) {

				 long  lRXSamples = vReadout.at(0).lSamples;	// Take the first in the vector
				 double dDwellTime = vReadout.at(0).fDwellTime;	// Take the first in the vector
			 
				 // Check for space
				 if (m_lCurrentRXSampleLength + lRXSamples > m_lRXSampleLength)
					 mexErrMsgOriginAndTxt(ptModule, "Cannot append RX sampling times. Allocated array is too small.");

				 // Fill in array
				 for (long t = 0; t < lRXSamples; t++) {
					 double dStartime = (10.0 *m_lCurrentGCSampleLength ) + dDelay; //10 us per GRT tick
					 m_adRXTimes[m_lCurrentRXSampleLength + t] = (dStartime + dDwellTime*(double)t)*1e-6; // Final unit is seconds
				 }

				 // Update counter
				 m_lCurrentRXSampleLength += lRXSamples;
				 m_uiRXEventLength[m_lCurrentRXNumber] = m_lCurrentRXSampleLength;
				 m_lCurrentRXNumber++;

			 }

			 //---------------------------------------------------------------------
			 // Erase the first element and relocate all the other elements                                                      
			 //---------------------------------------------------------------------
			 vReadout.erase(vReadout.begin());
		 }

		 //---------------------------------------------------------------------
		 // TX trigger                                                       
		 //---------------------------------------------------------------------
		 else if (sTrigger.compare("TX") == 0) {
	 
			 if (vTx.size() == 0) {
				 mexErrMsgOriginAndTxt(ptModule, "DSP::processTrigger: Size of vTX is zero.");
			 }

			 // Get RX trigger delay
			 double dDelay = 0.0;
			 if (pugi::xml_attribute attr = node.attribute("Delay")) // attribute really exists
			 {
				 dDelay = attr.as_double();
			 }

			 //---------------------------------------------------------------------
			 // Is it an imaging scan?                                                   
			 //---------------------------------------------------------------------
			 long  lTXSamples  = vTx.at(0).lSamples;	// Take the first in the vector
			 double dDwellTime = vTx.at(0).fDwellTime;	// Take the first in the vector
			 double dAsymmetry = vTx.at(0).fAsymmetry;	// Take the first in the vector

		 
			 // Check for space
			 if (m_lCurrentTXNumber + 1 > m_lTXEvents)
				 mexErrMsgOriginAndTxt(ptModule, "Cannot append TX sampling times. Allocated array is too small.");

			 // Fill in array
		 
			 m_adTXCenterTimes[m_lCurrentTXNumber] = 10.0*m_lCurrentGCSampleLength + dDwellTime*lTXSamples*dAsymmetry; //10 us per GRT tick 
			 m_adTXCenterTimes[m_lCurrentTXNumber] *= 1e-6; // Convert to seconds

			 //PRINT("dDwellTime %f, lTXSamples %d, dAsymmetry %f, m_adTXCenterTimes[%d] %f\n", dDwellTime, lTXSamples, dAsymmetry, m_lCurrentTXNumber, m_adTXCenterTimes[m_lCurrentTXNumber]);

			 // Update counter
			 m_lCurrentTXNumber++;
		 

			 //---------------------------------------------------------------------
			 // Erase the first element and relocate all the other elements                                                      
			 //---------------------------------------------------------------------
			 vTx.erase(vTx.begin());
	 
		 }

		 //---------------------------------------------------------------------
		 // FreqPhase trigger                                                       
		 //---------------------------------------------------------------------
		 else if (sTrigger.compare("FreqPhase") == 0) {}

		 //---------------------------------------------------------------------
		 // RX trigger                                                       
		 //---------------------------------------------------------------------
		 else { 
			 mexErrMsgOriginAndTxt(ptModule, "Unknown trigger value: ", data.value());
		 }


	 }

	 //*******************************************************************
	 // PROCESSCONTROL()                                                  
	 //******************************************************************* 

	 // Process constrol information
	 void DSP::processControl(pugi::xml_node &node, vector<Readout> &vReadout, vector<TX> &vTx) {

		 static const char *ptModule = { "DSP::processControl(): " };

		 // Get attribute 'Ticks'
		 if (pugi::xml_attribute attr = node.attribute("Ticks")) // attribute really exists
		 {
			 long lTicks = attr.as_int();

			 // Loop through children
			 for (pugi::xml_node child = node.first_child(); child; child = child.next_sibling())
			 {

				 // Convert name to string
				 std::string sChildName(child.name());

				 //---------------------------------------------------------------------
				 // Sync                                                        
				 //---------------------------------------------------------------------
				 if (sChildName.compare("Sync") == 0) {

					 // Get attribute 'Value'
					 if (pugi::xml_attribute attr = child.attribute("Ext")) // attribute really exists
					 {
						 m_adTrigTimes[m_lCurrentTrigNumber] = sfGRT*(double)m_lCurrentGCSampleLength;
						 m_lCurrentTrigNumber++;
					 }
				 
				 }
				 //---------------------------------------------------------------------
				 // Trigger                                                         
				 //---------------------------------------------------------------------
				 else if (sChildName.compare("Trigger") == 0) {
					 processTrigger(child, vReadout, vTx);
				 }
				 //---------------------------------------------------------------------
				 // Axis                                                        
				 //---------------------------------------------------------------------
				 else if (sChildName.compare("Axis") == 0) {

					 //---------------------------------------------------------------------
					 // Init axis object                                                        
					 //---------------------------------------------------------------------
					 Axis axis[2] = {};

					 //---------------------------------------------------------------------
					 // Process axis information                                                         
					 //---------------------------------------------------------------------
					 processAxis(child, axis);
				 
					 //---------------------------------------------------------------------
					 // Check for memory                                                        
					 //---------------------------------------------------------------------
					 if (m_lCurrentGCSampleLength + lTicks > m_lGradientShapeLength) {
						 mexErrMsgOriginAndTxt(ptModule, "Cannot append gradient shape. Allocated memory is insufficient.");
					 }

					 //---------------------------------------------------------------------
					 // Apply axis information                                                         
					 //---------------------------------------------------------------------
					 applyAxis(axis, &m_afMULTI_PURPOSEX[m_lCurrentGCSampleLength], &m_afMULTI_PURPOSEY[m_lCurrentGCSampleLength], &m_afMULTI_PURPOSEZ[m_lCurrentGCSampleLength], lTicks);
				 }
				 else {
					 mexErrMsgOriginAndTxt(ptModule, "Unknown Control instruction: ", child.name());
				 }
			 }

			 // Increment current shape length
			 m_lCurrentGCSampleLength += attr.as_int();

		 }

		 else {
			 mexErrMsgOriginAndTxt(ptModule, "Could not find attribute 'Ticks' in 'Control'.");
		 }

	 }

	 void DSP::applyPhaseModulation(complex<float> *data, unsigned short numberOfSamples, unsigned int ulScanCounter) {

		// Scancounter should start at zero

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
				 sprintf_s(tLine, "Number of RX samples (%d) for current scan counter (%d) is not equal to expected number (%d).\n", numberOfSamples, ulScanCounter + 1, lExpectedNumberOfSamples);  //Scan counter (Siemens raw data) starts orginally at 1. Therefor we add here 1 again.
				 mexErrMsgOriginAndTxt(ptModule, tLine);
			 }
         
			 long lStartIndex = m_uiRXEventLength[ulScanCounter];

			 // Apply phase modulation	
			 for (long i = 0; i < numberOfSamples; i++) {
				  //data[i] =  data[i] * polar((float)1.0, (float)m_afMULTI_PURPOSE_INTERPX[lStartIndex + i]);
				  data[i] = complex<float>((float)m_afMULTI_PURPOSE_INTERPX[lStartIndex + i],0.0);  //For testing 
			 }
        
		 }
		 
	 }



	 double DSP::getADCStartTime(unsigned int ulScanCounter) {
		 static const char *ptModule = { "DSP::getADCStartTime(): " };

		 throw "Not implemented.";
		 
			//if (m_bB0DataAvailable) {
			//
			//	 // Checks
			//	 if (ulScanCounter == m_lRXEvents) {
			//		 PRINT(ptModule, "ACQEND will not be processed.\n");
			//		 // ACQEND is not processed by DSP.cpp because it appears in the XML file after the GC-<Halt> instruction.
			//		 return -1;
			//	 }
			//	 else if (ulScanCounter > m_lRXEvents - 1) {
			//		 mexErrMsgOriginAndTxt(ptModule, "Scan counter exceeds number of stored ADC events.\n");
			//	 }
			//	 else if (ulScanCounter < 0) {
			//		 mexErrMsgOriginAndTxt(ptModule, "Scan counter cannot be smaller than zero.\n");
			//	 }
			//
			//	 return m_adRXTimes[m_vStartIndexRXEvents.at(ulScanCounter)];
			//
			//}
			//return -1;
			//
		 
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
			PRINT("#         DSP READER (VB17)         #\n");
			PRINT("#              and                  #\n");
			PRINT("#         ECC calculator            #\n");
			PRINT("#       Version 1.1 (2018)          #\n");
			PRINT("#####################################\n\n");
		}

		//---------------------------------------------------------------------
		// Open file                                                         
		//---------------------------------------------------------------------
		this->openFile();

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
		// Run DSP instructions                                                   
		//---------------------------------------------------------------------
		this->runInstructions();

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

			//---------------------------------------------------------------------
			// Calculate phase by trapezoidal integration
			//---------------------------------------------------------------------
			// Do not null phase at center of TX pulse
			calculateIntegral(m_afMULTI_PURPOSEX, m_afMULTI_PURPOSEX, false);

			for (long t = 0; t < m_lConvolutionlLength; t++) {
				m_afMULTI_PURPOSEX[t] *= 2.0*sfPI*sfGAMMA_1H;
			}
						
		}

		//---------------------------------------------------------------------
		// Interpolate to RX events                      
		//---------------------------------------------------------------------
		if (m_lOutputMode == OutputMode::INTERPOLATED_TO_RX) {
			this->interpolateData(m_afMULTI_PURPOSEX, m_afMULTI_PURPOSE_INTERPX);

			if (m_lDataType != DataType::EDDYPHASE) {
				this->interpolateData(m_afMULTI_PURPOSEY, m_afMULTI_PURPOSE_INTERPY);
				this->interpolateData(m_afMULTI_PURPOSEZ, m_afMULTI_PURPOSE_INTERPZ);
			}

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

		if (m_lOutputMode == OutputMode::INTERPOLATED_TO_RX && m_lDataType == DataType::EDDYPHASE) {
			m_bEccCompensationAvailable = true;
		}
	}
}
