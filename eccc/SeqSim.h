// 	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
//	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
// 	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
// 	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
// 	OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
//  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//	  Project:	ECC simulator
//	     File:	SeqSim.h
//	  Version:	2.0
//	     Date:	26.08.2020
//	     Lang:	C++
//
//Description:	This class can read simulated gradient shapes and calculate the eddy current compensation
// 				applied by the scanner based on the provided amplitude and time constants of the eddy currents.
// 				
//				Aside from the B0 eddy current phase, this class can also output the following types of data:
//				 	- 	gradients
//					-	k-space
// 					- 	slew rate
// 					-	B0 eddy currents
//				
// 				The class is intended to be compiled as a dll and integrated into the SiemensToMRD converter. 
//
// 				A wrapper for Matlab (SeqSim_mex.cpp) and a console application (SeqSimRun.cpp) are also available. 
//				The compilation instructions can be found within the source files. 
//
//	 Options:   setDataType defines the type to be returned by the class: 
//					- GRADIENT: gradient
//					- KSPACE		: k-space	
//					- SLEWRATE		: slew rate
//					- EDDYCURRENT	: B0 eddy currents
//					- EDDYPHASE		: B0 eddy current phase
//
//				setOutputMode:
//					FULL 				: return the full simulated data shape on the gradient raster time (GRT)
//					INTERPOLATED_TO_RX 	: return the interpolated data shape at the RX sample locations
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

#include <stdio.h>
#include <string>        
#include <math.h>  
#include <vector>  
#include <ctime>
#include <complex>


#ifdef MATLAB_MEX_FILE 
	#include "mex.h"
#endif

#pragma once  

using namespace std;

struct FileParts
{
	std::string path; // containing folder, if provided, including trailing slash
	std::string name; // base file name, without extension
	std::string ext;  // extension, including '.'
};

namespace SEQSIM
{
	//*******************************************************************
	// STRUCTURES                                                        
	//******************************************************************* 
	// Coefficients describing the decay of the B0 eddy currents
	#ifndef ECC_COEFF
		#define ECC_COEFF
		struct ECC_Coeff 
		{
			std::vector<double> vfTau;	// decay constant [s]
			std::vector<double> vfAmp;	// decay amplitude
		};
	#endif
	
	//*******************************************************************
	// ENUMS                                                             
	//******************************************************************* 
	// Data type returned by class
	enum DataType {
		GRADIENT,
		KSPACE,
		SLEWRATE,
		EDDYCURRENT,
		EDDYPHASE
	};

	// Data type returned by class
	enum OutputMode {
		FULL,
		INTERPOLATED_TO_RX
	};

	enum Verbose {
		DISPLAY_NONE = 0b000000,
		DISPLAY_BASIC = 0b000001,
		DISPLAY_ADVANCED = 0b000011,
		DISPLAY_ROTMAT = 0b000100,
		DISPLAY_INCR_OFFSET = 0b001000,
		DISPLAY_ECC_COEFFS = 0b010000,
		DISPLAY_DSP_INFO = 0b100000,
		DISPLAY_ALL = 0b111111
	};	

	//*******************************************************************
	// DSP CLASS                                                         
	//******************************************************************* 

	class DSP {
	public:

		// Constructor
		DSP();

		// Destructor
		~DSP();

		// Set file name of xml file to be read
		void setFileName(const char *pFileName);

		// Set path to dsp folder
		void setDSVFolderPath(const char *pFolderPath);

		// Read xml file and simulate the DSP instructions 
		#ifdef MATLAB_MEX_FILE 
			void run(mxArray *plhs[]);
		#else
			void run();
		#endif

		// Set the amount of information printed to the console (see enum Verbose in Type.h)
		void setVerboseMode(int verbose);

		// Set the type of data return by the class
		void setDataType(DataType type);

		// Set the output mode
		void setOutputMode(OutputMode mode);


		// Set the B0 coefficients Coeff.Tau and Coeff.Amp
		void setB0CorrCoeff(ECC_Coeff &CoeffX, ECC_Coeff &CoeffY, ECC_Coeff &CoeffZ);

		void setB0CorrCoeff(vector<double> &vfCoeffXamp, vector<double> &vfCoeffXtau,
									vector<double> &vfCoeffYamp, vector<double> &vfCoeffYtau,
									vector<double> &vfCoeffZamp, vector<double> &vfCoeffZtau);

		void applyPhaseModulation(complex<float> *data, unsigned short numberOfSamples, unsigned int ulScanCounter);

		double getADCStartTime(unsigned int ulScanCounter);

		void writeToFile();
		
		void setDSVFileNamePrefix(const char * pDSVFileNamePrefix);

	protected:

		//*******************************************************************
		// FUNCTIONS                                                         
		//*******************************************************************

		// Open the text file
		void openFile();

		// Loops through the text file and counts the number of gradient and receiver events
		// Thereafter the text file is rewound till the first Halt gradient instruction
		void calcMemoryRequirement();

		// Calculate the exponential decay curves
		void computeExponentials();

		// Calculates the derivatives of the provided gradient
		// If the second argument is a Null pointer, the first input array will used as output array 
		void calculateDerivative(double *adData);

		// Calculates the derivatives for the gradient on each axis
		void calculateDerivative();

		void calculateIntegral();
		void calculateIntegral(double *afG, double *afS, bool bNullAtTXCenter = true);

		// Determines the longest decay time on all axes
		void determineLongestTimeConstant();

		// Computes the ECC via convolution for the specifed axis
		// If afS is a Null pointer, the first input array will used as input and output array 
		// If afS is a not a Null pointer, afS will used as input and afECC as output array 
		void computeECC(double *afG, double *afExponential);

		// Computes the ECC on each axis
		void computeECC();

		// Allocates memory. In debug mode all arrays are handed over to MATLAB at the end of 
		// the calculation
		#ifdef MATLAB_MEX_FILE 
			void allocateMemory(mxArray *plhs[]);
		#else
			void allocateMemory();
		#endif

		// Interpolate the phase at the positions of the readout samples
		void interpolateData(double *afData, double *afDataInterp);
		

		// Read arbitrary gradient shapes
		void readGCShapes();

		// Get DSV data
		bool getDSV(ifstream &in, double *y);

		// Get dsp parameters
		void getDSVParams(ifstream &in, long &lSamples, double &dHoriDelta, double &dVertFactor);

		// Get RX parameters
		void getInfoParams(ifstream &in, long &lRXSampleLength, long &lRXEvents, long &lTXEvents);

		// Get RX parameters
		void getInfo(ifstream &in);

		//*******************************************************************
		// PROPERTIES                                                        
		//*******************************************************************

		// Text file name
		FileParts m_fpXMLFile;

		// Path to DSV folder
		string m_sDSVFolderPath;
		
		string m_sDSVFileNamePrefix;

		//---------------------------------------------------------------------
		// Dimensions & Counters
		//---------------------------------------------------------------------

		// Length of gradient shape in gradient raster times
		long m_lGradientShapeLength;

		// Length of exponential decay in gradient raster times
		long m_lExponentialLength;

		// Length of convolution = m_lGradientShapeLength + m_lExponentialLength - 1
		long m_lConvolutionlLength;

		// Number of total RX samples
		long m_lRXSampleLength;

		// Number of RX events
		long m_lRXEvents;

		// Number of RX events
		long m_lTXEvents;

		// Number of Trigger events
		long m_lTrigEvents;

		// Current number of RX samples stored
		long m_lCurrentRXSampleLength;

		// Current number of gradient samples stores
		long m_lCurrentGCSampleLength;

		// Current number of stored TX events
		long m_lCurrentTXNumber;

		// Current number of stored RX events
		long m_lCurrentRXNumber;

		// Current number of stored trigger events
		long m_lCurrentTrigNumber;

		// Settings
		DataType m_lDataType;
		OutputMode m_lOutputMode;

		//---------------------------------------------------------------------
		// Matrices
		//---------------------------------------------------------------------

		// Logical matrix A
		double **m_aiMatrixA;

		// Logical matrix B
		double **m_aiMatrixB;

		//---------------------------------------------------------------------
		// Arrays
		//---------------------------------------------------------------------

		// In normal mode, no intermediate results are stored, but the array m_afMULTI_PURPOSE
		// is overwritten several times during computation: GX->SX->ECC->Phase
		// In debug mode, it holds the X gradient shape.
		double *m_afMULTI_PURPOSEX;
		double *m_afMULTI_PURPOSEY;
		double *m_afMULTI_PURPOSEZ;

		double *m_afMULTI_PURPOSE_INTERPX;
		double *m_afMULTI_PURPOSE_INTERPY;
		double *m_afMULTI_PURPOSE_INTERPZ;

		// Stores the exponential decay curve for X [on the gradient raster time (10 us)]
		double *m_afExponentialX;

		// Stores the exponential decay curve for Y [on the gradient raster time (10 us)]
		double *m_afExponentialY;

		// Stores the exponential decay curve for Z [on the gradient raster time (10 us)]
		double *m_afExponentialZ;

		// Stores the sampling times of the RX events [s]
		// Here we need double-precision because very large (hundreds of seconds) 
		// and very small numbers (hundreds of ns) are added
		double *m_adRXTimes;

		// TX center times
		double *m_adTXCenterTimes;

		// Length of RX events
		uint32_t *m_uiRXEventLength;

		// Trigger times
		double *m_adTrigTimes;
		
		//---------------------------------------------------------------------
		// Vectors
		//---------------------------------------------------------------------

		// Vector holding the gradient vector shapes
		vector< vector<double> > m_vGCShapes;
		
		//---------------------------------------------------------------------
		// Coefficients and values
		//---------------------------------------------------------------------
		// Coefficients describing the decay of the B0 eddy currents on the X axis
		ECC_Coeff m_ECC_Coeff_X;

		// Coefficients describing the decay of the B0 eddy currents on the Y axis
		ECC_Coeff m_ECC_Coeff_Y;

		// Coefficients describing the decay of the B0 eddy currents on the Z axis
		ECC_Coeff m_ECC_Coeff_Z;

		// Longest decay time on all axes
		double m_fLargestTau;

		//---------------------------------------------------------------------
		// Debug and output
		//---------------------------------------------------------------------
		// Set the amount of information printed to the console (see enum Verbose in Type.h)
		unsigned int m_uiVerbose;

		// Save intermediate results to output
		bool m_bDebugMode;

	private:
	
		// Converts seconds to a time string: X hours, Y minutes and Z seconds
		void getDurationString(long duration, char *str);

		// Start time of calculation
		time_t m_tstart;

		// End time of calculation
		time_t m_tend;

		// Check if B0 compensation is available
		bool m_bEccCompensationAvailable;

	};

}
