#include <stdio.h>
#include <string>        
#include <math.h>  
#include <vector>  
#include <ctime>
#include <complex>
#include "pugixml.hpp" 


#ifdef MATLAB_MEX_FILE 
	#include "mex.h"
#endif

#pragma once  

using namespace std;

struct FileParts
{
	std::string path; //!< containing folder, if provided, including trailing slash
	std::string name; //!< base file name, without extension
	std::string ext;  //!< extension, including '.'
};


namespace SEQSIM
{
	//*******************************************************************
	// STRUCTURES                                                        
	//******************************************************************* 
	// Coefficients defining the evolution of the gradient waveform
	struct GradDef {
		double flIncr;		// gradient increment [mT/m]
		double flOffset;		// gradient offset [mT/m]
	};

	// Coefficients describing the decay of the B0 eddy currents
	#ifndef ECC_COEFF
	#define ECC_COEFF
	struct ECC_Coeff 
	{
		std::vector<double> vfTau;	// decay constant [s]
		std::vector<double> vfAmp;	// decay amplitude
	};
	#endif
	// Gradient instructions for linear ramp
	struct Increment {
		double fValue;
		double fOffset;
	};

	// Store gradient shape
	struct Shape {
		int ID;
	};

	//*******************************************************************
	// ENUMS                                                             
	//******************************************************************* 
	// Two gradients are allowed to overlap. Therefore two logical axes are needed.
	enum LogicalAxis { 
		LOG_A,
		LOG_B
	};

	// The three logical gradient axes
	enum GradientAxis {
		GRAD_PE,
		GRAD_RO,
		GRAD_SL
	};

	// TX Type
	enum TxType {
		Undefined,
		Excitation,
		Inversion
	};

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
	
	enum MdhBitField {
		MDH_FLAG_ACQEND				= 1,
		MDH_FLAG_RTFEEDBACK			= 1 << 1,
		MDH_FLAG_HPFEEDBACK			= 1 << 2,
		MDH_FLAG_ONLINE				= 1 << 3,
		MDH_FLAG_OFFLINE			= 1 << 4,
		MDH_FLAG_SYNCDATA			= 1 << 5,       // readout contains synchroneous data
		MDH_FLAG_LASTSCANINCONCAT	= 1 << 8,       // Flag for last scan in concatination

		MDH_FLAG_RAWDATACORRECTION	= 1 << 10,      // Correct the rawadata with the rawdata correction factor
		MDH_FLAG_LASTSCANINMEAS		= 1 << 11,      // Flag for last scan in measurement
		MDH_FLAG_SCANSCALEFACTOR	= 1 << 12,      // Flag for scan specific additional scale factor
		MDH_FLAG_2NDHADAMARPULSE	= 1 << 13,      // 2nd RF exitation of HADAMAR
		MDH_FLAG_REFPHASESTABSCAN	= 1 << 14,      // reference phase stabilization scan         
		MDH_FLAG_PHASESTABSCAN		= 1 << 15,      // phase stabilization scan
		MDH_FLAG_D3FFT				= 1 << 16,      // execute 3D FFT         
		MDH_FLAG_SIGNREV			= 1 << 17,      // sign reversal
		MDH_FLAG_PHASEFFT			= 1 << 18,      // execute phase fft     
		MDH_FLAG_SWAPPED			= 1 << 19,      // swapped phase/readout direction
		MDH_FLAG_POSTSHAREDLINE		= 1 << 20,      // shared line               
		MDH_FLAG_PHASCOR			= 1 << 21,      // phase correction data    
		MDH_FLAG_PATREFSCAN			= 1 << 22,      // additonal scan for PAT reference line/partition
		MDH_FLAG_PATREFANDIMASCAN	= 1 << 23,      // additonal scan for PAT reference line/partition that is also used as image scan
		MDH_FLAG_REFLECT			= 1 << 24,      // reflect line              
		MDH_FLAG_NOISEADJSCAN		= 1 << 25,      // noise adjust scan --> Not used in NUM4        
		MDH_FLAG_SHARENOW			= 1 << 26,      // all lines are acquired from the actual and previous e.g. phases
		MDH_FLAG_LASTMEASUREDLINE	= 1 << 27,      // indicates that the current line is the last measured line of all succeeding e.g. phases
		MDH_FLAG_FIRSTSCANINSLICE	= 1 << 28,      // indicates first scan in slice (needed for time stamps)
		MDH_FLAG_LASTSCANINSLICE	= 1 << 29       // indicates  last scan in slice (needed for time stamps)

		
		//MDH_TREFFECTIVEBEGIN			= 1 << 30,      // indicates the begin time stamp for TReff (triggered measurement)
		//MDH_TREFFECTIVEEND				= 1 << 31,      // indicates the   end time stamp for TReff (triggered measurement)
		//MDH_MDS_REF_POSITION			= 1 << 32,      // indicates the reference position for move during scan images (must be set once per slice/partition in MDS mode)
		//MDH_SLC_AVERAGED				= 1 << 33,      // indicates avveraged slice for slice partial averaging scheme
        //
		//MDH_FIRST_SCAN_IN_BLADE			= 1 << 40,		// Marks the first line of a blade
		//MDH_LAST_SCAN_IN_BLADE			= 1 << 41,		// Marks the last line of a blade
		//MDH_LAST_BLADE_IN_TR			= 1 << 42,		// Set for all lines of the last BLADE in each TR interval
        //
		//MDH_RETRO_LASTPHASE				= 1 << 45,		// Marks the last phase in a heartbeat
		//MDH_RETRO_ENDOFMEAS				= 1 << 46,		// Marks an ADC at the end of the measurement
		//MDH_RETRO_REPEATTHISHEARTBEAT	= 1 << 47,		// Repeat the current heartbeat when this bit is found
		//MDH_RETRO_REPEATPREVHEARTBEAT	= 1 << 48,		// Repeat the previous heartbeat when this bit is found
		//MDH_RETRO_ABORTSCANNOW			= 1 << 49,		// Just abort everything
		//MDH_RETRO_LASTHEARTBEAT			= 1 << 50,		// This adc is from the last heartbeat (a dummy)
		//MDH_RETRO_DUMMYSCAN				= 1 << 51,		// This adc is just a dummy scan, throw it away
		//MDH_RETRO_ARRDETDISABLED		= 1 << 52		// Disable all arrhythmia detection when this bit is found
		
	};
	

	//*******************************************************************
	// HELPER CLASSES                                                    
	//******************************************************************* 

	// Axis properties in event block
	class Axis {
	public:
		Axis();

		Increment Increment_PE;
		Increment Increment_RO;
		Increment Increment_SL;

		Shape Shape_PE;
		Shape Shape_RO;
		Shape Shape_SL;
	};

	// Readout properties in event block
	class Readout {
	public:
		Readout();
		bool bIsValidScan;
		double fDwellTime;
		long lSamples;
	};

	// TX properties in event block
	class TX {
	public:
		TX();
		TxType txType;
		double fAsymmetry;
		double fDwellTime;
		long lSamples;
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

		void applyPhaseModulation(complex<double> *data, unsigned short numberOfSamples, unsigned int ulScanCounter);

		double getADCStartTime(unsigned int ulScanCounter);

		void writeToFile();

	protected:

		//*******************************************************************
		// FUNCTIONS                                                         
		//*******************************************************************

		// Open the text file
		void openFile();

		// Loops through the text file and counts the number of gradient and receiver events
		// Thereafter the text file is rewound till the first Halt gradient instruction
		void calcMemoryRequirement();

		// Runs all the gradient and receiver instructions 
		void runInstructions();

		// Calculate the exponential decay curves
		void computeExponentials();

		// Calculates the derivatives of the provided gradient
		// If the second argument is a Null pointer, the first input array will used as output array 
		void calculateDerivative(double *adData);

		// Calculates the derivatives for the gradient on each axis
		void calculateDerivative();

		void calculateIntegral();
		void calculateIntegral(double *afG, double *afS);

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

		// Process Control node
		void processControl(pugi::xml_node &node, vector<Readout> &vReadout, vector<TX> &vTx);

		// Process trigger node
		void processTrigger(pugi::xml_node &node, vector<Readout> &vReadout, vector<TX> &vTx);

		// Process rotation node
		void processRotation(pugi::xml_node &node);

		// Process Axis node
		void processAxis(pugi::xml_node &node, Axis *axis);

		// Process shape node
		void processShape(pugi::xml_node &node, Axis *axis);

		// Process Increment node
		void processIncrement(pugi::xml_node &node, Axis *axis);

		// Apply axis instructions
		void applyAxis(Axis *axis, double *GX, double *GY, double *GZ, long lTicks);

		// Read arbitrary gradient shapes
		void readGCShapes();


		//*******************************************************************
		// PROPERTIES                                                        
		//*******************************************************************

		// Text file name
		FileParts m_fpXMLFile;


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

		// In normal mode, no intermediate resuts are stored, but the array m_afMULTI_PURPOSE
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

		// XML document
		pugi::xml_document m_doc;

		bool m_bB0DataAvailable;

	};

}
