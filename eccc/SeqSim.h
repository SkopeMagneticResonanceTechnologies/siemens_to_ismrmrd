#include <vector>  
#include <complex>

#ifdef MATLAB_MEX_FILE 
	#include "mex.h"
#endif

#pragma once  

using namespace std;

namespace SEQSIM
{
	#ifndef ECC_COEFF
	#define ECC_COEFF
		struct ECC_Coeff
		{
			std::vector<double> vfTau;
			std::vector<double> vfAmp;	
		};
	#endif

	enum DataType {
		GRADIENT,
		KSPACE,
		SLEWRATE,
		EDDYCURRENT,
		EDDYPHASE
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

	enum OutputMode {
		FULL,
		INTERPOLATED_TO_RX
	};

	class DSP {
	public:

		DSP();
		~DSP();
		void setFileName(const char *pFileName);
		void setDSVFolderPath(const char *pFolderPath);
		#ifdef MATLAB_MEX_FILE 
			void run(mxArray *plhs[]);
		#else
			void run();
		#endif
		void setVerboseMode(int verbose);
		void setDataType(DataType type);
		void setOutputMode(OutputMode mode);
		void setB0CorrCoeff(ECC_Coeff &CoeffX, ECC_Coeff &CoeffY, ECC_Coeff &CoeffZ);
		void setB0CorrCoeff(vector<double> &vfCoeffXamp, vector<double> &vfCoeffXtau,
									vector<double> &vfCoeffYamp, vector<double> &vfCoeffYtau,
									vector<double> &vfCoeffZamp, vector<double> &vfCoeffZtau);
		void applyPhaseModulation(complex<float> *data, unsigned short numberOfSamples, unsigned int ulScanCounter);
		void writeToFile();
		void setDSVFileNamePrefix(const char * pDSVFileNamePrefix);

	};

}