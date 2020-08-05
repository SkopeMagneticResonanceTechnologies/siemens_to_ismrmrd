#include "SeqSim.h"		

using namespace std;

namespace SEQSIM
{
	DSP::DSP(){}
	DSP::~DSP(){}
	void DSP::setFileName(const char * pFileName){}
	void DSP::setDSVFolderPath(const char * pFolderPath) {}
	void DSP::setVerboseMode(int verbose) {}
	void DSP::setDataType(DataType type) {}
	void DSP::setOutputMode(OutputMode mode) {}
	void DSP::setB0CorrCoeff(ECC_Coeff &CoeffX, ECC_Coeff &CoeffY, ECC_Coeff &CoeffZ) {}
	void DSP::setB0CorrCoeff(	vector<double> &vfCoeffXamp, vector<double> &vfCoeffXtau,
								vector<double> &vfCoeffYamp, vector<double> &vfCoeffYtau,
								vector<double> &vfCoeffZamp, vector<double> &vfCoeffZtau     ) {
	}
	void DSP::writeToFile() {}
	void DSP::applyPhaseModulation(complex<float> *data, unsigned short numberOfSamples, unsigned int ulScanCounter) {}

	#ifdef MATLAB_MEX_FILE
	 void DSP::run(mxArray *plhs[])
	#else
	 void DSP::run() 
	#endif
	{}

}
