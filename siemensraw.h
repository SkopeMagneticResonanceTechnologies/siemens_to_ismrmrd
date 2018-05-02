#ifndef SIEMENSRAW_H_
#define SIEMENSRAW_H_

#include <stdint.h>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#define MDH_NUMBEROFEVALINFOMASK   2

#define MDH_NUMBEROFICEPROGRAMPARA_VB 4
#define MDH_NUMBEROFICEPROGRAMPARA_VD 24

#define MDH_FREEHDRPARA_VB 4

#define MDH_DMA_LENGTH_MASK   (0x01FFFFFFL)
#define MDH_PACK_BIT_MASK     (0x02000000L)
#define MDH_ENABLE_FLAGS_MASK (0xFC000000L)

#define MDH_ACQEND				(1)
#define MDH_RTFEEDBACK			(1 << 1)
#define MDH_HPFEEDBACK			(1 << 2)
#define MDH_ONLINE				(1 << 3)
#define MDH_OFFLINE				(1 << 4)
#define MDH_SYNCDATA			(1 << 5)       // readout contains synchronous data
#define MDH_LASTSCANINCONCAT	(1 << 8)       // Flag for last scan in concatenation

#define MDH_RAWDATACORRECTION	(1 << 10)      // Correct the rawdata with the rawdata correction factor
#define MDH_LASTSCANINMEAS		(1 << 11)      // Flag for last scan in measurement
#define MDH_SCANSCALEFACTOR		(1 << 12)      // Flag for scan specific additional scale factor
#define MDH_2NDHADAMARPULSE		(1 << 13)      // 2nd RF excitation of HADAMAR
#define MDH_REFPHASESTABSCAN	(1 << 14)      // reference phase stabilization scan         
#define MDH_PHASESTABSCAN		(1 << 15)      // phase stabilization scan
#define MDH_D3FFT				(1 << 16)      // execute 3D FFT         
#define MDH_SIGNREV				(1 << 17)      // sign reversal
#define MDH_PHASEFFT			(1 << 18)      // execute phase fft     
#define MDH_SWAPPED				(1 << 19)      // swapped phase/readout direction
#define MDH_POSTSHAREDLINE		(1 << 20)      // shared line               
#define MDH_PHASCOR				(1 << 21)      // phase correction data    
#define MDH_PATREFSCAN			(1 << 22)      // additional scan for PAT reference line/partition
#define MDH_PATREFANDIMASCAN	(1 << 23)      // additional scan for PAT reference line/partition that is also used as image scan
#define MDH_REFLECT				(1 << 24)      // reflect line              
#define MDH_NOISEADJSCAN		(1 << 25)      // noise adjust scan --> Not used in NUM4        
#define MDH_SHARENOW			(1 << 26)      // all lines are acquired from the actual and previous e.g. phases
#define MDH_LASTMEASUREDLINE	(1 << 27)      // indicates that the current line is the last measured line of all succeeding e.g. phases
#define MDH_FIRSTSCANINSLICE	(1 << 28)      // indicates first scan in slice (needed for time stamps)
#define MDH_LASTSCANINSLICE		(1 << 29)      // indicates  last scan in slice (needed for time stamps)

enum class PMU_Type {
	END =  0x01FF0000,
	ECG1 = 0x01010000,
	ECG2 = 0x01020000,
	ECG3 = 0x01030000,
	ECG4 = 0x01040000,
	PULS = 0x01050000,
	RESP = 0x01060000,
	EXT1 = 0x01070000,
	EXT2 = 0x01080000

};

enum class Trajectory {
  TRAJECTORY_CARTESIAN = 0x01,
  TRAJECTORY_RADIAL    = 0x02,
  TRAJECTORY_SPIRAL    = 0x04,
  TRAJECTORY_BLADE     = 0x08
};

struct mdhLC {
  uint16_t ushLine;
  uint16_t ushAcquisition;
  uint16_t ushSlice;
  uint16_t ushPartition;
  uint16_t ushEcho;
  uint16_t ushPhase;
  uint16_t ushRepetition;
  uint16_t ushSet;
  uint16_t ushSeg;
  uint16_t ushIda;
  uint16_t ushIdb;
  uint16_t ushIdc;
  uint16_t ushIdd;
  uint16_t ushIde;
};

struct mdhCutOff {
  uint16_t ushPre;
  uint16_t ushPost;
};

struct mdhSlicePosVec {
  float flSag;
  float flCor;
  float flTra;
};

/*
struct mdhSD {
  mdhSlicePosVec sSlicePosVec;
};

*/
struct PMUdata {
	uint16_t data;
	uint16_t trigger;
};

struct mdhSliceData {
  mdhSlicePosVec sSlicePosVec;
  float aflQuaternion[4];
};

/* This is the VB line header */
struct sMDH
{
  uint32_t ulFlagsAndDMALength;
  int32_t lMeasUID;
  uint32_t ulScanCounter;
  uint32_t ulTimeStamp;
  uint32_t ulPMUTimeStamp;
  uint32_t aulEvalInfoMask[2];
  uint16_t ushSamplesInScan;
  uint16_t ushUsedChannels;
  mdhLC sLC;
  mdhCutOff sCutOff;

  uint16_t ushKSpaceCentreColumn;
  uint16_t ushCoilSelect;
  float fReadOutOffcentre;
  uint32_t ulTimeSinceLastRF;
  uint16_t ushKSpaceCentreLineNo;
  uint16_t ushKSpaceCentrePartitionNo;
  uint16_t aushIceProgramPara[4];
  uint16_t aushFreePara[4];

  mdhSliceData sSliceData;

  uint16_t ushChannelId;
  uint16_t ushPTABPosNeg;
};

/* This is the VD line header */
struct sScanHeader
{
	  uint32_t ulFlagsAndDMALength;
	  int32_t lMeasUID;
	  uint32_t ulScanCounter;
	  uint32_t ulTimeStamp;
	  uint32_t ulPMUTimeStamp;
	  uint16_t ushSystemType;
	  uint16_t ulPTABPosDelay;
	  int32_t lPTABPosX;
	  int32_t lPTABPosY;
	  int32_t lPTABPosZ;
	  uint32_t ulReserved1;
	  uint32_t aulEvalInfoMask[2];
	  uint16_t ushSamplesInScan;
	  uint16_t ushUsedChannels;
	  mdhLC sLC;
	  mdhCutOff sCutOff;
	  uint16_t ushKSpaceCentreColumn;
	  uint16_t ushCoilSelect;
	  float fReadOutOffcentre;
	  uint32_t ulTimeSinceLastRF;
	  uint16_t ushKSpaceCentreLineNo;
	  uint16_t ushKSpaceCentrePartitionNo;
	  mdhSliceData sSliceData;
	  uint16_t aushIceProgramPara[24];
	  uint16_t aushReservedPara[4];
	  uint16_t ushApplicationCounter;
	  uint16_t ushApplicationMask;
	  uint32_t ulCRC;
};

typedef struct sChannelHeader
{
  uint32_t ulTypeAndChannelLength;
  int32_t lMeasUID;
  uint32_t ulScanCounter;
  uint32_t ulReserved1;
  uint32_t ulSequenceTime;
  uint32_t ulUnused2;
  uint16_t ulChannelId;
  uint16_t ulUnused3;
  uint32_t ulCRC;
} sChannelHeader;

struct MrParcRaidFileEntry
{
  uint32_t measId_;
  uint32_t fileId_;
  uint64_t off_;
  uint64_t len_;
  char patName_[64];
  char protName_[64];
};

struct MrParcRaidFileHeader
{
  uint32_t hdSize_;
  uint32_t count_;
};


typedef struct
{
  sMDH mdh;
  void* previous;
  void* next;
  float* data;
} SiemensMdhNode;

typedef struct
{
  uint32_t matrix_size[3];
  uint32_t pat_matrix_size[3];
  uint32_t base_resolution;
  uint32_t phase_encoding_lines;
  uint32_t partitions;
  uint32_t dimensions;
  float    phase_resolution;
  float    slice_resolution;
  float    dwell_time_us;
  uint32_t acceleration_factor_pe;
  uint32_t acceleration_factor_3d;
} SiemensBaseParameters;

class SiemensRawData
{
public:
  SiemensRawData();
  ~SiemensRawData();

  int ReadRawFile(char* filename);
  long GetNumberOfNodes();
  SiemensMdhNode* GetFirstNode();
  sMDH* GetMinValues();
  sMDH* GetMaxValues();

  SiemensBaseParameters GetBaseParameters() {return m_base_parameters;}
  
  int GetMeasYapsParameter(std::string parameter_name, std::string& value);

  const std::string& GetParameterBuffer(std::string name);

protected:
  int ReadMdhNode(std::ifstream* f);
  int DeleteNode(SiemensMdhNode* node);
  int DeleteNodeList();
  int UpdateMinMax();
  int ParseMeasYaps();

  SiemensMdhNode* m_first_node;
  SiemensMdhNode* m_last_node;
  long            m_number_of_nodes;
  sMDH            m_mdh_min;
  sMDH            m_mdh_max;
  bool            m_min_max_is_valid;

  std::map< std::string, std::string >  m_parameter_buffers;
  std::map< std::string, std::string >  m_meas_yaps;
  SiemensBaseParameters  m_base_parameters;
};

#endif //SIEMENSRAW_H_
