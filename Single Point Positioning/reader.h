//***********************************************************************************
//Filename: reader.h
//
//Date				  Author				Revison History
//05/08/2017			  Guled Sheikh				version 1
//
//***********************************************************************************
//
//Copyright  2017~2018. All Rights Reserved
//
//***********************************************************************************

#ifndef READER_H
#define READER_H
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>
#include <cstddef>
#include <cstring>
#include <array>
#include <map>
#include <utility>
#include <iomanip>      // setprecision
#include <math.h>       /* pow */
#include <Eigen/Dense>
#include <Eigen/Core>

#define endd '\n'
#define MAX_STRLEN 100
#define MAX_BYTE 1024
#define PI 3.14159265358979323846264
using namespace std;
using namespace Eigen;

class reader {
public:

	//Contructors
	reader(const char *obsFilename, const char *navFilename);
	
	//File pointers for observation and navigation file respectively
	FILE * pFile;  // File pointer for observation
	FILE* NavpFile;  // File pointer for navigation
	
	//Flags
	bool obsHeader = true;
	bool navHeader = true;

	int SENTINEL = 0;
	
	//Observation file structs
	//Holds the observation file header information
	struct Header {
		string version;
		string markername;
		int nobstype;
		vector< string > obstype;
		double interval;
		vector< double > antpos;
	};
	Header head;
	
	//Holds observation file records
	struct RECORDS {
		vector< vector<double> > record;
		double time;
		double YEAR;
		double MONTH;
		double DAY;
		double HOUR;
		double MIN;
		double Indivi_Sec;
		double SEC;
		vector<int>  prn;
		int NumOfSat;
	};
	RECORDS obs;

	//Navigation file struct
	//Holds the navigation file records
	struct NAVRECORDS {
		double IODE, Crs, DeltaN, Mo, Cuc, Eccentricity, Cus, Sqrta, Toe, Cic,
			OMEGA, Cis, Io, Crc, omega, OMEGADOT, IDOT, L2CodesChannel, GPSWeek, L2PDataFlag, SVAccuracy, SVHealth, TGD, IODC,
			TransmissionTime, FitInterval, SVClockBias, SVClockDrift, SVClockDriftRate, PRN, SEC, gpstime;
	};
	vector<NAVRECORDS> SV_data;


	//+++++++FUNCTIONS+++++++++++++++//
	double GPSTime(vector<double> time);    //Function for converting time into gps seconds
	double scientificTodouble(string str);  //Function for converting scientific string notation to double
	void rinexOBSReader(Header &head, RECORDS &obs_RECORDS, int &SENTINEL);
	void rinexNAVReader(vector<NAVRECORDS>& nav_RECORDS);
	void IMUReader(IMU_RECORDS &imu_obs, int &IMUSENTINEL);
};

#endif

