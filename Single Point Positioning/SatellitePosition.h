//***********************************************************************************
//Filename: SatellitePosition.h
//
//Date					  Author					Revison History
//05/08/2017			  Guled Sheikh				version 1
//
//***********************************************************************************
//
//Copyright  2017~2018. All Rights Reserved
//
//***********************************************************************************
#ifndef SATELLITEPOSITION_H
#define SATELLITEPOSITION_H
#define _CRT_SECURE_NO_DEPRECATE

#include "reader.h"
#define FILE_SIZE 200

class SatellitePosition {
public:


	reader * RinexFile; // CREATE A READER POINTER
	
	//Observation filename holder
	char obsFile[FILE_SIZE];
	//Navigation filename holder
	char navFile[FILE_SIZE];
	
	double epoch_count = 0.0;

	///ofstream myfile;
	FILE * Env_File;  // File pointer for observation

	string Filename1_Obs;
	string Filename1_Nav;
	string Filename2_Obs;
	string Filename2_Nav;
	string Filename3_Obs;

	//+++++++ENVI_FILE STRUCTS+++++++++++++++//
	struct ENVI_FILE {
		string version;
		vector< double > base_final_coord;
	};
	ENVI_FILE envi_records;
	void OpenEnvFile(ENVI_FILE& envi_records, SatellitePosition& object1, SatellitePosition & object2);

	struct SATPOS_struct {

		//Variables for Least Squares
		vector< vector<double> > Final_Epoch_XYZt;		// For epoch by epoch
		vector<double> Epoch_XYZt;				// For epoch by epoch
		vector<double> Epoch_XYZt_ave;				// Average coordinates of each epoch
		vector<double> Sat_elev;

		// Variables for Satellite Positioning function
		vector< vector<double> > pos_rot_epoch;
		vector< vector<double> > Initial_SatPos;
		vector <double> sat_offset_epoch;
		vector <double> pr_epoch_iono;
		vector <double> pr_epoch;
		vector <double> P1_epoch;
		vector <double> P2_epoch;
		vector <double> C1_epoch;
		vector <double> L1_epoch;
		vector <double> L2_epoch;
		double TROPO;
	};
	SATPOS_struct epoch_data;

	/* FUNCTIONS*/
	void SetFilenames(const char* ObsFile, const char* NavFile);
	double EucldNorm(double x, double y, double z);   // EUCLIDEAN DISTANCE
	vector <double>  EarthRot_Corr(double time, vector <double> pos);   // EARTH ROTATION CORRECTION
	double Tropo_Corr(vector <double> pos, vector <double> ini_pos);    // TROPOSPHERE CORRECTION
	double iono_free(double count);
	double klobuchar(double fi, double lambda, double elev, double azimuth, double tow);
	double gpstimeCalc(double year, double month, double day, double hour, double min, double sec);    // TROPOSPHERE CORRECTION
	double NMF_hyd(double sat_elev);    // TROPOSPHERE CORRECTION NIELL HYDROSTATIC
	void Sat_Pos(SATPOS_struct &EPOCH_DATA);  //SATELLITE POSITION
};




#endif
