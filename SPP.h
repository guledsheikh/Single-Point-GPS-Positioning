#ifndef SPP_H
#define SPP_H
#define _CRT_SECURE_NO_DEPRECATE

#include "SatellitePosition.h"

class SPP {
public:

	// Create an object of satelliteposition 
	SatellitePosition CSatellitePosition;
	
	 // File pointer to an environment file
	FILE * Env_File; 
	
	// Filename holders
	string Filename1_Obs;
	string Filename1_Nav;
	string Filename2_Obs;
	string Filename2_Nav;

	// Flags
	bool epochflag = true;
	bool initial_flag = true;

	// Variables
	double SUM_X = 0;
	double SUM_Y = 0;
	double SUM_Z = 0;
	double AVE_X = 0;
	double AVE_Y = 0;
	double AVE_Z = 0;
	double epoch = 0;
	double SIGMA_SQRD;
	double SIGMA_X;
	double SIGMA_Y;
	double SIGMA_Z;
	double SIGMA_t;
	double SIGMA_SQRD_epoch;
	VectorXf Residual;
	VectorXf Residual_std;
	MatrixXf Prev_D_xx;
	
	
	//FUNCTIONS
	void MSGFILE(char* outfilename, char* residual, SatellitePosition& object); // Basic output file function
	vector <double> Initialize_Rec_pos(SatellitePosition& object); // to initialize the receiver position
	void OpenEnvFile(); //Opens the observation and navigation files
	void Synthesize(SatellitePosition& object);  // Function that performs the single point positioning static
};




#endif
