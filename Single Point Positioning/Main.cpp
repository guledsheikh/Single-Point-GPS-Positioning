//***********************************************************************************
//Filename: Main.CPP
//
//Date			          Author			        Revison History
//01/01/2018			  Guled Sheikh				version 1
//
//***********************************************************************************
//
//Copyright  2017~2018. All Rights Reserved
//
//***********************************************************************************

#include "reader.h"
#include "SatellitePosition.h"
#include "SPP.h"

int main(int argc, const char * argv[]) {
	int start_s = clock();

	// Create object of each module
	SPP CSPP; //Single point positioning static

	// Key words
	int arg = 1;

	//Single point positioning static
	if (arg == 1) {
		CSPP.Synthesize(CSPP.CSatellitePosition);
	}
	
	//YOU CAN ADD OTHER POSITIONING MODES HERE!!!!!
	// Single point positioning kinematic
	if (arg == 2) {}

	// Relative positioning static
	if (arg == 3) {}

	// Relative positioning kinematic
	if (arg == 4) {}

	// GPS and IMU integrated positioning
	if (arg == 5) {}
	
	int stop_s = clock();
	cout << "time: " << (stop_s - start_s) / double(CLOCKS_PER_SEC) << endd;
	cout << "Complete!!!" << endd;
	cin.get();
}


