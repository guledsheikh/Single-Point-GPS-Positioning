//***********************************************************************************
//Filename: reader.cpp
//
//Date				  Author				Revison History
//05/08/2017			  Guled Sheikh				version 1
//
//***********************************************************************************
//
//Copyright  2017~2018. All Rights Reserved
//
//***********************************************************************************

#include "reader.h"

//Reader constructor
reader::reader(const char * obsFilename, const char * navFilename)
{
	pFile = fopen(obsFilename, "rt");
	if (pFile == NULL) perror("Error opening file");
	NavpFile = fopen(navFilename, "rt");
	if (NavpFile == NULL) perror("Error opening file");
}

//This fucntions computes the gps time in seconds
double reader::GPSTime(vector<double> time)
{
	double julian_date, JDholder1, JDholder2, gpsweekholder, y, m, gpsweek, gpstime;

	double UT = time[3] + (time[4] / 60) + (time[5] / 3600);
	if (time[1] > 2) {
		y = time[0] + 2000;
		m = time[1];
	}
	else {
		y = time[0] + 2000 - 1;
		m = time[1] + 12;
	}
	JDholder1 = 365.25*y;
	JDholder2 = 30.6001*(m + 1);
	if (JDholder1 > 0 && JDholder2 > 0) {
		julian_date = floor(JDholder1) + floor(JDholder2) + time[2] + (UT / 24) + 1720981.5;
	}
	else {
		julian_date = ceil(JDholder1) + ceil(JDholder2) + time[2] + (UT / 24) + 1720981.5;
	}
	gpsweekholder = (julian_date - 2444244.5) / 7;

	if (gpsweekholder > 0) {
		gpsweek = floor(gpsweekholder);
	}
	else {
		gpsweek = ceil(gpsweekholder);
	}
	gpstime = round((((julian_date - 2444244.5) / 7) - gpsweek) * 7 * 24 * 3600);

	return gpstime;
}

//This functions converts scientific notations to double
double reader::scientificTodouble(string str)
{
	istringstream os(str.replace(15, 1, "e")); //replace D with e in rinex Nav file
	double str1;
	os >> str1;
	return str1;
}

//RINEX OBSERVATION FILE READER
void reader::rinexOBSReader(Header &head, RECORDS &obs_RECORDS, int &SENTINEL)
{
	//Character line for storing line by line characters
	char line[MAX_STRLEN];
	
	//Read the header.......

	if (obsHeader == true) {
		fgets(line, MAX_STRLEN, pFile);
		while (!feof(pFile) && strstr(line, "END OF HEADER") == NULL) {
			if (sizeof(line) >= 61) {
				if (strstr(line, "RINEX VERSION")) {
					string sline(line);  // convert the line (char) to string class type sLine
					string token = sline.substr(5, 3);   // from 5th character and five spaces
					head.version = token;
				}
				if (strstr(line, "MARKER NAME")) {
					string sline(line);  // convert the line (char) to string class type sLine
					string token = sline.substr(0, 10);   // from 5th character and five spaces
					head.markername = token;
				}
				if (strstr(line, "INTERVAL")) {
					string sline(line);  // convert the line (char) to string class type sLine
					string token = sline.substr(5, 5);   // from 5th character and five spaces
					head.interval = stod(token);
				}
				if (strstr(line, "APPROX POSITION XYZ")) {
					string sline(line);  // convert the line (char) to string class type sLine
					string token = sline.substr(3, 56);   // from 3rd character upto 56
					char cline[MAX_BYTE];     // cline is a char line
					strncpy(cline, token.c_str(), sizeof(cline));
					cline[sizeof(cline) - 1] = 0;
					for (char * r = strtok(cline, " "); r; r = strtok(NULL, " "))
					{
						double rr = atof(r);
						head.antpos.push_back(rr);
					}
				}
				if (strstr(line, "# / TYPES OF OBSERV")) {
					string sline_new(line);  // convert the line (char) to string class type sLine
					string token_new = sline_new.substr(5, 2);   // from 5th character and two spaces
					head.nobstype = stoi(token_new);
					if (head.nobstype <= 9)
					{
						string sline(line);  // convert the line (char) to string class type sLine
						string token = sline.substr(10, 50);   // from 10th character upto #
						char cline[MAX_BYTE];     // cline is a char line
						strncpy(cline, token.c_str(), sizeof(cline));
						cline[sizeof(cline) - 1] = 0;

						for (char * p = strtok(cline, " "); p; p = strtok(NULL, " "))
						{
							head.obstype.push_back(p);
						}
					}
					else {
						string sline(line);  // convert the line (char) to string class type sLine
						string token = sline.substr(10, 50);   // from 10th character upto #
						char cline[MAX_BYTE];     // cline is a char line
						strncpy(cline, token.c_str(), sizeof(cline));
						cline[sizeof(cline) - 1] = 0;

						for (char * p = strtok(cline, " "); p; p = strtok(NULL, " "))
						{
							head.obstype.push_back(p);
						}
						//......... Get a new line.......................
						fgets(line, MAX_STRLEN, pFile);
						string sline2(line);  // convert the line (char) to string class type sLine
						string token2 = sline.substr(10, 50);   // from 10th character upto #
						char cline2[MAX_BYTE];     // cline is a char line
						strncpy(cline2, token2.c_str(), sizeof(cline2));
						cline2[sizeof(cline2) - 1] = 0;
						for (char * pp = strtok(cline2, " "); pp; pp = strtok(NULL, " "))
						{
							head.obstype.push_back(pp);
						}
					}
				}
			}
			//......... Get a new line.......................
			fgets(line, MAX_STRLEN, pFile);
		}
	}
	obsHeader = false;


	//Read the observation file......
	int Nobs = head.nobstype;
	int Num_Sat;
	
	//---------------------------------------------------------------------------------//
	while (!feof(pFile) && SENTINEL != -1)
	{
		fgets(line, MAX_STRLEN, pFile);  // Get a newline

		if (strstr(line, "G") && !strstr(line, "COMMENT"))
		{
			double Y, M, D, H, Min, Sec, timeStamp, timeSEC;
			RECORDS obs_temp;
			std::vector< int > SV;
			std::vector< string > tempEpoch;

			string sline(line);  // convert the line (char) to string class type sLine
			string token = sline.substr(0, 32);   // first line of each epoch from start upto G
			char epoch[MAX_BYTE];     // cline is a char line rep the starting line of each epoch
			strncpy(epoch, token.c_str(), sizeof(epoch));
			epoch[sizeof(epoch) - 1] = 0;

			for (char * p = strtok(epoch, " "); p; p = strtok(NULL, " "))
			{
				tempEpoch.push_back(p);
			}

			//Get the time from the file.....
			Y = stod(tempEpoch[0]);
			M = stod(tempEpoch[1]);
			D = stod(tempEpoch[2]);
			H = stod(tempEpoch[3]);
			Min = stod(tempEpoch[4]);
			Sec = stod(tempEpoch[5]);


			vector<double> timeHolder;

			timeHolder.push_back(Y);
			timeHolder.push_back(M);
			timeHolder.push_back(D);
			timeHolder.push_back(H);
			timeHolder.push_back(Min);
			timeHolder.push_back(Sec);

			timeStamp = GPSTime(timeHolder);
			timeSEC = (H * 3600) + (Min * 60) + (Sec);

			int Num_Sat_i, L_To_SN, S_Of_Str;
			Num_Sat = 0;                    //Initialize a variable to hold the number of sat at each epoch
			Num_Sat_i = stoi(tempEpoch[7]); //initial num of satellite
			L_To_SN = 33;                   //Length to satellite number
			S_Of_Str = 0;					//Start of string

			//..................................................................................//
			//     Loop for getting the satellite vehicles in each epoch and num of sat
			//.................................................................................//

			string sline1(line);  // convert the line (char) to string class type sLine
			string token1 = sline1.substr(L_To_SN, sline1.length());

			//Prn storage
			for (int r = 0; r < Num_Sat_i; r++) {
				string token2 = token1.substr(S_Of_Str, 2);
				int p = stoi(token2);
				SV.push_back(p);
				Num_Sat = Num_Sat + 1;
				S_Of_Str = S_Of_Str + 3;
			}
			obs_RECORDS.NumOfSat = Num_Sat;
			obs_RECORDS.prn = SV;
			obs_RECORDS.time = timeStamp;
			obs_RECORDS.YEAR = Y;
			obs_RECORDS.MONTH = M;
			obs_RECORDS.DAY = D;
			obs_RECORDS.HOUR = H;
			obs_RECORDS.MIN = Min;
			obs_RECORDS.Indivi_Sec = Sec;    // individual second
			obs_RECORDS.SEC = timeSEC;       // Total time in seconds
			//------------------------------------------------------------------------------------------------------------------------------------------------------

			double con = floor(Nobs / 5);   // Condition for looping according to continuation lines if past 5 obs in a line

			vector< vector<double> > array;

			for (int i = 0; i < Num_Sat; i++)
			{
				int curObsType = 1;		// current observation type
				size_t obsSpace = 14;		// space of an observation in a line
				std::vector< double > recVec;	// Get a newline

				fgets(line, MAX_STRLEN, pFile);

				// Each observation line contains a max of 5 records.... //So if less or equal(we have less or equal to 5 obs types)

				if (Nobs <= 5)
				{
					for (int j = 0; j < Nobs + 1; j++) {
						string sline2(line);	// convert the line (char) to string class type sLine
						size_t S_Of_Rec = 0;
						if (strlen(line) > 1) {
							string token3 = sline2.substr(S_Of_Rec, 14);
							double p = stod(token3);
							recVec.push_back(p);
							S_Of_Rec = S_Of_Rec + 16;
						}
						else {
							//(An empty line)
							double p = NAN;
							recVec.push_back(p);
							S_Of_Rec = S_Of_Rec + 16;
						}
					}
				}
				else
				{
					for (int k = 0; k < con; k++)
					{
						string sline2(line);	// convert the line (char) to string class type sLine
						size_t S_Of_Rec = 0;
						double count = 0;
						double strlen_sline2 = sline2.length();
						for (int l = 0; l < 5; l++) {  //MAX NUMBER OF OBS PER LINE = 5	
							count = count + 1;
							//****************************ADDED 8/26/2018***************************************//
							if (strlen(line) > 1) {
								if (count * 14 > strlen_sline2) {
									double p = 0.0;
									recVec.push_back(p);
									S_Of_Rec = S_Of_Rec + 16;
								}
								else {
									string token3 = sline2.substr(S_Of_Rec, 14);
									//****************************ADDED 8/26/2018***************************************//
									double p;
									bool whiteSpacesOnly = std::all_of(token3.begin(), token3.end(), isspace);
									if (whiteSpacesOnly) {
										p = 0.0;
									}
									else {
										p = stod(token3);
									}
									//****************************CHANGED**************************************//							
									recVec.push_back(p);
									S_Of_Rec = S_Of_Rec + 16;
								}
							}
							else {
								//(An empty line)
								double p = 0.0;
								recVec.push_back(p);
								S_Of_Rec = S_Of_Rec + 16;
							}
						}

						//..........Continuation line...........//
						fgets(line, MAX_STRLEN, pFile);
					}
					string sline3(line);	// convert the line (char) to string class type sLine
					size_t S_Of_Rec = 0;
					double count = 0;
					double strlen_sline3 = sline3.length();

					for (int o = 0; o < (Nobs - 5 * con); o++)   //MAX NUMBER OF OBS PER LINE = 5	
					{
						count = count + 1;
						if (strlen(line) > 1)
						{
							if (count * 14 > strlen_sline3) {
								double p = 0.0;
								recVec.push_back(p);
								S_Of_Rec = S_Of_Rec + 16;
							}
							else {
								string token4 = sline3.substr(S_Of_Rec, 14);
								double p;
								bool whiteSpacesOnly = std::all_of(token4.begin(), token4.end(), isspace);
								if (whiteSpacesOnly) {
									p = 0.0;
								}
								else {
									p = stod(token4);
								}
								recVec.push_back(p);
								S_Of_Rec = S_Of_Rec + 16;
							}
						}
						else
						{  //An empty line
							double p = 0.0;
							recVec.push_back(p);
							S_Of_Rec = S_Of_Rec + 16;
						}
					}
					// Put recVec in array double vector
					array.push_back(recVec);
					obs_RECORDS.record = array;
				}
			}
			array.clear();
			SENTINEL = -1;
		}
		else {
			//Do nothing
		}
	}
}

//RINEX NAVIGATION FILE READER
void reader::rinexNAVReader(vector<NAVRECORDS>& nav_RECORDS)
{
	char aline[MAX_STRLEN];

	// SKIP THE HEADER......
	if (navHeader == true) {
		fgets(aline, MAX_STRLEN, NavpFile);
		while (!feof(NavpFile) && strstr(aline, "END OF HEADER") == NULL) {
			if (strstr(aline, "END OF HEADER")) {
				break;
			}
			fgets(aline, MAX_STRLEN, NavpFile);
		}
	}
	navHeader = false;

	// Read the navigation file...
	//Get a new line......
	fgets(aline, MAX_STRLEN, NavpFile);

	while (!feof(NavpFile))
	{
		double prn, Y, M, D, H, Min, Sec, timeStamp, timeSEC;
		NAVRECORDS navrecord_temp;

		// First line of the nav records at each epoch
		string sline(aline); // convert the line (char) to string class type sLine
		
		//Store the time in the file
		string prn_s = sline.substr(0, 2);     prn = stod(prn_s);
		string Y_s = sline.substr(3, 2);       Y = stod(Y_s);
		string M_s = sline.substr(6, 2);	   M = stod(M_s);
		string D_s = sline.substr(9, 2);	   D = stod(D_s);
		string H_s = sline.substr(12, 2);	   H = stod(H_s);
		string Min_s = sline.substr(15, 2);    Min = stod(Min_s);
		string Sec_s = sline.substr(18, 4);    Sec = stod(Sec_s);

		vector<double> timeHolder;

		timeHolder.push_back(Y);
		timeHolder.push_back(M);
		timeHolder.push_back(D);
		timeHolder.push_back(H);
		timeHolder.push_back(Min);
		timeHolder.push_back(Sec);
		
		//Get gps time in seconds....
		timeStamp = GPSTime(timeHolder);
		timeSEC = (H * 3600) + (Min * 60) + (Sec);

		string SVCB_s = sline.substr(22, 19);
		double SVCB = scientificTodouble(SVCB_s);  // Convert scientific string in Nav to double (and replace D with e)

		string SVCD_s = sline.substr(41, 19);
		double SVCD = scientificTodouble(SVCD_s);

		string SVCDR_s = sline.substr(60, 19);
		double SVCDR = scientificTodouble(SVCDR_s);
		
		//Get a new line.......
		fgets(aline, MAX_STRLEN, NavpFile);
		vector<double> navRecHolder;

		for (int i = 0; i < 6; i++) {
			string sline(aline);
			size_t S_Of_Rec = 3;
			for (int j = 0; j < 4; j++) {
				string p_s = sline.substr(S_Of_Rec, 19);
				double p = scientificTodouble(p_s);
				navRecHolder.push_back(p);
				S_Of_Rec = S_Of_Rec + 19;
			}
			fgets(aline, MAX_STRLEN, NavpFile);
		}
		string sline6(aline);
		if (sline6.length() <= 23) {
			string TransmissionTime_s = sline6.substr(3, 19);
			double TransmissionTime = scientificTodouble(TransmissionTime_s);
			navrecord_temp.TransmissionTime = TransmissionTime;
			navrecord_temp.FitInterval = NAN;
		}
		else {
			string TransmissionTime_s = sline6.substr(3, 19);
			double TransmissionTime = scientificTodouble(TransmissionTime_s);
			string FitInterval_s = sline6.substr(3, 19);
			double FitInterval = scientificTodouble(FitInterval_s);

			navrecord_temp.TransmissionTime = TransmissionTime;
			navrecord_temp.FitInterval = FitInterval;
		}

		// Store all the nav obs in the local struct
		navrecord_temp.PRN = prn;
		navrecord_temp.SVClockBias = SVCB;
		navrecord_temp.SVClockDrift = SVCD;
		navrecord_temp.SVClockDriftRate = SVCDR;
		navrecord_temp.SEC = timeSEC;
		navrecord_temp.gpstime = timeStamp;
		navrecord_temp.IODE = navRecHolder[0];
		navrecord_temp.Crs = navRecHolder[1];
		navrecord_temp.DeltaN = navRecHolder[2];
		navrecord_temp.Mo = navRecHolder[3];
		navrecord_temp.Cuc = navRecHolder[4];
		navrecord_temp.Eccentricity = navRecHolder[5];
		navrecord_temp.Cus = navRecHolder[6];
		navrecord_temp.Sqrta = navRecHolder[7];
		navrecord_temp.Toe = navRecHolder[8];
		navrecord_temp.Cic = navRecHolder[9];
		navrecord_temp.OMEGA = navRecHolder[10];
		navrecord_temp.Cis = navRecHolder[11];
		navrecord_temp.Io = navRecHolder[12];
		navrecord_temp.Crc = navRecHolder[13];
		navrecord_temp.omega = navRecHolder[14];
		navrecord_temp.OMEGADOT = navRecHolder[15];
		navrecord_temp.IDOT = navRecHolder[16];
		navrecord_temp.L2CodesChannel = navRecHolder[17];
		navrecord_temp.GPSWeek = navRecHolder[18];
		navrecord_temp.L2PDataFlag = navRecHolder[19];
		navrecord_temp.SVAccuracy = navRecHolder[20];
		navrecord_temp.SVHealth = navRecHolder[21];
		navrecord_temp.TGD = navRecHolder[22];
		navrecord_temp.IODC = navRecHolder[23];

		// Store everything in the vector of global struct
		nav_RECORDS.push_back(navrecord_temp);
		
		//Get a new line
		fgets(aline, MAX_STRLEN, NavpFile);
	}
}
