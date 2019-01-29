//***********************************************************************************
//Filename: SatellitePosition.cpp
//
//Date					  Author					Revison History
//05/08/2017			  Guled Sheikh				version 1
//
//***********************************************************************************
//
//Copyright  2017~2018. All Rights Reserved
//
//***********************************************************************************
#include "SatellitePosition.h"

//Computes the Euclidean distance
double SatellitePosition::EucldNorm(double x, double y, double z)
{
	double norm = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	return norm;
}

// Fucntion that returns satellite positions corrected for Earth rotation errors
vector<double> SatellitePosition::EarthRot_Corr(double time, vector<double> pos)
{
	/*...Returns rotated satellite ECEF coordinates due to Earth rotation during signal travel time
	/*....time represents the travel time...*/
	/*....WGS 84 value of the earth's rotation rate (rad/sec)...*/
	double OMEGA_e = 7.2921151467e-5;
	double OMEGA_tau = OMEGA_e*time;
	Eigen::Matrix3d R3;
	Eigen::Vector3d POS_ROT;
	Eigen::Vector3d POS(pos[0], pos[1], pos[2]);
	vector <double> POSITION;

	double R11, R12, R13, R21, R22, R23, R31, R32, R33;
	R11 = cos(OMEGA_tau);   R21 = -sin(OMEGA_tau);  R31 = 0;
	R12 = sin(OMEGA_tau);   R22 = cos(OMEGA_tau);	R32 = 0;
	R13 = 0;				R23 = 0;				R33 = 1;

	R3 << R11, R12, R13,
		R21, R22, R23,
		R31, R32, R33;

	POS_ROT = R3*POS;

	POSITION.push_back(POS_ROT[0]);
	POSITION.push_back(POS_ROT[1]);
	POSITION.push_back(POS_ROT[2]);

	return POSITION;
}

// Function that corrects for tropospheric effects
double SatellitePosition::Tropo_Corr(vector<double> pos, vector<double> ini_pos)
{
	Eigen::Vector3d POS(pos[0], pos[1], pos[2]);
	Eigen::Vector3d INI_POS(ini_pos[0], ini_pos[1], ini_pos[2]);

	double dtrop;

	double tot_p = 1013; // Mbar
	double pwv = 16;    //Mbar
	double pd = tot_p - pwv; // Mbar

	double T = 300; // Kelvins
	double R = 6371e3;

	double sat_pos_norm = EucldNorm(pos[0], pos[1], pos[2]);
	double rec_pos_norm = EucldNorm(ini_pos[0], ini_pos[1], ini_pos[2]);


	double g = acos((POS.dot(INI_POS)) / (sat_pos_norm*rec_pos_norm));
	double e = atan((sat_pos_norm*cos(g) - R) / (sat_pos_norm*sin(g)));
	//cout << (e*180)/ PI << endd;
	//cin.get();
	double z = (PI / 2) - e;

	if (z == 0) {
		dtrop = (0.002277 / cos(z))*(pd + (((1255 / T) + 0.05)*pwv));
	}
	else {
		dtrop = (0.002277 / cos(z))*(pd + (((1255 / T) + 0.05)*pwv) - (1.16*(pow(tan(z), 2))));
	}
	return dtrop;
}

// This function removes effects of ionosphere based on the iono-free model
double SatellitePosition::iono_free(double count)
{
	double gamma = pow((1575.42 / 1227.6), 2);
	double P1_temp = 0.0, P2_temp = 0.0, PR_temp = 0.0, C1_temp = 0.0;
	// Finding the pseudoranges
	if (find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "P1") - RinexFile->head.obstype.begin() != 0) {
		ptrdiff_t pos_P1 = find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "P1") - RinexFile->head.obstype.begin();
		if (pos_P1 == RinexFile->head.obstype.size()) {

		}
		else {
			P1_temp = RinexFile->obs.record[count][pos_P1];
		}
	}

	if (find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "P2") - RinexFile->head.obstype.begin() != 0) {
		ptrdiff_t pos_P2 = find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "P2") - RinexFile->head.obstype.begin();
		if (pos_P2 == RinexFile->head.obstype.size()) {

		}
		else {
			P2_temp = RinexFile->obs.record[count][pos_P2];
		}
	}

	if (P2_temp != 0 && P1_temp != 0) {
		PR_temp = (P2_temp - (gamma*P1_temp)) / (1 - gamma);
	}
	else {
		ptrdiff_t pos_PR = find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "C1") - RinexFile->head.obstype.begin();
		PR_temp = RinexFile->obs.record[count][pos_PR];
	}
	return PR_temp;
}

// This function removes effects of ionosphere based on the Klobuchar model
double SatellitePosition::klobuchar(double fi, double lambda, double elev, double azimuth, double tow)
{
	VectorXd alfa;
	VectorXd beta;

	alfa << 9.3132e-09, 1.4901e-08, -5.9605e-08, -1.1921e-07;
	beta << 8.8064e+04, 4.9152e+04, -1.3107e+05, -3.2768e+05;

	double c = 299792458, deg2semi = 1.0 / 180.0, semi2rad = PI, deg2rad = PI / 180.0, dIon1;

	double a = azimuth*deg2rad;                  // asimuth in radians
	double e = elev*deg2semi;                    // elevation angle in  semicircles
	double psi = 0.0137 / (e + 0.11) - 0.022;           //Earth Centered angle
	double lat_i = fi*deg2semi + psi*cos(a);     // Subionospheric lat

	if (lat_i > 0.416) { lat_i = 0.416; }
	else if (lat_i < -0.416) { lat_i = -0.416; }

	//Subionospheric long
	double long_i = lambda*deg2semi + (psi*sin(a) / cos(lat_i*semi2rad));

	// Geomagnetic latitude
	double lat_m = lat_i + 0.064*cos((long_i - 1.617)*semi2rad);

	double t = 4.32e4*long_i + tow;
	t = fmod(t, 86400.0);         // Seconds of day

	if (t > 86400.0) { t = t - 86400.0; }
	if (t < 86400.0) { t = t + 86400.0; }

	double sF = 1.0 + 16.0* pow((0.53 - e), 3);    // Slant factor 

	// Period of model
	double PER = beta(0) + beta(1)*lat_m + beta(2)* pow(lat_m, 2) + beta(3)* pow(lat_m, 3);

	if (PER < 72000.0) { PER = 72000.0; }
	double x = 2.0*PI*(t - 50400.0) / PER;    // Phase of the model

	double AMP = alfa(0) + alfa(1) *lat_m + alfa(2) * pow(lat_m, 2) + alfa(3) * pow(lat_m, 3); // Amplitud of the model
	if (AMP < 0.0) { AMP = 0.0; }

	if (abs(x) > 1.57) { dIon1 = sF * (5.0e-9); }
	else { dIon1 = sF * (5.0e-9 + AMP*(1.0 - x*x / 2.0 + x*x*x*x / 24.0)); }

	return dIon1;
}

// This function computes the gps time in seconds
double SatellitePosition::gpstimeCalc(double year, double month, double day, double hour, double min, double sec)
{
	double y, m, JD, gps_week, gps_seconds;
	double secs_per_week = 604800;     // Seconds in one week
	
	//Converts the two digit year to a four digit year.
    //Two digit year represents a year in the range 1980 - 2079.
	if (year >= 80 && year <= 99) { year = 1900 + year; }
	if (year >= 0 && year <= 79) { year = 2000 + year; }

	//Calculates the 'm' term used below from the given calendar month.
	if (month <= 2) { y = year - 1; m = month + 12; }
	if (month > 2) { y = year; m = month; }

	//Computes the Julian date corresponding to the given calendar date.
	JD = floor((365.25 * y)) + floor((30.6001 * (m + 1))) +
		day + ((hour + min / 60 + sec / 3600) / 24) + 1720981.5;

	//Computes the GPS week corresponding to the given calendar date.
	gps_week = floor((JD - 2444244.5) / 7);

	//Computes the GPS seconds corresponding to the given calendar date.
	gps_seconds = round(((((JD - 2444244.5) / 7) - gps_week)*secs_per_week) / 0.5)*0.5;
	return gps_seconds;
}

//This functions corrects for tropospheric effects based on initial values given in the navigation file
double SatellitePosition::NMF_hyd(double sat_elev)
{
	double A = 33.748;
	double B = -0.782;
	double A1 = 38.079;
	double B1 = -0.8088;
	double Y1 = A*pow(sat_elev, B);
	double Y2 = A1*pow(sat_elev, B1);
	double Y = Y1 + Y2;
	return Y;
}

// This function initializes the observation and navigation filenames
void SatellitePosition::SetFilenames(const char * ObsFile, const char * NavFile)
{
	size_t obs_size = strlen(ObsFile);
	size_t nav_size = strlen(NavFile);

	if (ObsFile != NULL) {
		strcpy(obsFile, ObsFile);
	}
	if (NavFile != NULL) {
		strcpy(navFile, NavFile);
	}

	obsFile[obs_size] = '\0';
	navFile[nav_size] = '\0';

	RinexFile = new reader(SatellitePosition::obsFile, SatellitePosition::navFile);
	RinexFile->rinexNAVReader(RinexFile->SV_data);
}

//This function opens the observation and navigation files from the environment file
void SatellitePosition::OpenEnvFile(ENVI_FILE & envi_records, SatellitePosition & object1, SatellitePosition & object2)
{
	Env_File = fopen("EnvFileGMIS.txt", "rt");
	if (Env_File == NULL) perror("Error opening file");

	char line[MAX_STRLEN];
	fgets(line, MAX_STRLEN, Env_File);

	while (!feof(Env_File)) {
		if (sizeof(line) >= 61) {
			if (strstr(line, "VERSION / TYPE")) {
				string sline(line);  // convert the line (char) to string class type sLine
				string token = sline.substr(4, 4);   // from 4th character and 4 spaces
				envi_records.version = token;
			}
			if (strstr(line, "BASE FILE(OBS/NAV)")) {
				std::istringstream iss(line);
				std::vector<std::string> results((std::istream_iterator<std::string>(iss)),
					std::istream_iterator<std::string>());
				Filename1_Obs = results[0]; Filename1_Nav = results[1];
				object1.SetFilenames(Filename1_Obs.c_str(), Filename1_Nav.c_str());
			}
			if (strstr(line, "ROVER FILE(OBS/NAV)")) {
				std::istringstream iss2(line);
				std::vector<std::string> results2((std::istream_iterator<std::string>(iss2)),
					std::istream_iterator<std::string>());
				Filename2_Obs = results2[0]; Filename2_Nav = results2[1];
				object2.SetFilenames(Filename2_Obs.c_str(), Filename2_Nav.c_str());
			}
			if (strstr(line, "BASE FINAL COORDINATES")) {
				string sline(line);  // convert the line (char) to string class type sLine
				string token = sline.substr(4, 56);   // from 4th character upto 56
				char cline[MAX_BYTE];     // cline is a char line
				strncpy(cline, token.c_str(), sizeof(cline));
				cline[sizeof(cline) - 1] = 0;
				for (char * r = strtok(cline, " "); r; r = strtok(NULL, " "))
				{
					double rr = atof(r);
					envi_records.base_final_coord.push_back(rr);
				}
			}
		}
		//......... Get a new line.......................
		fgets(line, MAX_STRLEN, Env_File);
	}
}

// This function computes the satellite positions and stores them in a global file
void SatellitePosition::Sat_Pos(SATPOS_struct &EPOCH_DATA)
{
	RinexFile->rinexOBSReader(RinexFile->head, RinexFile->obs, RinexFile->SENTINEL);  // Read the observation file
	epoch_count = epoch_count + 1;
	vector< vector<double> > pos_rot_epoch_temp;
	vector< vector<double> > initial_satpos_temp;
	vector <double> pr_epoch_temp;
	vector <double> P1_epoch_temp;
	vector <double> P2_epoch_temp;
	vector <double> C1_epoch_temp;
	vector <double> L1_epoch_temp;
	vector <double> L2_epoch_temp;
	vector <double> pr_epoch_Iono_temp;
	vector <double> sat_offset_epoch_temp;
	vector<double> sat_elev_temp;
	double Tropo;
	reader::NAVRECORDS storer;

	if (RinexFile->obs.NumOfSat != 0) {
		double count = 0;
		for (int j = 0; j < RinexFile->obs.NumOfSat; j++) { // Each satellite
			double gamma = pow((1575.42 / 1227.6), 2);
			double c = 299792458;
			double P1_temp = 0.0, P2_temp = 0.0, PR_temp = 0.0, L1_temp = 0.0, L2_temp = 0.0, C1_temp = 0.0, TIME_temp;
			//reader::NAVRECORDS storer;

			int curPrn = RinexFile->obs.prn[count];        // Current observation satelite vehicle
			double time_LOGGER = std::numeric_limits<double>::max();

			//search through the navigation file
			for (size_t k = 0; k < RinexFile->SV_data.size(); k++) {
				if (RinexFile->SV_data[k].PRN == curPrn) { //Compare obs prn with nav prn
					if (time_LOGGER > abs(RinexFile->SV_data[k].gpstime - RinexFile->obs.time)) {
						time_LOGGER = abs(RinexFile->SV_data[k].gpstime - RinexFile->obs.time);
						storer = RinexFile->SV_data[k];
					}
				}
			}
			// Finding the pseudoranges P1 code
			if (find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "P1") - RinexFile->head.obstype.begin() != RinexFile->head.nobstype) {

				ptrdiff_t pos_P1 = find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "P1") - RinexFile->head.obstype.begin();
				if (pos_P1 == RinexFile->head.obstype.size()) {

				}
				else {
					P1_temp = RinexFile->obs.record[j][pos_P1];
				}
			}
			// Finding the pseudoranges P2 code
			if (find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "P1") - RinexFile->head.obstype.begin() != RinexFile->head.nobstype) {

				ptrdiff_t pos_P2 = find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "P1") - RinexFile->head.obstype.begin();
				if (pos_P2 == RinexFile->head.obstype.size()) {

				}
				else {
					P2_temp = RinexFile->obs.record[j][pos_P2];
				}
			}
			// Finding the pseudoranges L1 carrier
			if (find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "L1") - RinexFile->head.obstype.begin() != RinexFile->head.nobstype) {

				ptrdiff_t pos_L1 = find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "L1") - RinexFile->head.obstype.begin();

				if (pos_L1 == RinexFile->head.obstype.size()) {

				}
				else {
					L1_temp = RinexFile->obs.record[j][pos_L1];
				}
			}
			// Finding the pseudoranges L2 carrier
			if (find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "L2") - RinexFile->head.obstype.begin() != RinexFile->head.nobstype) {

				ptrdiff_t pos_L2 = find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "L2") - RinexFile->head.obstype.begin();

				if (pos_L2 == RinexFile->head.obstype.size()) {

				}
				else {
					L2_temp = RinexFile->obs.record[j][pos_L2];
				}
			}
			// Finding the pseudoranges C1 code
			if (find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "C1") - RinexFile->head.obstype.begin() != RinexFile->head.nobstype) {

				ptrdiff_t pos_C1 = find(RinexFile->head.obstype.begin(), RinexFile->head.obstype.end(), "C1") - RinexFile->head.obstype.begin();
				if (pos_C1 == RinexFile->head.obstype.size()) {

				}
				else {
					C1_temp = RinexFile->obs.record[j][pos_C1];
				}
			}
			PR_temp = iono_free(j);
			TIME_temp = RinexFile->obs.time - (PR_temp / c);  //GPS system time at time of transmission corrected for transit time


		    // SATELLITE POSITIONING
			// Constants
			//WGS 84 value of the earth's gravitational constant for GPS user
			double mu = 3.986005e14;

			//WGS 84 value of the earth's rotation rate (rad/sec)
			double OMEGA_e = 7.2921151467e-5;

			// Initialize necessary variables from NAV
			double A = pow(storer.Sqrta, 2);   //Semi-major axis
			double toe = storer.Toe;           //toe = Time of Ephemeris
			double deltaN = storer.DeltaN;     //Delta n
			double Mo = storer.Mo;		       // Mo angle
			double e = storer.Eccentricity;    //Orbit Eccentricity
			double omega = storer.omega;       //omega Angle

			// Second Harmonic Perturbations
			double Cus = storer.Cus;          //Latitude Correction Sinus Component
			double Cuc = storer.Cuc;		  //Latitude Correction Cosinus Component
			double Crs = storer.Crs;		  //Radius Correction Sinus Component
			double Crc = storer.Crc;		  //Radius Correction Cosinus Component
			double Cis = storer.Cis;          //Angular Velocity
			double Cic = storer.Cic;		  //Inclination Correction Cosinus Component
			double IDOT = storer.IDOT;        //Inclination Rate
			double Io = storer.Io;            //Initial Inclination

			double OMEGA = storer.OMEGA;      //OMEGA Angle
			double OMEGADOT = storer.OMEGADOT; //Angular Velocity

			double SVCB = storer.SVClockBias;  //SVClockBias (sec)
			double SVCD = storer.SVClockDrift;  //SVClockDrift (sec)
			double SVCDR = storer.SVClockDriftRate;  //SVClockDriftRate (sec)
			double SAT_TIME = storer.gpstime;
			//double SAT_TIME = storer.SEC;

			double tk = TIME_temp - toe;       //Time from ephemeris reference epoch

			// Correct tk within range of -302400 to 302400
			if (tk > 302400) {
				tk = tk - 604800;
			}
			else if (tk < -302400) {
				tk = tk + 604800;
			}

			// Computed mean motion (rad/sec)
			double n0 = sqrt(mu / (pow(A, 3)));

			// Corrected mean motion
			double n = n0 + deltaN;

			//Mean anomaly
			double Mk = Mo + n*tk;

			// Keep Mk with in 360 degrees ....// Mk = rem(Mk+2*pi,2*pi);

			// Kepler's Equation for Eccentric Anomaly (may be solved by iteration)(radians)
			double Ek = Mk;
			double Ek0 = Mk + e*sin(Ek);
			double limit = 1e-11;   // limit

			while (abs(Ek - Ek0) > limit) {
				Ek0 = Ek;
				Ek = Mk + e*sin(Ek0);
			}

			// Keep Ek with in 360 degrees.....// Ek = rem(Ek + 2 * pi, 2 * pi);

			// True Anomaly
			double Y = (sqrt(1 - pow(e, 2))*sin(Ek)) / (1 - e*cos(Ek));
			double X = (cos(Ek) - e) / (1 - e*cos(Ek));
			double vk = atan2(Y, X);

			//Argument of Latitude
			double Phi_k = vk + omega;

			double del_uk = Cus*sin(2 * Phi_k) + Cuc*cos(2 * Phi_k);      // Argument of Latitude Correction
			double del_rk = Crs*sin(2 * Phi_k) + Crc*cos(2 * Phi_k);      //Radius Correction
			double del_ik = Cis*sin(2 * Phi_k) + Cic*cos(2 * Phi_k);      //Inclination Correction

			//Corrected Argument of Latitude
			double uk = Phi_k + del_uk;

			//Corrected Radius
			double rk = A*(1 - e*cos(Ek)) + del_rk;

			//Corrected Inclination
			double ik = Io + del_ik + IDOT*tk;

			//Positions in orbital plane.
			double xk_prime = rk*cos(uk);
			double yk_prime = rk*sin(uk);

			//Corrected longitude of ascending node
			double OMEGA_k = OMEGA + (OMEGADOT - OMEGA_e)*tk - OMEGA_e*toe;

			//Earth-fixed coordinates
			double xk = xk_prime*cos(OMEGA_k) - yk_prime*cos(ik)*sin(OMEGA_k);
			double yk = xk_prime*sin(OMEGA_k) + yk_prime*cos(ik)*cos(OMEGA_k);
			double zk = yk_prime*sin(ik);

			// Satellite position without tropo, iono and Earth rotation corrections
			vector <double> POSITION;

			//Store the Earth-fixed coordinates (satellite position xk,yk,zk of satellite k)
			POSITION.push_back(xk);
			POSITION.push_back(yk);
			POSITION.push_back(zk);

			// Satellite elevation angle calculation
			//--------------------------------------------------------------------------------------------------------------//
			Eigen::Vector3d POS(POSITION[0], POSITION[1], POSITION[2]);
			Eigen::Vector3d INI_POS(RinexFile->head.antpos[0], RinexFile->head.antpos[1], RinexFile->head.antpos[2]);
			double R = 6371e3;
			double sat_pos_norm = EucldNorm(POSITION[0], POSITION[1], POSITION[2]);
			double rec_pos_norm = EucldNorm(RinexFile->head.antpos[0], RinexFile->head.antpos[1], RinexFile->head.antpos[2]);
			double g = acos((POS.dot(INI_POS)) / (sat_pos_norm*rec_pos_norm));
			double sat_elev = (atan((sat_pos_norm*cos(g) - R) / (sat_pos_norm*sin(g))))*(180 / PI); // in degrees
			sat_elev_temp.push_back(sat_elev); // store the current satellite elevation
			//----------------------------------------------------------------------------------------------------------------//

			// CUT_OFF ANGLE OF 15 DEGREES
			if (sat_elev >= 15) {
				// Realativistic Correction Term (sec)
				double F = -4.442807633e-10;
				double Del_tr = F*e*sqrt(A)*sin(Ek);   // Realativistic Correction Term(sec)

				// Satellite Clock Errors (sec)
				double Del_tsv = SVCB + SVCD*(SAT_TIME - toe) + SVCDR*(pow((SAT_TIME - toe), 2)) + Del_tr;

				// CORRECTIONS
				// Range before earth rotation correction
				double geom_range = EucldNorm((POSITION[0] - RinexFile->head.antpos[0]), (POSITION[1] - RinexFile->head.antpos[1]), (POSITION[2] - RinexFile->head.antpos[2]));
				double t = (geom_range) / c;

				// CORRECTIONS
				// Earth rotation correction
				vector <double> POSITION_ROT = EarthRot_Corr(t, POSITION);
				double geom_range_corr = EucldNorm((POSITION_ROT[0] - RinexFile->head.antpos[0]), (POSITION_ROT[1] - RinexFile->head.antpos[1]), (POSITION_ROT[2] - RinexFile->head.antpos[2]));
				pos_rot_epoch_temp.push_back(POSITION_ROT);

				double check_pseudo = (PR_temp - geom_range_corr) + c*(Del_tsv);

				//	Tropo correction
				//double dtropo = NMF_hyd(sat_elev);   // Niell Mapping Function
				double dtropo = Tropo_Corr(POSITION_ROT, RinexFile->head.antpos);   // SAASTAMOINEN
				//Tropo = dtropo;

				// Corrected pseudorange for iono only
				pr_epoch_Iono_temp.push_back(PR_temp);

				// Corrected pseudorange for iono and tropo
				//double Corr_Pseudo = PR_temp - dtropo; doing tropo correc in kalman!!!!!
				double Corr_Pseudo = PR_temp - dtropo;
				pr_epoch_temp.push_back(Corr_Pseudo);
				P1_epoch_temp.push_back(P1_temp);
				P2_epoch_temp.push_back(P2_temp);
				C1_epoch_temp.push_back(C1_temp);
				L1_epoch_temp.push_back(L1_temp);
				L2_epoch_temp.push_back(L2_temp);
				sat_offset_epoch_temp.push_back(Del_tsv);
				initial_satpos_temp.push_back(POSITION);
				count = count + 1;
			}
			else {
				//remove the satellite from the prn vector
				count = count;
				RinexFile->obs.prn.erase(RinexFile->obs.prn.begin() + count);
			}
		}
		EPOCH_DATA.Sat_elev = sat_elev_temp; // Satellite elevations for the epoch
		EPOCH_DATA.pos_rot_epoch = pos_rot_epoch_temp;     // corrected satellite positions for the epoch
		EPOCH_DATA.Initial_SatPos = initial_satpos_temp;
		EPOCH_DATA.pr_epoch = pr_epoch_temp;                // Corrected pseudorange
		EPOCH_DATA.P1_epoch = P1_epoch_temp;
		EPOCH_DATA.P2_epoch = P2_epoch_temp;
		EPOCH_DATA.C1_epoch = C1_epoch_temp;
		EPOCH_DATA.L1_epoch = L1_epoch_temp;
		EPOCH_DATA.L2_epoch = L2_epoch_temp;
		EPOCH_DATA.pr_epoch_iono = pr_epoch_Iono_temp;
		EPOCH_DATA.sat_offset_epoch = sat_offset_epoch_temp;  // satellite clock offset
	}
}

