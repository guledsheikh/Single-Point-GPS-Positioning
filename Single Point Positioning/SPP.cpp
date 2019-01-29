//***********************************************************************************
//Filename: "SPP.cpp
//
//Date				  Author				Revison History
//05/08/2017			  Guled Sheikh				version 1
//
//***********************************************************************************
//
//Copyright  2017~2018. All Rights Reserved
//
//***********************************************************************************

#include "SPP.h"

// Basic output file
void SPP::MSGFILE(char* outfilename, char* residual, SatellitePosition& object)
{
	ofstream fout;
	ofstream fout2;

	if (epoch == 1) {
		fout.open(outfilename); fout.close();
		fout.open(outfilename, ios::app);
		fout << "--------------------" << endd;
		fout << left << setw(20) << "GPS  OUTPUT  REPORT" << endd;
		fout << "-------------------" << endd;
		fout << left << setw(20) << "APROX. ANTENNA POSITION" << endd;
		fout << left << "X: " << object.RinexFile->head.antpos[0] << endd;
		fout << left << "Y: " << object.RinexFile->head.antpos[1] << endd;
		fout << left << "Z: " << object.RinexFile->head.antpos[2] << endd;
		fout << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endd;
		fout << left << setw(7) << "HOUR" << setw(7) << "MIN" << setw(10) << "SEC" << setw(16) << "X[m]"
			<< setw(16) << "Y[m]" << setw(16) << "Z[m]" << setw(15) << "Time[Sec]"
			<< setw(15) << "Sig_X[m]" << setw(15) << "Sig_Y[m]" << setw(15) << "Sig_Z[m]"
			<< setw(18) << "Sig_t[sec]" << setw(18) << "Sig_sqr[m2]"
			<< setw(18) << "LS_RclockE[sec]"
			<< endd;
		fout << "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endd;
	}
	else {
		fout.open(outfilename, ios::app);
	}

	bool curr_rec = true;
	if (object.RinexFile->obs.NumOfSat < 4) {
		curr_rec = false; //SIGMA_X
	}
	// For epochwise	
	if (curr_rec == true) {
		fout << std::fixed << setprecision(0) << left << setw(7) << object.RinexFile->obs.HOUR
			<< setw(7) << object.RinexFile->obs.MIN
			<< setprecision(4) << setw(10) << object.RinexFile->obs.Indivi_Sec
			<< setw(16) << object.epoch_data.Epoch_XYZt[0]
			<< setw(16) << object.epoch_data.Epoch_XYZt[1]
			<< setw(16) << object.epoch_data.Epoch_XYZt[2]
			<< setw(15) << setprecision(0) << object.RinexFile->obs.SEC
			<< setprecision(4) << setw(15) << SIGMA_X << setw(15) << SIGMA_Y << setw(15) << SIGMA_Z
			<< setw(18) << SIGMA_SQRD_epoch
			<< endd;
	}

	fout.close();
	fout2.close();
}

//This function initializes the receiver position
vector<double> SPP::Initialize_Rec_pos(SatellitePosition& object)
{
	vector<double> X, Y, Z, Pu, Ps, pos;
	double R = 6371000;
	MatrixXd A(object.RinexFile->obs.prn.size(), 4);
	MatrixXd A2(object.RinexFile->obs.prn.size(), 4);
	VectorXd L(object.RinexFile->obs.prn.size(), 1), xhat;
	MatrixXd Norm_Mat;

	// Pseudorange for receiver from first epoch
	for (size_t j = 0; j < object.RinexFile->obs.prn.size(); j++) {
		Pu.push_back(object.epoch_data.pr_epoch[j]);
	}

	// Satellite positions from first epoch
	for (size_t j = 0; j < object.RinexFile->obs.prn.size(); j++) {
		X.push_back(object.epoch_data.Initial_SatPos[j][0]);
		Y.push_back(object.epoch_data.Initial_SatPos[j][1]);
		Z.push_back(object.epoch_data.Initial_SatPos[j][2]);
	}

	// The satellite pseudorange
	for (size_t j = 0; j < object.RinexFile->obs.prn.size(); j++) {
		Ps.push_back(sqrt(pow(X[j], 2) + pow(Y[j], 2) + pow(Z[j], 2)));
	}

	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < object.RinexFile->obs.prn.size(); j++) {
			//First Loop
			if (i == 0)
			{
				//Design Matrix (Satellite Vehicle Position Matrix)
				A(j, 0) = -X[j];
				A(j, 1) = -Y[j];
				A(j, 2) = -Z[j];
				A(j, 3) = object.epoch_data.pr_epoch[j];

				//Observation Vector
				L(j) = 0.5*(pow(Pu[j], 2) - pow(Ps[j], 2) - pow(R, 2));
			}
			else
			{
				//Receiver Range from Earth Center
				R = sqrt(pow(xhat(0), 2) + pow(xhat(1), 2) + pow(xhat(2), 2));
				//Observations Vector
				L(j) = 0.5*(pow(Pu[j], 2) - pow(Ps[j], 2) - pow(R, 2));
				//Design Matrix
				A(j, 3) = object.epoch_data.pr_epoch[j] - (0.5*(xhat(3)));
			}

		}
		//Normal Equation Matrix
		Norm_Mat = ((A.transpose())*A).inverse();
		//Solution to unknowns using Least Square principle
		xhat = Norm_Mat*A.transpose()*L;
	}

	// Update the position vector with the final solution
	pos.push_back(xhat(0));
	pos.push_back(xhat(1));
	pos.push_back(xhat(2));
	pos.push_back(0);

	// Makes sure this function runs only once
	initial_flag = false;
	return pos;
}

//This function opens the observation and navigation files
void SPP::OpenEnvFile()
{
	Env_File = fopen("EnvFileSPP.txt", "rt");
	if (Env_File == NULL) perror("Error opening file");
	char line[MAX_STRLEN];

	fgets(line, MAX_STRLEN, Env_File);
	std::istringstream iss(line);
	std::vector<std::string> results((std::istream_iterator<std::string>(iss)),
		std::istream_iterator<std::string>());
	Filename1_Obs = results[0]; Filename1_Nav = results[1];
	CSatellitePosition.SetFilenames(Filename1_Obs.c_str(), Filename1_Nav.c_str());
}


//Epoch by Epoch Least Squares
void SPP::Synthesize(SatellitePosition& object)
{
	//Open the observation and navigation files for reading!!!
	OpenEnvFile();

	while (!feof(object.RinexFile->pFile)) {

		object.Sat_Pos(object.epoch_data);

		if (epochflag == true) {
			// First EPOCH
			//*........ LEAST SQUARES.................*//
			double x, y, z, t, c = 299792458;
			vector<double> X, Y, Z, ini_pos;
			vector< double> receiver_pos;
			vector< double> receiver_pos_ave;
			MatrixXf Norm_mat; // Normal eqn Matrix
			MatrixXf Norm_mat_temp; // Normal eqn Matrix
			MatrixXf Q_xx;
			MatrixXf D_xx;
			VectorXf residual(object.RinexFile->obs.prn.size());    // Residual
			epoch = epoch + 1;

			vector<double> sat_pos; //Satellite position definition
			vector<double> rec_pos; // Receiver position definition

			// Receiver position initialization
			ini_pos = Initialize_Rec_pos(object);
			x = ini_pos[0]; // 3d position
			y = ini_pos[1];
			z = ini_pos[2];
			t = 0;  // time
			
			// Satellite positions
			for (size_t j = 0; j <object.RinexFile->obs.prn.size(); j++) {
				X.push_back(object.epoch_data.pos_rot_epoch[j][0]);
				Y.push_back(object.epoch_data.pos_rot_epoch[j][1]);
				Z.push_back(object.epoch_data.pos_rot_epoch[j][2]);
			}

			for (int i = 0; i < 3; i++) {

				VectorXf delta(4);  delta(0) = 0;  delta(1) = 0;  delta(2) = 0;  delta(3) = 0;  // Initialize delta(4x1) to zero
				MatrixXf A_mat(object.RinexFile->obs.prn.size(), 4); // Design Matrix
				VectorXf W(object.RinexFile->obs.prn.size());    // Misclosure

				for (size_t j = 0; j < object.RinexFile->obs.prn.size(); j++) {
					//Geometric Range
					double dist_temp = sqrt(pow((X[j] - x), 2) + pow((Y[j] - y), 2) + pow((Z[j] - z), 2));

					//Design Matrix
					A_mat(j, 0) = -((X[j] - x) / dist_temp);
					A_mat(j, 1) = -((Y[j] - y) / dist_temp);
					A_mat(j, 2) = -((Z[j] - z) / dist_temp);
					A_mat(j, 3) = 1;

					// Misclosure
					//W(j) = PR_iono - dtropo - dist_temp - c*(t - CSatellitePosition.epoch_data.sat_offset_epoch[j]);
					W(j) = object.epoch_data.pr_epoch[j] - dist_temp - c*(t - object.epoch_data.sat_offset_epoch[j]);
				}

				//// Normal eqn Matrix
				Norm_mat_temp = (A_mat.transpose())*A_mat;
				Q_xx = Norm_mat_temp.inverse();

				// Solution to least squares
				delta = Q_xx*(A_mat.transpose() * W);
				
				//residual vector
				residual = A_mat*delta - W;

				x = x + delta(0);
				y = y + delta(1);
				z = z + delta(2);
				t = t + delta(3) / c;
			}
			
			// Store the residual vector in a global variable
			Residual = residual;
			double dof = object.epoch_data.pos_rot_epoch.size() - 4; //degrees of freedom
			double sig_sq = (((residual.transpose())*residual) / dof)[0]; //a-posteriori variance
			
			//Covariance matrix for the parameters
			D_xx = sig_sq*Q_xx;
			Prev_D_xx = D_xx;   // Creation of previous cov mat for the following epoch

			//Std deviations for X,Y and Z
			SIGMA_X = sqrt(D_xx(0, 0));
			SIGMA_Y = sqrt(D_xx(1, 1));
			SIGMA_Z = sqrt(D_xx(2, 2));
			SIGMA_t = sqrt(D_xx(3, 3));

			// Store the Receiver final position
			receiver_pos.push_back(x);
			receiver_pos.push_back(y);
			receiver_pos.push_back(z);
			receiver_pos.push_back(t);
			
			//Summing all the coordinates of each epoch
			SUM_X = SUM_X + x;
			SUM_Y = SUM_Y + y;
			SUM_Z = SUM_Z + z;
			
			//Taking the mean of the single point positions
			AVE_X = SUM_X / epoch;
			AVE_Y = SUM_Y / epoch;
			AVE_Z = SUM_Z / epoch;
			
			//Store the mean position
			receiver_pos_ave.push_back(AVE_X);
			receiver_pos_ave.push_back(AVE_Y);
			receiver_pos_ave.push_back(AVE_Z);

			// Storing the current position in Epoch_XYZt to be used in the next epoch initialization
			object.epoch_data.Epoch_XYZt = receiver_pos;
			// Storing the current position epoch in variable that stores all epoch solutions
			object.epoch_data.Epoch_XYZt_ave = receiver_pos_ave;
			epochflag = false;
		}
		else {

			// ALL NEXT EPOCHS
			double x, y, z, t, c = 299792458;
			vector<double> X, Y, Z, ini_pos;
			vector< double> receiver_pos;
			vector< double> receiver_pos_ave;
			MatrixXf Norm_mat; // Normal eqn Matrix
			MatrixXf Norm_mat_temp; // Normal eqn Matrix
			MatrixXf Q_xx;
			MatrixXf D_xx;
			VectorXf residual(object.RinexFile->obs.prn.size());    // Residual
			epoch = epoch + 1;

			// Receiver position initialization	
			x = object.epoch_data.Epoch_XYZt[0];
			y = object.epoch_data.Epoch_XYZt[1];
			z = object.epoch_data.Epoch_XYZt[2];
			t = 0;


			// Satellite positions
			for (size_t j = 0; j <object.RinexFile->obs.prn.size(); j++) {
				X.push_back(object.epoch_data.pos_rot_epoch[j][0]);
				Y.push_back(object.epoch_data.pos_rot_epoch[j][1]);
				Z.push_back(object.epoch_data.pos_rot_epoch[j][2]);
			}

			for (int i = 0; i < 3; i++) {

				VectorXf delta(4);  delta(0) = 0;  delta(1) = 0;  delta(2) = 0;  delta(3) = 0;  // Initialize delta(4x1) to zero
				MatrixXf A_mat(object.RinexFile->obs.prn.size(), 4); // Design Matrix
				VectorXf W(object.RinexFile->obs.prn.size());    // Misclosure

				for (size_t j = 0; j < object.RinexFile->obs.prn.size(); j++) {

					//Geometric Range
					double dist_temp = sqrt(pow((X[j] - x), 2) + pow((Y[j] - y), 2) + pow((Z[j] - z), 2));

					//Design Matrix
					A_mat(j, 0) = -((X[j] - x) / dist_temp);
					A_mat(j, 1) = -((Y[j] - y) / dist_temp);
					A_mat(j, 2) = -((Z[j] - z) / dist_temp);
					A_mat(j, 3) = 1;

					// Misclosure
					W(j) = object.epoch_data.pr_epoch[j] - dist_temp - c*(t - object.epoch_data.sat_offset_epoch[j]);
				}

				//// Normal eqnuation Matrix
				Norm_mat_temp = (A_mat.transpose())*A_mat;
				Q_xx = Norm_mat_temp.inverse();

				// Solution to least squares
				delta = Q_xx*(A_mat.transpose() * W);
				
				//residual vector
				residual = A_mat*delta - W;

				x = x + delta(0);
				y = y + delta(1);
				z = z + delta(2);
				//t = t + delta(3);   // rho and not delta t
				t = t + delta(3) / c;
			}
			
			// Store the residual vector in a global variable
			Residual = residual;

			double dof = object.RinexFile->obs.prn.size() - 4; //degrees of freedom
			double sig_sq; //a-posteriori variance

			if (dof <= 0.0) {
				sig_sq = 1.0;
			}
			else {
				sig_sq = (((residual.transpose())*residual) / dof)[0];
			}

			// Covariance matrix of the parameters
			D_xx = sig_sq*Q_xx;

			//Std deviations for X,Y and Z
			SIGMA_X = sqrt(D_xx(0, 0));
			SIGMA_Y = sqrt(D_xx(1, 1));
			SIGMA_Z = sqrt(D_xx(2, 2));
			SIGMA_t = sqrt(D_xx(3, 3));
			Prev_D_xx = D_xx;   // Creation of previous cov mat for the following epoch

			// Receiver final position
			receiver_pos.push_back(x);
			receiver_pos.push_back(y);
			receiver_pos.push_back(z);
			receiver_pos.push_back(t);

			//Summing all the coordinates of each epoch
			SUM_X = SUM_X + x;
			SUM_Y = SUM_Y + y;
			SUM_Z = SUM_Z + z;
			
			//Taking the mean of the single point positions
			AVE_X = SUM_X / epoch;
			AVE_Y = SUM_Y / epoch;
			AVE_Z = SUM_Z / epoch;
			
			//Store the mean position
			receiver_pos_ave.push_back(AVE_X);
			receiver_pos_ave.push_back(AVE_Y);
			receiver_pos_ave.push_back(AVE_Z);

			object.epoch_data.Epoch_XYZt = receiver_pos;
			object.epoch_data.Epoch_XYZt_ave = receiver_pos_ave;
		}
		
		// Output the results
		MSGFILE("Test\\SPP_rover.txt", "epoch_residual.txt", object);
		object.RinexFile->SENTINEL = 0;
	}
}
