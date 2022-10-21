// 2021/08/06 update
// Satoru OKUDA - okuda@mechgen.jp

#include <iostream>
#include <math.h>
#include <vector>
#include <random>
#include <fstream>
#include <unistd.h>

#define PATH_SIZE 512

using namespace std;

std::string workspace;

const int skip = 10000;

const double TimePeriod = 10.0;
const double dt = 0.000001;
const double dx = 0.025;
const double Lx = 10.0;

const double nAFinit = 1.0;

const double Mass = 0.01; // Unit
const double xi = 1.0; // Fix
const double chi = 1.0; // Unit
const double kAAoff = 1.0; // Unit
const double LAF = 1.0; // Unit
double rhoMY = 1.0;
double rhoAA = 1.0;
const double E = 1.0; // Fix
const double rateAF = 1.0; // Unit
const double T0 = 1.0;
const double EPS = 1.0E-12;
const double LIM = 1.0E+6;

std::vector<double> nAFCenter;
std::vector<double> nAFTotal;
std::vector<double> Pressure;
std::vector<double> GradPressure;
std::vector<double> SumAMForce;

std::vector<double> SingleAFVelocity;
std::vector<double> RelativeFrictionForce;

std::vector<double> Velocity;

int Nx = (int)(Lx/dx)+2;
int NxLAFHalf = (int)(0.5*LAF/dx);
int NxLAFTotal = 2*NxLAFHalf;

int step;

void Output(double t){
	double Vmax = 0.0;
	for(int i=0; i<Nx-1; i++){
		if(Vmax < fabs(Velocity.at(i))) Vmax = fabs(Velocity.at(i));
	}
	
	int num = 0;
	double AFNumAVE = 0.0;
	double AFNumSD = 0.0;
	double AFNumMIN = 1.0e+24;
	for(int i=NxLAFHalf+1; i<Nx-NxLAFHalf-1; i++){
		num++;
		AFNumAVE += nAFCenter.at(i);
		AFNumSD += nAFCenter.at(i)*nAFCenter.at(i);
		if(AFNumMIN > nAFCenter.at(i)) AFNumMIN = nAFCenter.at(i);
	}
	AFNumAVE /= (double)num;
	AFNumSD /= (double)num;
	AFNumSD -= AFNumAVE*AFNumAVE;
	AFNumSD = std::sqrt(AFNumSD);
	
	double RelForce = 0.0;
	for(int i=1; i<NxLAFHalf+1; i++){
		RelForce += -(SumAMForce.at(i) - RelativeFrictionForce.at(i));
	}
	
	ofstream FileAFCenter("op_AFCenterNumber.dat", ios::app);
	ofstream FileAFTotal("op_AFTotalNumber.dat", ios::app);
	ofstream FileTS("op_SFSumAMForce.dat", ios::app);
	ofstream FileVL("op_Velocity.dat", ios::app);
	ofstream FileAFT("op_SumAMForcePerAF.dat", ios::app);
	ofstream FileMV("op_MaxVelocity.dat", ios::app);
	ofstream FileEF("op_EndForce.dat", ios::app);
	ofstream FileAVE("op_AFNumAVE.dat", ios::app);
	ofstream FileSD("op_AFNumSD.dat", ios::app);
	ofstream FileMIN("op_AFNumMIN.dat", ios::app);
	
	FileAFCenter << t;
	FileAFTotal << t;
	FileTS << t;
	FileVL << t;
	FileAFT << t;
	FileMV << t;
	FileEF << t;
	FileAVE << t;
	FileSD << t;
	FileMIN << t;
	
	for(int i=1; i<Nx-1; i++){
		FileAFCenter << "\t" << nAFCenter.at(i);
		FileAFTotal << "\t" << nAFTotal.at(i);
	}
	for(int i=1; i<Nx-2; i++){
		FileTS << "\t" << SumAMForce.at(i);
		FileAFT << "\t" << SumAMForce.at(i)/nAFCenter.at(i);
		FileVL << "\t" << Velocity.at(i);
	}
	FileMV << "\t" << Vmax;
	FileEF << "\t" << RelForce;
	FileAVE << "\t" << AFNumAVE;
	FileSD << "\t" << AFNumSD;
	FileMIN << "\t" << AFNumMIN;
	
	FileAFCenter << std::endl;
	FileAFTotal << std::endl;
	FileTS << std::endl;
	FileVL << std::endl;
	FileAFT << std::endl;
	FileMV << std::endl;
	FileEF << std::endl;
	FileAVE << std::endl;
	FileSD << std::endl;
	FileMIN << std::endl;
}

void CalculateSumAMForce(){
	for(int i=0; i<Nx; i++){
		nAFTotal.at(i) = 0.0;
		for(int j=i-NxLAFHalf; j<i+NxLAFHalf; j++){
		if( (j >= 0) && (j < Nx) ){
			nAFTotal.at(i) += nAFCenter.at(j) * dx;
		}}
	}
	
	for(int i=0; i<Nx; i++){
		Pressure.at(i) = - nAFTotal.at(i) * rhoMY * chi * LAF;
	}
	
	for(int i=1; i<Nx-2; i++){
		GradPressure.at(i) = (Pressure.at(i+1) - Pressure.at(i)) / dx;
	}
	
	for(int i=1; i<Nx-2; i++){
		SumAMForce.at(i) = 0.0;
		for(int j=i-NxLAFHalf; j<i+NxLAFHalf; j++){
		if( (j >= 1) && (j < Nx-2) ){
			SumAMForce.at(i) += GradPressure.at(j) / nAFTotal.at(j) * dx;
		}}
		SumAMForce.at(i) *= nAFCenter.at(i);
		
		if(std::isnan(SumAMForce[i])||std::isinf(SumAMForce[i])||std::fabs(SumAMForce[i])>LIM){
			std::cout << "i = " << i << std::endl;
			std::cout << "SumAMForce.at(i) = " << SumAMForce.at(i) << std::endl;
			exit(1);
		}
	}
}

void CalculateRelativeFrictionForce(){
	for(int i=1; i<Nx-2; i++){
		double tmpN = 0.0;
		RelativeFrictionForce.at(i) = 0.0;
		for(int j=i-2*NxLAFHalf; j<i+2*NxLAFHalf; j++){
		if( (j >= 1) && (j < Nx-2) ){
			double Vi = 0.5 * (Velocity.at(i) + Velocity.at(i+1));
			double Vj = 0.5 * (Velocity.at(j) + Velocity.at(j+1));
			double gamma = rhoAA * E / kAAoff * (Vi - Vj) * (LAF - dx * (double)abs(i-j));
			
			double Ni = 0.5 * (nAFCenter.at(i) + nAFCenter.at(i+1));
			double Nj = 0.5 * (nAFCenter.at(j) + nAFCenter.at(j+1));
			RelativeFrictionForce.at(i) += Ni * Nj * gamma * dx;
		}}
		RelativeFrictionForce.at(i) *= rhoAA * E / kAAoff;
		
		if(std::isnan(RelativeFrictionForce.at(i))||std::isinf(RelativeFrictionForce.at(i))||std::fabs(RelativeFrictionForce.at(i))>LIM){
			std::cout << "i = " << i << std::endl;
			std::cout << "RelativeFrictionForce.at(i) = " << RelativeFrictionForce.at(i) << std::endl;
			exit(1);
		}
	}
}

void CalculateTimeDevelopment(){
	for(int i=1; i<Nx-1; i++){
		//if(SumAMForce.at(i) / (nAFCenter.at(i)+EPS) / T0 > 1.0){
		//	nAFCenter[i] += nAFCenter.at(i) * rateAF * dt;
		//}
		
		nAFCenter[i] += (nAFinit - nAFCenter.at(i)) * rateAF * dt;
		if(nAFCenter[i] < 0.0) nAFCenter[i] = 0.0;
	}
	
	for(int i=1; i<Nx-2; i++){
		double vel = Velocity.at(i);
		if(vel > 0.0){
			double dn = nAFCenter.at(i) * vel * dt / dx;
			if(dn > nAFCenter[i]) dn = nAFCenter[i];
			nAFCenter[i]   -= dn;
			nAFCenter[i+1] += dn;
		}
		else if(vel < 0.0){
			double dn = - nAFCenter.at(i+1) * vel * dt / dx;
			if(dn > nAFCenter[i+1]) dn = nAFCenter[i+1];
			nAFCenter[i]   += dn;
			nAFCenter[i+1] -= dn;
		}
	}
	
	nAFCenter[0] = 0.0;
	for(int i=1; i<NxLAFHalf+1; i++){
		nAFCenter[i] = nAFinit;
	}
	for(int i=Nx-NxLAFHalf-1; i<Nx-1; i++){
		nAFCenter[i] = nAFinit;
	}
	nAFCenter[Nx-1] = 0.0;
	
	for(int i=NxLAFHalf+1; i<Nx-NxLAFHalf-2; i++){
		double tmpN = 0.5 * (nAFCenter[i] + nAFCenter[i+1]);
		double A = (SumAMForce.at(i)
						- tmpN * xi * LAF * Velocity.at(i)
				 		- RelativeFrictionForce.at(i)
				 ) / (tmpN * Mass);
		Velocity[i] += A * dt;
		if(std::isnan(Velocity[i])||std::isinf(Velocity[i])||std::fabs(Velocity[i])>LIM){
			std::cout << "step = " << step << std::endl;
			std::cout << "i = " << i << std::endl;
			std::cout << "Velocity[i] = " << Velocity[i] << std::endl;
			std::cout << "SumAMForce.at(i+1) = " << SumAMForce.at(i+1) << std::endl;
			std::cout << "SumAMForce.at(i-1) = " << SumAMForce.at(i-1) << std::endl;
			exit(1);
		}
	}
	Velocity[0] = 0.0;
	for(int i=1; i<NxLAFHalf+1; i++){
		Velocity[i] = 0.0;
	}
	for(int i=Nx-NxLAFHalf-2; i<Nx-2; i++){
		Velocity[i] = 0.0;
	}
	Velocity[Nx-2] = 0.0;
}

int main(int argc, char *argv[]){
	if(argc < 3){
		std::cout << "Error command line value" << std::endl;
		std::cout << "argv[1]: directory" << std::endl;
		std::cout << "argv[2]: rhoMY" << std::endl;
		std::cout << "argv[3]: rhoAA" << std::endl;
		exit(1);
	}
	else {
		workspace = std::string(argv[1]);
		rhoMY = std::atof(argv[2]);
		rhoAA = std::atof(argv[3]);
		
		std::cout << "Workspace: " << workspace << std::endl;
		std::cout << "rhoMY: " << rhoMY << std::endl;
		std::cout << "rhoAA: " << rhoAA << std::endl;
	}
	
	char topspace[260];
	#ifdef _MSC_VER
	_getcwd(topspace,260);
	short sw_chdir = _chdir(workspace.c_str());
	#else
	getcwd(topspace,260);
	short sw_chdir = chdir(workspace.c_str());
	#endif
	if(sw_chdir != 0) {
		std::cout << "Error: chdir to workspace" << endl;
		exit(1);
	}
	
	char curspace[260];
	#ifdef _MSC_VER
	_getcwd(curspace,260);
	#else
	getcwd(curspace,260);
	#endif
	std::cout << "Cur dir: " << curspace << std::endl;
	
	
	std::random_device rnd;
	std::default_random_engine engine(rnd());
	uniform_real_distribution<> dist(0.0,1.0);
	
	ofstream FileP("op_parameter.dat");
	FileP << "xi = " << xi << "\n";
	FileP << "chi = " << chi << "\n";
	FileP << "kAAoff = " << kAAoff << "\n";
	FileP << "LAF = " << LAF << "\n";
	FileP << "rhoMY = " << rhoMY << "\n";
	FileP << "rhoAA = " << rhoAA << "\n";
	FileP << "E = " << E << "\n";
	FileP << "rateAF = " << rateAF << "\n";
	FileP << "T0 = " << T0;
	FileP << std::endl;
	
	nAFCenter.resize(Nx);
	nAFTotal.resize(Nx);
	Pressure.resize(Nx);
	for (int i = 0; i < Nx; i++){
		nAFCenter[i] = nAFinit + 0.1*(2.0*dist(engine)-1.0);
		nAFTotal[i] = 0.0;
		Pressure[i] = 0.0;
	}
	nAFCenter[0] = 0.0;
	for(int i=1; i<NxLAFHalf+1; i++){
		nAFCenter[i] = nAFinit;
	}
	for(int i=Nx-NxLAFHalf-1; i<Nx-1; i++){
		nAFCenter[i] = nAFinit;
	}
	nAFCenter[Nx-1] = 0.0;
	
	GradPressure.resize(Nx-1);
	SumAMForce.resize(Nx-1);
	Velocity.resize(Nx-1);
	RelativeFrictionForce.resize(Nx-1);
	for (int i = 0; i < Nx-1; i++){
		GradPressure[i] = 0.0;
		SumAMForce[i] = 0.0;
		Velocity[i] = 0.0;
		RelativeFrictionForce[i] = 0.0;
	}
	
	std::cout << "Initial step" << std::endl;
	{
		step = 0;
		double t = (double)step * dt;
		
		CalculateSumAMForce();
		CalculateRelativeFrictionForce();
		Output(t);
	}
	
	std::cout << "Calculation start" << std::endl;
	for (step = 1; step <= TimePeriod/dt; step++){
		double t = (double)step * dt;
		
		if(step%skip==0) std::cout << "step = " << step << " t = " << t << std::endl;
		
		CalculateSumAMForce();
		CalculateRelativeFrictionForce();
		
		CalculateTimeDevelopment();
		
		if(step%skip==0){
			CalculateSumAMForce();
			CalculateRelativeFrictionForce();
			Output(t);
		}
	}
	
	return 0;
}

