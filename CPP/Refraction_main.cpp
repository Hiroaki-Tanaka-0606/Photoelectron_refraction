// Photoelectron refraction
// main program

#include <H5Cpp.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <omp.h>
#include <iostream>
#include <chrono>
#include <random>

#include "HDF5_tools.hpp"

using namespace std;
using namespace H5;

// Load a line from input
void readLine(char* line, int buffer_length, FILE* input){
	char* fgets_status=fgets(line, buffer_length, input);
	if(fgets_status==NULL){
		printf("Error in loading input\n");
		exit(0);
	}
}

// gauss function
double gauss(double x, double s){
	return 1.0/(sqrt(2*M_PI)*s)*exp(-x*x/(s*s*2));
}

// inner product
double inProd(double* a, double* b){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

// Calculate refracted wavevector K
bool calcK(double* k, double epk, double* n, double* K){
	double KLength=sqrt(2*epk);

	double kPerpLength=inProd(k, n);
	double kPara[3];
	int i;
	for(i=0; i<3; i++){
		kPara[i]=k[i]-kPerpLength*n[i];
	}

	double KPerpLength2=KLength*KLength-inProd(kPara, kPara);
	if(KPerpLength2<0){
		return false;
	}
	double KPerpLength=sqrt(KPerpLength2);
	if(kPerpLength<0){
		KPerpLength*=-1;
	}
	for(i=0; i<3; i++){
		K[i]=kPara[i]+KPerpLength*n[i];
	}
	return true;
}

int main(int argc, char** argv){
	cout << "Photoelectron refraction calculation" << endl;
	if(argc<2){
		printf("Usage: %s input\n", argv[0]);
		return 0;
	}

	// constatnts
	/// 1 Hartree (eV)
	double Eh_eV=27.2114;
	/// 1 Bohr (Ang)
	double Bohr_ang=0.529177;
	/// Calculation range of the gauss function
	double sigmaMax=5.0;

	// load input and convert to atomic unit
	FILE* input=fopen(argv[1], "r");
	// variables
	double W;                  // work function (Eh)
	double V0;                 // inner potential (Eh)
	double k0[3];              // parabola origin
	double a;                  // parabola paramter
	double V1;                 // parabola top potision from EF (Eh)
	bool kFlat;                // k plane is flat (true, 1) or curved (false, 0)
	double kFlat_kz;           // to specify flat plane
	double kCurved_k;          // to specify curved plane
	bool surfaceConst;         // surface is constant (true, 1) or random (false, 0)
	double surfaceConst_theta; // to specify a constant surface
	double surfaceConst_phi;
	int surfaceRandom_samples; // to specify random surfaces
	double kxMin;              // x axis
	double kxMax;
	int kxCount;
	double dkx;
	double kyMin;              // y axis
	double kyMax;
	int kyCount;
	double dky;
	double eMin;               // E axis
	double eMax;
	int eCount;
	double de;
	double sigmak;             // k broadening
	double sigmae;             // e broadening
	char outputName[256];

	int buffer_length=1024;
	char line[buffer_length+1];
	int i, j, k;
	chrono::system_clock::time_point start;
	chrono::system_clock::time_point end;
	double duration;
	// line 1: W
	readLine((char*)line, buffer_length, input);
	double W_eV;
  sscanf(line, "%lf", &W_eV);
	W=W_eV/Eh_eV;
	printf("%32s = %.2f eV = %.2f Eh\n", "Work function", W_eV, W);
	// line 2: V0
	readLine((char*)line, buffer_length, input);
	double V0_eV;
  sscanf(line, "%lf", &V0_eV);
	V0=V0_eV/Eh_eV;
	printf("%32s = %.2f eV = %.2f Eh\n", "Inner potential", V0_eV, V0);
	// line 3: k0
	readLine((char*)line, buffer_length, input);
  sscanf(line, "%lf %lf %lf", &k0[0], &k0[1], &k0[2]);
	for(i=0; i<3; i++){
		k0[i]*=Bohr_ang;
	}
	// line 4: a
	readLine((char*)line, buffer_length, input);
  sscanf(line, "%lf", &a);
	// line 5: V1
	readLine((char*)line, buffer_length, input);
	double V1_eV;
  sscanf(line, "%lf", &V1_eV);
	V1=V1_eV/Eh_eV;
	printf("%32s = (k - (%.2f, %.2f, %.2f))^2 * %.2f /2 - %.2f + %.2f (Eh)\n", "Initial state dispersion", k0[0], k0[1], k0[2], a, W, V1);
	// line 6: kFlat (0->curved, 1->flat), parameter
	readLine((char*)line, buffer_length, input);
	int kFlat_i;
	double k_param;
  sscanf(line, "%d %lf", &kFlat_i, &k_param);
	if(kFlat_i==0){
		kFlat=false;
		kCurved_k=k_param*Bohr_ang;
		printf("%32s = curved, k = %.2f Ang^-1 = %.2f Bohr^-1\n", "k plane", k_param, kCurved_k);
	}else if(kFlat_i==1){
		kFlat=true;
		kFlat_kz=k_param*Bohr_ang;
		printf("%32s = flat, k = %.2f Ang^-1 = %.2f Bohr^-1\n", "k plane", k_param, kFlat_kz);
	}else{
		printf("Invalid value of kFlat: %d\n", kFlat_i);
		return 0;
	}
	// line 7: surfaceConst (0->random, 1->const), parameter
	readLine((char*)line, buffer_length, input);
	int surfaceConst_i;
	double surfaceParam1;
	double surfaceParam2;
	sscanf(line, "%d %lf %lf", &surfaceConst_i, &surfaceParam1, &surfaceParam2);
	if(surfaceConst_i==0){
		surfaceConst=false;
		surfaceRandom_samples=int(surfaceParam1);
		printf("%32s = random, samples = %d\n", "Surface orientation", surfaceRandom_samples);
	}else if(surfaceConst_i==1){
		surfaceConst=true;
		surfaceConst_theta=surfaceParam1*M_PI/180.0;
		surfaceConst_phi=surfaceParam2*M_PI/180.0;
		printf("%32s = constant, theta = %.1f deg = %.2f rad, phi = %.1f deg = %.2f rad\n", "Surface orientation", surfaceParam1, surfaceConst_theta, surfaceParam2, surfaceConst_phi);
	}
	// line 8: kx range
	readLine((char*)line, buffer_length, input);
	double kxMin_ang;
	double kxMax_ang;
	sscanf(line, "%lf %lf %d", &kxMin_ang, &kxMax_ang, &kxCount);
	if(kxCount<=1){
		printf("Invalid value of kxCount: %d\n", kxCount);
		return 0;
	}
	kxMin=kxMin_ang*Bohr_ang;
	kxMax=kxMax_ang*Bohr_ang;
	printf("%32s = %.2f to %.2f, %d points\n", "kx range (Ang^-1)", kxMin_ang, kxMax_ang, kxCount);
	printf("%32s = %.2f to %.2f, %d points\n", "kx range (Bohr^-1)", kxMin, kxMax, kxCount);
	// line 9: ky range
	readLine((char*)line, buffer_length, input);
	double kyMin_ang;
	double kyMax_ang;
	sscanf(line, "%lf %lf %d", &kyMin_ang, &kyMax_ang, &kyCount);
	if(kyCount<=1){
		printf("Invalid value of kyCount: %d\n", kyCount);
		return 0;
	}
	kyMin=kyMin_ang*Bohr_ang;
	kyMax=kyMax_ang*Bohr_ang;
	printf("%32s = %.2f to %.2f, %d points\n", "ky range (Ang^-1)", kyMin_ang, kyMax_ang, kyCount);
	printf("%32s = %.2f to %.2f, %d points\n", "ky range (Bohr^-1)", kyMin, kyMax, kyCount);
	// line 10: e range
	readLine((char*)line, buffer_length, input);
	double eMin_eV;
	double eMax_eV;
	sscanf(line, "%lf %lf %d", &eMin_eV, &eMax_eV, &eCount);
	if(eCount<=1){
		printf("Invalid value of eCount: %d\n", eCount);
		return 0;
	}
	eMin=eMin_eV/Eh_eV-W;
	eMax=eMax_eV/Eh_eV-W;
	printf("%32s = %.2f to %.2f, %d points\n", "E range from EF (eV)", eMin_eV, eMax_eV, eCount);
	printf("%32s = %.2f to %.2f, %d points\n", "E range from Vacuum (eV)", eMin*Eh_eV, eMax*Eh_eV, eCount);
	printf("%32s = %.2f to %.2f, %d points\n", "E range from Vacuum (Eh)", eMin, eMax, eCount);
	// line 11: sigma k
	readLine((char*)line, buffer_length, input);
	double sigmak_ang;
	sscanf(line, "%lf", &sigmak_ang);
	sigmak=sigmak_ang*Bohr_ang;
	printf("%32s = %.2f Ang^-1 = %.2f Bohr^-1\n", "k broadening", sigmak_ang, sigmak);
	// line 12: sigma e
	readLine((char*)line, buffer_length, input);
	double sigmae_eV;
	sscanf(line, "%lf", &sigmae_eV);
	sigmae=sigmae_eV/Eh_eV;
	printf("%32s = %.2f eV = %.2f Eh\n", "E broadening", sigmae_eV, sigmae);
	// line 13: output
	readLine((char*)line, buffer_length, input);
	sscanf(line, "%s", outputName);
	printf("%32s = %s\n", "Output file", outputName);
	
	fclose(input);

	cout << "Finished loading input" << endl;
	dkx=(kxMax-kxMin)/(kxCount-1);
	dky=(kyMax-kyMin)/(kyCount-1);
	de=(eMax-eMin)/(eCount-1);

	printf("%32s = %.3f Ang^-1, %.3f Ang^-1, %.3f eV\n", "dkx, dky, de", dkx/Bohr_ang, dky/Bohr_ang, de*Eh_eV);
	printf("%32s = %.3f Bohr^-1, %.3f Bohr^-1, %.3f Eh\n", "dkx, dky, de", dkx, dky, de);

	// profile cube
	cout << "Making profile cube" << endl;
  start=chrono::system_clock::now();
	int kxCenter=ceil(sigmak*sigmaMax/dkx);
	int kyCenter=ceil(sigmak*sigmaMax/dky);
	int eCenter=ceil(sigmae*sigmaMax/de);

	// Gaussian broadened cube
	double cube[kxCenter*2+1][kyCenter*2+1][eCenter*2+1];
#pragma omp parallel private(j, k) firstprivate(kyCenter, eCenter, dkx, dky, de, sigmak, sigmae)
#pragma omp for
	for(i=0; i<=kxCenter; i++){
		for(j=0; j<=kyCenter; j++){
			for(k=0; k<=eCenter; k++){
				double weight=gauss(i*dkx, sigmak)*gauss(j*dky, sigmak)*gauss(k*de, sigmae)*dkx*dky*de;
				cube[kxCenter+i][kyCenter+j][eCenter+k]=weight;
				cube[kxCenter+i][kyCenter+j][eCenter-k]=weight;
				cube[kxCenter+i][kyCenter-j][eCenter+k]=weight;
				cube[kxCenter+i][kyCenter-j][eCenter-k]=weight;
				cube[kxCenter-i][kyCenter+j][eCenter+k]=weight;
				cube[kxCenter-i][kyCenter+j][eCenter-k]=weight;
				cube[kxCenter-i][kyCenter-j][eCenter+k]=weight;
				cube[kxCenter-i][kyCenter-j][eCenter-k]=weight;
			}
		}
	}
	end=chrono::system_clock::now();
	duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	printf("Finished: %.3f [ms]\n", duration);

	/*
	double sum=0.0;
	for(i=0; i<=kxCenter*2; i++){
		for(j=0; j<=kyCenter*2; j++){
			for(k=0; k<=eCenter*2; k++){
				sum+=cube[i][j][k]*dkx*dky*de;
			}
		}
	}
	printf("Weight sum: %.4f\n", sum);*/

	// Original dispersion
	double dispCube1[kxCount][kyCount][eCount]={0};
	// Refracted dispersion
	double dispCube2[kxCount][kyCount][eCount]={0};

	// Generate vectors perpendicular to surface(s)
	int nCount= (surfaceConst==true) ? 1 : surfaceRandom_samples;
	double nList[nCount][3];
	if(surfaceConst==true){
		nList[0][0]=sin(surfaceConst_theta)*cos(surfaceConst_phi);
		nList[0][1]=sin(surfaceConst_theta)*sin(surfaceConst_phi);
		nList[0][2]=cos(surfaceConst_theta);
	}else{
		// random surface generation
		random_device seed_gen;
		mt19937 mt(seed_gen());
		uniform_real_distribution<> xDist(-1.0, 1.0);
		uniform_real_distribution<> yDist(-1.0, 1.0);
		uniform_real_distribution<> zDist(0.0, 1.0);
		double vec[3];
		for(i=0; i<surfaceRandom_samples; i++){
			while(true){
				vec[0]=xDist(mt);
				vec[1]=yDist(mt);
				vec[2]=zDist(mt);

				double vec_length=sqrt(inProd(vec, vec));
				if(0.1 < vec_length && vec_length < 1){
					for(j=0; j<3; j++){
						nList[i][j]=vec[j]/vec_length;
					}
					break;
				}
			}
		}
	}
	/*
	for(i=0; i<nCount; i++){
		printf("n[%4d] = (%6.2f, %6.2f, %6.2f)\n", i, nList[i][0], nList[i][1], nList[i][2]);
		}*/
	
	cout << "Calculating dispersion" << endl;
  start=chrono::system_clock::now();
  // Calculate dispersions
#pragma omp parallel private(j, k) firstprivate(kFlat_kz, kCurved_k, a, W, V1, V0, cube, kxCenter, kyCenter, eCenter, nCount, nList)
# pragma omp for
	for(i=0; i<kxCount; i++){
		double k1[3];
		double kdiff[3];
		k1[0]=kxMin+dkx*i;
		kdiff[0]=k1[0]-k0[0];
		for(j=0; j<kyCount; j++){
			k1[1]=kyMin+dky*j;
			kdiff[1]=k1[1]-k0[1];
			if(kFlat==true){
				k1[2]=kFlat_kz;
			}else{
				k1[2]=sqrt(kCurved_k*kCurved_k-k1[0]*k1[0]-k1[1]*k1[1]);
			}
			kdiff[2]=k1[2]-k0[2];

			double esk=inProd(kdiff, kdiff)*a/2.0-W+V1;
			double epk=inProd(k1, k1)/2.0-V0;
			int eK=round((esk-eMin)/de);

			int i2, j2, k2;
			int i3, j3, k3;
			// cube1 (original)
			for(i2=-kxCenter; i2<=kxCenter; i2++){
				i3=i+i2;
				for(j2=-kyCenter; j2<=kyCenter; j2++){
					j3=j+j2;
					for(k2=-eCenter; k2<=eCenter; k2++){
						k3=eK+k2;
						if(0<=i3 && i3<kxCount && 0<=j3 && j3<kyCount && 0<=k3 && k3<eCount){
							dispCube1[i3][j3][k3]+=cube[i2+kxCenter][j2+kyCenter][k2+eCenter];
						}
					}
				}
			}

			// cube2 (refracted)
			double K[3];
			int iK, jK;
			for(k=0; k<nCount; k++){
				// calcK returns false is full reflection (no refracted wave) happens
				if(calcK(k1, epk, nList[k], K)){
					iK=round((K[0]-kxMin)/dkx);
					jK=round((K[1]-kyMin)/dky);
					for(i2=-kxCenter; i2<=kxCenter; i2++){
						i3=iK+i2;
						for(j2=-kyCenter; j2<=kyCenter; j2++){
							j3=jK+j2;
							for(k2=-eCenter; k2<=eCenter; k2++){
								k3=eK+k2;
								if(0<=i3 && i3<kxCount && 0<=j3 && j3<kyCount && 0<=k3 && k3<eCount){
									dispCube2[i3][j3][k3]+=cube[i2+kxCenter][j2+kyCenter][k2+eCenter];
								}
							}
						}
					}
				}
			}
		}
	}
	
	end=chrono::system_clock::now();
	duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	printf("Finished: %.3f [ms]\n", duration);

	// output
	H5File output(outputName, H5F_ACC_TRUNC);

	Group rootG(output.openGroup("/"));

	w_data_3d(rootG, "Original", kxCount, kyCount, eCount, (double***) dispCube1);
	w_data_3d(rootG, "Refracted", kxCount, kyCount, eCount, (double***) dispCube2);

	double offset[3]={kxMin_ang, kyMin_ang, eMin_eV};
	double delta[3]={dkx/Bohr_ang, dky/Bohr_ang, de*Eh_eV};
	int size[3]={kxCount, kyCount, eCount};

	w_att_1d(rootG, "Offset", 3, offset);
	w_att_1d(rootG, "Delta", 3, delta);
	w_att_1i(rootG, "Size", 3, size);
	w_att_double(rootG, "W", W_eV);
	w_att_double(rootG, "V0", V0_eV);
	w_att_double(rootG, "V1", V1_eV);
	w_att_double(rootG, "a", a);
	w_att_double(rootG, "sigmak", sigmak_ang);
	w_att_double(rootG, "sigmae", sigmae_eV);

	double k0_ang[3];
	for(i=0; i<3; i++){
		k0_ang[i]=k0[i]/Bohr_ang;
	}
	w_att_1d(rootG, "k0", 3, k0_ang);

	if(kFlat==true){
		w_att_str(rootG, "kPlane", string("Flat"));
		w_att_double(rootG, "kPlane_kz", k_param);
	}else{
		w_att_str(rootG, "kPlane", string("Curved"));
		w_att_double(rootG, "kPlane_k", k_param);
	}

	if(surfaceConst==true){
		w_att_str(rootG, "Surface", string("Constant"));
		w_att_double(rootG, "Surface_theta", surfaceParam1);
		w_att_double(rootG, "Surface_phi", surfaceParam2);
	}else{
		w_att_str(rootG, "Surface", string("Random"));
		w_att_int(rootG, "Surface_samples", surfaceRandom_samples);
	}
	
	cout << "Output finished" << endl;
}
