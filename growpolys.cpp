
///////////////////////////////////////////////////////////////////////////////////////////////////////
//  Date:  18 July 2016 
//  Description:    rosenbluth sampling of growing ligands on a spheroid
//
//  Usage Syntax:
//
///////////////////////////////////////////////////////////////////////////////////////////////////////



#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include "genarray.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

// Global Variables
const double PI = 3.14159265359;
gsl_rng * mrRand;
double mySeed = 1;
double TEMP = 1;
double mass = 1.0;
int totalSteps = 10000000;
int writeEvery = 1000;
double xlo, xhi, ylo, yhi, zlo, zhi;
int ngraft;
int totalMono;
string spheroidFile;
double r_cut = 1;
double sigma = 0.5;
double epsilon = 10;
double theta_0 = PI;
double k_angle = 2;
//double elp_a = 4.189;
//double elp_c = 10.4725;
double elp_a = 7.5;
double elp_c = 3.0;

// Rosenbluth params:
int nTrial = 10;


class Mono {
    public:
	double x, y, z;
//    int next;			    // The next atomID in the polymer
};

class Poly {
    public: 
	int nMono;
	genarray < Mono > chain;
};

// Functions
void readConfig (string readName, genarray<double> &atomPositions);

void paramReader (string fileName);

void initPolys ( genarray<double> &spheroidPos, genarray < Poly > &brush, genarray < Poly > &trialMonos, genarray < Poly > &trialRose); 

double calcAddOneEn(Mono trialMono, int whichPoly, genarray< Poly > &brush); 

void genTrialPts( int whichPoly, genarray< Poly > &brush, genarray < Poly > &trialMonos, int previousPoly); 

//void genAllTrialPts( int whichPoly, genarray< Poly > &brush, genarray < Poly > &trialMonos); 

void calcRosenbluth (genarray< Poly > &brush, genarray < Poly > &trialMonos, genarray < Poly > &trialRose, int previousPoly); 

double addMonomer (genarray< Poly > &brush, genarray < Poly > &trialMonos, genarray < Poly > &trialRose, int &previousPoly); 

void writeBrush (string outName, genarray< Poly > &brush); 

void writeLammps (string outName, genarray< Poly > &brush);

void clearBrush (genarray< Poly > &brush);

void calcPolarAngle (genarray< Poly > &brush, genarray< double > &polarAngles); 

// Main

int main(int argv, char *argc[]) {
    // Read parameter file
    string paramFile(argc[1]);
    cout << "reading parameters from: " << paramFile << endl;
    paramReader (paramFile);

    // Seed!!!! arfgljkad;lkfjadsfjsda;'
    mrRand = gsl_rng_alloc (gsl_rng_taus);
    cout << "mySeed is: " << mySeed << endl;
    gsl_rng_set(mrRand, mySeed);

    ofstream myOut;
    myOut.open("brush.xyz", ios::out);
    myOut.close();
    myOut.open("brush.lammpstrj", ios::out);
    myOut.close();

    genarray<double> atomPositions;
    readConfig (spheroidFile, atomPositions);




    genarray< double > polarAngles;
    genarray< Poly > brush;
    genarray< Poly > trialMonos;
    genarray< Poly > trialRose;
    initPolys (atomPositions, brush, trialMonos, trialRose);
    cout << "finished initPolys " << "\n";
    


    double rose_w = 1;
    ofstream outer, outer1;
    outer.open("r_fac.dat", ios::out);
    outer1.open("anglevlen.dat", ios::out);
    int previous;
    for (int j = 0; j < 1000; j++)  {
	previous = -1;
	for (int i = 0; i < totalMono; i++) {

	    // generate all trial points
	    for (int j = 0; j < brush.length(); j++) {
//    	    cout << "\tgenerating for poly " << j << endl;
		genTrialPts (j, brush, trialMonos, previous);
	    }

	    // calc rosenbluth
//	    cout << "\tcalcRosenbluth!" << endl;
	    calcRosenbluth (brush, trialMonos, trialRose, previous);

	    // add monomer
//	    cout << "adding monomer: " << i << endl;
	    rose_w *= addMonomer (brush, trialMonos, trialRose, previous);

	    // write:
	    if (i%500 == 0)
		cout << "\tfinished addMonomer " << i << "\n";
	}

	if (j%1 ==0 ) {
	    cout << "finished brush: " << j << endl;
	}

	writeBrush ("brush.xyz", brush);
	writeLammps ("brush.lammpstrj", brush);
	outer << j << "\t" << rose_w << "\n";
	rose_w = 1;

	// Analyze:
	calcPolarAngle(brush, polarAngles);
	for (int i = 0; i < brush.length(); i++) {
	    outer1 << brush(i).nMono-1 << "\t" << polarAngles(i) << "\n";
	}

	// Reset brush:
	clearBrush(brush);
    }


    return 0;
}

void paramReader (string fileName)
{
// Read in Parameters
    ifstream in;
    in.open(fileName.c_str(), ios::in);
    string junk1;
    in >> junk1 >> totalMono;
//    in >> junk1 >> writeEvery;
    in >> junk1 >> mySeed;
    in >> junk1 >> TEMP;
    in >> junk1 >> ngraft;
    in >> junk1 >> spheroidFile;
    in >> junk1 >> nTrial;
    in >> junk1 >> elp_a;
    in >> junk1 >> elp_c;
    in.close();

}


void readConfig (string readName, genarray<double> &atomPositions)
{
    int ctr = 0;
    int atomctr = 0;
    double x,y,z;
    double theta;
    int type;
    cout << "reading file " << readName << " for restart." << endl;
    ifstream in(readName.c_str());
    string junk, myLine;
    atomPositions.resize(ngraft*3);

    while (getline(in, myLine)) {
	if (myLine != "") {
	    stringstream ss(myLine);
	    if (ctr == 5) {
		ss >> xlo >> xhi;
	    }
	    if (ctr == 6) {
		ss >> ylo >> yhi;
	    }
	    if (ctr == 7) {
		ss >> zlo >> zhi;
	    }
	    if (ctr > 8) {
		ss >> junk >> type >> x >> y >> z;
		// Each frame has natoms*3 elements.
		if (type == 2) {
		    atomPositions(atomctr*3) = x*(xhi-xlo) + xlo;
		    atomPositions(atomctr*3+1) = y*(yhi-ylo) + ylo;
		    atomPositions(atomctr*3+2) = z*(zhi-zlo) + zlo;
		    atomctr++;
		}
	    }
	    ctr++;
	}
    }


}

void initPolys ( genarray<double> &spheroidPos, genarray < Poly > &brush, genarray < Poly > &trialMonos, genarray < Poly > &trialRose) { 

    int maxMono = 2000;

    // Init monomer
    Mono dummy;
    dummy.x = 0;
    dummy.y = 0;
    dummy.z = 0;
//    dummy.next = -1;
    
    // Init single polymer
    Poly tempPoly;
    tempPoly.nMono = 0;
    tempPoly.chain.resize(maxMono);
    cout << "0 \n";
    for (int i = 0; i < maxMono; i++) {
	tempPoly.chain(i) = dummy;
    }
    cout << "1 \n";

    // Init all polymers, set first monomer to the spheroid monomers.
    cout << "ngraft is: " << ngraft << endl;
    int ic;
    brush.resize(ngraft);
    for (int i = 0; i < ngraft; i++) {
	brush(i) = tempPoly;
	brush(i).nMono = 1;

	ic = i*3;
	dummy.x = spheroidPos(ic);
	dummy.y = spheroidPos(ic+1);
	dummy.z = spheroidPos(ic+2);
	brush(i).chain(0) = dummy;
    }
    cout << "2 \n";

    // Init trialMonos, trialRose:
    dummy.x = 0;
    dummy.y = 0;
    dummy.z = 0;
    tempPoly.chain.resize(nTrial);
    for (int i = 0; i < nTrial; i++) {
	tempPoly.chain(i) = dummy;
    }
    cout << "3 \n";

    tempPoly.nMono = nTrial;
    trialMonos.resize(ngraft);
    trialRose.resize(ngraft);
    cout << "4 \n";
    for (int i = 0; i < ngraft; i++) {
	trialMonos(i) = tempPoly;
	trialRose(i) = tempPoly;
    }

}


double calcEnergy(genarray< Poly > &brush) {
// Soft potential - V (r) = epsilon*(1 - r/sigma)^2
// Harmonic angles - V_a (r) = k*(theta - theta_0)^2

    double r_cutsq = r_cut*r_cut;
    double e = 0;
    int startMono;
    int ic, jc, nMono1, nMono2;
    double dx, dy, dz, drsq, dr;
    double e_one;
    Mono m1, m2, m3;

    // Outer loop - monomer 1
    for (int i = 0; i < brush.length(); i++) {
	nMono1 = brush(i).nMono;
	for (int j = 0; j < nMono1; j++) {
	    m1 = brush(i).chain(j);

	    // Inner loop - monomer 2:
	    for (int k = i; k < brush.length(); k++) {
		if (k == i) 
		    startMono = j + 1;
		else
		    startMono = 0;
		nMono2 = brush(k).chain.length();
		for (int l = startMono; l < nMono2; l++) {
		    m2 = brush(k).chain(l);

		    // Pair energy: 
		    dx = m1.x - m2.x;
		    dy = m1.y - m2.y;
		    dz = m1.z - m2.z;
		    drsq = dx*dx+dy*dy+dz*dz;
		    if (drsq < r_cutsq){
			dr = sqrt(dr);
			e_one = (1-dr/sigma);
			e_one = e_one*e_one;
			e += epsilon*e_one;
		    }

		}
	    }
	}
    }

    // Angles:
    double dx1, dy1, dz1, dx2, dy2, dz2, rsq1, rsq2, r1, r2;
    double c, s, dtheta;
    int jc1, jc2, jc3;

    for (int i = 0; i < brush.length(); i++) {
	nMono1 = brush(i).nMono;
	for (int j = 0; j < (nMono1-2); j++) {
	    m1 = brush(i).chain(j);
	    m2 = brush(i).chain(j+1);
	    m3 = brush(i).chain(j+2);

	    dx1 = m1.x - m2.x;
	    dy1 = m1.y - m2.y;
	    dz1 = m1.z - m2.z;
	    rsq1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
	    r1 = sqrt(r1);

	    dx2 = m2.x - m3.x;
	    dy2 = m2.y - m3.y;
	    dz2 = m2.z - m3.z;
	    rsq2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
	    r2 = sqrt(r2);

	    c = dx1*dx2 + dy1*dy2 + dz1*dz2;
	    c /= r1*r2;
	    if (c > 1.0) c = 1.0;
	    if (c < -1.0) c = -1.0;
	    dtheta = acos(c) - theta_0;
	    e += k_angle * dtheta*dtheta;
	    
	}
    }

    return e;
}

double calcAddOneEn(Mono trialMono, int whichPoly, genarray< Poly > &brush) {
// Calculate pair energies and angle of adding trialMono to polymer whichPoly. 
//    cout << "\t in calcAddOneEn! " << endl;
    double asqinv = 1/(elp_a*elp_a);
    double csqinv = 1/(elp_c*elp_c);

    double r_cutsq = r_cut*r_cut;
    double e = 0;
    int startMono;
    int ic, jc, nMono1, nMono2;
    double dx, dy, dz, drsq, dr;
    double e_one;
    Mono m1, m2, m3;

    m2 = trialMono;
    // temp
    double elp_rad = (m2.x*m2.x + m2.y*m2.y)*asqinv + m2.z*m2.z*csqinv;
    if (elp_rad < 1) {
	e = 10000;
    }
    
    else {
	Mono tempm = brush(whichPoly).chain( brush(whichPoly).nMono-1 );
    //    cout << "\t\t\t last Mono: " << tempm.x << " " << tempm.y << " " << tempm.z << endl;
    //    cout << "\t\t\t trialMono: " << m2.x << " " << m2.y << " " << m2.z << endl;
	for (int i = 0; i < brush.length(); i++) {
	    nMono1 = brush(i).nMono;
	    for (int j = 0; j < nMono1; j++) {
		m1 = brush(i).chain(j);

		// Pair energy: 
		dx = m1.x - m2.x;
		dy = m1.y - m2.y;
		dz = m1.z - m2.z;
		drsq = dx*dx+dy*dy+dz*dz;
		if (drsq < r_cutsq) {
    //		cout << "\t\t\t: overlap! x,y,z: " << m1.x << " " << m1.y << " " << m1.z << endl;
		    dr = sqrt(drsq);
		    e_one = (1-dr/sigma);
		    e_one = e_one*e_one;
		    if (isnan(e_one)) {
			cout << "BAD PAIR! " << m1.x << " " << m1.y << " " << m1.z << endl;
			cout << "dr: " << dr << endl;
			cout << "i, j: " << i << " " << j << endl;
			cout << "nMono1: " << nMono1 << endl;
		    }
		    e += epsilon*e_one;
		}
	    }
	}

	// Angles:
	int len = brush(whichPoly).nMono;
	double e_angle = 0;
	if (len > 1) {
	    double dx1, dy1, dz1, dx2, dy2, dz2, rsq1, rsq2, r1, r2;
	    double c, s, dtheta;
	    int jc1, jc2, jc3;

	    m1 = brush(whichPoly).chain(len-2);
	    m2 = brush(whichPoly).chain(len-1);
	    m3 = trialMono;

	    dx1 = m1.x - m2.x;
	    dy1 = m1.y - m2.y;
	    dz1 = m1.z - m2.z;
	    rsq1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
	    r1 = sqrt(rsq1);

	    dx2 = m3.x - m2.x;
	    dy2 = m3.y - m2.y;
	    dz2 = m3.z - m2.z;
	    rsq2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
	    r2 = sqrt(rsq2);

	    c = dx1*dx2 + dy1*dy2 + dz1*dz2;
	    c /= r1*r2;
	    if (c > 1.0) c = 1.0;
	    if (c < -1.0) c = -1.0;
	    dtheta = acos(c) - theta_0;
	    e_angle = k_angle * dtheta*dtheta;
	    if (e_angle > 40) {
		cout << "WEIRD ANGLE!!!!" << endl;
		cout << "len: " << len << endl;
		cout << "r1, r2, c: " << r1 << " " << r2 << " " << c << endl;
		cout << "\tacos(c): " << acos(c) << endl;
		cout << "\tdtheta: " << dtheta << endl;
		cout << "\te_angle: " << e_angle << endl;
		cout << "\t mono1: " << m1.x <<  " " << m1.y << " " << m1.z << endl;
		cout << "\t mono2: " << m2.x <<  " " << m2.y << " " << m2.z << endl;
		cout << "\t mono3: " << m3.x <<  " " << m3.y << " " << m3.z << endl;
		cout << "d1s: " << dx1 << " " << dy1 << " " << dz1 << endl;
		cout << "d2s: " << dx2 << " " << dy2 << " " << dz2 << endl;

	    }
	    e += e_angle;
	}
    }

    return e;
}

void genTrialPts( int whichPoly, genarray< Poly > &brush, genarray < Poly > &trialMonos, int previousPoly) {
// Generates trial additions to polymer whichPoly
//    cout << "\t in genTrialPts \n";

    double dx, dy, dz, mag;
    Mono m1, m2, mt;

//    cout << "whichPoly: " << whichPoly << endl;
    int nMono = brush(whichPoly).nMono;
//    cout << "nMono-1: " << (nMono-1) << endl;
    if (nMono > 1) {
	m1 = brush(whichPoly).chain(nMono-2);
    }
    else {
	m1.x = 0;
	m1.y = 0;
	m1.z = 0;
    }
    m2 = brush(whichPoly).chain(nMono-1);
    Mono cool;

    // Check if distance between the trialPoly's tail and previousPoly tail monomers is close enough to warrant a new trial addition:
    double drsq = 0;
    double r_cutsq = 10;
    if (previousPoly > -1) {
	Mono mprev = brush(previousPoly).chain(brush(previousPoly).nMono-1);
	dx = m2.x - mprev.x;
	dy = m2.y - mprev.y;
	dz = m2.z - mprev.z;
	drsq = dx*dx+dy*dy+dz*dz;
    }
    else {
	drsq = 0;
    }
    if ( drsq < r_cutsq ) {

	// First trial move is 180 degrees (straight)
	dx = m2.x - m1.x;
	dy = m2.y - m1.y;
	dz = m2.z - m1.z;
	mag = sqrt(dx*dx + dy*dy + dz*dz);

	mt.x = m2.x + dx/mag;
	mt.y = m2.y + dy/mag;
	mt.z = m2.z + dz/mag;

	trialMonos(whichPoly).chain(0) = mt;
	
	// Generate random points in the hemisphere above the last monomer's bond-vector:
	double myangle = 0;
	double x, y, z, u, v, theta, phi;
    //    cout << "ntrial: " << nTrial << endl;

	double dbgmag;
	double anglecut;
	double elp_rad;
	double asqinv = 1/(elp_a*elp_a);
	double csqinv = 1/(elp_c*elp_c);
	double cc;
	for (int i = 1; i < nTrial; i++) {
	    u = gsl_rng_uniform(mrRand);
	    v = gsl_rng_uniform(mrRand);
	    theta = 2*PI*u;
	    phi = acos(2*v-1.0);
	    x = cos(theta)*sin(phi);
	    y = sin(theta)*sin(phi);
	    z = cos(phi);

    //	myangle = acos( (-dx*x - dy*y - dz*z)/mag );		// The negative signs account for the direction of the bond vectors
	    cc = (-dx*x - dy*y - dz*z)/mag;
	    if (cc > 1.0) cc = 1.0;
	    if (cc < -1.0) cc = -1.0;
	    myangle = acos( cc );

	    
	    mt.x = m2.x + x;
	    mt.y = m2.y + y;
	    mt.z = m2.z + z;


	    anglecut = PI/3.0;
	    if ( fabs(myangle - theta_0) < anglecut) {
		trialMonos(whichPoly).chain(i) = mt;

		// Check angle:
		double dx1, dy1, dz1, dx2, dy2, dz2, c, r1, r2, rsq1, rsq2;
		dx1 = m2.x - m1.x;
		dy1 = m2.y - m1.y;
		dz1 = m2.z - m1.z;
		rsq1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
		r1 = sqrt(rsq1);

		dx2 = m2.x - mt.x;
		dy2 = m2.y - mt.y;
		dz2 = m2.z - mt.z;
		rsq2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
		r2 = sqrt(rsq2);

		c = dx1*dx2 + dy1*dy2 + dz1*dz2;
		c /= r1*r2;
		if (c > 1.0) c = 1.0;
		if (c < -1.0) c = -1.0;
		if (fabs(myangle - acos(c)) > 0.0001) {
		    cout << "angle, anglecheck: " << myangle << " " << acos(c) << endl;
		    cout << "\tcc: " << cc << endl;
		    cout << "\td: " << dx << " " << dy << " " << dz << endl;
		    cout << "\td/mag: " << dx/mag << " " << dy/mag << " " << dz/mag << endl;
		    cout << "\ttrial: " << x << " " << y << " " << z << endl;
		    cout << "\td1s: " << dx1 << " " << dy1 << " " << dz1 << endl;
		    cout << "\td2s: " << dx2 << " " << dy2 << " " << dz2 << endl;
		    cout << "\tr1, r2, c: " << r1 << " " << r2 << " " << c << endl;
		}
	    }
	    else {
		i--;
	    }
	}
    }

//    cout << "finished genTrialPts\n";

}


void calcRosenbluth (genarray< Poly > &brush, genarray < Poly > &trialMonos, genarray < Poly > &trialRose, int previousPoly) {
// Loop: calculate the energy of each trialMono:
//  - use a Mono object for rosenbluth.  Mono.x - energy of the trial addition.  Mono.y - rosenbluth factor of the same. 

//    cout << "in calcRosenbluth! " << endl;

    double en, rose;
    double dx, dy, dz;
    double drsq;
    double r_cutsq = 10;
    double beta = -1/TEMP;
    Mono store;		
    store.z = 0;
    double temp = 0;
    //
    Mono dbg;
    //
    for (int i = 0; i < brush.length(); i++) {
	if (previousPoly > -1) {
	    Mono mprev = brush(previousPoly).chain(brush(previousPoly).nMono-1);
	    Mono mtail = brush(i).chain(brush(i).nMono-1);
	    dx = mtail.x - mprev.x;
	    dy = mtail.y - mprev.y;
	    dz = mtail.z - mprev.z;
	    drsq = dx*dx+dy*dy+dz*dz;
	}
	else {
	    drsq = 0;
	}
	if (drsq < r_cutsq) {
	    for (int j = 0; j < nTrial; j++) {
		dbg = trialMonos(i).chain(j);
		en = calcAddOneEn(trialMonos(i).chain(j), i, brush);
		if (isnan(en)) {
		    cout << " BAD ENERGY! trialmono: " << dbg.x << " " << dbg.y << " " << dbg.z << endl;
		}
		temp = beta*en;
		rose = exp(temp);
    //	    cout << "\t\t  en, rose: " << en << " " <<  rose << " " << endl;
		store.x = en;
		store.y = rose;
		trialRose(i).chain(j) = store;
	    }
	}
    }
//    cout << "finished calcRosenbluth" << endl;

}

double addMonomer (genarray< Poly > &brush, genarray < Poly > &trialMonos, genarray < Poly > &trialRose, int &previousPoly) {
// Probabilities (rosenbluth factors) are stored in trialRose.  Select one of the trialMonomers from brush.

    double w_move = 0;					// rosenbluth factor for the (entire) trial move
    for (int i = 0; i < trialRose.length(); i++) {
	for (int j = 0; j < nTrial; j++) {
	    w_move = w_move + trialRose(i).chain(j).y;

	}
    }
    if (w_move < 10) {
	cout << "bad w_move!!!: \n";
	for (int i = 0; i < trialRose.length(); i++) {
	    cout << "poly " << i << ": " << endl;
	    for (int j = 0; j < nTrial; j++) {
		cout << "\t trial, en, rose: " << j << ", " << trialRose(i).chain(j).x << ", " <<  trialRose(i).chain(j).y << endl;
	    }
	}
    }

//    cout << "w_move is: " << w_move << endl << endl;
    
//    int is, js;		YOU HAVE A BUG!!! SOME HOW w_move is zero!!!
    int is = 0;
    int js = 0; 
    // SELECT (Frenkel and Smit Algorithm 41):
    double ws = gsl_rng_uniform(mrRand)*w_move;
    double cumw = trialRose(0).chain(0).y;
    int start;
//    cout << "ws: " << ws << endl;
    for (int i = 0; i < trialRose.length(); i++) {
	if (i == 0) 
	    start = 1;
	else
	    start = 0;
	for (int j = start; j < nTrial; j++) {
//	    cout << "\tcumw is: " << cumw << endl;
	    if (cumw > ws) {
		is = i;
		js = j;
//		if (js > 10) {
//		    cout << "js, nTrial, len: " << js << " " << nTrial << " " << brush(i).chain.length() << endl;
//
//		}
//		if (is > 1211) {
//		    cout << "is, trialRose.len: " << is << " " << trialRose.length() << endl;
//
//		}
		i = trialRose.length()+1;
		j = nTrial+1;
	    }
	    else
		cumw = cumw + trialRose(i).chain(j).y;
	}
    }


//    cout << "is: " << is << endl;
    int end = brush(is).nMono;
//    cout << "is, js, end: " << is << " " << js << " " << end << endl;
    brush(is).chain(end) = trialMonos(is).chain(js);
    brush(is).nMono = end + 1;

    previousPoly = is;
    Mono myMono = trialMonos(is).chain(js);
    Mono prev = brush(is).chain(end-1);
    double mag1 = prev.x*prev.x + prev.y*prev.y + prev.z*prev.z;
    mag1 = sqrt(mag1);
    double mag2 = myMono.x*myMono.x + myMono.y*myMono.y + myMono.z*myMono.z;
    mag2 = sqrt(mag2);
//    cout << "\tpicking i, j: " << is << " " << js << endl;
//    cout << "\ttail mono: " << prev.x << " " << prev.y << " " << prev.z << ", " << mag1 << endl;
//    cout << "\tadding mono: " << myMono.x << " " << myMono.y << " " << myMono.z << ", " << mag2 << endl;
//    cout << "\n";

    return w_move;
}


void writeBrush (string outName, genarray< Poly > &brush) {

    int nxyz = 2*(ngraft+totalMono);
    ofstream myOut;
    myOut.open(outName.c_str(), ios::app);
    myOut << nxyz << "\n\n";
    Mono m;
    int ctr = 0;
    for (int i = 0; i < brush.length(); i++) {
	for (int j = 0; j < brush(i).nMono; j++) {
	    m = brush(i).chain(j);
	    myOut << "O " << m.x << " " << m.y << " " << m.z << "\n";
	    myOut << "N 0 0 0 \n";
	    ctr++;
	}
    }
    int dummies = (ngraft+totalMono) - ctr;
    for (int i = 0; i < dummies; i++) {
	myOut << "O 0 0 0 \n";
	myOut << "N 0 0 0 \n";
    }
//    myOut << "\n";
//    myOut << "ITEM: TIMESTEP\n" << "0" << "\nITEM: NUMBER OF ATOMS\n";
//    myOut << natoms << "\n" << "ITEM: BOX BOUNDS pp pp pp\n";
//    myOut << -boxrad << " " << boxrad << "\n" << -boxrad << " " << boxrad << "\n" << -boxrad << " " << boxrad << "\n";
//    myOut << "ITEM: ATOMS id type xs ys zs q\n";
//    for (int i = 0; i < natoms; i++) {
//	myOut << i+1 << " 1 " << (atomPositions(i*3)+boxrad)*length_inv << " " << (atomPositions(i*3+1)+boxrad)*length_inv;
//	myOut << " 0 " << atomPositions(i*3+2) << "\n";	    // z-coordinate and angle
//    }


}


void clearBrush (genarray< Poly > &brush) {
    for (int i = 0; i < brush.length(); i++) {
	brush(i).nMono = 1;
    }
}

void calcPolarAngle (genarray< Poly > &brush, genarray< double > &polarAngles) {
    
    double c, r2;
    Mono m1, m2;
    m1.x = 0;
    m1.y = 0;
    m1.z = 1;

    polarAngles.resize(brush.length());

    for (int i = 0; i < brush.length(); i++) {
	m2 = brush(i).chain(0);
	r2 = m2.x*m2.x + m2.y*m2.y + m2.z*m2.z;
	r2 = sqrt(r2);
	c = m2.z/r2;
	if (c > 1.0) c = 1.0;
	if (c < -1.0) c = -1.0;
	polarAngles(i) = acos(c);
    }
}


void writeLammps (string outName, genarray< Poly > &brush) {

    int nxyz = 2*(ngraft+totalMono);
    ofstream myOut;
    myOut.open(outName.c_str(), ios::app);
    Mono m;
    int ctr = 0;

    double boxrad = 20;
    double length_inv = 1/40.0;
//    myOut << "\n";
    myOut << "ITEM: TIMESTEP\n" << "0" << "\nITEM: NUMBER OF ATOMS\n";
    myOut << brush.length() << "\n" << "ITEM: BOX BOUNDS pp pp pp\n";
    myOut << -boxrad << " " << boxrad << "\n" << -boxrad << " " << boxrad << "\n" << -boxrad << " " << boxrad << "\n";
    myOut << "ITEM: ATOMS id type xs ys zs vz\n";
    for (int i = 0; i < brush.length(); i++) {
	m = brush(i).chain(0);
	myOut << i+1 << " 1 " << (m.x+boxrad)*length_inv << " " << (m.y+boxrad)*length_inv << " " << (m.z+boxrad)*length_inv;
	myOut << " " << brush(i).nMono << "\n";

    }
    myOut.close();
}




