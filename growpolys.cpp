
///////////////////////////////////////////////////////////////////////////////////////////////////////
//  Date:  18 July 2016 
//	    29 Dec 2016
//  Description:    rosenbluth sampling of growing ligands on a planar surface (PBCs in x-y)
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
//#include <gsl_rng.h>
//#include <gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

const double smalleps = 1e-6;

// Global Variables
const double PI = 3.14159265359;
gsl_rng * mrRand;
double mySeed = 1;
double TEMP = 1;
double mass = 1.0;
int totalSteps = 10000000;
int writeEvery = 1000;
double xlo, xhi, ylo, yhi, zlo, zhi;
double xhalf, yhalf, xbox, ybox;
int ngraft;
int totalMono;
//string spheroidFile;
string thomsonFile;
double r_cut = 2.1;
double r_cutsq = r_cut*r_cut;
double sigma = 1.0;
double sigma_sq = sigma*sigma;
double epsilon = 10;
double theta_0 = PI;
double k_angle = 2;
//double elp_a = 4.189;
//double elp_c = 10.4725;
double elp_a = 7.5;
double elp_c = 3.0;
double iterations = 1;
double lattice_size = 1.0;
double box_len = 10;
int writeGrowth = 1000;

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

double calcAddOneEn(Mono &trialMono, int whichPoly, genarray< Poly > &brush); 

void genTrialPts( int whichPoly, genarray< Poly > &brush, genarray < Poly > &trialMonos, int previousPoly, 
	genarray< Poly > &trialThomson); 

void calcRosenbluth (genarray< Poly > &brush, genarray < Poly > &trialMonos, genarray < Poly > &trialRose, int previousPoly); 

double addMonomer (genarray< Poly > &brush, genarray < Poly > &trialMonos, genarray < Poly > &trialRose, int &previousPoly); 

void writeBrush (string outName, genarray< Poly > &brush); 

void writeLammps (string outName, genarray< Poly > &brush);

void clearBrush (genarray< Poly > &brush);

void calcPolarAngle (genarray< Poly > &brush, genarray< double > &polarAngles); 

void readThomson (string inName, genarray< Poly > &trialThomson);

void makeSquareLattice (genarray <double> &atomPositions);

double monoDistSq (Mono m1, Mono m2);

void rotate_pt (double &x1, double &y1, double &z1, double theta, double ax, double ay, double az); 

void quat_mult (double &a1, double &b1, double &c1, double &d1, double a2, double b2, double c2, double d2);


void genTrialPts2( int whichPoly, genarray< Poly > &brush, genarray < Poly > &trialMonos, int previousPoly, 
	genarray< Poly > &trialThomson, genarray< Poly > &trialRose); 

void genTrialPts2dbg( int whichPoly, genarray< Poly > &brush, genarray < Poly > &trialMonos, int previousPoly, 
	genarray< Poly > &trialThomson, genarray< Poly > &trialRose, int &dbgctr); 

void writeLengths (string outName, genarray< Poly > &brush, int growthStep); 

// Main

int main(int argv, char *argc[]) {
    cout << "Command line arguments: paramFile thomsonFile" << endl;
    // Read parameter file
    string paramFile(argc[1]);
    cout << "reading parameters from: " << paramFile << endl;
    paramReader (paramFile);

    mrRand = gsl_rng_alloc (gsl_rng_taus);
    cout << "mySeed is: " << mySeed << endl;
    gsl_rng_set(mrRand, mySeed);

    ofstream myOut;
    myOut.open("brush.xyz", ios::out);
    myOut.close();
    myOut.open("brush.lammpstrj", ios::out);
    myOut.close();

    genarray< Poly > trialThomson;
    readThomson(thomsonFile, trialThomson);
    Poly d;
    for (int i = 0; i < trialThomson.length(); i++) {
	cout << "sphere " << i << endl;
	d = trialThomson(i);
	for (int j = 0; j < d.nMono; j++) {
	    cout << "\t" << d.chain(j).x << " " << d.chain(j).y << " " << d.chain(j).z << endl;
	}
    }


    genarray< double > atomPositions;
    makeSquareLattice (atomPositions);

    genarray< double > polarAngles;
    genarray< Poly > brush;
    genarray< Poly > trialMonos;
    genarray< Poly > trialRose;
    initPolys (atomPositions, brush, trialMonos, trialRose);
    cout << "finished initPolys " << "\n";
    
    genarray<double> hist(1000);
    genarray<double> histctr(1000);
    for (int i = 0; i < hist.length(); i++) {
	hist(i) = 0;
    }


    double rose_w = 1;
    ofstream outer, outer1, outer2;
    outer.open("r_fac.dat", ios::out);
    outer1.open("anglevlen.dat", ios::out);
    outer2.open("all_lens.dat", ios::out);
    int previous;
    int dbgctr = 0;
    for (int j = 0; j < iterations; j++)  {
	dbgctr = 0;
	previous = -1;
	for (int i = 0; i < totalMono; i++) {

	    // generate all trial points
	    for (int k = 0; k < brush.length(); k++) {
//    	    cout << "\tgenerating for poly " << j << endl;
//		genTrialPts2 (j, brush, trialMonos, previous, trialThomson, trialRose);
		genTrialPts2dbg (k, brush, trialMonos, previous, trialThomson, trialRose, dbgctr);
	    }

	    // calc rosenbluth
//	    cout << "\tcalcRosenbluth!" << endl;
//	    calcRosenbluth (brush, trialMonos, trialRose, previous);

	    // add monomer
//	    cout << "adding monomer: " << i << endl;
	    rose_w *= addMonomer (brush, trialMonos, trialRose, previous);

	    // write:
	    if (i%1000 == 0) {
//		cout << "\tfinished addMonomer " << i << "\n";
//		cout << "\t\tcomputed " << dbgctr << " monomers: " << endl;
		cout << "addMonomer, computed: " << i << " " << dbgctr << endl;
		dbgctr = 0;
	    }
	    if (i%writeGrowth == 0 && (i > 0)) {
		writeLengths("all_lens.dat", brush, i);
	    }	
		
	}

	if (j%1 ==0 ) {
	    cout << "finished brush: " << j << endl;
	}

	writeBrush ("brush.xyz", brush);
	writeLammps ("brush.lammpstrj", brush);
	outer << j << "\t" << rose_w << "\n";
	rose_w = 1;

	// Analyze:
//	calcPolarAngle(brush, polarAngles);
	for (int k = 0; k < brush.length(); k++) {
//	    outer1 << brush(i).nMono-1 << "\t" << polarAngles(i) << "\n";
	    hist(brush(k).nMono) = hist(brush(k).nMono) + 1;
	}
	

	// Reset brush:
	clearBrush(brush);
    }
    ofstream outman;
    outman.open("planar_hist.dat", ios::out);
    for (int i = 0; i < hist.length(); i++) {
	outman << i << "\t" << hist(i)/double(iterations) << endl;
    }

    return 0;
    outer1.close();
    outer2.close();
}

void paramReader (string fileName)
{
// Read in Parameters
    ifstream in;
    in.open(fileName.c_str(), ios::in);
    string junk1;
    in >> junk1 >> thomsonFile;
    in >> junk1 >> nTrial;
    in >> junk1 >> iterations;
    in >> junk1 >> totalMono;
//    in >> junk1 >> writeEvery;
    in >> junk1 >> mySeed;
    in >> junk1 >> TEMP;
    in >> junk1 >> box_len;
    in >> junk1 >> lattice_size;
    in.close();
    xhi = 0.5*box_len;
    yhi = 0.5*box_len;
    zhi = 100000;
    xlo = -xhi;
    xhalf = xhi;
    yhalf = yhi;
    xbox = box_len;
    ybox = box_len;
    ylo = -yhi;
    zlo = 0;

    int ngrafter = box_len/lattice_size;
    writeGrowth = (ngrafter*ngrafter)*0.5;

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

//    int maxMono = 50000;
    int maxMono = 20000;

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
    for (int i = 0; i < maxMono; i++) {
	tempPoly.chain(i) = dummy;
    }

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

double calcAddOneEn(Mono &trialMono, int whichPoly, genarray< Poly > &brush) {
// Calculate pair energies and angle of adding trialMono to polymer whichPoly. 
// New - use step potential - e = epsilon at overlap.
//    cout << "\t in calcAddOneEn! " << endl;
//    double asqinv = 1/(elp_a*elp_a);
//    double csqinv = 1/(elp_c*elp_c);

//    double r_cutsq = r_cut*r_cut;
//    double sigma_sq = sigma*sigma;
    double e = 0;
    int startMono;
    int ic, jc, nMono1, nMono2;
    double dx, dy, dz, drsq, dr;
    double e_one;
//    Mono m1, m2, m3;

//    m2 = trialMono;
    // temp
//    double elp_rad = (m2.x*m2.x + m2.y*m2.y)*asqinv + m2.z*m2.z*csqinv;
//    if (elp_rad < 1) {
//	e = 10000;
//    }
    if (trialMono.z < zlo) {
	e = 10000000;
//	cout << "m2.z less than zlo" << endl;
    }
    
    else {
//	Mono tempm = brush(whichPoly).chain( brush(whichPoly).nMono-1 );
    //    cout << "\t\t\t last Mono: " << tempm.x << " " << tempm.y << " " << tempm.z << endl;
    //    cout << "\t\t\t trialMono: " << m2.x << " " << m2.y << " " << m2.z << endl;
	for (int i = 0; i < brush.length(); i++) {
	    nMono1 = brush(i).nMono;
	    for (int j = 0; j < nMono1; j++) {
//		m1 = brush(i).chain(j);
		dx = brush(i).chain(j).x - trialMono.x;
		dy = brush(i).chain(j).y - trialMono.y;
		dz = brush(i).chain(j).z - trialMono.z;
		if (dx > xhalf) dx -= xbox;
		if (dy > yhalf) dy -= ybox;
		if (dx < -xhalf) dx += xbox;
		if (dy < -yhalf) dy += ybox;
		drsq = dx*dx+dy*dy+dz*dz;

		// Pair energy: 
//		dx = m1.x - m2.x;
//		dy = m1.y - m2.y;
//		dz = m1.z - m2.z;
//		drsq = dx*dx+dy*dy+dz*dz;
//		drsq = monoDistSq(m1, m2);
		if (drsq < sigma_sq) {
//		    dr = sqrt(drsq);
//		    e_one = (1-dr/sigma);
//		    e_one = e_one*e_one;
//		    if (isnan(e_one)) {
//			cout << "BAD PAIR! " << m1.x << " " << m1.y << " " << m1.z << endl;
//			cout << "dr: " << dr << endl;
//			cout << "i, j: " << i << " " << j << endl;
//			cout << "nMono1: " << nMono1 << endl;
//		    }
		    e += epsilon;
		    // What is temperature in a hard sphere system with this Rosenbluth stuff?
		}
	    }
	}
    }

//	// Angles:
//	int len = brush(whichPoly).nMono;
//	double e_angle = 0;
//	if (len > 1) {
//	    double dx1, dy1, dz1, dx2, dy2, dz2, rsq1, rsq2, r1, r2;
//	    double c, s, dtheta;
//	    int jc1, jc2, jc3;
//
//	    m1 = brush(whichPoly).chain(len-2);
//	    m2 = brush(whichPoly).chain(len-1);
//	    m3 = trialMono;
//
//	    dx1 = m1.x - m2.x;
//	    dy1 = m1.y - m2.y;
//	    dz1 = m1.z - m2.z;
//	    rsq1 = dx1*dx1 + dy1*dy1 + dz1*dz1;
//	    r1 = sqrt(rsq1);
//
//	    dx2 = m3.x - m2.x;
//	    dy2 = m3.y - m2.y;
//	    dz2 = m3.z - m2.z;
//	    rsq2 = dx2*dx2 + dy2*dy2 + dz2*dz2;
//	    r2 = sqrt(rsq2);
//
//	    c = dx1*dx2 + dy1*dy2 + dz1*dz2;
//	    c /= r1*r2;
//	    if (c > 1.0) c = 1.0;
//	    if (c < -1.0) c = -1.0;
//	    dtheta = acos(c) - theta_0;
//	    e_angle = k_angle * dtheta*dtheta;
//	    if (e_angle > 40) {
//		cout << "WEIRD ANGLE!!!!" << endl;
//		cout << "len: " << len << endl;
//		cout << "r1, r2, c: " << r1 << " " << r2 << " " << c << endl;
//		cout << "\tacos(c): " << acos(c) << endl;
//		cout << "\tdtheta: " << dtheta << endl;
//		cout << "\te_angle: " << e_angle << endl;
//		cout << "\t mono1: " << m1.x <<  " " << m1.y << " " << m1.z << endl;
//		cout << "\t mono2: " << m2.x <<  " " << m2.y << " " << m2.z << endl;
//		cout << "\t mono3: " << m3.x <<  " " << m3.y << " " << m3.z << endl;
//		cout << "d1s: " << dx1 << " " << dy1 << " " << dz1 << endl;
//		cout << "d2s: " << dx2 << " " << dy2 << " " << dz2 << endl;
//
//	    }
//	    e += e_angle;
//	}
//    }

    return e;
}

void genTrialPts( int whichPoly, genarray< Poly > &brush, genarray < Poly > &trialMonos, int previousPoly, 
	genarray< Poly > &trialThomson) { 

// Generates trial additions to polymer whichPoly
//    cout << "\t in genTrialPts \n";

    double r_cutsq = r_cut*r_cut;
    double dx, dy, dz, mag;
    Mono m2, mt;

//    cout << "whichPoly: " << whichPoly << endl;
    int nMono = brush(whichPoly).nMono;
//    cout << "nMono-1: " << (nMono-1) << endl;

    m2 = brush(whichPoly).chain(nMono-1);	// trialPoly's tail monomer.

    // Check if distance between the trialPoly's tail and previousPoly tail monomers is close enough to warrant a new trial addition:
    double drsq = 0;
    if (previousPoly > -1) {
	Mono mprev = brush(previousPoly).chain(brush(previousPoly).nMono-1);
	drsq = monoDistSq(m2, mprev);
    }
    else {
	drsq = 0;
    }
    double xo, yo, zo;
    if ( drsq < r_cutsq ) {
	// Use Thomson pts to generate trial pts around the tail monomer:
	int rand = gsl_rng_uniform_int (mrRand, trialThomson.length());
	Poly trialSphere = trialThomson(rand);
	Mono ms;
	for (int i = 0; i < trialSphere.nMono; i++) {
	    ms = trialSphere.chain(i);
	    ms.x = ms.x + m2.x;
	    ms.y = ms.y + m2.y;
	    ms.z = ms.z + m2.z;
	    if(ms.x < xlo) ms.x = ms.x + xbox;
	    if(ms.y < ylo) ms.y = ms.y + ybox;
	    if(ms.x > xhi) ms.x = ms.x - xbox;
	    if(ms.y > yhi) ms.y = ms.y - ybox;
	    trialMonos(whichPoly).chain(i) = ms;
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
    double r_cutsq = r_cut*r_cut;
    double beta = -1/TEMP;
    Mono store;		
    store.z = 0;
    double temp = 0;
//    double sigma_sq = sigma*sigma;
    //
    Mono dbg;
    //
    for (int i = 0; i < brush.length(); i++) {
	if (previousPoly > -1) {
	    Mono mprev = brush(previousPoly).chain(brush(previousPoly).nMono-1);
	    Mono mtail = brush(i).chain(brush(i).nMono-1);
//	    dx = mtail.x - mprev.x;
//	    dy = mtail.y - mprev.y;
//	    dz = mtail.z - mprev.z;
//	    drsq = dx*dx+dy*dy+dz*dz;
	    drsq = monoDistSq(mtail, mprev);
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
		if (en > 1000)
		    rose = 0;
		else
		    rose = exp(temp);
		    
//		rose = exp(temp);
//		cout << "\t\t  en, rose: " << en << " " <<  rose << " " << endl;
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

//    cout << "w_move is: " << w_move << endl << endl;
    
    int is = -1;
    int js = -1; 
    // SELECT (Frenkel and Smit Algorithm 41):
    double ws = gsl_rng_uniform(mrRand)*w_move;
    double cumw = trialRose(0).chain(0).y;
    int start;
//    cout << "ws, cumw: " << ws << " " << cumw << endl;
    for (int i = 0; i < trialRose.length(); i++) {
	if (i == 0) 
	    start = 1;
	else
	    start = 0;
	for (int j = start; j < nTrial; j++) {
//	    cout << "\ti, j, cumw is: " << i << " " << j << " " << cumw << endl;
	    cumw = cumw + trialRose(i).chain(j).y;
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
	}
    }
    if (is < 0 && js < 0) {
	cout << "stuck!  can't add monomer!" << endl;
	cout << "w_move is: " << w_move << endl;
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
    if (myMono.z < zlo || myMono.x > xhi || myMono.x < xlo || myMono.y > yhi || myMono.y < ylo) {
    cout << "bad mono! is, js: " << is << " " << js << endl;
    cout << "trialRose en, rose: " << trialRose(is).chain(js).x << " " << trialRose(is).chain(js).y << endl;
    cout << "\ttail mono: " << prev.x << " " << prev.y << " " << prev.z << ", " << mag1 << endl;
    cout << "\tadding mono: " << myMono.x << " " << myMono.y << " " << myMono.z << ", " << mag2 << endl;
    }
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
	myOut << "O 0 0 -2 \n";
	myOut << "N 0 0 -2 \n";
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

void readThomson (string inName, genarray< Poly > &trialThomson) {
// Read in a XYZ file with points on a sphere (thomson problem).
// At the end, trialThomson has 10 different (rotated) configs of points on sphere.
    int ctr = 0;
    int atomctr = 0;
    double x,y,z;
    double theta;
    int type;
    cout << "reading file " << inName << " for Thomson points." << endl;
    ifstream in(inName.c_str());
    string junk, myLine;
    Poly mySphere;
    int nPts;

    trialThomson.resize(10);

    while (getline(in, myLine)) {
	if (myLine != "") {
	    stringstream ss(myLine);
	    if (ctr == 0) {
		ss >> nPts;
		mySphere.chain.resize(nPts);
		mySphere.nMono = nPts;
	    }
	    if (ctr > 1) {
		ss >> junk >> x >> y >> z;
		mySphere.chain(ctr-2).x = x;
		mySphere.chain(ctr-2).y = y;
		mySphere.chain(ctr-2).z = z;
	    }
	}
	ctr++;
    }
    trialThomson(0) = mySphere;
    
// Fill in trialThomson with rotations of the sphere (write quaternion crap laterz)
    double ax, ay, az, at;
    double latThomson = sqrt(nPts/(4*PI));
    double angleThomson = 2*PI/latThomson;
    double at3 = angleThomson/3.0;
    at = at3;
    Mono myM;
    Poly newSphere = mySphere;
    for (int i = 1; i < 10; i++) {
//	trialThomson(i) = mySphere;
	if (i < 4) {ax = 1; ay = 0; az = 0;}
	if (i >= 4 && i < 7) {ax = 0; ay = 1; az = 0;}
	if (i >= 7 && i < 10) {ax = 0; ay = 0; az = 1;}
	for (int j = 0; j < mySphere.nMono; j++) {
	    x = mySphere.chain(j).x;
	    y = mySphere.chain(j).y;
	    z = mySphere.chain(j).z;
	    rotate_pt(x, y, z, at, ax, ay, az);
	    myM.x = x;
	    myM.y = y;
	    myM.z = z;
	    newSphere.chain(j) = myM;
	}
	at += at3; 
	if (at > 3*at3)
	    at = at3;
	trialThomson(i) = newSphere;
    }
}

void makeSquareLattice (genarray <double> &atomPositions) {

    int ngraft_sqrt = double(box_len)/lattice_size;
    ngraft = ngraft_sqrt*ngraft_sqrt;
    atomPositions.resize(ngraft*3);
    int ic;
    for (int i = 0; i < ngraft_sqrt; i++) {
	for (int j = 0; j < ngraft_sqrt; j++) {
	    ic = (i*ngraft_sqrt + j)*3;
	    atomPositions(ic) = xlo + lattice_size*i;
	    atomPositions(ic+1) = ylo + lattice_size*j;
	    atomPositions(ic+2) = zlo;
	}
    }
}

double monoDistSq (Mono m1, Mono m2) {
    double dx = m1.x - m2.x;
    double dy = m1.y - m2.y;
    double dz = m1.z - m2.z;
    if (dx > xhalf)
	dx -= xbox;
    if (dx < -xhalf)
	dx += xbox;
    if (dy > yhalf)
	dy -= ybox;
    if (dy < -yhalf)
	dy += ybox;

    return (dx*dx + dy*dy + dz*dz);
}

void rotate_pt (double &x1, double &y1, double &z1, double theta, double ax, double ay, double az) {
// Rotates point (x,y,z) about axis (ax, ay, az) by angle 'theta'
    double cc = cos(theta*0.5);
    double ss = sin(theta*0.5);
    double bq = ss*ax;
    double cq = ss*ay;
    double dq = ss*az;
    double bqi = -bq;
    double cqi = -cq;
    double dqi = -dq;
    quat_mult(cc, bq, cq, dq, 0, x1, y1, z1);
    quat_mult(cc, bq, cq, dq, cos(theta*0.5), bqi, cqi, dqi);
    x1 = bq;
    y1 = cq;
    z1 = dq;

}

void quat_mult (double &a1, double &b1, double &c1, double &d1, double a2, double b2, double c2, double d2) {
    double aa = a1*a2 - b1*b2 - c1*c2 - d1*d2;
    double bb = a1*b2 + b1*a2 + c1*d2 - d1*c2;
    double cc = a1*c2 - b1*d2 + c1*a2 + d1*b2;
    double dd = a1*d2 + b1*c2 - c1*b2 + d1*a2;
    a1 = aa;
    b1 = bb;
    c1 = cc;
    d1 = dd;
}

//void genTrialPts2( int whichPoly, genarray< Poly > &brush, genarray < Poly > &trialMonos, int previousPoly, 
//	genarray< Poly > &trialThomson, genarray< Poly > &trialRose) { 
//
void genTrialPts2dbg( int whichPoly, genarray< Poly > &brush, genarray < Poly > &trialMonos, int previousPoly, 
	genarray< Poly > &trialThomson, genarray< Poly > &trialRose, int &dbgctr) { 
// Combines gentrialpts and calcrosenbluth:

    double beta = -1/TEMP;
    double r_cutsq = r_cut*r_cut;
    double dx, dy, dz, mag;
    Mono m2, mt, store;
    double en, rose, myexp;

    store.x = 0;
    store.y = 0;
    int nMono = brush(whichPoly).nMono;
//    m2 = brush(whichPoly).chain(nMono-1);	// trialPoly's tail monomer.

    // Check if distance between the trialPoly's tail and previousPoly tail monomers is close enough to warrant a new trial addition:
    double drsq = 0;
    if (previousPoly > -1) {
//	Mono mprev = brush(previousPoly).chain(brush(previousPoly).nMono-1);
//	drsq = monoDistSq(m2, mprev);

	dx = brush(whichPoly).chain(nMono-1).x - brush(previousPoly).chain(brush(previousPoly).nMono-1).x;
	dy = brush(whichPoly).chain(nMono-1).y - brush(previousPoly).chain(brush(previousPoly).nMono-1).y;
	dz = brush(whichPoly).chain(nMono-1).z - brush(previousPoly).chain(brush(previousPoly).nMono-1).z;
	if (dx > xhalf)
	    dx -= xbox;
	if (dx < -xhalf)
	    dx += xbox;
	if (dy > yhalf)
	    dy -= ybox;
	if (dy < -yhalf)
	    dy += ybox;
	drsq = dx*dx + dy*dy + dz*dz;
    }
    else {
	drsq = 0;
    }
    double xo, yo, zo;
    int trialnMono;
    if ( drsq < r_cutsq ) {
	dbgctr++;
//	cout << "\t\t" << whichPoly << endl;
	// Use Thomson pts to generate trial pts around the tail monomer:
	int rand = gsl_rng_uniform_int (mrRand, trialThomson.length());
	Poly trialSphere = trialThomson(rand);
	Mono ms;
	trialnMono = trialThomson(rand).nMono;
	for (int i = 0; i < trialSphere.nMono; i++) {
//	for (int i = 0; i < trialnMono; i++) {
	    // start old stuff
	    ms = trialSphere.chain(i);
//	    ms.x = ms.x + m2.x;
//	    ms.y = ms.y + m2.y;
//	    ms.z = ms.z + m2.z;
	    ms.x = ms.x + brush(whichPoly).chain(nMono-1).x;
	    ms.y = ms.y + brush(whichPoly).chain(nMono-1).y;
	    ms.z = ms.z + brush(whichPoly).chain(nMono-1).z;

	    if(ms.x < xlo) ms.x = ms.x + xbox;
	    if(ms.y < ylo) ms.y = ms.y + ybox;
	    if(ms.x > xhi) ms.x = ms.x - xbox;
	    if(ms.y > yhi) ms.y = ms.y - ybox;
	    trialMonos(whichPoly).chain(i) = ms;
	    // end old stuff


	    en = calcAddOneEn(ms, i, brush);
//	    en = calcAddOneEn(trialMonos(whichPoly).chain(i), i, brush);
//	    if (isnan(en)) {
//		cout << " BAD ENERGY! trialmono: " << dbg.x << " " << dbg.y << " " << dbg.z << endl;
//	    }
	    myexp = beta*en;
	    if (en > 1000)
		rose = 0;
	    else
		rose = exp(myexp);
		
	    store.x = en;
	    store.y = rose;
	    trialRose(whichPoly).chain(i) = store;
	}
    }


}

void writeLengths (string outName, genarray< Poly > &brush, int growthStep) {

    int ctr = 0;

    ofstream myOut;
    myOut.open(outName.c_str(), ios::app);

    for (int i = 0; i < brush.length(); i++) {
	myOut << growthStep << "\t" << brush(i).nMono << "\n";
    }
    myOut.close();


}
