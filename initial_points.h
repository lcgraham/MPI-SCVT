/*
 * initial_routines.h
 *
 *  Created on: Jan 13, 2017
 *      Author: huanhuan
 */

#ifndef INITIALIZE_ROUTINES_H_
#define INITIALIZE_ROUTINES_H_

#include <iostream>
#include <fstream>
#include "triangulation.h"

#define SEED	3729

using namespace std;
double density(const pnt &p);
double pop_lowres_density(const pnt &p);
double pop_highres_density(const pnt &p);
double ellipse_density(const pnt &p, double lat_c, double lon_c, double lat_width, double lon_width);


/* ***** Point Init Routines ***** {{{*/
void readPoints(int& n, vector<pnt>& points){/*{{{*/
	//Read in initial point set from SaveVertices
	ifstream points_in("SaveVertices");
	pnt p;
	int i;

	i = 0;
	while(!points_in.eof()){
		points_in >> p;
		points_in.ignore(10000,'\n');
		p.idx = i;
		p.isBdry = 0;
		p.normalize();

		if(points_in.good()){
			points.push_back(p);
		}
		i++;
	}

	n = points.size();

	cout << "Read in " << n << " points from SaveVertices" << endl;

}/*}}}*/


void makeMCPoints(int n, vector<pnt>& points){
	//Create Monte Carlo random point set
	int i, j;
	srand48(time(NULL));
	double dlon, dlat;
	double dtr;
	double lat, lon;
	double x, y, z;
	double dens_check, dens_comp;
	pnt p;

	dtr = M_PI/180.0;

	// Uniform Spherical Points
	for(i = 0; i < n; i++){
		x = drand48()*2.0-1.0;
		y = drand48()*2.0-1.0;
		z = drand48()*2.0-1.0;

		p = pnt(x,y,z);
		p.idx = i;
		p.isBdry = 0;
		p.normalize();


		points.push_back(p);
	}
	cout << "Created " << points.size() << " points using monte carlo." << endl;
}

void makeMCPoints_rejection(int n, vector<pnt>& points){
	//Create Monte Carlo random point set
	int i, j;
	srand48(time(NULL));
	double dlon, dlat;
	double dtr;
	double lat, lon;
	double x, y, z;
	double dens_check, dens_comp;
	pnt p;

	dtr = M_PI/180.0;

	// Non-Uniform Spherical Points w/ rejection sampling
	for(i = 0; i < n; i++){
		x = drand48()*2.0-1.0;
		y = drand48()*2.0-1.0;
		z = drand48()*2.0-1.0;

		p = pnt(x,y,z);
		p.idx = i;
		p.isBdry = 0;
		p.normalize();

		dens_check = density(p);

		if drand48() <= dens_check:
			points.push_back(p);
	}
	cout << "Created " << points.size() << " points using monte carlo." << endl;
}


void makeGeneralizedSpiralPoints(int n, vector<pnt>& points){/*{{{*/
	//Create Generalize Spiral point set
	int i;
	int idx;
	pnt p;
	double phi_curr, h, theta, aa, bb, cc;
	double gsC = 3.809;
	double twopi_dp;

	twopi_dp = 2.0*M_PI;

	p = pnt(0.0,0.0,0.0,0,0);

	//	first pt, loop primer
	i = 0;
	h = -1.0;
	theta = acos(h);
	phi_curr = 0;

	p.x = cos( phi_curr ) * sin( theta );
	p.y = sin( phi_curr ) * sin( theta );
	p.z = cos( theta );
	p.idx = i;
	p.isBdry = 0;
	points.push_back(p);

	for(i = 1; i < n-1; i++){
		h = -1.0 + (2.0*(double)(i))/(double)(n-1);
		theta = acos(h);
		aa = phi_curr;
		bb = gsC / sqrt((double)n);
		cc = 1.0 / sqrt(1.0-h*h);
		phi_curr =  fmod(aa + bb * cc,twopi_dp);

		p.x = cos( phi_curr ) * sin( theta );
		p.y = sin( phi_curr ) * sin( theta );
		p.z = cos( theta );
		p.idx = i;
		p.isBdry = 0;
		points.push_back(p);
	}

	p.x = 0.0;
	p.y = 0.0;
	p.z = 1.0;
	p.idx = n-1;
	p.isBdry = 0;
	points.push_back(p);

	cout << "Created " << points.size() << " points using generalized spiral." << endl;
}/*}}}*/

void makeGeneralizedSpiralPoints_rejection(int n, vector<pnt>& points){/*{{{*/
	//Create Generalize Spiral point set
	int i;
	int idx;
	pnt p;
	double phi_curr, h, theta, aa, bb, cc;
	double gsC = 3.809;
	double twopi_dp;
	double dens_check;

	twopi_dp = 2.0*M_PI;

	p = pnt(0.0,0.0,0.0,0,0);

	//	first pt, loop primer
	i = 0;
	h = -1.0;
	theta = acos(h);
	phi_curr = 0;

	p.x = cos( phi_curr ) * sin( theta );
	p.y = sin( phi_curr ) * sin( theta );
	p.z = cos( theta );
	p.idx = i;
	p.isBdry = 0;
	points.push_back(p);

	for(i = 1; i < n-1; i++){
		h = -1.0 + (2.0*(double)(i))/(double)(n-1);
		theta = acos(h);
		aa = phi_curr;
		bb = gsC / sqrt((double)n);
		cc = 1.0 / sqrt(1.0-h*h);
		phi_curr =  fmod(aa + bb * cc,twopi_dp);

		p.x = cos( phi_curr ) * sin( theta );
		p.y = sin( phi_curr ) * sin( theta );
		p.z = cos( theta );
		p.idx = i;
		p.isBdry = 0;

		dens_check = density(p);
		if drand48() <= dens_check:
			points.push_back(p);
	}

	p.x = 0.0;
	p.y = 0.0;
	p.z = 1.0;
	p.idx = n-1;
	p.isBdry = 0;
	points.push_back(p);

	cout << "Created " << points.size() << " points using generalized spiral." << endl;
}/*}}}*/



void makeFibonacciGridPoints(int n, vector<pnt>& points){/*{{{*/
	const double g_ratio = (1.0+sqrt(5)) / 2.0;
	double lambda, phi, sinphi, x, y, z;
	int i, j, m;
	pnt p;

	m = n/2;

	for (i = 0; i < n; i++){
		j = i - m;

		lambda = i / g_ratio;
		sinphi = (2.0*j) / (n + 1);
		phi = asin(sinphi);
		x = cos(phi) * sin(lambda);
		y = cos(phi) * cos(lambda);
		z = sin(phi);
		p = pntFromLatLon(phi, lambda);
		p.idx = i;
		p.isBdry = 0;
		p.normalize();
		points.push_back(p);
	}
}/*}}}*/
/*}}}*/


void makeFibonacciGridPoints_rejction(int n, vector<pnt>& points){/*{{{*/
	const double g_ratio = (1.0+sqrt(5)) / 2.0;
	double lambda, phi, sinphi, x, y, z;
	int i, j, m;
	pnt p;
	double dens_check;

	m = n/2;

	for (i = 0; i < n; i++){
		j = i - m;

		lambda = i / g_ratio;
		sinphi = (2.0*j) / (n + 1);
		phi = asin(sinphi);
		x = cos(phi) * sin(lambda);
		y = cos(phi) * cos(lambda);
		z = sin(phi);
		p = pntFromLatLon(phi, lambda);
		p.idx = i;
		p.isBdry = 0;
		p.normalize();
		dens_check = density(p);
		if drand48() <= dens_check:
			points.push_back(p);
	}
}/*}}}*/
/*}}}*/



double density(const pnt &p){/*{{{*/
	//density returns the value of the density function at point p
	//return 1.0; // Uniform density

	/* Density function for Shallow Water Test Case 5 
	pnt cent;
	double r;
	double norm;
	double density;
	double min_val;
	double width, trans_cent;

	cent = pnt(0.0, -0.866025403784439, 0.5);
	cent.normalize();

	width = 0.15;
	min_val = 1.0/8.0;
	min_val = pow(min_val,4);
	trans_cent = 30.0*M_PI/180.0;
	norm = 1.0/(1.0-min_val);

	r = cent.dotForAngle(p);

	density = ((tanh((trans_cent-r)*(1.0/width))+1.0)/2.0)/norm + min_val;

	return density;
	// */
	
	/* Ellipse density function.
	
	return ellipse_density(p, 40.0, 0.0, 1.0, 0.5);
	// */
    
    /* Pop low resolution density function.
    return pop_lowres_density(p);
    // */
    
    // /* Pop high resolution density function. 
    return pop_highres_density(p);
    // */

}/*}}}*/
double ellipse_density(const pnt &p, double lat_c, double lon_c, double lat_width, double lon_width){/*{{{*/
	//density returns the value of the density function at point p
	//	return 1.0; // Uniform density
	//	/* Ellipse Density function
	pnt work;
	double r1, r2, r;
	double dtr;
	double width, trans_cent, min_val, norm, density;

	dtr = M_PI/180.0;

	work = pntFromLatLon(p.getLat(), lon_c*dtr);
	r1 = work.dotForAngle(p);

	work = pntFromLatLon(lat_c*dtr, p.getLon());
	r2 = work.dotForAngle(p);

	r1 = r1/lon_width;
	r2 = r2/lat_width;
	r = sqrt (r1*r1 + r2*r2);

	width = 0.15;
	trans_cent = 30.0 * dtr;
	min_val = 1.0/12.0;
	min_val = pow(min_val,4);
	norm = 1.0/(1.0-min_val);
	density = ((tanh((trans_cent-r)*(1.0/width))+1.0)/2)/norm + min_val;

	return density;
	// */

}/*}}}*/
double pop_lowres_density(const pnt &p){/*{{{*/
    double dtr;
    double density, gamma, lat, lat_cent, width;

	dtr = M_PI/180.0;

    lat_cent = 25.0 * dtr;
    width = 6.0 * dtr;
    gamma = powf(1.0 / 1.95, 4.0);
    lat = p.getLat();

    density = ((1.0-gamma) * (1.0 + tanh( (lat_cent - fabs(lat)) / width)) / 2.0) + gamma;

    return density;
}/*}}}*/
double pop_highres_density(const pnt &p){/*{{{*/
    double dtr;
    double density, lat, gamma;

    dtr = M_PI/180.0;

    gamma = pow(1.0 / 3.0, 4);
    lat = p.getLat();

    density = (1-gamma) * pow( sin(lat), 4) + gamma;

    return density;
}/*}}}*/
/*}}}*/

#endif /* INITIALIZE_ROUTINES_H_ */
