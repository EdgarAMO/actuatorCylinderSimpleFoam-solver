/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    actuatorCylinder

Description
    Add volume force in an actuator cylinder.

    Use with actuatorCylinderSimpleFoam

    Modified by Edgar Martínez, September 2020.

    The actuator disk can be defined by adding the following 
    lines in fvSolutions:

    actuatorCylinder
    {
   
    }
\*---------------------------------------------------------------------------*/

#include "actuatorCylinder.H"

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "faceAreaPairGAMGAgglomeration.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// Pi
scalar PI = 3.14159265359;

// Offset factor
scalar F = 1.01;

// Convergence tolerance
scalar ERR = 0.01;

// Relaxation factor
scalar BETA = 0.35;

// Kinematic viscosity
scalar NU = 1.5e-5;

// Slice width in radians
scalar DT = 2.0 * PI / N;

// Sublice width in radians
scalar DP = DT / NS;

// Conversion to degrees factor
scalar TODEG = 180.0 / PI;

// Turbine's NACA digit
label NACA = 0;

// Blade's pitch angle
scalar BPA = 0.;

// Turbine's number of blades
label NB = 0;

// Turbine's chord
scalar CH = 0.;

// Turbines radius
scalar R = 0.;

// Turbine's RPM
scalar RPM = 0.;

// Turbine's tip-speed ratio
scalar LAMBDA = 0.;

// Freestream velocity
scalar U0 = 0.;

namespace Foam {

    defineTypeNameAndDebug(actuatorCylinder, 0);



    // Default constructor
    actuatorCylinder::actuatorCylinder() {
       
        // Initialize cylinder's centerline
        center = vector::zero;

        // Cylinder's internal and external radii
        intRadius = 0.0;
        extRadius = 0.0;

        // Rotation vector
        spin = 0;

    }



    // Destructor
    actuatorCylinder::~actuatorCylinder() {}



    // Convert ID to string
    // *** PUBLIC ***
    void actuatorCylinder::setID(label number){ ID = std::to_string(number); }



    // Read the cylinder's axis points
    // *** PUBLIC ***
    void actuatorCylinder::readGeometry(const fvMesh &mesh) {

        if(debug >= 1) { Info << "Reading actuator cylinder geometry.\n"; }

        // Get each turbine's characteristics
        std::string object;
        object.append("actuatorCylinder");
        object.append(ID);  

        Istream& is1 = mesh.solutionDict().subDict(object).lookup("center");
        is1.format(IOstream::ASCII);
        is1 >> center; 

        Istream& is2 = mesh.solutionDict().subDict(object).lookup("intRadius");
        is2.format(IOstream::ASCII);
        is2 >> intRadius;

        Istream& is3 = mesh.solutionDict().subDict(object).lookup("extRadius");
        is3.format(IOstream::ASCII);
        is3 >> extRadius;

        Istream& is4 = mesh.solutionDict().subDict(object).lookup("spin");
        is4.format(IOstream::ASCII);
        is4 >> spin;

        // Get the general geometrical and operational parameters
        Istream& gis1 = mesh.solutionDict().subDict("general").lookup("NACA");
        gis1.format(IOstream::ASCII);
        gis1 >> NACA;

        Istream& gis2 = mesh.solutionDict().subDict("general").lookup("BPA");
        gis2.format(IOstream::ASCII);
        gis2 >> BPA;

        Istream& gis3 = mesh.solutionDict().subDict("general").lookup("NB");
        gis3.format(IOstream::ASCII);
        gis3 >> NB;

        Istream& gis4 = mesh.solutionDict().subDict("general").lookup("CH");
        gis4.format(IOstream::ASCII);
        gis4 >> CH;

        Istream& gis5 = mesh.solutionDict().subDict("general").lookup("R");
        gis5.format(IOstream::ASCII);
        gis5 >> R;

        Istream& gis6 = mesh.solutionDict().subDict("general").lookup("RPM");
        gis6.format(IOstream::ASCII);
        gis6 >> RPM;

        Istream& gis7 = mesh.solutionDict().subDict("general").lookup("LAMBDA");
        gis7.format(IOstream::ASCII);
        gis7 >> LAMBDA;

        Istream& gis8 = mesh.solutionDict().subDict("general").lookup("U0");
        gis8.format(IOstream::ASCII);
        gis8 >> U0;
 
        Info << "\n";
        Info << "Center .........." << center << "\n";
        Info << "Ri .............." << intRadius << "\n";
        Info << "Ro .............." << extRadius << "\n";
        Info << "Spin ............" << spin << "\n";
        Info << "NACA00 .........." << NACA << "\n";
        Info << "BPA ............." << BPA << "\n";
        Info << "NB .............." << NB << "\n";
        Info << "CH .............." << CH << "\n";
        Info << "R ..............." << R << "\n";
        Info << "RPM ............." << RPM << "\n";
        Info << "LAMBDA .........." << LAMBDA << "\n";
        Info << "U0 .............." << U0 << "\n\n";

        BPA *= PI / 180.0;

        // Find the cells inside the hollow cylinder
        collectCells(mesh);

        Info << "Cylinder cells have been collected.\n";

    }



    // Initialize arrays
    // *** PUBLIC ***
    void actuatorCylinder::initializeArrays() {

    	 // Fill AA[]
        std::ifstream aaFile;
        aaFile.open("nacaaa.csv");

        label n = 0;
        std::string line;

        while(aaFile.good()){
            std::getline(aaFile, line, '\n');
            AA[n] = atof(line.c_str());
            ++n;
        }

        aaFile.close();

        // Fill RN[]
        std::ifstream reFile;
        reFile.open("nacare.csv");

        n = 0;

        while(reFile.good()){
            std::getline(reFile, line, '\n');
            RN[n] = atof(line.c_str());
            ++n;
        }

        reFile.close();

        // Choose file names for coefficients
        std::string clname = " ";
        std::string cdname = " ";

        if (NACA == 12){clname = "naca0012cl.csv"; cdname = "naca0012cd.csv";}
        if (NACA == 15){clname = "naca0015cl.csv"; cdname = "naca0015cd.csv";}
        if (NACA == 18){clname = "naca0018cl.csv"; cdname = "naca0018cd.csv";}
        if (NACA == 21){clname = "naca0021cl.csv"; cdname = "naca0021cd.csv";}
        if (NACA == 20){clname = "du06w200cl.csv"; cdname = "du06w200cd.csv";}

        // Initialize coefficients look-up tables CL[][] & CD[][]
        initializeLookUpTables(CL);
        initializeLookUpTables(CD);

        //Parse coefficients files
        parse(clname, CL);
        parse(cdname, CD);

        // Initialize azimuthal angle array
        for (label j = 0; j < N; j++) { t[j] = (j + 1.0 / 2.0) * DT; }

        Info << "Arrays were initialized.\n";

    }



    // Add the volumetric force in the appropriate cells
    // *** PUBLIC ***
    void actuatorCylinder::addForce(const fvMesh &mesh, vectorField &volumeForce, const vectorField &U) {

        if(debug >= 1) { Info << "Calculating volume force from actuator cylinder.\n"; }

        vector volForce(0, 0, 0);

        vector y_axis(0, 1, 0);

        collectVelocities(mesh, U);
        //Info << "Velocities in cylinder have been collected.\n";

       	computeVolForce();
        //Info << "Actuator Cylinder has been run.\n";

        //volumeForce[i] += volForce;
       	
        // Loop over all cells in the hollow cylinder
        forAll (indices, ix) {

        	label index = indices[ix];

            // current point being inscpected
            vector point = mesh.C()[index];

            // vector from center to current point
            vector centerToPoint(point - center);

            centerToPoint /= mag(centerToPoint);

            // angle between cylinder's y-axis and point vector
            scalar angle = std::acos( (centerToPoint & y_axis) / ( mag(centerToPoint) * mag(y_axis) ) );

            if ( point.x() > center.x() ) {
                angle = 2 * PI - angle;
            }

        	volForce.x() = interpolateForces(angle, volFx);
        	volForce.y() = interpolateForces(angle, volFy);
        	volumeForce[index] += volForce;
            	
        }

    }



    // *** PRIVATE ***
    bool actuatorCylinder::pointIsInside(const vector &point) {

        // start point to test point vector "r = point - center"
        vector centerToPoint(point - center);

	   // squared magnitude of "r"
        scalar radialDist = magSqr(centerToPoint);

        // check magSqr(r) < extRadius^2 && mag(r) > intRadius^2
        scalar r2 = extRadius * extRadius;
        scalar r1 = intRadius * intRadius;
        
        return (radialDist <= r2 && radialDist >= r1);

    }



    // *** PRIVATE ***
    void actuatorCylinder::collectCells(const fvMesh &mesh) {

        for (label i = 0; i < mesh.C().size(); i++) {

            if ( pointIsInside(mesh.C()[i]) ) {

                indices.append(i);

            }
        }

    }



    // *** PRIVATE ***
    void  actuatorCylinder::collectVelocities(const fvMesh &mesh, const vectorField &U) {

        vector y_axis(0, 1, 0);

        // Check for every sector
        for (label k = 0; k < N; k++){

            label count = 0;
            Vx[k] = 0;
            Vy[k] = 0;

            // Loop over all cells and check if the cell center is in the actuator disc region
            forAll (indices, ix) {

            	label index = indices[ix];

                // current point being inscpected
                vector point = mesh.C()[index];

                // vector from center to current point
                vector centerToPoint(point - center);

                centerToPoint /= mag(centerToPoint);

                // angle between cylinder's y-axis and point vector
                scalar angle = std::acos( (centerToPoint & y_axis) / ( mag(centerToPoint) * mag(y_axis) ) );

                if ( point.x() > center.x() ) {
                    angle = 2 * PI - angle;
                }
                
                // check whether the angle is inside sector
                if ( (angle > (t[k] - DT/2)) && (angle < (t[k] + DT/2)) ){

                    Vx[k] += U[index].x();
                    Vy[k] += U[index].y();
                    count++;

                } 

            }

            // Normalized average
            Vx[k] /= (count * U0);
            Vy[k] /= (count * U0);

        }

    }



    // *** PRIVATE ***
    void actuatorCylinder::computeVolForce() {

    	// Solidity
		scalar SIGMA = (NB * CH) / (2.0 * R);

		// Global Reynolds number
		scalar RE = (RPM * R * CH * PI) / (30 * NU);

        // Control point
        label j;

        // Evaluation point
        label i;

        // Coordinates:

        // Slightly offsetted 'x' coordinate
        scalar x[N];

        // Slightly offsetted 'y' coordinate
        scalar y[N];

        for (j = 0; j < N; j++){

            x[j] = -F * std::sin(t[j]);
            y[j] = +F * std::cos(t[j]);

        }

        // Influence coefficients:

        // Influence coefficient for Wx
        scalar I1[N][N];

        // Influence coefficient for Wy
        scalar I2[N][N];

        // Dummy angle
        scalar p;

        // Auxiliary numerator
        scalar nom;

        // Auxiliary denominator
        scalar den;

        for (j = 0; j < N; j++){

            for (i = 0; i < N; i++){

                I1[j][i] = 0;
                I2[j][i] = 0;

                for (label k = 0; k <= NS; k++){

                    p = t[i] - (1.0 / 2.0) * DT + k * DP;
                    nom = -(x[j] + std::sin(p))*std::sin(p) + (y[j] - std::cos(p))*std::cos(p);
                    den = std::pow(x[j] + std::sin(p), 2) + std::pow(y[j] - std::cos(p), 2);
                    I1[j][i] += nom / den;
                    nom = -(x[j] + std::sin(p))*std::cos(p) - (y[j] - std::cos(p))*std::sin(p);
                    I2[j][i] += nom / den;

                }

                I1[j][i] = I1[j][i] * DP / (-2.0 * PI);
                I2[j][i] = I2[j][i] * DP / (-2.0 * PI);

            }

        }
        
        // Normal loads wake matrix:
        scalar WXn[N][N];
        for (j = 0; j < N; j++){

            for (i = 0; i < N; i++){

                WXn[j][i] = 0.0;

            }

        }

        for (j = N/2; j < N; j++){

            WXn[j][N-j-1] = -1.0;
            WXn[j][j] = 1.0;

        }

        // Tangential loads wake matrix:
        scalar WXt[N][N];
        for (j = 0; j < N; j++){

            for (i = 0; i < N; i++){

                WXt[j][i] = 0.0;

            }

        }

        for (j = N/2 + 1; j < N - 1; j++){

            // Excluded points 185° & 365°
            WXt[j][N-j-1] = -y[j] / std::sqrt(1 - std::pow(y[j], 2));
            // Avoid division by zero
            WXt[j][j] = -y[j] / std::sqrt(1 - std::pow(y[j], 2));

        }

        // Old values of induced velocities:

        // 'x' component
        scalar wx[N];
        // 'y' component
        scalar wy[N];

        for (i = 0; i < N; i++){ wx[i] = 0.0; wy[i] = 0.0; }

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

        // Aerodynamics:
        scalar wx_new[N];   // New induced velocity in 'x'
        scalar wy_new[N];   // New induced velocity in 'y'
        scalar vx[N];       // Non-dimensional velocity in the 'x' direction
        scalar vy[N];       // Non-dimensional velocity in the 'y' direction
        scalar vn;          // Non-dimensional normal velocity
        scalar vt;          // Non-dimensional tangential velocity
        scalar vr[N];       // Non-dimensional relative velocity
        scalar aoa[N];      // Angle of attack
        scalar cl[N];       // Lift coefficient
        scalar cd[N];       // Drag coefficient
        scalar cn[N];       // Normal force coefficient
        scalar ct[N];       // Tangential force coefficient
        scalar Qn[N];       // Non-dimensional normal load
        scalar Qt[N];       // Non-dimensional tangential load
        scalar mod_lin;     // Correction factor for the induced velocities
        scalar re;          // Local Reynolds number
        scalar indFac;      // Induction factor
        scalar thrCoe;      // Thrust coefficient
        bool   fail = true; // Convergence flag
        label    loops = 0;   // Number of iterations
        label    sign;        // Auxiliary sign to modify equations according to spin

        if (spin == 1) {sign = -1;}
        if (spin == 0) {sign =  1;}
        
        // Convergence loop:
        while(fail){

            for (i = 0; i < N; i++){

                // Inherit dimensionless velocities from the cells
                vx[i] = Vx[i];
                vy[i] = Vy[i];

                // Update the deficit velocities
                wx[i] = vx[i] - 1.0;
                wy[i] = vy[i];

                vn = vx[i] * std::sin(t[i]) - vy[i] * std::cos(t[i]);
                vt = -sign * vx[i] * std::cos(t[i]) - sign * vy[i] * std::sin(t[i]) + LAMBDA;
                vr[i] = std::sqrt(std::pow(vn, 2) + std::pow(vt, 2));
                aoa[i] = std::atan(vn / vt);
                re = (RE * vr[i]) / (LAMBDA * 1e6);
                cl[i] = interp2D(aoa[i], re, CL);
                cd[i] = interp2D(aoa[i], re, CD);
                cn[i] = cl[i] * std::cos(aoa[i]) + cd[i] * std::sin(aoa[i]);
                ct[i] = cl[i] * std::sin(aoa[i]) - cd[i] * std::cos(aoa[i]);
                Qn[i] =  (SIGMA / (2.0 * PI)) * std::pow(vr[i], 2) * ( cn[i] * std::cos(BPA) - ct[i] * std::sin(BPA) );
                Qt[i] = sign * (SIGMA / (2.0 * PI)) * std::pow(vr[i], 2) * ( cn[i] * std::sin(BPA) + ct[i] * std::cos(BPA) );

                // Volumetric forces
                scalar del = (extRadius - intRadius);

                volFx[i] = (-1) * Qn[i] * std::sin(t[i]) * (U0 * U0) / del;
                volFy[i] = ( 1) * Qn[i] * std::cos(t[i]) * (U0 * U0) / del;

            }

            mod_lin = correctVelocities(Qn, Qt, indFac, thrCoe);

            // New induced velocities:
            for (j = 0; j < N; j++){

                // Sum reset for new induced velocity 'x'
                wx_new[j] = 0.0;
                // Sum reset for new induced velocity 'y'
                wy_new[j] = 0.0;

                for (i = 0; i < N; i++){

                    wx_new[j] += (I1[j][i] + WXn[j][i])*Qn[i] + (I2[j][i] + WXt[j][i])*Qt[i];
                    wy_new[j] += I2[j][i]*Qn[i] - I1[j][i]*Qt[i];

                }

                wx_new[j] = wx_new[j] * mod_lin;
                wy_new[j] = wy_new[j] * mod_lin;

            }

            // Avoid instability by swapping signs
            if (wx_new[N / 4] > 0.0){

                for (label k = 0; k < N; k++){

                    wx_new[k] = wx_new[k] * (-1);
                    wy_new[k] = wy_new[k] * (-1);

                }

            }

            fail = checkConvergence(wx, wx_new, wy, wy_new);

            if (fail){

                for (j = 0; j < N; j++){

                    wx[j] = BETA * wx_new[j] + (1.0 - BETA) * wx[j];
                    wy[j] = BETA * wy_new[j] + (1.0 - BETA) * wy[j];

                }

            }

            ++loops;

            if (loops >= 100){

                break;
            }

        }

        Info << "Power coefficient[" << ID << "] ... " << getCP(Qt) << "\n\n";
		
		Cp = getCP(Qt);

    }



    // *** PRIVATE ***
    scalar actuatorCylinder::correctVelocities(scalar Qn[N], scalar Qt[N], scalar &indFac, scalar &thrCoe) {

        scalar factor;
        scalar function[N];

        const scalar K3 = 0.0892;
        const scalar K2 = 0.0544;
        const scalar K1 = 0.2511;
        const scalar K0 = -0.0017;

        //CALL INTEGRATION FUNCTION:
        for (label k = 0; k < N; k++){

            function[k] = Qn[k] * std::sin(t[k]) + Qt[k] * std::cos(t[k]);

        }

        thrCoe = simpsons(function);

        indFac = K3 * std::pow(thrCoe, 3) + K2 * std::pow(thrCoe, 2) + K1 * std::pow(thrCoe, 1) + K0;

        //HEAVY LOADED ROTOR CORRECTION:
        if (indFac <= 0.15){

            factor = 1.0 / (1.0 - indFac);

        }

        else factor = (1.0 / (1.0 - indFac)) * (0.65 + 0.35 * std::exp(-4.5 * (indFac - 0.15)));

        return factor;
    }

    

    // *** PRIVATE ***
    bool actuatorCylinder::checkConvergence(scalar oldx[N], scalar newx[N], scalar oldy[N], scalar newy[N]) {

        bool loop_again = false;
        scalar condition1;
        scalar condition2;

        for (label j = 0; j < N; j++){

            condition1 = std::fabs((oldx[j] - newx[j]) / newx[j]);
            condition2 = std::fabs((oldy[j] - newy[j]) / newy[j]);

            if (condition1 > ERR || condition2 > ERR){

                loop_again = true;
                break;

            }
        }

        return loop_again;

    }

    

    // *** PRIVATE ***
    scalar actuatorCylinder::getCP(scalar Qt[N]) {

    	label sign = 0;
        if (spin == 1) {sign = -1;}
        if (spin == 0) {sign =  1;}

        scalar cq;

        cq = simpsons(Qt);

        return sign * LAMBDA * cq;

    }

    
    // *** PRIVATE ***
    scalar actuatorCylinder::simpsons(scalar f[N]) {

        const label out_mask[11] = {1, 4, 2, 4, 2, 4, 2, 4, 2, 4, 1};
        const label in_mask[7] = {1, 4, 2, 4, 2, 4, 1};

        scalar sum1 = 0;
        scalar sum2 = 0;
        scalar sum3 = 0;
        scalar sum4 = 0;

        label k;

        for (k = 0; k < 11; k++){
            sum1 += out_mask[k] * f[k];
        }
        for (k = 11; k < 18; k++){
            sum2 += in_mask[k-11] * f[k];
        }
        for (k = 18; k < 25; k++){
            sum3 += in_mask[k-18] * f[k];
        }
        for (k = 25; k < N; k++){
            sum4 += out_mask[k-25] * f[k];
        }
        return (DT / 3.0) * (sum1 + sum2 + sum3 + sum4);

    }

    

    // *** PRIVATE ***
    scalar actuatorCylinder::interp2D(scalar aa, scalar re, scalar C[ROWS][COLS]) {

        // Convert radians to degrees
        aa *= TODEG;

        label i0 = 0;
        label i1 = 0;
        label j0 = 0;
        label j1 = 0;

        // AA and RE locii
        scalar aLoc = 0.;
        scalar rLoc = 0.;

        // Coefficient at aa
        scalar cRej0 = 0.;
        scalar cRej1 = 0.;

        // Find i0 and i1
        for (label k = 0; k < ROWS; k++){

            if (aa >= AA[k] && aa<= AA[k + 1]) {
                i0 = k;
                i1 = k + 1;
                aLoc = (aa - AA[k]) / (AA[k + 1] - AA[k]);
                break;
            }
        }

        // Find j0 and j1
        for (label k = 0; k < COLS; k++){

            if (re >= RN[k] && re<= RN[k + 1]) {
                j0 = k;
                j1 = k + 1;
                rLoc = (re - RN[k]) / (RN[k + 1] - RN[k]);
                break;
            }
        }

        // Coefficients at aa
        cRej0 = C[i0][j0] + aLoc * (C[i1][j0] - C[i0][j0]);
        cRej1 = C[i0][j1] + aLoc * (C[i1][j1] - C[i0][j1]);

        // Coefficient at re
        return cRej0 + rLoc * (cRej1 - cRej0);

    }



    // *** PRIVATE ***
    void actuatorCylinder::parse(std::string cname, scalar C[ROWS][COLS]) {
    	
        // cname: ...... coefficients file name
        // C[][]: ...... coefficients arrays

        std::ifstream cFile;
        cFile.open(cname.c_str());

        std::string line, field, element;

        std::vector< std::vector<std::string> > array;
        std::vector< std::string > v;

        while ( std::getline(cFile, line) )

        {
            v.clear();
            std::stringstream ss(line);

            while ( std::getline(ss, field, ',') )
            {
                v.push_back(field);
            }

            array.push_back(v);
        }

        for (std::size_t i = 0; i < array.size(); ++i)
        {
            for (std::size_t j = 0; j < array[i].size(); ++j)
            {
                element = array[i][j];
                C[i][j] = atof(element.c_str());
            }
        }

        cFile.close();

    }



    // *** PRIVATE ***
    scalar actuatorCylinder::interpolateForces(scalar angle, scalar volF[N]) {

    	label i;
    	scalar loc = 0.;

    	if (angle < t[0] || angle > t[N - 1]) {

    		loc = (angle - t[N - 1]) / (t[0] - t[N - 1]);
    		return volF[N - 1] + loc * (volF[0] - volF[N - 1]);

    	}

    	for (i = 0; i < N; i++) {

    		if (angle >= t[i] && angle <= t[i + 1]) {

    			loc = (angle - t[i]) / (t[i + 1] - t[i]);
    			break;
    			
    		}

    	}

    	return volF[i] + loc * (volF[i + 1] - volF[i]);

    }



    // *** PRIVATE ***
    void actuatorCylinder::initializeLookUpTables(scalar C[ROWS][COLS]) {

    	for (label i = 0; i < ROWS; i++) {

    		for (label j = 0; j < COLS; j++) {

    			C[i][j] = 0;
    		}
    	}

    }

} 
