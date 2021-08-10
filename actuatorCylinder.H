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
    Add volume force in an actuator cylinder region

    Use with actuatorCylinderSimpleFoam

    Modified by Edgar Martínez, september 2020.
\*---------------------------------------------------------------------------*/

#ifndef ACTUATORCYLINDER_H_
#define ACTUATORCYLINDER_H_

#include "fvCFD.H"
#include <string>

#define ROWS    117
#define COLS    11
#define N       36
#define NS      100

namespace Foam {

    class actuatorCylinder {

        public:

        // Runtime type information 
        TypeName("actuatorCylinder");

        // Constructor
        actuatorCylinder();

        // Destructor
        ~actuatorCylinder();

        // Set cylinder ID
        void setID(label number);

        // Read turbine axis
        void readGeometry(const fvMesh &mesh);

        // Initializae global arrays
        void initializeArrays();

        // 'x' and 'y' force acting on the shell of the cylinder
        void addForce(const fvMesh &mesh, vectorField &volumeForce, const vectorField &U);

        // Cylinder's ID
        std::string ID;

        // Power coefficient
        scalar Cp;

        private:

    	// Center vector
    	vector center;

        // Internal radius
        scalar intRadius;

        // External radius
        scalar extRadius;

        // Rotation vector
        label spin;

        // Lift coefficient array
        scalar CL[ROWS][COLS];

        // Drag coefficient array
        scalar CD[ROWS][COLS];

        // Angle of attack array
        scalar AA[ROWS];

        // Reynolds number array
        scalar RN[COLS];

        // Local azimuthal angle
        scalar t[N];

        // Control points x-velocity
        scalar Vx[N];

        // Control points y-velocity
        scalar Vy[N];

        // Volumetric force x-direction
        scalar volFx[N];

        // Volumetric force y-direction
        scalar volFy[N];

        // List of cell indices that belong to the hollow cylinder
        DynamicList<label> indices;

        // Check if a point is inside the ring
        bool pointIsInside(const vector &point);

        // Collect cells that belong to the hollow cylinder
        void collectCells(const fvMesh &mesh);

        // Collect velocities on the control points located inside the ring
        void collectVelocities(const fvMesh &mesh, const vectorField &U);

        // Compute volume forces on 'x' and 'y' directions
        void computeVolForce();

        // Induced velocities correction factor
        scalar correctVelocities(scalar Qn[N], scalar Qt[N], scalar &indFac, scalar &thrCoe);

        // Compare old and new values of induced velocities
        bool checkConvergence(scalar oldx[N], scalar newx[N], scalar oldy[N], scalar newy[N]);

        // Power coefficient
        scalar getCP(scalar Qt[N]);

        // Simpsons rule
        scalar simpsons(scalar f[N]);

        // Bilinear interpolation for lift or drag coefficient
        scalar interp2D(scalar aa, scalar re, scalar C[ROWS][COLS]);

        // Two-dimensional array parser
        void parse(std::string cname, scalar C[ROWS][COLS]);

        // Linear interpolation for volumetric forces
        scalar interpolateForces(scalar angle, scalar volF[N]);

        // Two-dimensional array initializaer
        void initializeLookUpTables(scalar C[ROWS][COLS]);

  };
  
}

#endif
