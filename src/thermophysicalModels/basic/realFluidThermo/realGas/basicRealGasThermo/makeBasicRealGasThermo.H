/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

InClass
    Foam::basicRealFluidThermo

Description
    Macros for creating 'basic' compresibility-based thermo packages

Modified by
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany


\*---------------------------------------------------------------------------*/

#ifndef makeBasicRealGasThermo_H
#define makeBasicRealGasThermo_H

#include "basicRealGasThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeBasicRealGasThermo(Cthermo,Mixture,Transport,SpecieThermo,Thermo,EqnOfState)       \
                                                                              \
typedef Cthermo<Mixture<Transport<SpecieThermo<Thermo<EqnOfState> > > > >     \
    Cthermo##Mixture##Transport##SpecieThermo##Thermo##EqnOfState;                          \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
(                                                                             \
    Cthermo##Mixture##Transport##SpecieThermo##Thermo##EqnOfState,                          \
    #Cthermo                                                                  \
        "<"#Mixture"<"#Transport"<"#SpecieThermo"<"#Thermo"<"#EqnOfState">>>>>", \
    0                                                                         \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    basicRealGasThermo,                                                           \
    Cthermo##Mixture##Transport##SpecieThermo##Thermo##EqnOfState,                          \
    fvMesh                                                                    \
)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
