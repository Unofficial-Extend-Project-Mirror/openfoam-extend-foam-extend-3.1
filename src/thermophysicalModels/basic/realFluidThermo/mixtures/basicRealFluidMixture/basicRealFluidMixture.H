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

Class
    Foam::basicRealFluidMixture

Description
    Foam::basicRealFluidMixture

SourceFiles
    basicRealFluidMixture.C

Modified by
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany
\*---------------------------------------------------------------------------*/

#ifndef basicRealFluidMixture_H
#define basicRealFluidMixture_H

#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class dictionary;

/*---------------------------------------------------------------------------*\
                        Class basicRealFluidMixture Declaration
\*---------------------------------------------------------------------------*/

class basicRealFluidMixture
{

public:

    // Runtime type information
    TypeName("basicRealFluidMixture");


    // Constructors

        //- Construct from dictionary and mesh
        basicRealFluidMixture(const dictionary&, const fvMesh&);


    //- Destructor
    virtual ~basicRealFluidMixture();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
