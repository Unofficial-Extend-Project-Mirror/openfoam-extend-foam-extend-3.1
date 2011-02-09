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
    Foam::pureRealFluidMixture

Description
    Foam::pureRealFluidMixture

SourceFiles
    pureRealFluidMixture.C

Modified by
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#ifndef pureRealFluidMixture_H
#define pureRealFluidMixture_H

#include "basicRealFluidMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pureMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class pureRealFluidMixture
:
    public basicRealFluidMixture
{
    // Private data

        ThermoType mixture_;

        //- Construct as copy (not implemented)
        pureRealFluidMixture(const pureRealFluidMixture<ThermoType>&);


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    //- Runtime type information
    TypeName("pureRealFluidMixture");


    // Constructors

        //- Construct from dictionary and mesh
        pureRealFluidMixture(const dictionary&, const fvMesh&);


    // Destructor

        ~pureRealFluidMixture();


    // Member functions

        const ThermoType& cellMixture(const label) const
        {
            return mixture_;
        }

        const ThermoType& patchFaceMixture
        (
            const label,
            const label
        ) const
        {
            return mixture_;
        }

        //- Read dictionary
        void read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
#   include "pureRealFluidMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
