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
    Foam::basicRealGasThermo

Description
    Basic thermodynamic properties class for real gas fluids
SourceFiles
    basicRealGasThermo.C
    newBasicRealGasThermo.C

Modified by
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany


\*---------------------------------------------------------------------------*/

#ifndef basicRealGasThermo_H
#define basicRealGasThermo_H


#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "typeInfo.H"
#include "basicThermo.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class basicRealFluidThermo Declaration
\*---------------------------------------------------------------------------*/

class basicRealGasThermo
:
    public basicThermo
{

protected:

    // Protected member functions

        //- Construct as copy (not implemented)
        basicRealGasThermo(const basicRealGasThermo&);

            //- Return the enthalpy field boundary types by interrogating the
            //  temperature field boundary types
            wordList hRealBoundaryTypes();

            //- Correct the enthalpy field boundaries
            void hRealBoundaryCorrection(volScalarField& h);


        // Internal energy

            //- Return the internal energy  field boundary types by
            //  interrogating the temperature field boundary types
            wordList eRealBoundaryTypes();

            //- Correct the internal energy field boundaries
            void eRealBoundaryCorrection(volScalarField& e);


public:

    //- Runtime type information
    TypeName("basicRealGasThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicRealGasThermo,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );

    // Constructors

        //- Construct from mesh
        basicRealGasThermo(const fvMesh&);


    //- Selector
    static autoPtr<basicRealGasThermo> New(const fvMesh&);

    //- Destructor
    virtual ~basicRealGasThermo();


    // Member functions

        //- Update properties
        virtual void correct() = 0;

    // Access to thermodynamic state variables for the Boundary Functions

            //-  enthalpy for cell-set [J/kg]
            //function needed for BC
            virtual tmp<scalarField> hBC
            (
                const scalarField& rho,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- enthalpy for patch [J/kg]
	    //function needed for BC
            virtual tmp<scalarField> hBC
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //-  internal energy for cell-set [J/kg]
            //function needed for BC
            virtual tmp<scalarField> eBC
            (
                const scalarField& rho,
                const scalarField& T,
                const labelList& cells
            ) const;

            //-  internal energy for patch [J/kg]
	    //function needed for BC
            virtual tmp<scalarField> eBC
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for cell-set [J/kg/K]
            //function needed for BC
            virtual tmp<scalarField> CpBC
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            //function needed for BC
            virtual tmp<scalarField> CpBC
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant spec.Volume for cell-set [J/kg/K]
            //function needed for BC
            virtual tmp<scalarField> CvBC
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Heat capacity at constant spec.Volume for patch [J/kg/K]
            //function needed for BC
            virtual tmp<scalarField> CvBC
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

        //- Read thermophysicalProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
