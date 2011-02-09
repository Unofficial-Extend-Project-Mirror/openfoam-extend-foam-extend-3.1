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

Modified by
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "basicRealGasThermo.H"
#include "fvMesh.H"
#include "HashTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedEnthalpyRealFluids.H"
#include "gradientEnthalpyRealFluids.H"
#include "mixedEnthalpyRealFluids.H"
#include "fixedInternalEnergyRealFluids.H"
#include "gradientInternalEnergyRealFluids.H"
#include "mixedInternalEnergyRealFluids.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicRealGasThermo, 0);
    defineRunTimeSelectionTable(basicRealGasThermo, fvMesh);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::basicRealGasThermo::hRealBoundaryTypes()
{
    const volScalarField::GeometricBoundaryField& tbf = T_.boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = fixedEnthalpyRealFluids::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            hbt[patchi] = gradientEnthalpyRealFluids::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            hbt[patchi] = mixedEnthalpyRealFluids::typeName;
        }
    }

    return hbt;
}


void Foam::basicRealGasThermo::hRealBoundaryCorrection(volScalarField& h)
{
    volScalarField::GeometricBoundaryField& hbf = h.boundaryField();

    forAll(hbf, patchi)
    {
        if (isA<gradientEnthalpyRealFluids>(hbf[patchi]))
        {
            refCast<gradientEnthalpyRealFluids>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnthalpyRealFluids>(hbf[patchi]))
        {
            refCast<mixedEnthalpyRealFluids>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


Foam::wordList Foam::basicRealGasThermo::eRealBoundaryTypes()
{
    const volScalarField::GeometricBoundaryField& tbf = T_.boundaryField();

    wordList ebt = tbf.types();

    forAll(tbf, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
        {
            ebt[patchi] = fixedInternalEnergyRealFluids::typeName;
        }
        else if
        (
            isA<zeroGradientFvPatchScalarField>(tbf[patchi])
         || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
        )
        {
            ebt[patchi] = gradientInternalEnergyRealFluids::typeName;
        }
        else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
        {
            ebt[patchi] = mixedInternalEnergyRealFluids::typeName;
        }
    }

    return ebt;
}


void Foam::basicRealGasThermo::eRealBoundaryCorrection(volScalarField& e)
{
    volScalarField::GeometricBoundaryField& ebf = e.boundaryField();

    forAll(ebf, patchi)
    {
        if (isA<gradientInternalEnergyRealFluids>(ebf[patchi]))
        {
            refCast<gradientInternalEnergyRealFluids>(ebf[patchi])
                .gradient() = ebf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedInternalEnergyRealFluids>(ebf[patchi]))
        {
            refCast<mixedInternalEnergyRealFluids>(ebf[patchi])
                .refGrad() = ebf[patchi].fvPatchField::snGrad();
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicRealGasThermo::basicRealGasThermo(const fvMesh& mesh)
:
    basicThermo(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicRealGasThermo::~basicRealGasThermo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



Foam::tmp<Foam::scalarField> Foam::basicRealGasThermo::hBC
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    notImplemented
    (
        "basicRealGasThermo::hBC"
        "(const scalarField& p,const scalarField& T, const label& patchi) const"
    );
    return tmp<scalarField>(NULL);
  
}

Foam::tmp<Foam::scalarField> Foam::basicRealGasThermo::hBC
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
       notImplemented
    (
        "basicRealGasThermo::hBC"
        "(const scalarField& p, const scalarField& T, const labelList& cells) const"
    );
    return tmp<scalarField>(NULL);
}

Foam::tmp<Foam::scalarField> Foam::basicRealGasThermo::eBC
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
           notImplemented
    (
        "basicRealGasThermo::eBC"
        "(const scalarField& p, const scalarField& T, const label& patchi ) const"
    );
    return tmp<scalarField>(NULL);
}

Foam::tmp<Foam::scalarField> Foam::basicRealGasThermo::eBC
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
           notImplemented
    (
        "basicRealGasThermo::eBC"
        "(const scalarField& p, const scalarField& T, const labelList& cells ) const"
    );
    return tmp<scalarField>(NULL);
}

Foam::tmp<Foam::scalarField> Foam::basicRealGasThermo::CpBC
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
           notImplemented
    (
        "basicRealGasThermo::CpBC"
        "(const scalarField& p, const scalarField& T, const label& patchi ) const"
    );
    return tmp<scalarField>(NULL);
}

Foam::tmp<Foam::scalarField> Foam::basicRealGasThermo::CpBC
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    notImplemented
    (
        "basicRealGasThermo::CpBC"
        "(const scalarField& p, const scalarField& T, const labelList& cells ) const"
    );
    return tmp<scalarField>(NULL);
}

Foam::tmp<Foam::scalarField> Foam::basicRealGasThermo::CvBC
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
           notImplemented
    (
        "basicRealGasThermo::CvBC"
        "(const scalarField& p, const scalarField& T, const label& patchi ) const"
    );
    return tmp<scalarField>(NULL);
}

Foam::tmp<Foam::scalarField> Foam::basicRealGasThermo::CvBC
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    notImplemented
    (
        "basicRealGasThermo::CvBC"
        "(const scalarField& p, const scalarField& T, const labelList& cells ) const"
    );
    return tmp<scalarField>(NULL);
}

bool Foam::basicRealGasThermo::read()
{
    return regIOobject::read();
}
// ************************************************************************* //
