/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "allToCell.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(allToCell, 0);

addToRunTimeSelectionTable(topoSetSource, allToCell, word);

addToRunTimeSelectionTable(topoSetSource, allToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::allToCell::usage_
(
    allToCell::typeName,
    "\n    Usage: allToCell \n\n"
    "    Select all available cells\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::allToCell::combine(topoSet& set, const bool add) const
{
    // We simply add everything to the set 
    const cellList& cells = mesh_.cells();

    forAll(cells, cellI)
    {
      addOrDelete(set, cellI, add);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::allToCell::allToCell
(
    const polyMesh& mesh
)
:
    topoSetSource(mesh)
{}


// Construct from dictionary
Foam::allToCell::allToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh)
{}


// Construct from Istream
Foam::allToCell::allToCell
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::allToCell::~allToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::allToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all cells "  << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all cells " << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
