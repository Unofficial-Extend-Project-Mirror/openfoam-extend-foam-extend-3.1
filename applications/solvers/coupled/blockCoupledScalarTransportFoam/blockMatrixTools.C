/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace blockMatrixTools
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class BlockType>
void blockInsert
(
    const direction dir,
    const scalarField& x,
    Field<BlockType>& blockX
)
{
    // Reverse pointer loop below is an optimised implementation of:
    //forAll (x, i)
    //{
    //    blockX[i](dir) = x[i];
    //}
    
    // check the two fields
    checkFields
    (
        blockX,
        x,
        "blockDiag.component (s) = diag"
    );
    
    // set access to f1 and f2 at end of each field
    List_ACCESS(BlockType, blockX, f1P);
    List_CONST_ACCESS(scalar, x, f2P);

    // loop through fields performing operations
    List_FOR_ALL(x, i)
        List_ELEM(blockX, f1P, i) .component((dir))
        = List_ELEM(x, f2P, i);
    List_END_FOR_ALL
}


template<class BlockType>
void blockInsert
(
    const direction dirI,
    const direction dirJ,
    const scalarField& x,
    Field<BlockType>& blockX
)
{
    const direction dir(dirI * BlockType::rowLength + dirJ);
    blockInsert(dir, x, blockX);
}


template<class BlockType>
void blockIncrement
(
    const direction dir,
    const scalarField& x,
    Field<BlockType>& blockX
)
{
    // Reverse pointer loop below is an optimised implementation of:
    //forAll (x, i)
    //{
    //    blockX[i](dir) += x[i];
    //}

    // check the two fields
    checkFields
    (
        blockX,
        x,
        "blockDiag.component (s) += diag"
    );
    
    // set access to f1 and f2 at end of each field
    List_ACCESS(BlockType, blockX, f1P);
    List_CONST_ACCESS(scalar, x, f2P);

    // loop through fields performing operations
    List_FOR_ALL(x, i)
        List_ELEM(blockX, f1P, i) .component((dir))
        += List_ELEM(x, f2P, i);
    List_END_FOR_ALL
}


template<class BlockType>
void blockIncrement
(
    const direction dirI,
    const direction dirJ,
    const scalarField& x,
    Field<BlockType>& blockX
)
{
    const direction dir(dirI * BlockType::rowLength + dirJ);
    blockIncrement(dir, x, blockX);
}


template<class BlockType>
void blockRetrieve
(
    const direction dir,
    scalarField& x,
    const Field<BlockType>& blockX
)
{
    // Reverse pointer loop below is an optimised implementation of:
    //forAll (x, i)
    //{
    //    x[i] = blockX[i](dir);
    //}
    
    TFOR_ALL_F_OP_F_FUNC_S
    (
        scalar, x, =, BlockType, blockX, .component, const direction, dir
    )
}


template<class BlockType>
void insertDiagSource
(
    const direction dir,
    const fvScalarMatrix& m,
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& blockB
)
{
    // Prepare the diagonal and source

    scalarField diag = m.diag();
    scalarField source = m.source();

    // Add boundary source contribution
    m.addBoundaryDiag(diag, 0);
    m.addBoundarySource(source, false);

    switch (blockM.diag().activeType())
    {
        case blockCoeffBase::UNALLOCATED:
        {
            blockM.diag().asScalar() = diag;
            break;
        }
        case blockCoeffBase::SCALAR:
        case blockCoeffBase::LINEAR:
        {
            typename CoeffField<BlockType>::linearTypeField& blockDiag =
                blockM.diag().asLinear();

            // Reverse pointer loop below is an optimised implementation of:
            //forAll (diag, i)
            //{
            //    blockDiag[i](dir) = diag[i];
            //}
            
            typedef typename CoeffField<BlockType>::linearType linearType;

            // check the two fields
            checkFields
            (
                blockDiag,
                diag,
                "blockDiag.component (s) = diag"
            );
            
            // set access to f1 and f2 at end of each field
            List_ACCESS(linearType, blockDiag, f1P);
            List_CONST_ACCESS(scalar, diag, f2P);

            // loop through fields performing operations
            List_FOR_ALL(diag, i)
                List_ELEM(blockDiag, f1P, i) .component((dir))
                = List_ELEM(diag, f2P, i);
            List_END_FOR_ALL
            break;
        }
        case blockCoeffBase::SQUARE:
        {
            typename CoeffField<BlockType>::squareTypeField& blockDiag =
                blockM.diag().asSquare();

            // Reverse pointer loop below is an optimised implementation of:
            //forAll (diag, i)
            //{
            //    blockDiag[i](dir, dir) = diag[i];
            //}

            typedef typename CoeffField<BlockType>::squareType squareType;

            // check the two fields
            checkFields
            (
                blockDiag,
                diag,
                "blockDiag.component (s) = diag"
            );
            
            // set access to f1 and f2 at end of each field
            List_ACCESS(squareType, blockDiag, f1P);
            List_CONST_ACCESS(scalar, diag, f2P);

            // loop through fields performing operations
            List_FOR_ALL(diag, i)
                List_ELEM(blockDiag, f1P, i) ((dir), (dir))
                = List_ELEM(diag, f2P, i);
            List_END_FOR_ALL
        }
        default:
        {
            // do nothing
        }
    }

    blockInsert(dir, source, blockB);
}


template<class BlockType>
void addDiagSource
(
    const direction dir,
    const fvScalarMatrix& m,
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& blockB
)
{
    // Prepare the diagonal and source

    scalarField diag = m.diag();
    scalarField source = m.source();

    // Add boundary source contribution
    m.addBoundaryDiag(diag, 0);
    m.addBoundarySource(source, false);

    switch (blockM.diag().activeType())
    {
        case blockCoeffBase::UNALLOCATED:
        {
            blockM.diag().asScalar() = diag;
            break;
        }
        case blockCoeffBase::SCALAR:
        case blockCoeffBase::LINEAR:
        {
            typename CoeffField<BlockType>::linearTypeField& blockDiag =
                blockM.diag().asLinear();

            // Reverse pointer loop below is an optimised implementation of:
            //forAll (diag, i)
            //{
            //    blockDiag[i](dir) += diag[i];
            //}
            
            typedef typename CoeffField<BlockType>::linearType linearType;

            // check the two fields
            checkFields
            (
                blockDiag,
                diag,
                "blockDiag.component (s) += diag"
            );
            
            // set access to f1 and f2 at end of each field
            List_ACCESS(linearType, blockDiag, f1P);
            List_CONST_ACCESS(scalar, diag, f2P);

            // loop through fields performing operations
            List_FOR_ALL(diag, i)
                List_ELEM(blockDiag, f1P, i) .component((dir))
                += List_ELEM(diag, f2P, i);
            List_END_FOR_ALL
            break;
        }
        case blockCoeffBase::SQUARE:
        {
            typename CoeffField<BlockType>::squareTypeField& blockDiag =
                blockM.diag().asSquare();

            // Reverse pointer loop below is an optimised implementation of:
            //forAll (diag, i)
            //{
            //    blockDiag[i](dir, dir) += diag[i];
            //}

            typedef typename CoeffField<BlockType>::squareType squareType;

            // check the two fields
            checkFields
            (
                blockDiag,
                diag,
                "blockDiag.component (s) += diag"
            );
            
            // set access to f1 and f2 at end of each field
            List_ACCESS(squareType, blockDiag, f1P);
            List_CONST_ACCESS(scalar, diag, f2P);

            // loop through fields performing operations
            List_FOR_ALL(diag, i)
                List_ELEM(blockDiag, f1P, i) ((dir), (dir))
                += List_ELEM(diag, f2P, i);
            List_END_FOR_ALL
        }
        default:
        {
            // do nothing
        }
    }

    blockInsert(dir, source, blockB);
}


template<class BlockType>
void insertUpperLower
(
    const direction dir,
    const fvScalarMatrix& m,
    BlockLduMatrix<BlockType>& blockM
)
{
    if (m.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
    }

    if (m.hasUpper())
    {
        const scalarField& upper = m.upper();

        switch (blockM.upper().activeType())
        {
            case blockCoeffBase::UNALLOCATED:
            {
                blockM.upper().asScalar() = upper;
                break;
            }
            case blockCoeffBase::SCALAR:
            case blockCoeffBase::LINEAR:
            {
                typename CoeffField<BlockType>::linearTypeField& blockUpper =
                    blockM.upper().asLinear();

                // Reverse pointer loop below is an optimised implementation
                // of:
                //forAll (upper, i)
                //{
                //    blockUpper[i](dir) = upper[i];
                //}

                typedef typename CoeffField<BlockType>::linearType linearType;

                // check the two fields
                checkFields
                (
                    blockUpper,
                    upper,
                    "blockUpper.component (s) = upper"
                );
                
                // set access to f1 and f2 at end of each field
                List_ACCESS(linearType, blockUpper, f1P);
                List_CONST_ACCESS(scalar, upper, f2P);

                // loop through fields performing operations
                List_FOR_ALL(upper, i)
                    List_ELEM(blockUpper, f1P, i) .component((dir))
                    = List_ELEM(upper, f2P, i);
                List_END_FOR_ALL
                break;
            }
            case blockCoeffBase::SQUARE:
            {
                typename CoeffField<BlockType>::squareTypeField& blockUpper =
                    blockM.upper().asSquare();

                // Reverse pointer loop below is an optimised implementation
                // of:
                //forAll (upper, i)
                //{
                //    blockUpper[i](dir, dir) = upper[i];
                //}

                typedef typename CoeffField<BlockType>::squareType squareType;

                // check the two fields
                checkFields
                (
                    blockUpper,
                    upper,
                    "blockUpper.component (s) = upper"
                );
                
                // set access to f1 and f2 at end of each field
                List_ACCESS(squareType, blockUpper, f1P);
                List_CONST_ACCESS(scalar, upper, f2P);

                // loop through fields performing operations
                List_FOR_ALL(upper, i)
                    List_ELEM(blockUpper, f1P, i) ((dir), (dir))
                    = List_ELEM(upper, f2P, i);
                List_END_FOR_ALL
                break;
            }
            default:
            {
                // do nothing
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void insertUpperLower\n"
            "(\n"
            "    const direction dir,\n"
            "    const fvScalarMatrix& m,\n"
            "    BlockLduMatrix<BlockType>& blockM\n"
            ")"
        )   << "Error in matrix insertion: problem with block structure"
            << abort(FatalError);
    }

    if (!m.symmetric() || !blockM.symmetric())
    {
        // Either scalar or block matrix is asymmetric: insert lower triangle
        const scalarField& lower = m.lower();

        switch (blockM.lower().activeType())
        {
            case blockCoeffBase::UNALLOCATED:
            {
                blockM.lower().asScalar() = lower;
                break;
            }
            case blockCoeffBase::SCALAR:
            case blockCoeffBase::LINEAR:
            {
                typename CoeffField<BlockType>::linearTypeField& blockLower =
                    blockM.lower().asLinear();

                // Reverse pointer loop below is an optimised implementation
                // of:
                //forAll (lower, i)
                //{
                //    blockLower[i](dir) = lower[i];
                //}

                typedef typename CoeffField<BlockType>::linearType linearType;

                // check the two fields
                checkFields
                (
                    blockLower,
                    lower,
                    "blockLower.component (s) = lower"
                );
                
                // set access to f1 and f2 at end of each field
                List_ACCESS(linearType, blockLower, f1P);
                List_CONST_ACCESS(scalar, lower, f2P);

                // loop through fields performing operations
                List_FOR_ALL(lower, i)
                    List_ELEM(blockLower, f1P, i) .component((dir))
                    = List_ELEM(lower, f2P, i);
                List_END_FOR_ALL
                break;
            }
            case blockCoeffBase::SQUARE:
            {
                typename CoeffField<BlockType>::squareTypeField& blockLower =
                    blockM.lower().asSquare();

                // Reverse pointer loop below is an optimised implementation
                // of:
                //forAll (lower, i)
                //{
                //    blockLower[i](dir, dir) = lower[i];
                //}

                typedef typename CoeffField<BlockType>::squareType squareType;

                // check the two fields
                checkFields
                (
                    blockLower,
                    lower,
                    "blockLower.component (s) = lower"
                );
                
                // set access to f1 and f2 at end of each field
                List_ACCESS(squareType, blockLower, f1P);
                List_CONST_ACCESS(scalar, lower, f2P);

                // loop through fields performing operations
                List_FOR_ALL(lower, i)
                    List_ELEM(blockLower, f1P, i) ((dir), (dir))
                    = List_ELEM(lower, f2P, i);
                List_END_FOR_ALL
                break;
            }
            default:
            {
                // do nothing
            }
        }
    }
}


template<class BlockType>
void addUpperLower
(
    const direction dir,
    const fvScalarMatrix& m,
    BlockLduMatrix<BlockType>& blockM
)
{
    if (m.diagonal())
    {
        // Matrix for insertion is diagonal-only: nothing to do
        return;
    }

    if (m.hasUpper())
    {
        const scalarField& upper = m.upper();

        switch (blockM.upper().activeType())
        {
            case blockCoeffBase::UNALLOCATED:
            {
                blockM.upper().asScalar() = upper;
                break;
            }
            case blockCoeffBase::SCALAR:
            case blockCoeffBase::LINEAR:
            {
                typename CoeffField<BlockType>::linearTypeField& blockUpper =
                    blockM.upper().asLinear();

                // Reverse pointer loop below is an optimised implementation
                // of:
                //forAll (upper, i)
                //{
                //    blockUpper[i](dir) += upper[i];
                //}

                typedef typename CoeffField<BlockType>::linearType linearType;

                // check the two fields
                checkFields
                (
                    blockUpper,
                    upper,
                    "blockUpper.component (s) += upper"
                );
                
                // set access to f1 and f2 at end of each field
                List_ACCESS(linearType, blockUpper, f1P);
                List_CONST_ACCESS(scalar, upper, f2P);

                // loop through fields performing operations
                List_FOR_ALL(upper, i)
                    List_ELEM(blockUpper, f1P, i) .component((dir))
                    += List_ELEM(upper, f2P, i);
                List_END_FOR_ALL
                break;
            }
            case blockCoeffBase::SQUARE:
            {
                typename CoeffField<BlockType>::squareTypeField& blockUpper =
                    blockM.upper().asSquare();

                // Reverse pointer loop below is an optimised implementation
                // of:
                //forAll (upper, i)
                //{
                //    blockUpper[i](dir, dir) += upper[i];
                //}

                typedef typename CoeffField<BlockType>::squareType squareType;

                // check the two fields
                checkFields
                (
                    blockUpper,
                    upper,
                    "blockUpper.component (s) += upper"
                );
                
                // set access to f1 and f2 at end of each field
                List_ACCESS(squareType, blockUpper, f1P);
                List_CONST_ACCESS(scalar, upper, f2P);

                // loop through fields performing operations
                List_FOR_ALL(upper, i)
                    List_ELEM(blockUpper, f1P, i) ((dir), (dir))
                    += List_ELEM(upper, f2P, i);
                List_END_FOR_ALL
                break;
            }
            default:
            {
                // do nothing
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void insertUpperLower\n"
            "(\n"
            "    const direction dir,\n"
            "    const fvScalarMatrix& m,\n"
            "    BlockLduMatrix<BlockType>& blockM\n"
            ")"
        )   << "Error in matrix insertion: problem with block structure"
            << abort(FatalError);
    }

    if (!m.symmetric() || !blockM.symmetric())
    {
        // Either scalar or block matrix is asymmetric: insert lower triangle
        const scalarField& lower = m.lower();

        switch (blockM.lower().activeType())
        {
            case blockCoeffBase::UNALLOCATED:
            {
                blockM.lower().asScalar() = lower;
                break;
            }
            case blockCoeffBase::SCALAR:
            case blockCoeffBase::LINEAR:
            {
                typename CoeffField<BlockType>::linearTypeField& blockLower =
                    blockM.lower().asLinear();

                // Reverse pointer loop below is an optimised implementation
                // of:
                //forAll (lower, i)
                //{
                //    blockLower[i](dir) += lower[i];
                //}

                typedef typename CoeffField<BlockType>::linearType linearType;

                // check the two fields
                checkFields
                (
                    blockLower,
                    lower,
                    "blockLower.component (s) += lower"
                );
                
                // set access to f1 and f2 at end of each field
                List_ACCESS(linearType, blockLower, f1P);
                List_CONST_ACCESS(scalar, lower, f2P);

                // loop through fields performing operations
                List_FOR_ALL(lower, i)
                    List_ELEM(blockLower, f1P, i) .component((dir))
                    += List_ELEM(lower, f2P, i);
                List_END_FOR_ALL
                break;
            }
            case blockCoeffBase::SQUARE:
            {
                typename CoeffField<BlockType>::squareTypeField& blockLower =
                    blockM.lower().asSquare();

                // Reverse pointer loop below is an optimised implementation
                // of:
                //forAll (lower, i)
                //{
                //    blockLower[i](dir, dir) += lower[i];
                //}

                typedef typename CoeffField<BlockType>::squareType squareType;

                // check the two fields
                checkFields
                (
                    blockLower,
                    lower,
                    "blockLower.component (s) += lower"
                );
                
                // set access to f1 and f2 at end of each field
                List_ACCESS(squareType, blockLower, f1P);
                List_CONST_ACCESS(scalar, lower, f2P);

                // loop through fields performing operations
                List_FOR_ALL(lower, i)
                    List_ELEM(blockLower, f1P, i) ((dir), (dir))
                    += List_ELEM(lower, f2P, i);
                List_END_FOR_ALL
                break;
            }
            default:
            {
                // do nothing
            }
        }
    }
}


template<class BlockType>
void insertEquation
(
    const direction dir,
    const fvScalarMatrix& m,
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& blockX,
    Field<BlockType>& blockB
)
{
    insertDiagSource(dir, m, blockM, blockB);
    insertUpperLower(dir, m, blockM);
    blockInsert(dir, m.psi(), blockX);
}


template<class BlockType>
void addEquation
(
    const direction dir,
    const fvScalarMatrix& m,
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& blockX,
    Field<BlockType>& blockB
)
{
    addDiagSource(dir, m, blockM, blockB);
    addUpperLower(dir, m, blockM);
    blockIncrement(dir, m.psi(), blockX);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockMatrixTools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
