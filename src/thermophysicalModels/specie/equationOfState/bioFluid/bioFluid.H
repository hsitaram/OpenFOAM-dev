/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 OpenFOAM Foundation
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

Class
    Foam::bioFluid

Description
    bioFluid gas equation of state.

SourceFiles
    bioFluidI.H
    bioFluid.C

\*---------------------------------------------------------------------------*/

#ifndef bioFluid_H
#define bioFluid_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class bioFluid;

template<class Specie>
inline bioFluid<Specie> operator+
(
    const bioFluid<Specie>&,
    const bioFluid<Specie>&
);

template<class Specie>
inline bioFluid<Specie> operator*
(
    const scalar,
    const bioFluid<Specie>&
);

template<class Specie>
inline bioFluid<Specie> operator==
(
    const bioFluid<Specie>&,
    const bioFluid<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const bioFluid<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class bioFluid Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class bioFluid
:
    public Specie
{
    // Private data

        scalar rhoref_;

        scalar rhomax_;

        scalar pref_;

        scalar Xi_;

public:

    // Constructors

        //- Construct from components
        inline bioFluid
        (
            const Specie& sp,
            const scalar& rhoref,
            const scalar& rhomax,
            const scalar& pref,
            const scalar& Xi
        );

        //- Construct from dictionary
        bioFluid(const dictionary& dict);

        //- Construct as named copy
        inline bioFluid(const word& name, const bioFluid&);

        //- Construct and return a clone
        inline autoPtr<bioFluid> clone() const;

        // Selector from dictionary
        inline static autoPtr<bioFluid> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "bioFluid<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


    // IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator=(const bioFluid&);
        inline void operator+=(const bioFluid&);
        inline void operator*=(const scalar);


    // Friend operators

        friend bioFluid operator+ <Specie>
        (
            const bioFluid&,
            const bioFluid&
        );

        friend bioFluid operator* <Specie>
        (
            const scalar s,
            const bioFluid&
        );

        friend bioFluid operator== <Specie>
        (
            const bioFluid&,
            const bioFluid&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const bioFluid&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bioFluidI.H"

#ifdef NoRepository
    #include "bioFluid.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
