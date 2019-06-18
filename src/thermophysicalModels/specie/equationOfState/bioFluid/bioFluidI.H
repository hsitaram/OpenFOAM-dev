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

\*---------------------------------------------------------------------------*/

#include "bioFluid.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::bioFluid<Specie>::bioFluid
(
    const Specie& sp,
    const scalar& rhoref,
    const scalar& rhomax,
    const scalar& pref,
    const scalar& Xi
)
:
    Specie(sp),
    rhoref_(rhoref),
    rhomax_(rhomax),
    pref_(pref),
    Xi_(Xi)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::bioFluid<Specie>::bioFluid
(
    const word& name,
    const bioFluid& bfl
)
:
    Specie(name, bfl),
    rhoref_(bfl.rhoref_),
    rhomax_(bfl.rhomax_),
    pref_(bfl.pref_),
    Xi_(bfl.Xi_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::bioFluid <Specie>>
Foam::bioFluid<Specie>::clone() const
{
    return autoPtr<bioFluid<Specie>>
    (
        new bioFluid<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::bioFluid<Specie>>
Foam::bioFluid<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<bioFluid<Specie>>
    (
        new bioFluid<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::bioFluid<Specie>::rho
(
    scalar p,
    scalar T
) const
{

    const scalar expnt = scalar(1.0)-pow((p/pref_), scalar(1.0)-Xi_);
    const scalar rho   = rhoref_*pow((rhomax_/rhoref_), expnt);
    //Info <<"expnt,rho:"<<expnt<<"\t"<<rho<<"\n";
    return rho;
}


template<class Specie>
inline Foam::scalar Foam::bioFluid<Specie>::H(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::bioFluid<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::bioFluid<Specie>::S
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::bioFluid<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    const scalar expnt = scalar(1.0)-pow((p/pref_), scalar(1.0)-Xi_);
    const scalar rho   = rhoref_*pow((rhomax_/rhoref_), expnt);

    const scalar compressibility = (rho/pref_)*log(rhomax_/rhoref_)*(Xi_-scalar(1.0))*pow((p/pref_),-Xi_);
    return compressibility;
    //return (0);
}


template<class Specie>
inline Foam::scalar Foam::bioFluid<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    const scalar expnt = scalar(1.0)-pow((p/pref_), scalar(1.0)-Xi_);
    const scalar rho   = rhoref_*pow((rhomax_/rhoref_), expnt);
    return (p/(rho*this->R()*T));
    //return (0);
}


template<class Specie>
inline Foam::scalar Foam::bioFluid<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

//Don't know if I need equal operator (Hari S, Sep 17, 2018)

template<class Specie>
inline void Foam::bioFluid<Specie>::operator=
(
    const bioFluid<Specie>& bfl
)
{
    Specie::operator=(bfl);

    rhoref_  = bfl.rhoref_;
    rhomax_  = bfl.rhomax_;
    pref_    = bfl.pref_;
    Xi_      = bfl.Xi_;
}

template<class Specie>
inline void Foam::bioFluid<Specie>::operator+=
(
    const bioFluid<Specie>& bfl
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(bfl);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = bfl.Y()/this->Y();

        rhoref_ = Y1*rhoref_ + Y2*bfl.rhoref_;
        rhomax_ = Y1*rhomax_ + Y2*bfl.rhomax_;
        pref_   = Y1*pref_   + Y2*bfl.pref_;
        Xi_     = Y1*Xi_     + Y2*bfl.Xi_;
    }
}


template<class Specie>
inline void Foam::bioFluid<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::bioFluid<Specie> Foam::operator+
(
    const bioFluid<Specie>& bfl1,
    const bioFluid<Specie>& bfl2
)
{
    Specie sp
    (
        static_cast<const Specie&>(bfl1)
      + static_cast<const Specie&>(bfl2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return bioFluid<Specie>
        (
            sp,
            bfl1.rhoref_,
            bfl1.rhomax_,
            bfl1.pref_,
            bfl1.Xi_
        );
    }
    else
    {
        const scalar Y1 = bfl1.Y()/sp.Y();
        const scalar Y2 = bfl2.Y()/sp.Y();

        const scalar rhoref = Y1*bfl1.rhoref_ + Y2*bfl2.rhoref_;
        const scalar rhomax = Y1*bfl1.rhomax_ + Y2*bfl2.rhomax_;
        const scalar pref   = Y1*bfl1.pref_   + Y2*bfl2.pref_;
        const scalar Xi     = Y1*bfl1.Xi_     + Y2*bfl2.Xi_;

        return bioFluid<Specie>
        (
            sp,
            rhoref,
            rhomax,
            pref,
            Xi
        );
    }
}


template<class Specie>
Foam::bioFluid<Specie> Foam::operator*
(
    const scalar s,
    const bioFluid<Specie>& bfl
)
{
    return bioFluid<Specie>
    (
        s*static_cast<const Specie&>(bfl),
        bfl.rhoref_,
        bfl.rhomax_,
        bfl.pref_,
        bfl.Xi_
    );
}


template<class Specie>
Foam::bioFluid<Specie> Foam::operator==
(
    const bioFluid<Specie>& bfl1,
    const bioFluid<Specie>& bfl2
)
{
    Specie sp
    (
        static_cast<const Specie&>(bfl1)
     == static_cast<const Specie&>(bfl2)
    );

    const scalar Y1 = bfl1.Y()/sp.Y();
    const scalar Y2 = bfl2.Y()/sp.Y();

    const scalar rhoref = Y2*bfl2.rhoref_ - Y1*bfl1.rhoref_;
    const scalar rhomax = Y2*bfl2.rhomax_ - Y1*bfl1.rhomax_;
    const scalar pref   = Y2*bfl2.pref_   - Y1*bfl1.pref_;
    const scalar Xi     = Y2*bfl2.Xi_     - Y1*bfl1.Xi_;

    return bioFluid<Specie>
    (
        sp,
        rhoref,
        rhomax,
        pref,
        Xi
    );
}


// ************************************************************************* //
