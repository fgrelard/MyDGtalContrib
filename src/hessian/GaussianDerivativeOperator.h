/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file GaussianDerivativeOperator.h
 * @author Florent Grelard (\c florent.grelard@u-bordeaux.fr )
 * LaBRI, Bordeaux University
 *
 * @date 2017/03/27
 *
 *
 * This file is part of the DGtal library.
 */


#if defined(GaussianDerivativeOperator_RECURSES)
#error Recursive header files inclusion detected in GaussianDerivativeOperator.h
#else // defined(GaussianDerivativeOperator_RECURSES)
/** Prevents recursive inclusion of headers. */
#define GaussianDerivativeOperator_RECURSES

#if !defined GaussianDerivativeOperator_h
/** Prevents repeated inclusion of headers. */
#define GaussianDerivativeOperator_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <cmath>
#include <numeric>
#include <limits>
#include "DerivativeOperator.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal {

    /////////////////////////////////////////////////////////////////////////////
    // template class GaussianDerivativeOperator
    /**
     * Description of template class 'GaussianDerivativeOperator' <p>
     * \brief Aim:
     */
    template <typename TValue>
    class GaussianDerivativeOperator {

    public:
        typedef TValue Scalar;
        typedef DerivativeOperator<Scalar> Derivative;
        typedef typename Derivative::CoefficientVector CoefficientVector;

        // ----------------------- Standard services ------------------------------
    public:

        /**
         * Default constructor.
         */
        GaussianDerivativeOperator() {
            myNormalizeAcrossScale = true;
            myOrder = 1;
            myVariance = 1.0;
        }

        /**
         * Destructor.
         */
        ~GaussianDerivativeOperator() = default;

        GaussianDerivativeOperator(unsigned int anOrder, bool aNormalizeAcrossScale, double aVariance = 1.0) : myOrder(
                anOrder), myNormalizeAcrossScale(aNormalizeAcrossScale), myVariance(aVariance) {}

        /**
         * Copy constructor.
         * @param other the object to clone.
         */
        GaussianDerivativeOperator(const GaussianDerivativeOperator &other) = default;

        /**
         * Move constructor.
         * @param other the object to move.
         */
        GaussianDerivativeOperator(GaussianDerivativeOperator &&other) = delete;

        /**
         * Copy assignment operator.
         * @param other the object to copy.
         * @return a reference on 'this'.
         */
        GaussianDerivativeOperator &operator=(const GaussianDerivativeOperator &other) = delete;

        /**
         * Move assignment operator.
         * @param other the object to move.
         * @return a reference on 'this'.
         */
        GaussianDerivativeOperator &operator=(GaussianDerivativeOperator &&other) = delete;

        // ----------------------- Interface --------------------------------------
    public:

        CoefficientVector computeCoefficients();

        /**
         * Writes/Displays the object on an output stream.
         * @param out the output stream where the object is written.
         */
        void selfDisplay(std::ostream &out) const;

        /**
         * Checks the validity/consistency of the object.
         * @return 'true' if the object is valid, 'false' otherwise.
         */
        bool isValid() const;

        void setVariance(double variance);

        double getVariance() const;

        void setOrder(unsigned int anOrder);

        void setNormalizeAcrossScale(bool normalize);

        // ------------------------- Private Methods ------------------------------
    private:


        // ------------------------- Protected Datas ------------------------------
    protected:


        // ------------------------- Private Datas --------------------------------
    private:
        unsigned int myOrder;
        bool myNormalizeAcrossScale;
        double myVariance;

        // ------------------------- Hidden services ------------------------------
    protected:
        double modifiedBesselI0(double) const;

        double modifiedBesselI1(double) const;

        double modifiedBesselI(int, double) const;

        CoefficientVector computeGaussianCoefficients() const;



        // ------------------------- Internals ------------------------------------
    private:

    };// end of class GaussianDerivativeOperator


    /**
     * Overloads 'operator<<' for displaying objects of class 'GaussianDerivativeOperator'.
     * @param out the output stream where the object is written.
     * @param object the object of class 'GaussianDerivativeOperator' to write.
     * @return the output stream after the writing.
     */
    template <typename TValue>
    std::ostream &
    operator<<(std::ostream &out, const GaussianDerivativeOperator<TValue> &object);


} // namespace DGtal

template <typename TValue>
typename DGtal::GaussianDerivativeOperator<TValue>::CoefficientVector
DGtal::GaussianDerivativeOperator<TValue>::
computeCoefficients() {

    // compute gaussian kernel of 0-order
    CoefficientVector coeff = this->computeGaussianCoefficients();

    if (myOrder == 0) {
        return coeff;
    }


    // Calculate scale-space normalization factor for derivatives
    double norm;
    if (myNormalizeAcrossScale && myOrder) {
        norm = std::pow(myVariance, myOrder / 2.0);
    } else {
        norm = 1.0;
    }


    Derivative derivOp(myOrder);
    CoefficientVector coeffDeriv = derivOp.computeCoefficients();

    // The input gaussian kernel needs to be padded with a clamped
    // boundary condition. If N is the radius of the derivative
    // operator, then the output kernel needs to be padded by N-1. For
    // these values to be computed the input kernel needs to be padded
    // by 2N-1 on both sides.
    unsigned int N = (coeffDeriv.size() - 1) / 2;

    // copy the gaussian operator adding clamped boundary condition
    CoefficientVector paddedCoeff(coeff.size() + 4 * N - 2);

    // copy the whole gaussuan operator in coeff to paddedCoef
    // starting after the padding
    std::copy(coeff.begin(), coeff.end(), paddedCoeff.begin() + 2 * N - 1);

    // padd paddedCoeff with 2*N-1 number of boundary conditions
    std::fill(paddedCoeff.begin(), paddedCoeff.begin() + 2 * N, coeff.front());
    std::fill(paddedCoeff.rbegin(), paddedCoeff.rbegin() + 2 * N, coeff.back());

    // clear for output kernel/coeffs
    coeff = CoefficientVector();

    // Now perform convolution between derivative operators and padded gaussian
    for (unsigned int i = N; i < paddedCoeff.size() - N; ++i) {
        double sum = 0.0;

        // current index in derivative op
        for (unsigned int j = 0; j < coeffDeriv.size(); ++j) {
            unsigned int k = i + j - coeffDeriv.size() / 2;
            sum += paddedCoeff[k] * coeffDeriv[coeffDeriv.size() - 1 - j];
        }

        // normalize for scale-space and spacing
        coeff.push_back(norm * sum);
    }

    return coeff;
}


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
template <typename TValue>
typename DGtal::GaussianDerivativeOperator<TValue>::CoefficientVector
DGtal::GaussianDerivativeOperator<TValue>::
computeGaussianCoefficients() const {
    CoefficientVector coeff;

    // Use image spacing to modify variance
    const double pixelVariance = myVariance;

    // Now create coefficients as if they were zero order coeffs
    const double et = std::exp(-pixelVariance);
    const double cap = 0.995;
    double sum = 0.0;

    // Create the kernel coefficients as a std::vector
    coeff.push_back(et * modifiedBesselI0(pixelVariance));
    sum += coeff[0];
    coeff.push_back(et * modifiedBesselI1(pixelVariance));
    sum += coeff[1] * 2.0;

    for (int i = 2; sum < cap; i++) {
        coeff.push_back(et * modifiedBesselI(i, pixelVariance));
        sum += coeff[i] * 2.0;
        if (coeff[i] < sum * std::numeric_limits<double>::epsilon()) {
            break;
        }
    }

    // re-accumulate from smallest number to largest for maximum precision
    sum = std::accumulate(coeff.rbegin(), coeff.rend() - 1, 0.0);
    sum *= 2.0;
    sum += coeff[0]; // the first is only needed once

    // Normalize the coefficients so they sum one
    for (typename CoefficientVector::iterator it = coeff.begin(); it != coeff.end(); ++it) {
        *it /= sum;
    }

    // Make symmetric
    size_t s = coeff.size() - 1;
    coeff.insert(coeff.begin(), s, 0);
    std::copy(coeff.rbegin(), coeff.rbegin() + s, coeff.begin());

    return coeff;
}


template <typename TValue>
void
DGtal::GaussianDerivativeOperator<TValue>::
setVariance(double variance) {
    myVariance = variance;
}

template <typename TValue>
double
DGtal::GaussianDerivativeOperator<TValue>::
getVariance() const {
    return myVariance;
}

template <typename TValue>
void
DGtal::GaussianDerivativeOperator<TValue>::
setOrder(unsigned int anOrder) {
    myOrder = anOrder;
}

template <typename TValue>
void
DGtal::GaussianDerivativeOperator<TValue>::
setNormalizeAcrossScale(bool normalize) {
    myNormalizeAcrossScale = normalize;
}


template <typename TValue>
double
DGtal::GaussianDerivativeOperator<TValue>::
modifiedBesselI0(double y) const {
    double d, accumulator;
    double m;

    if ((d = std::fabs(y)) < 3.75) {
        m = y / 3.75;
        m *= m;
        accumulator = 1.0 + m * (3.5156229 + m * (3.0899424 + m * (1.2067492
                                                                   + m
                                                                     * (0.2659732 +
                                                                        m * (0.360768e-1 + m * 0.45813e-2)))));
    } else {
        m = 3.75 / d;
        accumulator = (std::exp(d) / std::sqrt(d)) * (0.39894228 + m * (0.1328592e-1
                                                                        + m
                                                                          * (0.225319e-2 + m
                                                                                           * (-0.157565e-2 +
                                                                                              m * (0.916281e-2
                                                                                                   +
                                                                                                   m
                                                                                                   * (-0.2057706e-1
                                                                                                      + m *
                                                                                                        (0.2635537e-1 +
                                                                                                         m *
                                                                                                         (-0.1647633e-1
                                                                                                          +
                                                                                                          m
                                                                                                          *
                                                                                                          0.392377e-2))))))));
    }
    return accumulator;
}

template <typename TValue>
double
DGtal::GaussianDerivativeOperator<TValue>::
modifiedBesselI1(double y) const {
    double d, accumulator;
    double m;

    if ((d = std::fabs(y)) < 3.75) {
        m = y / 3.75;
        m *= m;
        accumulator = d * (0.5 + m * (0.87890594 + m * (0.51498869 + m * (0.15084934
                                                                          + m
                                                                            * (0.2658733e-1 + m
                                                                                              * (0.301532e-2 +
                                                                                                 m * 0.32411e-3))))));
    } else {
        m = 3.75 / d;
        accumulator = 0.2282967e-1 + m * (-0.2895312e-1 + m * (0.1787654e-1
                                                               - m * 0.420059e-2));
        accumulator = 0.39894228 + m * (-0.3988024e-1 + m * (-0.362018e-2
                                                             + m *
                                                               (0.163801e-2 + m * (-0.1031555e-1 + m * accumulator))));

        accumulator *= (std::exp(d) / std::sqrt(d));
    }

    if (y < 0.0) { return -accumulator; }
    else { return accumulator; }
}

template <typename TValue>
double
DGtal::GaussianDerivativeOperator<TValue>::
modifiedBesselI(int n, double y) const {
    const double DIGITS = 10.0;
    int j;
    double qim, qi, qip, toy;
    double accumulator;

    if (y == 0.0) { return 0.0; }
    else {
        toy = 2.0 / std::fabs(y);
        qip = accumulator = 0.0;
        qi = 1.0;
        for (j = 2 * (n + (int) (DIGITS * std::sqrt((double) n))); j > 0; j--) {
            qim = qip + j * toy * qi;
            qip = qi;
            qi = qim;
            if (std::fabs(qi) > 1.0e10) {
                accumulator *= 1.0e-10;
                qi *= 1.0e-10;
                qip *= 1.0e-10;
            }
            if (j == n) { accumulator = qip; }
        }
        accumulator *= modifiedBesselI0(y) / qi;
        if (y < 0.0 && (n & 1)) { return -accumulator; }
        else { return accumulator; }
    }
}


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined GaussianDerivativeOperator_h

#undef GaussianDerivativeOperator_RECURSES
#endif // else defined(GaussianDerivativeOperator_RECURSES)