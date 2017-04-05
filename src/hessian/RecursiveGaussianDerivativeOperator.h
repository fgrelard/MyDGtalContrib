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
 * @file RecursiveGaussianDerivativeOperator.h
 * @author Florent Grelard (\c florent.grelard@u-bordeaux.fr )
 * LaBRI, Bordeaux University
 *
 * @date 2017/04/04
 *
 *
 * This file is part of the DGtal library.
 */

#if defined(RecursiveGaussianDerivativeOperator_RECURSES)
#error Recursive header files inclusion detected in RecursiveGaussianDerivativeOperator.h
#else // defined(RecursiveGaussianDerivativeOperator_RECURSES)
/** Prevents recursive inclusion of headers. */
#define RecursiveGaussianDerivativeOperator_RECURSES

#if !defined RecursiveGaussianDerivativeOperator_h
/** Prevents repeated inclusion of headers. */
#define RecursiveGaussianDerivativeOperator_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <cmath>
#include <vector>
#include <DGtal/images/ImageContainerBySTLVector.h>
//////////////////////////////////////////////////////////////////////////////

namespace DGtal {

    /////////////////////////////////////////////////////////////////////////////
    // template class RecursiveGaussianDerivativeOperator
    /**
     * Description of template class 'RecursiveGaussianDerivativeOperator' <p>
     * \brief Aim:
     */
    template <typename TInputImage>
    class RecursiveGaussianDerivativeOperator {
    public:
        typedef TInputImage InputImage;
        typedef typename InputImage::Domain Domain;
        typedef typename Domain::Point Point;
        typedef double Scalar;
        typedef InputImage OutputImage;


        // ----------------------- Standard services ------------------------------
    public:

        /**
         * Default constructor.
         */
        RecursiveGaussianDerivativeOperator() = delete;

        /**
         * Destructor.
         */
        ~RecursiveGaussianDerivativeOperator() = default;

        RecursiveGaussianDerivativeOperator(const InputImage &anImage) : myImage(anImage), mySigma(1.0), myOrder(0),
                                                                         myNormalizeAcrossScale(true) {}

        RecursiveGaussianDerivativeOperator(const InputImage &anImage, double sigma, unsigned int anOrder,
                                            bool normalize = true) : myImage(anImage), mySigma(sigma), myOrder(anOrder),
                                                                     myNormalizeAcrossScale(normalize) {}


        /**
         * Copy constructor.
         * @param other the object to clone.
         */
        RecursiveGaussianDerivativeOperator(const RecursiveGaussianDerivativeOperator &other) = default;

        /**
         * Move constructor.
         * @param other the object to move.
         */
        RecursiveGaussianDerivativeOperator(RecursiveGaussianDerivativeOperator &&other) = delete;

        /**
         * Copy assignment operator.
         * @param other the object to copy.
         * @return a reference on 'this'.
         */
        RecursiveGaussianDerivativeOperator &operator=(const RecursiveGaussianDerivativeOperator &other) = delete;

        /**
         * Move assignment operator.
         * @param other the object to move.
         * @return a reference on 'this'.
         */
        RecursiveGaussianDerivativeOperator &operator=(RecursiveGaussianDerivativeOperator &&other) = delete;

        // ----------------------- Interface --------------------------------------
    public:

        OutputImage gaussianImage();

        void computeCoefficients();

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

        void setSigma(double sigma);

        double getSigma() const;

        void setOrder(unsigned int anOrder);

        void setNormalizeAcrossScale(bool normalize);


        // ------------------------- Protected Datas ------------------------------
    protected:

        // ------------------------- Private Datas --------------------------------
    private:
        InputImage myImage;
        double mySigma;
        unsigned int myOrder;
        bool myNormalizeAcrossScale;
        // ------------------------- Hidden services ------------------------------
    protected:
        void filterDataArray(Scalar *outs, const Scalar *data, Scalar *scratch, size_t ln);

        void computeNCoefficients(Scalar sigmad,
                                  Scalar A1, Scalar B1, Scalar W1, Scalar L1,
                                  Scalar A2, Scalar B2, Scalar W2, Scalar L2,
                                  Scalar &N0, Scalar &N1,
                                  Scalar &N2, Scalar &N3,
                                  Scalar &SN, Scalar &DN, Scalar &EN);


        /** Compute the D coefficients in the recursive filter. */
        void computeDCoefficients(Scalar sigmad,
                                  Scalar W1, Scalar L1, Scalar W2, Scalar L2,
                                  Scalar &SD, Scalar &DD, Scalar &ED);

        /** Compute the M coefficients and the boundary coefficients in the
   * recursive filter. */
        void computeRemainingCoefficients(bool symmetric);

        template <typename T1, typename T2>
        inline void linearCombination(T1 &out,
                                      const T1 &a1, const T2 &b1,
                                      const T1 &a2, const T2 &b2,
                                      const T1 &a3, const T2 &b3,
                                      const T1 &a4, const T2 &b4) {
            out = a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4;
        }


        template <typename T1, typename T2>
        inline void linearCombination(std::vector<T1> &out,
                                      const std::vector<T1> &a1, const T2 &b1,
                                      const std::vector<T1> &a2, const T2 &b2,
                                      const std::vector<T1> &a3, const T2 &b3,
                                      const std::vector<T1> &a4, const T2 &b4) {
            const unsigned int sz = a1.size();
            if (sz != out.size()) {
                out.resize(sz);
            }
            for (unsigned int i = 0; i < sz; ++i) {
                out[i] = a1[i] * b1 + a2[i] * b2 + a3[i] * b3 + a4[i] * b4;
            }
        }

        template <typename T1, typename T2>
        inline void subLinearCombination(T1 &out,
                                         const T1 &a1, const T2 &b1,
                                         const T1 &a2, const T2 &b2,
                                         const T1 &a3, const T2 &b3,
                                         const T1 &a4, const T2 &b4) {
            out -= a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4;
        }

        template <typename T1, typename T2>
        inline void subLinearCombination(std::vector<T1> &out,
                                         const std::vector<T1> &a1, const T2 &b1,
                                         const std::vector<T1> &a2, const T2 &b2,
                                         const std::vector<T1> &a3, const T2 &b3,
                                         const std::vector<T1> &a4, const T2 &b4) {
            const unsigned int sz = a1.size();
            if (sz != out.size()) {
                out.resize(sz);
            }
            for (unsigned int i = 0; i < sz; ++i) {
                out[i] -= a1[i] * b1 + a2[i] * b2 + a3[i] * b3 + a4[i] * b4;
            }
        }

        // ------------------------- Internals ------------------------------------
    private:
        /** Causal coefficients that multiply the input data. */
        Scalar myN0 = 1.0;
        Scalar myN1 = 1.0;
        Scalar myN2 = 1.0;
        Scalar myN3 = 1.0;

        /** Recursive coefficients that multiply previously computed values
   * at the output. These are the same for the causal and
   * anti-causal parts of the filter. */
        Scalar myD1 = 0.0;
        Scalar myD2 = 0.0;
        Scalar myD3 = 0.0;
        Scalar myD4 = 0.0;

        /** Anti-causal coefficients that multiply the input data. */
        Scalar myM1 = 0.0;
        Scalar myM2 = 0.0;
        Scalar myM3 = 0.0;
        Scalar myM4 = 0.0;

        /** Recursive coefficients to be used at the boundaries to simulate
   * edge extension boundary conditions. */
        Scalar myBN1 = 0.0;
        Scalar myBN2 = 0.0;
        Scalar myBN3 = 0.0;
        Scalar myBN4 = 0.0;

        Scalar myBM1 = 0.0;
        Scalar myBM2 = 0.0;
        Scalar myBM3 = 0.0;
        Scalar myBM4 = 0.0;

    }; // end of class RecursiveGaussianDerivativeOperator


    /**
     * Overloads 'operator<<' for displaying objects of class 'RecursiveGaussianDerivativeOperator'.
     * @param out the output stream where the object is written.
     * @param object the object of class 'RecursiveGaussianDerivativeOperator' to write.
     * @return the output stream after the writing.
     */
    template <typename TInputImage>
    std::ostream &
    operator<<(std::ostream &out, const RecursiveGaussianDerivativeOperator<TInputImage> &object) {
        object.selfDisplay(out);
        return out;
    }

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
template <typename TInputImage>
void
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>::
selfDisplay(std::ostream &out) const {
    out << "[RecursiveGaussianDerivativeOperator] sigma=" << mySigma << ", order=" << myOrder << ", normalize="
        << (myNormalizeAcrossScale ? "true" : "false");
}

template <typename TInputImage>
bool
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>::
isValid() const {
    return true;
}

template <typename TInputImage>
typename DGtal::RecursiveGaussianDerivativeOperator<TInputImage>::OutputImage
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>::
gaussianImage() {

    size_t width = myImage.domain().upperBound()[0] - myImage.domain().lowerBound()[0] + 1;
    Domain domain = myImage.domain();
    OutputImage out(domain);

    computeCoefficients();

    Scalar *inps = new Scalar[width];
    Scalar *outs = new Scalar[width];
    Scalar *scratch = new Scalar[width];

    Point previous = *domain.begin();
    int i = 0;
    for (auto it = domain.begin(), ite = domain.end(); it != ite; ++it) {
        Point p = *it;
        if (p[0] - previous[0] < 0) {
            filterDataArray(outs, inps, scratch, width);
            Point current = previous;
            for (int j = 0; j < width; j++) {
                current[0] = j;
                out.setValue(current, outs[j]);
            }
            i = 0;
            inps[i++] = myImage(p);
            previous = p;
        } else {
            inps[i++] = myImage(p);
            previous = p;
            continue;
        }
    }
    return out;
}


template <typename TInputImage>
void
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>::
filterDataArray(Scalar *outs, const Scalar *data,
                Scalar *scratch, size_t ln) {

    Scalar *scratch1 = outs;
    Scalar *scratch2 = scratch;
    /**
   * Causal direction pass
   */

    // this value is assumed to exist from the border to infinity.
    const Scalar &outV1 = data[0];

    /**
   * Initialize borders
   */

    linearCombination(scratch1[0], outV1, myN0, outV1, myN1, outV1, myN2, outV1, myN3);
    linearCombination(scratch1[1], data[1], myN0, outV1, myN1, outV1, myN2, outV1, myN3);
    linearCombination(scratch1[2], data[2], myN0, data[1], myN1, outV1, myN2, outV1, myN3);
    linearCombination(scratch1[3], data[3], myN0, data[2], myN1, data[1], myN2, outV1, myN3);

    // note that the outV1 value is multiplied by the Boundary coefficients myBNi
    subLinearCombination(scratch1[0], outV1, myBN1, outV1, myBN2, outV1, myBN3, outV1, myBN4);
    subLinearCombination(scratch1[1], scratch1[0], myD1, outV1, myBN2, outV1, myBN3, outV1, myBN4);
    subLinearCombination(scratch1[2], scratch1[1], myD1, scratch1[0], myD2, outV1, myBN3, outV1, myBN4);
    subLinearCombination(scratch1[3], scratch1[2], myD1, scratch1[1], myD2, scratch1[0], myD3, outV1, myBN4);

    /**
   * Recursively filter the rest
   */
    for (unsigned int i = 4; i < ln; i++) {
        linearCombination(scratch1[i], data[i], myN0, data[i - 1], myN1, data[i - 2], myN2, data[i - 3], myN3);
        subLinearCombination(scratch1[i], scratch1[i - 1], myD1, scratch1[i - 2], myD2, scratch1[i - 3], myD3,
                             scratch1[i - 4], myD4);
    }

    /**
   * Store the causal result: outs = scratch already done via alias
   *
   */

    /**
   * AntiCausal direction pass
   */

    // this value is assumed to exist from the border to infinity.
    const Scalar &outV2 = data[ln - 1];

    /**
   * Initialize borders
   */
    linearCombination(scratch2[ln - 1], outV2, myM1, outV2, myM2, outV2, myM3, outV2, myM4);
    linearCombination(scratch2[ln - 2], data[ln - 1], myM1, outV2, myM2, outV2, myM3, outV2, myM4);
    linearCombination(scratch2[ln - 3], data[ln - 2], myM1, data[ln - 1], myM2, outV2, myM3, outV2, myM4);
    linearCombination(scratch2[ln - 4], data[ln - 3], myM1, data[ln - 2], myM2, data[ln - 1], myM3, outV2, myM4);

    // note that the outV2value is multiplied by the Boundary coefficients myBMi
    subLinearCombination(scratch2[ln - 1], outV2, myBM1, outV2, myBM2, outV2, myBM3, outV2, myBM4);
    subLinearCombination(scratch2[ln - 2], scratch2[ln - 1], myD1, outV2, myBM2, outV2, myBM3, outV2, myBM4);
    subLinearCombination(scratch2[ln - 3], scratch2[ln - 2], myD1, scratch2[ln - 1], myD2, outV2, myBM3, outV2, myBM4);
    subLinearCombination(scratch2[ln - 4], scratch2[ln - 3], myD1, scratch2[ln - 2], myD2, scratch2[ln - 1], myD3,
                         outV2, myBM4);

    /**
   * Recursively filter the rest
   */
    for (unsigned int i = ln - 4; i > 0; i--) {
        linearCombination(scratch2[i - 1], data[i], myM1, data[i + 1], myM2, data[i + 2], myM3, data[i + 3], myM4);
        subLinearCombination(scratch2[i - 1], scratch2[i], myD1, scratch2[i + 1], myD2, scratch2[i + 2], myD3,
                             scratch2[i + 3], myD4);
    }

    /**
   * Roll the antiCausal part into the output
   */
    for (unsigned int i = 0; i < ln; i++) {
        outs[i] += scratch2[i];
    }
}

template <typename TInputImage>
void
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>::
computeCoefficients() {
    /**  Parameters of exponential series. */
    Scalar A1[3];
    Scalar B1[3];
    Scalar W1;
    Scalar L1;
    Scalar A2[3];
    Scalar B2[3];
    Scalar W2;
    Scalar L2;

    Scalar direction = 1.0;


    const Scalar sigmad = mySigma;
    Scalar across_scale_normalization = 1.0;


    A1[0] = static_cast< Scalar >(  1.3530 );
    B1[0] = static_cast< Scalar >(  1.8151 );
    W1 = static_cast< Scalar >(  0.6681 );
    L1 = static_cast< Scalar >( -1.3932 );
    A2[0] = static_cast< Scalar >( -0.3531 );
    B2[0] = static_cast< Scalar >(  0.0902 );
    W2 = static_cast< Scalar >(  2.0787 );
    L2 = static_cast< Scalar >( -1.3732 );

    A1[1] = static_cast< Scalar >( -0.6724 );
    B1[1] = static_cast< Scalar >( -3.4327 );
    A2[1] = static_cast< Scalar >(  0.6724 );
    B2[1] = static_cast< Scalar >(  0.6100 );

    A1[2] = static_cast< Scalar >( -1.3563 );
    B1[2] = static_cast< Scalar >(  5.2318 );
    A2[2] = static_cast< Scalar >(  0.3446 );
    B2[2] = static_cast< Scalar >( -2.2355 );

    Scalar SD, DD, ED;
    computeDCoefficients(sigmad, W1, L1, W2, L2, SD, DD, ED);
    Scalar SN, DN, EN;

    switch (myOrder) {
        case 0: {
            // Approximation of convolution with a gaussian.
            computeNCoefficients(sigmad,
                                 A1[0], B1[0], W1, L1,
                                 A2[0], B2[0], W2, L2,
                                 myN0,
                                 myN1,
                                 myN2,
                                 myN3,
                                 SN, DN, EN);

            Scalar alpha0 = 2 * SN / SD - myN0;
            myN0 *= across_scale_normalization / alpha0;
            myN1 *= across_scale_normalization / alpha0;
            myN2 *= across_scale_normalization / alpha0;
            myN3 *= across_scale_normalization / alpha0;
            const bool symmetric = true;
            computeRemainingCoefficients(symmetric);
            break;
        }
        case 1: {
            if (myNormalizeAcrossScale) {
                across_scale_normalization = mySigma;
            }
            // Approximation of convolution with the first derivative of a  Gaussian
            computeNCoefficients(sigmad,
                                 A1[1], B1[1], W1, L1,
                                 A2[1], B2[1], W2, L2,
                                 myN0, myN1, myN2, myN3,
                                 SN, DN, EN);

            Scalar alpha1 = 2 * (SN * DD - DN * SD) / (SD * SD);
            // If negative spacing, negate the first derivative response.
            alpha1 *= direction;
            myN0 *= across_scale_normalization / alpha1;
            myN1 *= across_scale_normalization / alpha1;
            myN2 *= across_scale_normalization / alpha1;
            myN3 *= across_scale_normalization / alpha1;
            const bool symmetric = false;
            computeRemainingCoefficients(symmetric);
            break;
        }
        case 2: {
            if (myNormalizeAcrossScale) {
                across_scale_normalization = std::sqrt(mySigma);
            }
            // Approximation of convolution with the second derivative of a
            // Gaussian.
            Scalar N0_0, N1_0, N2_0, N3_0;
            Scalar N0_2, N1_2, N2_2, N3_2;
            Scalar SN0, DN0, EN0;
            Scalar SN2, DN2, EN2;
            computeNCoefficients(sigmad,
                                 A1[0], B1[0], W1, L1,
                                 A2[0], B2[0], W2, L2,
                                 N0_0, N1_0, N2_0, N3_0,
                                 SN0, DN0, EN0);
            computeNCoefficients(sigmad,
                                 A1[2], B1[2], W1, L1,
                                 A2[2], B2[2], W2, L2,
                                 N0_2, N1_2, N2_2, N3_2,
                                 SN2, DN2, EN2);

            Scalar beta = -(2 * SN2 - SD * N0_2) / (2 * SN0 - SD * N0_0);
            myN0 = N0_2 + beta * N0_0;
            myN1 = N1_2 + beta * N1_0;
            myN2 = N2_2 + beta * N2_0;
            myN3 = N3_2 + beta * N3_0;
            SN = SN2 + beta * SN0;
            DN = DN2 + beta * DN0;
            EN = EN2 + beta * EN0;

            Scalar alpha2;
            alpha2 = EN * SD * SD - ED * SN * SD - 2 * DN * DD * SD + 2 * DD * DD * SN;
            alpha2 /= SD * SD * SD;

            myN0 *= across_scale_normalization / alpha2;
            myN1 *= across_scale_normalization / alpha2;
            myN2 *= across_scale_normalization / alpha2;
            myN3 *= across_scale_normalization / alpha2;

            const bool symmetric = true;
            computeRemainingCoefficients(symmetric);
            break;
        }
    }
}


/**
 * Compute the N coefficients.
 */
template <typename TInputImage>
void
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>
::computeNCoefficients(Scalar sigmad,
                       Scalar A1, Scalar B1, Scalar W1, Scalar L1,
                       Scalar A2, Scalar B2, Scalar W2, Scalar L2,
                       Scalar &N0, Scalar &N1, Scalar &N2, Scalar &N3,
                       Scalar &SN, Scalar &DN, Scalar &EN) {
    Scalar Sin1 = std::sin(W1 / sigmad);
    Scalar Sin2 = std::sin(W2 / sigmad);
    Scalar Cos1 = std::cos(W1 / sigmad);
    Scalar Cos2 = std::cos(W2 / sigmad);
    Scalar Exp1 = std::exp(L1 / sigmad);
    Scalar Exp2 = std::exp(L2 / sigmad);

    N0 = A1 + A2;
    N1 = Exp2 * (B2 * Sin2 - (A2 + 2 * A1) * Cos2);
    N1 += Exp1 * (B1 * Sin1 - (A1 + 2 * A2) * Cos1);
    N2 = (A1 + A2) * Cos2 * Cos1;
    N2 -= B1 * Cos2 * Sin1 + B2 * Cos1 * Sin2;
    N2 *= 2 * Exp1 * Exp2;
    N2 += A2 * Exp1 * Exp1 + A1 * Exp2 * Exp2;
    N3 = Exp2 * Exp1 * Exp1 * (B2 * Sin2 - A2 * Cos2);
    N3 += Exp1 * Exp2 * Exp2 * (B1 * Sin1 - A1 * Cos1);

    SN = N0 + N1 + N2 + N3;
    DN = N1 + 2 * N2 + 3 * N3;
    EN = N1 + 4 * N2 + 9 * N3;
}

/**
 * Compute the D coefficients.
 */
template <typename TInputImage>
void
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>
::computeDCoefficients(Scalar sigmad,
                       Scalar W1, Scalar L1, Scalar W2, Scalar L2,
                       Scalar &SD, Scalar &DD, Scalar &ED) {
    //  const Scalar Sin1 = std::sin(W1 / sigmad);
    //  const Scalar Sin2 = std::sin(W2 / sigmad);
    const Scalar Cos1 = std::cos(W1 / sigmad);
    const Scalar Cos2 = std::cos(W2 / sigmad);
    const Scalar Exp1 = std::exp(L1 / sigmad);
    const Scalar Exp2 = std::exp(L2 / sigmad);

    myD4 = Exp1 * Exp1 * Exp2 * Exp2;
    myD3 = -2 * Cos1 * Exp1 * Exp2 * Exp2;
    myD3 += -2 * Cos2 * Exp2 * Exp1 * Exp1;
    myD2 = 4 * Cos2 * Cos1 * Exp1 * Exp2;
    myD2 += Exp1 * Exp1 + Exp2 * Exp2;
    myD1 = -2 * (Exp2 * Cos2 + Exp1 * Cos1);

    SD = 1.0 + myD1 + myD2 + myD3 + myD4;
    DD = myD1 + 2 * myD2 + 3 * myD3 + 4 * myD4;
    ED = myD1 + 4 * myD2 + 9 * myD3 + 16 * myD4;
}

/**
 * Compute the M coefficients and the boundary coefficients.
 */
template <typename TInputImage>
void
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>
::computeRemainingCoefficients(bool symmetric) {
    if (symmetric) {
        myM1 = myN1 - myD1 * myN0;
        myM2 = myN2 - myD2 * myN0;
        myM3 = myN3 - myD3 * myN0;
        myM4 = -myD4 * myN0;
    } else {
        myM1 = -(myN1 - myD1 * myN0);
        myM2 = -(myN2 - myD2 * myN0);
        myM3 = -(myN3 - myD3 * myN0);
        myM4 = myD4 * myN0;
    }

    // Compute coefficients to be used at the boundaries
    // in order to simulate edge extension boundary conditions.
    const Scalar SN = myN0 + myN1 + myN2 + myN3;
    const Scalar SM = myM1 + myM2 + myM3 + myM4;
    const Scalar SD = 1.0 + myD1 + myD2 + myD3 + myD4;

    myBN1 = myD1 * SN / SD;
    myBN2 = myD2 * SN / SD;
    myBN3 = myD3 * SN / SD;
    myBN4 = myD4 * SN / SD;

    myBM1 = myD1 * SM / SD;
    myBM2 = myD2 * SM / SD;
    myBM3 = myD3 * SM / SD;
    myBM4 = myD4 * SM / SD;
}

template <typename TInputImage>
void
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>::
setSigma(double sigma) {
    mySigma = sigma;
}

template <typename TInputImage>
double
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>::
getSigma() const {
    return mySigma;
}

template <typename TInputImage>
void
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>::
setOrder(unsigned int anOrder) {
    myOrder = anOrder;
}

template <typename TInputImage>
void
DGtal::RecursiveGaussianDerivativeOperator<TInputImage>::
setNormalizeAcrossScale(bool normalize) {
    myNormalizeAcrossScale = normalize;
}


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined RecursiveRecursiveGaussianDerivativeOperator_h

#undef RecursiveRecursiveGaussianDerivativeOperator_RECURSES
#endif // else defined(RecursiveRecursiveGaussianDerivativeOperator_RECURSES)