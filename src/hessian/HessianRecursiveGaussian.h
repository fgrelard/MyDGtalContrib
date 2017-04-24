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
 * @file HessianRecursiveGaussian.h
 * @author Florent Grelard (\c florent.grelard@u-bordeaux.fr )
 * LaBRI, Bordeaux University
 *
 * @date 2017/04/05
 *
 *
 * This file is part of the DGtal library.
 */

#if defined(HessianRecursiveGaussian_RECURSES)
#error Recursive header files inclusion detected in HessianRecursiveGaussian.h
#else // defined(HessianRecursiveGaussian_RECURSES)
/** Prevents recursive inclusion of headers. */
#define HessianRecursiveGaussian_RECURSES

#if !defined HessianRecursiveGaussian_h
/** Prevents repeated inclusion of headers. */
#define HessianRecursiveGaussian_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <DGtal/math/linalg/SimpleMatrix.h>
#include "RecursiveGaussianDerivativeOperator.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal {

    /////////////////////////////////////////////////////////////////////////////
    // template class HessianRecursiveGaussian
    /**
     * Description of template class 'HessianRecursiveGaussian' <p>
     * \brief Aim:
     */
    template <typename TInputImage>
    class HessianRecursiveGaussian {
    public:
        typedef TInputImage InputImage;
        typedef typename InputImage::Domain Domain;
        typedef typename Domain::Point Point;
        typedef double Scalar;
        typedef RecursiveGaussianDerivativeOperator<InputImage> DerivativeFilterA;
        typedef typename DerivativeFilterA::OutputImage SmoothedImage;
        typedef RecursiveGaussianDerivativeOperator<SmoothedImage> GaussianFilter;
        typedef RecursiveGaussianDerivativeOperator<SmoothedImage> DerivativeFilterB;
        typedef std::vector<GaussianFilter> GaussianFilterArray;
        typedef DGtal::SimpleMatrix<Scalar, Domain::dimension,
                Domain::dimension> HessianMatrix;
        typedef ImageContainerBySTLVector<Domain, HessianMatrix> OutputImage;

        // ----------------------- Standard services ------------------------------
    public:

        /**
         * Default constructor.
         */
        HessianRecursiveGaussian() : myImage(Domain(Point::zero, Point::zero)), mySigma(1.0),
                                     myNormalizeAcrossScale(false) {
            initialize();
        }

        HessianRecursiveGaussian(const InputImage &anImage) : myImage(anImage), mySigma(1.0),
                                                              myNormalizeAcrossScale(true) {
            initialize();
        }

        HessianRecursiveGaussian(const InputImage &anImage, double sigma, bool normalize = true) : myImage(anImage),
                                                                                                   mySigma(sigma),
                                                                                                   myNormalizeAcrossScale(
                                                                                                           normalize) {
            initialize();
        }

        /**
         * Destructor.
         */
        ~HessianRecursiveGaussian() = default;

        /**
         * Copy constructor.
         * @param other the object to clone.
         */
        HessianRecursiveGaussian(const HessianRecursiveGaussian &other) : myImage(other.myImage),
                                                                          mySigma(other.mySigma),
                                                                          myNormalizeAcrossScale(
                                                                                  other.myNormalizeAcrossScale) {
            initialize();
        }

        /**
         * Move constructor.
         * @param other the object to move.
         */
        HessianRecursiveGaussian(HessianRecursiveGaussian &&other) = delete;

        /**
         * Copy assignment operator.
         * @param other the object to copy.
         * @return a reference on 'this'.
         */
        HessianRecursiveGaussian &operator=(const HessianRecursiveGaussian &other) = delete;

        /**
         * Move assignment operator.
         * @param other the object to move.
         * @return a reference on 'this'.
         */
        HessianRecursiveGaussian &operator=(HessianRecursiveGaussian &&other) = delete;

        // ----------------------- Interface --------------------------------------
    public:
        OutputImage computeHessian();

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

        void setNormalizeAcrossScale(bool normalize);


    public:
        const unsigned int NUMBER_OF_FILTERS = (Domain::dimension > 2 ? Domain::dimension - 2 : 0);

        // ------------------------- Protected Datas ------------------------------
    protected:

        // ------------------------- Private Datas --------------------------------
    private:
        InputImage myImage;
        double mySigma;
        bool myNormalizeAcrossScale;

        // ------------------------- Hidden services ------------------------------
    protected:
        void initialize();

        // ------------------------- Internals ------------------------------------
    private:
        GaussianFilterArray mySmoothingFilters;
        DerivativeFilterA myDerivativeFilterA;
        DerivativeFilterB myDerivativeFilterB;

    }; // end of class HessianRecursiveGaussian


    /**
     * Overloads 'operator<<' for displaying objects of class 'HessianRecursiveGaussian'.
     * @param out the output stream where the object is written.
     * @param object the object of class 'HessianRecursiveGaussian' to write.
     * @return the output stream after the writing.
     */
    template <typename TInputImage>
    std::ostream &
    operator<<(std::ostream &out, const HessianRecursiveGaussian<TInputImage> &object) {
        object.selfDisplay(out);
        return out;
    }


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
template <typename TInputImage>
void
DGtal::HessianRecursiveGaussian<TInputImage>::
selfDisplay(std::ostream &out) const {
    out << "[HessianRecursiveGaussian] sigma=" << mySigma << ", normalize="
        << (myNormalizeAcrossScale ? "true" : "false");
}

template <typename TInputImage>
bool
DGtal::HessianRecursiveGaussian<TInputImage>::
isValid() const {
    return true;
}

template <typename TInputImage>
typename DGtal::HessianRecursiveGaussian<TInputImage>::OutputImage
DGtal::HessianRecursiveGaussian<TInputImage>::
computeHessian() {

    OutputImage out(myImage.domain());

    unsigned int number = NUMBER_OF_FILTERS;

    for (unsigned int dima = 0; dima < Domain::dimension; dima++) {
        for (unsigned int dimb = dima; dimb < Domain::dimension; dimb++) {
            // Manage the diagonal in a different way in order to avoid
            // applying a double smoothing to this direction, and missing
            // to smooth one of the other directions.
            bool inPlace = false;
            if (dimb == dima) {
                myDerivativeFilterA.setOrder(2);
                myDerivativeFilterB.setOrder(0);

                inPlace = true;

                unsigned int i = 0;
                unsigned int j = 0;
                // find the direction for the first filter.
                while (j < Domain::dimension) {
                    if (j != dima) {
                        myDerivativeFilterB.setDirection(j);
                        j++;
                        break;
                    }
                    j++;
                }
                // find the direction for all the other filters
                while (i < number) {
                    while (j < Domain::dimension) {
                        if (j != dima) {
                            mySmoothingFilters[i].setDirection(j);
                            j++;
                            break;
                        }
                        j++;
                    }
                    i++;
                }

                myDerivativeFilterA.setDirection(dima);
            } else {
                myDerivativeFilterA.setOrder(1);
                myDerivativeFilterB.setOrder(1);

                if (dimb < Domain::dimension - 1) {
                    // can reuse the output of m_DerivativeFilterA
                    inPlace = false;
                } else {
                    inPlace = true;
                }

                unsigned int i = 0;
                unsigned int j = 0;
                while (i < number) {
                    while (j < Domain::dimension) {
                        if (j != dima && j != dimb) {
                            mySmoothingFilters[i].setDirection(j);
                            j++;
                            break;
                        }
                        j++;
                    }
                    i++;
                }

                myDerivativeFilterA.setDirection(dima);
                myDerivativeFilterB.setDirection(dimb);
            }

            SmoothedImage derivativeImage(myImage.domain());

            // Deal with the 2D case.
            if (number > 0) {
                int temp_dim = static_cast< int >( Domain::dimension ) - 3;
                myDerivativeFilterB.setImage(myDerivativeFilterA.gaussianImage());
                GaussianFilter lastFilter = mySmoothingFilters[temp_dim];
                lastFilter.setImage(myDerivativeFilterB.gaussianImage());
                derivativeImage = lastFilter.gaussianImage();
            } else {
                if (inPlace) {
                    myDerivativeFilterB.setImage(myDerivativeFilterA.gaussianImage());
                }
                derivativeImage = myDerivativeFilterB.gaussianImage();
            }
            for (auto it = derivativeImage.domain().begin(), ite = derivativeImage.domain().end();
                 it != ite; ++it) {
                Point p = *it;
                HessianMatrix hessian = out(p);
                hessian(dima, dimb) = derivativeImage(p);
                hessian(dimb, dima) = derivativeImage(p);
                out.setValue(p, hessian);
            }
        }
    }
    return out;

}

template <typename TInputImage>
void
DGtal::HessianRecursiveGaussian<TInputImage>::
initialize() {
    unsigned int number = NUMBER_OF_FILTERS;
    mySmoothingFilters.resize(number);

    myDerivativeFilterA = DerivativeFilterA(myImage, mySigma, 1, myNormalizeAcrossScale);
    SmoothedImage smooth = myDerivativeFilterA.gaussianImage();

    myDerivativeFilterB = DerivativeFilterB(smooth, mySigma, 1, myNormalizeAcrossScale);

    if (number > 0) {
        mySmoothingFilters[0] = GaussianFilter(myDerivativeFilterB.gaussianImage(), mySigma, 0, false);
    }
    for (unsigned int i = 1; i < number; i++) {
        mySmoothingFilters[i] = GaussianFilter(mySmoothingFilters[i - 1].gaussianImage(), mySigma, 0, false);
    }
}

template <typename TInputImage>
void
DGtal::HessianRecursiveGaussian<TInputImage>::
setSigma(double sigma) {
    mySigma = sigma;
    unsigned int number = NUMBER_OF_FILTERS;

    for (unsigned int i = 0; i < number; i++) {
        mySmoothingFilters[i].setSigma(sigma);
    }
    myDerivativeFilterA.setSigma(sigma);
    myDerivativeFilterA.setSigma(sigma);
}

template <typename TInputImage>
double
DGtal::HessianRecursiveGaussian<TInputImage>::
getSigma() const {
    return mySigma;
}

template <typename TInputImage>
void
DGtal::HessianRecursiveGaussian<TInputImage>::
setNormalizeAcrossScale(bool normalize) {
    myNormalizeAcrossScale = normalize;
    unsigned int number = NUMBER_OF_FILTERS;


    myDerivativeFilterA.setNormalizeAcrossScale(normalize);
    myDerivativeFilterB.setNormalizeAcrossScale(normalize);
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined HessianRecursiveGaussian_h

#undef HessianRecursiveGaussian_RECURSES
#endif // else defined(HessianRecursiveGaussian_RECURSES)
