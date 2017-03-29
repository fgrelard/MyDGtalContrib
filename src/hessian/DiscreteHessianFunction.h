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
 * @file DiscreteHessianFunction.h
 * @author Florent Grelard (\c florent.grelard@u-bordeaux.fr )
 * LaBRI, Bordeaux University
 *
 * @date 2017/03/29
 *
 *
 * This file is part of the DGtal library.
 */

#if defined(DiscreteHessianFunction_RECURSES)
#error Recursive header files inclusion detected in DiscreteHessianFunction.h
#else // defined(DiscreteHessianFunction_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DiscreteHessianFunction_RECURSES

#if !defined DiscreteHessianFunction_h
/** Prevents repeated inclusion of headers. */
#define DiscreteHessianFunction_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include "GaussianDerivativeOperator.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal {

    /////////////////////////////////////////////////////////////////////////////
    // template class DiscreteHessianFunction
    /**
     * Description of template class 'DiscreteHessianFunction' <p>
     * \brief Aim:
     */
    template <typename TImage>
    class DiscreteHessianFunction {
    public:
        typedef TImage Image;
        typedef typename Image::Domain Domain;
        typedef typename Image::Dimension Dimension;
        typedef typename Image::Point Point;
        typedef double HessianValue;
        typedef GaussianDerivativeOperator<HessianValue> GaussianDerivative;
        typedef ImageContainerBySTLVector<HessianValue, Domain::dimension> OutputImage;

    public:
        static Domain emptyDomain(Point::zero, Point::zero);

        // ----------------------- Standard services ------------------------------
    public:

        /**
         * Default constructor.
         */
        DiscreteHessianFunction() = delete;

        DiscreteHessianFunction(const Image &image) : myImage(image), mySigma(1.0), myNormalizeAcrossScale(true) {}

        DiscreteHessianFunction(const Image &image, double sigma, bool normalize = true) : myImage(image),
                                                                                           mySigma(sigma),
                                                                                           myNormalizeAcrossScale(
                                                                                                   normalize) {}


        /**
         * Destructor.
         */
        ~DiscreteHessianFunction() = default;

        /**
         * Copy constructor.
         * @param other the object to clone.
         */
        DiscreteHessianFunction(const DiscreteHessianFunction &other) = default;

        /**
         * Move constructor.
         * @param other the object to move.
         */
        DiscreteHessianFunction(DiscreteHessianFunction &&other) = delete;

        /**
         * Copy assignment operator.
         * @param other the object to copy.
         * @return a reference on 'this'.
         */
        DiscreteHessianFunction &operator=(const DiscreteHessianFunction &other) = delete;

        /**
         * Move assignment operator.
         * @param other the object to move.
         * @return a reference on 'this'.
         */
        DiscreteHessianFunction &operator=(DiscreteHessianFunction &&other) = delete;

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

        bool getNormalizeAcrossScale() const;

        // ------------------------- Protected Datas ------------------------------
    protected:

        // ------------------------- Private Datas --------------------------------
    private:
        Image myImage;
        double mySigma;
        bool myNormalizeAcrossScale;



        // ------------------------- Hidden services ------------------------------
    protected:

        // ------------------------- Internals ------------------------------------
    private:

    }; // end of class DiscreteHessianFunction


    /**
     * Overloads 'operator<<' for displaying objects of class 'DiscreteHessianFunction'.
     * @param out the output stream where the object is written.
     * @param object the object of class 'DiscreteHessianFunction' to write.
     * @return the output stream after the writing.
     */
    template <typename T>
    std::ostream &
    operator<<(std::ostream &out, const DiscreteHessianFunction<T> &object) {
        object.selfDisplay(out);
        return out;
    }

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
template <typename TImage>
void
DGtal::DiscreteHessianFunction<TImage>::
selfDisplay(std::ostream &out) const {
    out << "[DiscreteHessianFunction] image=" << myImage << ", sigma=" << mySigma << ", normalize="
        << (myNormalizeAcrossScale ? "true" : "false");
}

/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
template <typename TImagiane>
bool
DGtal::DiscreteHessianFunction<TImage>::
isValid() const {
    return true;
}

template <typename TImage>
typename DGtal::DiscreteHessianFunction<TImage>::OutputImage
DGtal::DiscreteHessianFunction<TImage>::
computeHessian() {
    /* Create 3*N operators (N=ImageDimension) where the
 * first N are zero-order, the second N are first-order
 * and the third N are second order */
    typedef typename GaussianDerivative::CoefficientVector CoefficientVector;
    Dimension dimension = Domain::dimension;
    typedef unsigned int OrderArrayType[dimension];

    unsigned int idx;
    unsigned int maxRadius = 0;

    std::vector<CoefficientVector> operators(3 * dimension);
    for (unsigned int direction = 0; direction < dimension; direction++) {
        for (unsigned int order = 0; order <= 2; ++order) {
            idx = dimension * order + direction;
            GaussianDerivative g(mySigma, order, myNormalizeAcrossScale);
            CoefficientVector coeff = g.computeCoefficients();
            operators.push_back(coeff);
            if (coeff.size() > maxRadius)
                maxRadius = coeff.size();
        }
    }


// Now precompute the N-dimensional kernel. This fastest as we don't
// have to perform N convolutions for each point we calculate but
    // only one.

    Domain kernelDomain(Point::zero - Point::diagonal(2 * maxRadius),
                        Point::zero + Point::diagonal(2 * maxRadius));
    OutputImage kernelImage(kernelDomain);
    for (const Point &p : kernelDomain)
        kernelImage.setValue(p, DGtal::NumberTraits<HessianValue>::ZERO);


    OutputImage zeroKernel = kernelImage;
    Domain kernelRegion(Point::zero - Point::diagonal(maxRadius),
                        Point::zero + Point::diagonal(maxRadius));


    // Precalculate compound derivative kernels (n-dimensional)
    // The order of calculation in the 3D case is: dxx, dxy, dxz, dyy,
    // dyz, dzz
    unsigned int opidx; // current operator index in operators
    unsigned int kernelidx = 0;

    for (unsigned int i = 0; i < dimension; ++i) {
        for (unsigned int j = i; j < dimension; ++j) {
            OrderArrayType orderArray;
            for (unsigned int k = 0; k < dimension; k++) {
                orderArray[i] = 0;
            }
            ++orderArray[i];
            ++orderArray[j];

            // Reset kernel image
            kernelImage = zeroKernel;
            kernelImage.setValue(Point::zero, DGtal::NumberTraits<HessianValue>::ONE);

            for (unsigned int direction = 0; direction < dimension; ++direction) {
                opidx = dimension * orderArray[direction] + direction;
                convolutionFilter->SetInput(kernelImage);
                convolutionFilter->SetOperator(operators[opidx]);
                convolutionFilter->Update();
                kernelImage = convolutionFilter->GetOutput();

            }

            // Set the size of the current kernel
            m_KernelArray[kernelidx].SetRadius(maxRadius);

            // Copy kernel image to neighborhood. Do not copy boundaries.
            ImageRegionConstIterator <KernelImageType> it(kernelImage, kernelRegion);
            it.GoToBegin();
            idx = 0;

            while (!it.IsAtEnd()) {
                m_KernelArray[kernelidx][idx] = it.Get();
                ++idx;
                ++it;
            }
            kernelidx++;
        }
    }
}

}

template <typename TImage>
void
DGtal::DiscreteHessianFunction<TImage>::
setSigma(double sigma) {
    mySigma = sigma;
}

template <typename TImage>
double
DGtal::DiscreteHessianFunction<TImage>::
getSigma() const {
    return mySigma;
}

template <typename TImage>
void
DGtal::DiscreteHessianFunction<TImage>::
setNormalizeAcrossScale(bool normalize) {
    myNormalizeAcrossScale = normalize;
}

template <typename TImage>
bool
DGtal::DiscreteHessianFunction<TImage>::
getNormalizeAcrossScale() const {
    return myNormalizeAcrossScale;
}

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DiscreteHessianFunction_h

#undef DiscreteHessianFunction_RECURSES
#endif // else defined(DiscreteHessianFunction_RECURSES)