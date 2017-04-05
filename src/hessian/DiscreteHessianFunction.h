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
#include <DGtal/math/linalg/SimpleMatrix.h>
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
        typedef DGtal::SimpleMatrix<HessianValue, Domain::dimension,
                Domain::dimension> HessianMatrix;
//        typedef std::vector<HessianValue> HessianMatrix;

        typedef GaussianDerivativeOperator<HessianValue> GaussianDerivative;
        typedef typename GaussianDerivative::CoefficientVector CoefficientVector;

        typedef ImageContainerBySTLVector<Domain, HessianValue> KernelImage;
        typedef ImageContainerBySTLVector<Domain, HessianMatrix> OutputImage;


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

        HessianMatrix hessianAtPoint(const Point &p);


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

        // ------------------------- Hidden services ------------------------------
    protected:
        void computeKernels();

        HessianMatrix convertSymmetricToSquareMatrix(const std::vector<HessianValue> &values);

        KernelImage uniqueKernelFromMultipleDirections(const std::vector<CoefficientVector> &coeff, int size);

        HessianValue convolve(const Point &currentPoint,
                              const Domain &region,
                              const KernelImage &kernel);

        // ------------------------- Protected Datas ------------------------------
    protected:

        // ------------------------- Private Datas --------------------------------
    private:
        Image myImage;
        double mySigma;
        bool myNormalizeAcrossScale;


        // ------------------------- Internals ------------------------------------
    private:
        std::vector<KernelImage> myKernels;


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

    /*
     * Workaround for DGtal simple matrix
     */
    template <typename Scalar, DGtal::Dimension TM, DGtal::Dimension TN>
    bool operator!=(const DGtal::SimpleMatrix<Scalar, TM,
            TN> &m1,
                    const DGtal::SimpleMatrix<Scalar, TM, TN> &m2) { return !(m1 == m2); }

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
template <typename TImage>
bool
DGtal::DiscreteHessianFunction<TImage>::
isValid() const {
    return true;
}

template <typename TImage>
typename DGtal::DiscreteHessianFunction<TImage>::OutputImage
DGtal::DiscreteHessianFunction<TImage>::
computeHessian() {
    OutputImage out(myImage.domain());
    unsigned int size = out.domain().size();
    unsigned int i = 0;

#pragma omp parallel
#pragma omp single
    {
        for (auto it = out.domain().begin(), ite = out.domain().end(); it != ite; ++it) {
#pragma omp task firstprivate(it)
            {
                HessianMatrix value = hessianAtPoint(*it);
                out.setValue(*it, value);
            }
#pragma omp taskwait
        }
    }

    return out;
}

template <typename TImage>
typename DGtal::DiscreteHessianFunction<TImage>::HessianMatrix
DGtal::DiscreteHessianFunction<TImage>::
hessianAtPoint(const Point &p) {
    if (myKernels.size() == 0)
        computeKernels();

    HessianMatrix hessian;
    if (!myImage.domain().isInside(p))
        return hessian;

    std::vector<HessianValue> symmetricMatrix(Domain::dimension * (Domain::dimension + 1) / 2);
    for (unsigned int i = 0; i < myKernels.size(); i++) {
        KernelImage currentKernel = myKernels[i];
        Domain domain(p + currentKernel.domain().lowerBound(),
                      p + currentKernel.domain().upperBound());
        HessianValue value = convolve(p, domain, currentKernel);
        symmetricMatrix[i] = value;
    }

    hessian = convertSymmetricToSquareMatrix(symmetricMatrix);

    return hessian;
}

template <typename TImage>
typename DGtal::DiscreteHessianFunction<TImage>::HessianMatrix
DGtal::DiscreteHessianFunction<TImage>::
convertSymmetricToSquareMatrix(const std::vector<HessianValue> &values) {
    HessianMatrix m;
    for (int i = 0; i < Domain::dimension; i++) {
        for (int j = i; j < Domain::dimension; j++) {
            unsigned int index = j + (Domain::dimension - 1) * i - i * (i - 1) / 2;
            HessianValue value = values[index];
            m.setComponent(i, j, value);
            m.setComponent(j, i, value);
        }
    }
    return m;
}

template <typename TImage>
void
DGtal::DiscreteHessianFunction<TImage>::
computeKernels() {
    /* Create 3*N operators (N=ImageDimension) where the
 * first N are zero-order, the second N are first-order
 * and the third N are second order */
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
            operators[idx] = coeff;
            if (coeff.size() > maxRadius)
                maxRadius = coeff.size();
        }
    }

// Now precompute the N-dimensional kernel. This fastest as we don't
// have to perform N convolutions for each point we calculate but
    // only one.

    Domain kernelDomain(Point::zero - Point::diagonal(2 * maxRadius),
                        Point::zero + Point::diagonal(2 * maxRadius));



    // Precalculate compound derivative kernels (n-dimensional)
    // The order of calculation in the 3D case is: dxx, dxy, dxz, dyy,
    // dyz, dzz
    unsigned int opidx; // current operator index in operators
    unsigned int kernelidx = 0;

    for (unsigned int i = 0; i < dimension; ++i) {
        for (unsigned int j = i; j < dimension; ++j) {
            OrderArrayType orderArray;
            for (unsigned int k = 0; k < dimension; k++) {
                orderArray[k] = 0;
            }
            ++orderArray[i];
            ++orderArray[j];


            std::vector<CoefficientVector> coeffs;
            int maxIndex = 0;
            for (unsigned int direction = 0; direction < dimension; ++direction) {
                opidx = dimension * orderArray[direction] + direction;
                CoefficientVector coeff = operators[opidx];
                int indexFirst = coeff.size() / 2;
                if (indexFirst > maxIndex) {
                    maxIndex = indexFirst;
                }
                coeffs.push_back(coeff);
            }

            KernelImage kernelImage = uniqueKernelFromMultipleDirections(coeffs, maxIndex);
            myKernels.push_back(kernelImage);
        }
    }
}


template <typename TImage>
typename DGtal::DiscreteHessianFunction<TImage>::HessianValue
DGtal::DiscreteHessianFunction<TImage>::
convolve(const Point &currentPoint, const Domain &region, const KernelImage &kernel) {
    HessianValue innerProduct = DGtal::NumberTraits<HessianValue>::ZERO;
    for (const Point &p : region) {
        if (!myImage.domain().isInside(p)) continue;
        innerProduct += myImage(p) * kernel(p - currentPoint);
    }
    return innerProduct;
}


template <typename TImage>
typename DGtal::DiscreteHessianFunction<TImage>::KernelImage
DGtal::DiscreteHessianFunction<TImage>::
uniqueKernelFromMultipleDirections(const std::vector<CoefficientVector> &coeff, int radius) {

    Domain domain(Point::zero - Point::diagonal(radius),
                  Point::zero + Point::diagonal(radius));
    KernelImage out(domain);

    int previous = 1000;
    for (const Point &p : domain) {
        double value = 1.0;
        for (unsigned int i = 0; i < Domain::dimension; i++) {
            //flip axis
            value *= coeff[i][2 * radius - (p[i] + radius)];
        }
        if (std::abs(value - 0.0) < std::numeric_limits<double>::epsilon())
            value = 0.0;
        out.setValue(p, value);
    }
    return out;

}

template <typename TImage>
void
DGtal::DiscreteHessianFunction<TImage>::
setSigma(double sigma) {
    mySigma = sigma;
    myKernels.clear();
    computeKernels();
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
    myKernels.clear();
    computeKernels();
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