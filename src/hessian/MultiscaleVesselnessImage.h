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
 * @file MultiscaleVesselness.h
 * @author Florent Grelard (\c florent.grelard@u-bordeaux.fr )
 * LaBRI, Bordeaux University
 *
 * @date 2017/04/05
 *
 *
 * This file is part of the DGtal library.
 */

#if defined(MultiscaleVesselness_RECURSES)
#error Recursive header files inclusion detected in MultiscaleVesselness.h
#else // defined(MultiscaleVesselness_RECURSES)
/** Prevents recursive inclusion of headers. */
#define MultiscaleVesselness_RECURSES

#if !defined MultiscaleVesselness_h
/** Prevents repeated inclusion of headers. */
#define MultiscaleVesselness_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
//////////////////////////////////////////////////////////////////////////////

namespace DGtal {

    /////////////////////////////////////////////////////////////////////////////
    // template class MultiscaleVesselness
    /**
     * Description of template class 'MultiscaleVesselness' <p>
     * \brief Aim:
     */
    template <typename THessianMeasure, typename TVesselnessMeasure>
    class MultiscaleVesselness {
    public:
        typedef THessianMeasure HessianMeasure;
        typedef TVesselnessMeasure VesselnessMeasure;
        typedef typename VesselnessMeasure::HessianImage HessianImage;
        typedef typename HessianImage::Domain Domain;
        typedef typename Domain::Point Point;
        typedef typename VesselnessMeasure::OutputImage OutputImage;
        typedef typename OutputImage::Value Value;

        // ----------------------- Standard services ------------------------------
    public:

        /**
         * Default constructor.
         */
        MultiscaleVesselness() : mySigmaMin(1.0), mySigmaMax(1.0), myNumberSigmaSteps(1) {}

        MultiscaleVesselness(const HessianMeasure &hessianImage, const VesselnessMeasure &measure) : myHessianMeasure(
                hessianImage), myMeasure(measure), mySigmaMin(1.0), mySigmaMax(1.0), myNumberSigmaSteps(1) {}

        MultiscaleVesselness(const HessianMeasure &hessianImage, const VesselnessMeasure &measure, double sigmaMin,
                             double sigmaMax, unsigned int numberSigma) : myHessianMeasure(hessianImage),
                                                                          myMeasure(measure), mySigmaMin(sigmaMin),
                                                                          mySigmaMax(sigmaMax),
                                                                          myNumberSigmaSteps(numberSigma) {}

        /**
         * Destructor.
         */
        ~MultiscaleVesselness() = default;

        /**
         * Copy constructor.
         * @param other the object to clone.
         */
        MultiscaleVesselness(const MultiscaleVesselness &other) : myMeasure(other.myMeasure),
                                                                  mySigmaMin(other.mySigmaMin),
                                                                  mySigmaMax(other.mySigmaMax),
                                                                  myNumberSigmaSteps(other.myNumberSigmaSteps) {}

        /**
         * Move constructor.
         * @param other the object to move.
         */
        MultiscaleVesselness(MultiscaleVesselness &&other) = delete;

        /**
         * Copy assignment operator.
         * @param other the object to copy.
         * @return a reference on 'this'.
         */
        MultiscaleVesselness &operator=(const MultiscaleVesselness &other) = delete;

        /**
         * Move assignment operator.
         * @param other the object to move.
         * @return a reference on 'this'.
         */
        MultiscaleVesselness &operator=(MultiscaleVesselness &&other) = delete;

        // ----------------------- Interface --------------------------------------
    public:
        OutputImage computeMultiscaleVesselness();

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

        void setSigmaMin(double sigmaMin);

        double getSigmaMin() const;

        void setSigmaMax(double sigmaMax);

        double getSigmaMax() const;

        void setNumberSigmaSteps(double number);

        double getNumberSigmaSteps() const;

        // ------------------------- Protected Datas ------------------------------
    protected:

        // ------------------------- Private Datas --------------------------------
    private:
        HessianMeasure myHessianMeasure;
        VesselnessMeasure myMeasure;
        double mySigmaMin;
        double mySigmaMax;
        unsigned int myNumberSigmaSteps;


        // ------------------------- Hidden services ------------------------------
    protected:

        // ------------------------- Internals ------------------------------------
    private:

    }; // end of class MultiscaleVesselness


    /**
     * Overloads 'operator<<' for displaying objects of class 'MultiscaleVesselness'.
     * @param out the output stream where the object is written.
     * @param object the object of class 'MultiscaleVesselness' to write.
     * @return the output stream after the writing.
     */
    template <typename THessianMeasure, typename TVesselnessMeasure>
    std::ostream &
    operator<<(std::ostream &out, const MultiscaleVesselness<THessianMeasure, TVesselnessMeasure> &object) {
        object.selfDisplay(out);
        return out;
    }

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
template <typename THessianMeasure, typename TVesselnessMeasure>
void
DGtal::MultiscaleVesselness<THessianMeasure, TVesselnessMeasure>::
selfDisplay(std::ostream &out) const {
    out << "[MultiscaleVesselness] method=" << myMeasure << ", sigmaMin=" << mySigmaMin << ", sigmaMax=" << mySigmaMax
        << ", numberSigmaSteps="
        << myNumberSigmaSteps;
}

template <typename THessianMeasure, typename TVesselnessMeasure>
typename DGtal::MultiscaleVesselness<THessianMeasure, TVesselnessMeasure>::OutputImage
DGtal::MultiscaleVesselness<THessianMeasure, TVesselnessMeasure>::
computeMultiscaleVesselness() {
    myHessianMeasure.setSigma(mySigmaMin);
    HessianImage currentHessian = myHessianMeasure.computeHessian();

    OutputImage vesselnessImage(currentHessian.domain());
    for (const Point&  p : vesselnessImage.domain()) {
        vesselnessImage.setValue(p, DGtal::NumberTraits<Value>::ZERO);
    }

    myMeasure.setImage(currentHessian);
    OutputImage currentVesselnessImage = myMeasure.computeVesselness();
    for (const Point &p : currentVesselnessImage.domain()) {
        Value currentValue = currentVesselnessImage(p);
        if (currentValue > vesselnessImage(p)) {
            vesselnessImage.setValue(p, currentValue);
        }
    }
    for (int i = 1; i < myNumberSigmaSteps; i++) {
        double stepSize = (mySigmaMax - mySigmaMin) / (myNumberSigmaSteps - 1);
        double currentSigma = mySigmaMin + stepSize * i;
        myHessianMeasure.setSigma(currentSigma);
        currentHessian = myHessianMeasure.computeHessian();
        myMeasure.setImage(currentHessian);
        OutputImage currentVesselnessImage = myMeasure.computeVesselness();
        for (const Point &p : currentVesselnessImage.domain()) {
            Value currentValue = currentVesselnessImage(p);
            if (currentValue > vesselnessImage(p)) {
                vesselnessImage.setValue(p, currentValue);
            }
        }
    }
    return vesselnessImage;
}


/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
template <typename THessianMeasure, typename TVesselnessMeasure>
bool
DGtal::MultiscaleVesselness<THessianMeasure, TVesselnessMeasure>::
isValid() const {
    return true;
}

template <typename THessianMeasure, typename TVesselnessMeasure>
void
DGtal::MultiscaleVesselness<THessianMeasure, TVesselnessMeasure>::
setSigmaMin(double sigmaMin) {
    mySigmaMin = sigmaMin;
}

template <typename THessianMeasure, typename TVesselnessMeasure>
double
DGtal::MultiscaleVesselness<THessianMeasure, TVesselnessMeasure>::
getSigmaMin() const {
    return mySigmaMin;
}

template <typename THessianMeasure, typename TVesselnessMeasure>
void
DGtal::MultiscaleVesselness<THessianMeasure, TVesselnessMeasure>::
setSigmaMax(double sigmaMax) {
    mySigmaMax = sigmaMax;
}

template <typename THessianMeasure, typename TVesselnessMeasure>
double
DGtal::MultiscaleVesselness<THessianMeasure, TVesselnessMeasure>::
getSigmaMax() const {
    return mySigmaMax;
}

template <typename THessianMeasure, typename TVesselnessMeasure>
void
DGtal::MultiscaleVesselness<THessianMeasure, TVesselnessMeasure>::
setNumberSigmaSteps(double number) {
    myNumberSigmaSteps = number;
}

template <typename THessianMeasure, typename TVesselnessMeasure>
double
DGtal::MultiscaleVesselness<THessianMeasure, TVesselnessMeasure>::
getNumberSigmaSteps() const {
    return myNumberSigmaSteps;
}

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined MultiscaleVesselness_h

#undef MultiscaleVesselness_RECURSES
#endif // else defined(MultiscaleVesselness_RECURSES)