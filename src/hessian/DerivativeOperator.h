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
 * @file
 * @author Florent Grelard (\c florent.grelard@u-bordeaux.fr )
 * LaBRI, Bordeaux University
 *
 * @date 2017/03/28
 *
 *
 * This file is part of the DGtal library.
 */



#if defined(DerivativeOperator_RECURSES)
#error Recursive header files inclusion detected in DerivativeOperator.h
#else // defined(DerivativeOperator_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DerivativeOperator_RECURSES

#if !defined DerivativeOperator_h
/** Prevents repeated inclusion of headers. */
#define DerivativeOperator_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
//////////////////////////////////////////////////////////////////////////////

namespace DGtal {

    /////////////////////////////////////////////////////////////////////////////
    // template class DerivativeOperator
    /**
     * Description of template class 'DerivativeOperator' <p>
     * \brief Aim:
     */
    template <typename TValue>
    class DerivativeOperator {
    public:
        typedef TValue Value;
        typedef std::vector<Value> CoefficientVector;

        // ----------------------- Standard services ------------------------------
    public:

        /**
         * Default constructor.
         */
        DerivativeOperator() : myOrder(1) {}

        DerivativeOperator(unsigned int order) {
            myOrder = order;
        }

        /**
         * Destructor.
         */
        ~DerivativeOperator() = default;

        /**
         * Copy constructor.
         * @param other the object to clone.
         */
        DerivativeOperator(const DerivativeOperator &other) = default;

        /**
         * Move constructor.
         * @param other the object to move.
         */
        DerivativeOperator(DerivativeOperator &&other) = delete;

        /**
         * Copy assignment operator.
         * @param other the object to copy.
         * @return a reference on 'this'.
         */
        DerivativeOperator &operator=(const DerivativeOperator &other) = delete;

        /**
         * Move assignment operator.
         * @param other the object to move.
         * @return a reference on 'this'.
         */
        DerivativeOperator &operator=(DerivativeOperator &&other) = delete;

        // ----------------------- Interface --------------------------------------
    public:

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

        void setOrder(unsigned int order);

        unsigned int getOrder() const;

        CoefficientVector computeCoefficients();


        // ------------------------- Protected Datas ------------------------------
    protected:

        // ------------------------- Private Datas --------------------------------
    private:
        unsigned int myOrder;

        // ------------------------- Hidden services ------------------------------
    protected:

        // ------------------------- Internals ------------------------------------
    private:

    }; // end of class DerivativeOperator


    /**
     * Overloads 'operator<<' for displaying objects of class 'DerivativeOperator'.
     * @param out the output stream where the object is written.
     * @param object the object of class 'DerivativeOperator' to write.
     * @return the output stream after the writing.
     */
    template <typename TValue>
    std::ostream &
    operator<<(std::ostream &out, const DerivativeOperator<TValue> &object) {
        object.selfDisplay(out);
        return out;
    }

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

template <typename TValue>
void
DGtal::DerivativeOperator<TValue>::
selfDisplay(std::ostream &out) const {
    out << "[DerivativeOperator] order=" << this->getOrder();
}

template <typename TValue>
bool
DGtal::DerivativeOperator<TValue>::
isValid() const {
    return true;
}

template <typename TValue>
void
DGtal::DerivativeOperator<TValue>::
setOrder(unsigned int order) {
    myOrder = order;
}

template <typename TValue>
unsigned int
DGtal::DerivativeOperator<TValue>::
getOrder() const {
    return myOrder;
}

template <typename TValue>
typename DGtal::DerivativeOperator<TValue>::CoefficientVector
DGtal::DerivativeOperator<TValue>::
computeCoefficients() {
    unsigned int i;
    unsigned int j;
    Value previous;
    Value next;
    const unsigned int w = 2 * ((myOrder + 1) / 2) + 1;
    CoefficientVector coeff(w);

    coeff[w / 2] = 1.0;
    for (i = 0; i < myOrder / 2; i++) {
        previous = coeff[1] - 2 * coeff[0];
        for (j = 1; j < w - 1; j++) {
            next = coeff[j - 1] + coeff[j + 1] - 2 * coeff[j];
            coeff[j - 1] = previous;
            previous = next;
        }
        next = coeff[j - 1] - 2 * coeff[j];
        coeff[j - 1] = previous;
        coeff[j] = next;
    }
    for (i = 0; i < myOrder % 2; i++) {
        previous = 0.5 * coeff[1];
        for (j = 1; j < w - 1; j++) {
            next = -0.5 * coeff[j - 1] + 0.5 * coeff[j + 1];
            coeff[j - 1] = previous;
            previous = next;
        }
        next = -0.5 * coeff[j - 1];
        coeff[j - 1] = previous;
        coeff[j] = next;
    }

    return coeff;
}


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DerivativeOperator_h

#undef DerivativeOperator_RECURSES
#endif // else defined(DerivativeOperator_RECURSES)