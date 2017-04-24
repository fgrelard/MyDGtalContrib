#ifndef DELTA_VCM_H
#define DELTA_VCM_H

#include "DGtal/math/linalg/SimpleMatrix.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

namespace DGtal {
    typedef SimpleMatrix<double, 2, 2> MatrixDouble;

    bool operator!=(const MatrixDouble &m1, const MatrixDouble &m2) { return !(m1 == m2); }

    typedef SimpleMatrix<float, 2, 2> MatrixFloat;

    bool operator!=(const MatrixFloat &m1, const MatrixFloat &m2) { return !(m1 == m2); }

    namespace functors {
        bool operator==(Identity f1, Identity f2) { return true; }
    }
}


template<typename TDistanceLikeFunction>
class DeltaVCM {
public:
    typedef TDistanceLikeFunction DistanceLikeFunction;
    typedef typename DistanceLikeFunction::Value Value;
    typedef typename DistanceLikeFunction::Point Point;
    typedef typename DistanceLikeFunction::Domain Domain;
    typedef typename Domain::Space Space;
    typedef typename Domain::Dimension Dimension;
    typedef typename Space::Integer Integer;
    typedef typename Space::RealVector RealVector;
    typedef typename DistanceLikeFunction::Value Scalar;
    typedef DGtal::SimpleMatrix<Scalar,
            Space::dimension,
            Space::dimension> Matrix; ///< the type for nxn matrix of real numbers.
    // typedef typename Matrix::RowVector Vector;   ///< the type for N-vector of real numbers

    typedef DGtal::ImageContainerBySTLVector<Domain, Matrix> MatrixField;

public:
    DeltaVCM(const DistanceLikeFunction &delta, double R, double r)
            : myDelta(delta), myR(R), myr(r),
              myVCM(delta.domain()) {
        init();
    }

    void init() {
        Matrix M;
        for (typename Domain::ConstIterator it = myDelta.domain().begin(),
                     itE = myDelta.domain().end(); it != itE; ++it) {
            Point p = *it;
            // eliminates points too far away.
            if (myDelta(p) > myR) continue;
            RealVector n = myDelta.projection(p);
            Point q = Point((Integer) round(p[0] + n[0]),
                            (Integer) round(p[1] + n[1]));
            // eliminates projections going outside the domain.
            if (q != myDelta.box(q)) continue;
            for (Dimension i = 0; i < Space::dimension; ++i)
                for (Dimension j = 0; j < Space::dimension; ++j)
                    M.setComponent(i, j, n[i] * n[j]);
            myVCM.setValue(q, myVCM(q) + M); // add tensor n x n
        }
    }

public:
    inline Domain domain() const {
        return myDelta.domain();
    }

    /**
       Computes the Voronoi Covariance Measure of the function \a chi_r.

       @tparam Point2ScalarFunction the type of a functor
       Point->Scalar. For instance functors::HatPointFunction and
       functors::BallConstantPointFunction are models of this type.

       @param chi_r the kernel function whose support is included in
       the cube centered on the origin with edge size 2r.

       @param p the point where the kernel function is moved. It must lie within domain.
    */
    template<typename Point2ScalarFunction>
    Matrix measure(Point2ScalarFunction chi_r, Point p) const {
        Integer r = (Integer) ceil(myr);
        Point low = domain().lowerBound().sup(p - Point::diagonal(r));
        Point up = domain().upperBound().inf(p + Point::diagonal(r));
        //trace.info() << "r=" << r << " low=" << low << " up=" << up << std::endl;
        Domain local(low, up);
        Scalar mass = 0.0;
        Matrix M;
        for (typename Domain::ConstIterator it = local.begin(), itE = local.end();
             it != itE; ++it) {
            Point q = *it;
            Scalar chi = chi_r(q - p);
            if (chi <= 0.0) continue;
            // JOL: to check : I don't know if you should weight chi by the measure.
            // (0) no correction
            //chi *= myDelta.measure()( q );     // (1) more stable than (2) and (0)
            // chi *= myProjectedMeasure( q ); // (2)
            //trace.info() << "chi=" << chi << " VCM=" << myVCM( q ) << endl;
            M += ::operator*(chi, myVCM(q)); // workaround simplematrix bug in DGtal.
        }
        return M;
    }

    // chi_r is normalized to have mass 1.
    template<typename Point2ScalarFunction>
    Matrix measure1(Point2ScalarFunction chi_r, Point p) const {
        Integer r = (Integer) ceil(myr);
        Point low = domain().lowerBound().sup(p - Point::diagonal(r));
        Point up = domain().upperBound().inf(p + Point::diagonal(r));
        //trace.info() << "r=" << r << " low=" << low << " up=" << up << std::endl;
        Domain local(low, up);
        Scalar mass = 0.0;
        Matrix M;
        for (typename Domain::ConstIterator it = local.begin(), itE = local.end();
             it != itE; ++it) {
            Point q = *it;
            Scalar chi = chi_r(q - p);
            if (chi <= 0.0) continue;
            mass += chi;
            // JOL: to check : I don't know if you should weight chi by the measure.
            chi *= myDelta.measure()(q);
            //trace.info() << "chi=" << chi << " VCM=" << myVCM( q ) << endl;
            M += myVCM(q) * chi; // else ::operator*(chi, myVCM( q )); workaround simplematrix bug in DGtal.
        }
        return mass > 0.0 ? M / mass : M;
    }

private:
    DistanceLikeFunction myDelta;
    Scalar myR;
    Scalar myr;
    MatrixField myVCM;
};


#endif
