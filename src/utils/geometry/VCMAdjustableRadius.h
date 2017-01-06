#ifndef VORONOI_COVARIANCE_MEASURE_ADJUSTABLE_RADIUS_H
#define VORONOI_COVARIANCE_MEASURE_ADJUSTABLE_RADIUS_H

// Inclusions
#include <cmath>
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/geometry/volumes/estimation/VoronoiCovarianceMeasure.h"


template <typename TSpace, typename TSeparableMetric>
class VCMAdjustableRadius : public DGtal::VoronoiCovarianceMeasure<TSpace, TSeparableMetric>
{
public:
	typedef DGtal::VoronoiCovarianceMeasure<TSpace, TSeparableMetric> Base;
	typedef typename Base::MatrixNN MatrixNN;
	typedef typename Base::Point Point;
	typedef typename Base::Vector Vector;
	typedef typename Base::Point2ScalarFunction Point2ScalarFunction;
	typedef typename Base::Metric Metric;

	// ----------------------- Standard services ------------------------------
public:

	/**
	 * Constructor.
	 *
	 * @param _R the offset radius for the set of points. Voronoi cells
	 * are intersected with this offset. The unit corresponds to a step in the digital space.
	 *
	 * @param _r (an upper bound of) the radius of the support of
	 * forthcoming kernel functions (\f$ \chi_r \f$). The unit
	 * corresponds to a step in the digital space. This parameter is
	 * used for preparing the data structure that answers to proximity
	 * queries.
	 *
	 * @param aMetric an instance of the metric.
	 * @param verbose if 'true' displays information on ongoing computation.
	 */
	VCMAdjustableRadius( double _R, double _r, Metric aMetric = Metric(), bool verbose = false, bool inDomain = false );

	/**
	 * Destructor.
	 */
	~VCMAdjustableRadius();

	void setMySmallR(double r) {this->mySmallR = r;}

	template <typename Point2ScalarFunction>
	MatrixNN measure( const Point2ScalarFunction& chi_r, const Point& p) const;

	// template <typename Point2ScalarFunction>
	// MatrixNN measure( const std::vector<Point>& neighbors, Point2ScalarFunction chi_r, Point p ) const;

	// template <typename Point2ScalarFunction>
	// MatrixNN measureJunction( const Vector& dirVector, Point2ScalarFunction chi_r, Point p ) const;

	// ----------------------- Interface --------------------------------------

	// ------------------------- Protected Datas ------------------------------
	// ------------------------- Private Datas --------------------------------
private:


	// ------------------------- Hidden services ------------------------------
protected:

	/**
	 * Constructor.
	 * Forbidden by default (protected to avoid g++ warnings).
	 */
	VCMAdjustableRadius();

private:

	/**
	 * Copy constructor.
	 * @param other the object to clone.
	 * Forbidden by default.
	 */
	VCMAdjustableRadius ( const VCMAdjustableRadius & other );

	/**
	 * Assignment.
	 * @param other the object to copy.
	 * @return a reference on 'this'.
	 * Forbidden by default.
	 */
	VCMAdjustableRadius & operator= ( const VCMAdjustableRadius & other );

	// ------------------------- Internals ------------------------------------
private:

};

template <typename TSpace, typename TSeparableMetric>
inline
VCMAdjustableRadius<TSpace,TSeparableMetric>::
~VCMAdjustableRadius()
{
  this->clean();
}

template <typename TSpace, typename TSeparableMetric>
inline
VCMAdjustableRadius<TSpace,TSeparableMetric>::
VCMAdjustableRadius( double _R, double _r, Metric aMetric, bool verbose, bool isInDomain )
	:  Base::Base( _R, _r, aMetric, verbose, isInDomain) {}

template <typename TSpace, typename TSeparableMetric>
inline
VCMAdjustableRadius<TSpace,TSeparableMetric>::
VCMAdjustableRadius( const VCMAdjustableRadius& other )
	:  Base::Base( other ) {}


template <typename TSpace, typename TSeparableMetric>
inline
VCMAdjustableRadius<TSpace,TSeparableMetric>&
VCMAdjustableRadius<TSpace,TSeparableMetric>::
operator=( const VCMAdjustableRadius& other )
{
  Base::operator=( other );
  return *this;
}


template <typename TSpace, typename TSeparableMetric>
template <typename Point2ScalarFunction>
typename VCMAdjustableRadius<TSpace,TSeparableMetric>::MatrixNN
VCMAdjustableRadius<TSpace,TSeparableMetric>::measure( const VCMAdjustableRadius<TSpace,TSeparableMetric>::Point2ScalarFunction& chi_r,
													   const VCMAdjustableRadius<TSpace,TSeparableMetric>::Point& p )
{
  ASSERT( myProximityStructure != 0 );

  Point b = myProximityStructure->bin( p );
  Ball<Point> ball(p,mySmallR);
  std::vector<Point> neighbors = ball.intersection(myContainer);
  MatrixNN vcm;
  // std::cout << *it << " has " << neighbors.size() << " neighbors." << std::endl;
  for ( typename std::vector<Point>::const_iterator it_neighbors = neighbors.begin(),
          it_neighbors_end = neighbors.end(); it_neighbors != it_neighbors_end; ++it_neighbors )

  {
      Point q = *it_neighbors;
      Scalar coef = chi_r( q - p );
      if ( coef > 0.0 )
        {
          MatrixNN vcm_q = myVCM.at(q);
          vcm_q *= coef;
          vcm += vcm_q;
        }
    }
  return vcm;
}

// //-----------------------------------------------------------------------------
// template <typename TSpace, typename TSeparableMetric>
// template <typename Point2ScalarFunction>
// inline
// typename VCMAdjustableRadius<TSpace,TSeparableMetric>::MatrixNN
// VCMAdjustableRadius<TSpace,TSeparableMetric>::
// measure( const std::vector<typename VCMAdjustableRadius<TSpace,TSeparableMetric>::Point>& neighbors,
// 		 typename VCMAdjustableRadius<TSpace,TSeparableMetric>::Point2ScalarFunction chi_r,
// 		 typename VCMAdjustableRadius<TSpace,TSeparableMetric>::Point p ) const
// {
//   ASSERT( myProximityStructure != 0 );

//   Point b = myProximityStructure->bin( p );

//   MatrixNN vcm;
//   // std::cout << *it << " has " << neighbors.size() << " neighbors." << std::endl;
//   for ( typename std::vector<Point>::const_iterator it_neighbors = neighbors.begin(),
//           it_neighbors_end = neighbors.end(); it_neighbors != it_neighbors_end; ++it_neighbors )
//     {
//       Point q = *it_neighbors;
// 	   Scalar coef = chi_r( q - p );
//       // Scalar coef = 1.0;
//       if ( coef > 0.0 )
//         {
//           typename std::map<Point,MatrixNN>::const_iterator it = myVCM.find( q );
//           if ( it != myVCM.end() ) {
// 			  MatrixNN vcm_q = it->second;
// 			  vcm_q *= coef;
// 			  vcm += vcm_q;
// 		  }
//         }
//     }
//   return vcm;
// }


// template <typename TSpace, typename TSeparableMetric>
// template <typename Point2ScalarFunction>
// inline
// typename VCMAdjustableRadius<TSpace,TSeparableMetric>::MatrixNN
// VCMAdjustableRadius<TSpace,TSeparableMetric>::
// measureJunction( const Vector& dirVector, Point2ScalarFunction chi_r, Point p ) const
// {
//   ASSERT( myProximityStructure != 0 );
//   typedef DGtal::EigenDecomposition<3,double> LinearAlgebraTool;
//   typedef typename TSpace::RealVector Vector;

//   Ball<Point> ball(p,mySmallR);
//   std::vector<Point> neighbors = ball.pointsInHalfBall(dirVector);

//   MatrixNN vcm, evec, null;

//   Vector eval;
//   std::map<Point, Vector> mapPoint2Normal;
//   for ( typename std::vector<Point>::const_iterator it_neighbors = neighbors.begin(),
//           it_neighbors_end = neighbors.end(); it_neighbors != it_neighbors_end; ++it_neighbors )
//     {
//       Point q = *it_neighbors;
//       Scalar coef = chi_r( q - p );
//       if ( coef > 0.0 )
//         {
//           typename std::map<Point,MatrixNN>::const_iterator it = myVCM.find( q );
//           if ( it != myVCM.end() ) {
// 			  MatrixNN vcm_q = it->second;
//               vcm_q *= coef;
// 			  vcm += vcm_q;
// 		  }
//         }
//     }
//   return vcm;
// }



#endif
