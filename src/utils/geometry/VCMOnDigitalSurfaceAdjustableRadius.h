#ifndef VCM_ON_DIGITAL_SURFACE_ADJUSTABLE_RADIUS_H
#define VCM_ON_DIGITAL_SURFACE_ADJUSTABLE_RADIUS_H

#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/geometry/volumes/estimation/VoronoiCovarianceMeasure.h"
#include "DGtal/base/ConstAlias.h"
#include "geometry/VCMAdjustableRadius.h"
#include "geometry/MedialAxis.h"

template <typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
class VCMOnDigitalSurfaceAdjustableRadius :
	public DGtal::VoronoiCovarianceMeasureOnDigitalSurface<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>
{


	/////////////////////////////////////////////////////////////////////////////
	// template class VoronoiCovarianceMeasureOnDigitalSurface
	/**
	 * Description of template class
	 * 'VoronoiCovarianceMeasureOnDigitalSurface' <p> \brief Aim: This
	 * class specializes the Voronoi covariance measure for digital
	 * surfaces. It adds notably the embedding of surface elements, the
	 * diagonalisation of the VCM, and the orientation of the first VCM
	 * eigenvector toward the interior of the surface.
	 *
	 * @note Documentation in \ref moduleVCM_sec3_1.
	 *
	 * @see VoronoiCovarianceMeasure
	 *
	 * @tparam TDigitalSurfaceContainer the type of digital surface
	 * container (model of CDigitalSurfaceContainer).
	 *
	 * @tparam TSeparableMetric a model of CSeparableMetric used for
	 * computing the Voronoi map (e.g. Euclidean metric is
	 * DGtal::ExactPredicateLpSeparableMetric<TSpace, 2> )
	 *
	 * @tparam TKernelFunction the type of the kernel function chi_r used
	 * for integrating the VCM, a map: Point -> Scalar.
	 */

public:
	typedef DGtal::VoronoiCovarianceMeasureOnDigitalSurface<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction> Base;
	typedef typename Base::Metric                         Metric;  ///< the chosen metric
	typedef typename Base::KernelFunction                  KernelFunction;  ///< the kernel function
	typedef typename Base::Surface Surface;
	typedef typename Base::KSpace KSpace;
	typedef typename Base::Surfel2PointEmbedding Surfel2PointEmbedding;
	typedef typename Base::Space                    Space;  ///< the digital space
	typedef typename Base::Point                    Point;  ///< the digital points
	typedef typename Base::Scalar                     Scalar;  ///< the "real number" type
	typedef typename Base::Surfel Surfel;
	typedef typename Base::DigitalSurfaceContainer DigitalSurfaceContainer;
	typedef typename Base::EigenStructure EigenStructure;
	typedef typename Base::Normals Normals;
	typedef typename Base::ConstIterator ConstIterator;
	typedef typename Base::LinearAlgebraTool LinearAlgebraTool;
	typedef typename Base::VectorN                   VectorN;  ///< n-dimensional R-vector
	typedef typename Base::MatrixNN                 MatrixNN;  ///< nxn R-matrix
	typedef typename Base::Domain Domain;
	typedef typename DGtal::DigitalSetSelector< Domain, DGtal::BIG_DS + DGtal::HIGH_BEL_DS >::Type DigitalSet;
	typedef  VCMAdjustableRadius<Space, TSeparableMetric> VCM;

	typedef std::map<Point, double> Point2Radius;
	// ----------------------- Standard services ------------------------------
public:

	/**
	 * Destructor.
	 */
	~VCMOnDigitalSurfaceAdjustableRadius();

	/**
	 * Constructor. Computes the VCM of the given \a surface.
	 *
	 * @param _surface the digital surface that is aliased in this. The
	 * user can \b secure the aliasing by passing a
	 * CountedConstPtrOrConstPtr.
	 *
	 * @param _surfelEmbedding the chosen embedding for surfels.
	 *
	 * @param _R the offset radius for the set of points. Voronoi cells
	 * are intersected with this offset. The unit corresponds to a step in the digital space.
	 *
	 * @param _r (an upper bound of) the radius of the support of the
	 * kernel function \a chi_r (note \f$\chi_r\f$ in the VCM
	 * paper). The unit corresponds to a step in the digital
	 * space. This parameter is used for preparing the data structure
	 * that answers to proximity queries.
	 *
	 * @param chi_r the kernel function whose support has radius less
	 * or equal to \a r.
	 *
	 * @param t the radius for the trivial normal estimator, which is
	 * used for finding the correct orientation inside/outside for the
	 * VCM.
	 *
	 * @param aMetric an instance of the metric.
	 *
	 * @param verbose if 'true' displays information on ongoing computation.
	 */
	VCMOnDigitalSurfaceAdjustableRadius( DGtal::ConstAlias< Surface > _surface,
										 Surfel2PointEmbedding _surfelEmbedding,
										 Scalar _R, Scalar _r,
										 KernelFunction chi_r,
										 const MedialAxis<DigitalSet>& medialAxis,
										 Scalar t = 2.5, Metric aMetric = Metric(),
										 bool verbose = false );


protected:

	/**
	 * Constructor.
	 * Forbidden by default (protected to avoid g++ warnings).
	 */
	VCMOnDigitalSurfaceAdjustableRadius();

protected:
	MedialAxis<DigitalSet> myMedialAxis;
private:

	/**
	 * Copy constructor.
	 * @param other the object to clone.
	 * Forbidden by default.
	 */
	VCMOnDigitalSurfaceAdjustableRadius ( const VCMOnDigitalSurfaceAdjustableRadius & other );

	/**
	 * Assignment.
	 * @param other the object to copy.
	 * @return a reference on 'this'.
	 * Forbidden by default.
	 */
	VCMOnDigitalSurfaceAdjustableRadius & operator= ( const VCMOnDigitalSurfaceAdjustableRadius & other );


};

template <typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
~VCMOnDigitalSurfaceAdjustableRadius() {
	Base::~Base();
}


template <typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
VCMOnDigitalSurfaceAdjustableRadius( const VCMOnDigitalSurfaceAdjustableRadius& other ) : myMedialAxis(other.myMedialAxis) {
	Base::Base( other );
}

template <typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
VCMOnDigitalSurfaceAdjustableRadius() : Base::Base(), myMedialAxis(DigitalSet(Domain(Point(0,0,0), Point(0,0,0)))){

}

template <typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>&
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
operator=( const VCMOnDigitalSurfaceAdjustableRadius& other )
{
  Base::operator=( other );
  myMedialAxis = other.myMedialAxis;
  return *this;
}

template <typename TDigitalSurfaceContainer, typename TSeparableMetric, typename TKernelFunction>
VCMOnDigitalSurfaceAdjustableRadius<TDigitalSurfaceContainer, TSeparableMetric, TKernelFunction>::
VCMOnDigitalSurfaceAdjustableRadius(DGtal::ConstAlias< Surface > _surface,
									Surfel2PointEmbedding _surfelEmbedding,
									Scalar _R, Scalar _r,
									KernelFunction chi_r,const MedialAxis<DigitalSet>& medialAxis,
									Scalar t, Metric aMetric,
									bool verbose ) :
	Base::mySurface( _surface ), Base::mySurfelEmbedding( _surfelEmbedding ), Base::myChi( chi_r ),
    Base::myVCM( VCM(_R, _r, aMetric, verbose) ), Base::myRadiusTrivial( t ), myMedialAxis(medialAxis) {
	if ( verbose ) DGtal::trace.beginBlock( "Computing VCM on digital surface." );
	const KSpace & ks = this->mySurface->container().space();
	std::vector<Point> vectPoints;


// Get points.
	if ( verbose ) DGtal::trace.beginBlock( "Getting points." );
	std::set<Point> pointSet;
	for ( ConstIterator it = this->mySurface->begin(), itE = this->mySurface->end(); it != itE; ++it )
		getPoints( std::inserter( pointSet, pointSet.begin() ), *it );
	vectPoints.resize( pointSet.size() );
	std::copy( pointSet.begin(), pointSet.end(), vectPoints.begin() );
	if ( verbose ) DGtal::trace.endBlock();

// Compute Voronoi Covariance Matrix for all points.
	this->myVCM.init( vectPoints.begin(), vectPoints.end() );

// Compute VCM( chi_r ) for each point.
	if ( verbose ) DGtal::trace.beginBlock ( "Integrating VCM( chi_r(p) ) for each point." );
	int i = 0;
	DigitalSet medialAxisSet = medialAxis.pointSet();
// HatPointFunction< Point, Scalar > chi_r( 1.0, r );
	for (auto it = vectPoints.begin(), itE = vectPoints.end();
		  it != itE; ++it )
	{
		Point p = *it;
		if ( verbose ) DGtal::trace.progressBar( ++i, vectPoints.size() );

		Point closestPointToCurrent = *min_element(medialAxisSet.begin(), medialAxisSet.end(), [&](const Point& one, const Point& two) {
				return DGtal::Z3i::l2Metric(one, p) < DGtal::Z3i::l2Metric(two, p);
			});
		std::vector<DGtal::Z3i::RealPoint> normalsToPoint{DGtal::Z3i::RealPoint(0,0,1), DGtal::Z3i::RealPoint(0,1,0), DGtal::Z3i::RealPoint(1,0,0)};
		double radius = this->myMetric(closestPointToCurrent,p) + 2.0;
		this->myChi = KernelFunction( 1.0, radius);
		this->myVCM.setMySmallR(radius);


		MatrixNN measure = this->myVCM.measure( this->myChi, p );
// On diagonalise le résultat.
		EigenStructure & evcm = this->myPt2EigenStructure[ p ];
		LinearAlgebraTool::getEigenDecomposition( measure, evcm.vectors, evcm.values );
	}
	this->myVCM.clean(); // free some memory.
	if ( verbose ) DGtal::trace.endBlock();

	if ( verbose ) DGtal::trace.beginBlock ( "Computing average orientation for each surfel." );
	typedef DGtal::functors::HatFunction<Scalar> Functor;
	Functor fct( 1.0, this->myRadiusTrivial );
	typedef DGtal::functors::ElementaryConvolutionNormalVectorEstimator< Surfel, DGtal::CanonicSCellEmbedder<KSpace> >
		SurfelFunctor;
	typedef DGtal::LocalEstimatorFromSurfelFunctorAdapter< DigitalSurfaceContainer, Metric, SurfelFunctor, Functor>
		NormalEstimator;

	DGtal::CanonicSCellEmbedder<KSpace> canonic_embedder( ks );
	SurfelFunctor surfelFct( canonic_embedder, 1.0 );
	NormalEstimator estimator;
	estimator.attach( *(this->mySurface) );
	estimator.setParams( aMetric, surfelFct, fct , this->myRadiusTrivial);
	estimator.init( 1.0,  this->mySurface->begin(), this->mySurface->end());
	i = 0;
	std::vector<Point> pts;
	int surf_size = this->mySurface->size();
	for ( ConstIterator it = this->mySurface->begin(), itE = this->mySurface->end(); it != itE; ++it )
	{
		if ( verbose ) DGtal::trace.progressBar(++i, surf_size );
		Surfel s = *it;
		Normals & normals = this->mySurfel2Normals[ s ];
// get rough estimation of normal
		normals.trivialNormal = estimator.eval( it );
// get points associated with surfel s
		getPoints( std::back_inserter( pts ), s );
		for ( typename std::vector<Point>::const_iterator itPts = pts.begin(), itPtsE = pts.end();
			  itPts != itPtsE; ++itPts )
		{
			Point p = *itPts;
			const EigenStructure& evcm = this->myPt2EigenStructure[ p ];
			VectorN n = evcm.vectors.column( 2 );
			if ( n.dot( normals.trivialNormal ) < 0 ) normals.vcmNormal -= n;
			else                                      normals.vcmNormal += n;
		}
		if ( pts.size() > 1 ) normals.vcmNormal /= pts.size();
		pts.clear();
	}
	if ( verbose ) DGtal::trace.endBlock();

	if ( verbose ) DGtal::trace.endBlock();
}

#endif
