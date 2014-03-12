/*
 * Copyright (c) 2008-2009 Children's Hospital Boston.
 *
 * This software is licensed by the copyright holder under the terms of the
 * Open Software License version 3.0.
 * http://www.opensource.org/licenses/osl-3.0.php
 *
 * Attribution Notice.
 *
 * This research was carried out in the Computational Radiology Laboratory of
 * Children's Hospital, Boston and Harvard Medical School.
 * http://www.crl.med.harvard.edu
 * For more information contact: simon.warfield@childrens.harvard.edu
 *
 * This research work was made possible by Grant Number R01 RR021885 (Principal
 * Investigator: Simon K. Warfield, Ph.D.) to Children's Hospital, Boston
 * from the National Center for Research Resources (NCRR), a component of the
 * National Institutes of Health (NIH).
*/

#ifndef CRL_KRIGEAGE_INTERPOLATION_H
#define CRL_KRIGEAGE_INTERPOLATION_H

#include "itkSingleValuedCostFunction.h"
#include "itkNeighborhoodIterator.h"
#include "itkMatrix.h"


using namespace itk;

namespace crl{ 

	/**********************************************************************************************//**
	 * \class	KrigeageInterpolation
	 *
	 * \brief	Krigeage interpolation. 
	 *
	 * \author	Benoit Scherrer
	 * \date	February 2011
	*************************************************************************************************/
	template <typename TSpaceType, unsigned int SpaceDimension, typename TValueType, unsigned int ValueDimension> 
	class KrigeageInterpolation  
	{
	public:
		typedef vnl_vector_fixed<TSpaceType, SpaceDimension> SpacePointType;
		typedef vnl_vector_fixed<TValueType, ValueDimension> ValueType;
		typedef double TInternalMatrixType;

	public:

		/**********************************************************************************************//**
		 * \fn	KrigeageInterpolation()
		 *
		 * \brief	Default constructor. 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		*************************************************************************************************/
		KrigeageInterpolation()
		{
			m_InterpolatorEstimated=false;
		}

		virtual ~KrigeageInterpolation()
		{
		}

		/**********************************************************************************************//**
		 * \fn	void addObservation( const SpacePointType &point, TValueType val )
		 *
		 * \brief	Provided for convenience. Adds a scalar observation (ValueDimension should be 1). 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		 *
		 * \param	point	The point. 
		 * \param	val		The value. 
		*************************************************************************************************/
		void addObservation( const SpacePointType &point, TValueType val )
		{
			ValueType vals;
			vals[0]=val;
			addObservation(point, vals);
		}

		/**********************************************************************************************//**
		 * \fn	void addObservation( const SpacePointType &point, const ValueType &vals )
		 *
		 * \brief	Adds an observation. 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		 *
		 * \param	point	The point. 
		 * \param	vals	The vals. 
		*************************************************************************************************/
		void addObservation( const SpacePointType &point, const ValueType &vals )
		{
			m_Points.push_back(point);
			m_Values.push_back(vals);
		} ; // end of class

		/**********************************************************************************************//**
		 * \fn	int getNumberOfObservations() const
		 *
		 * \brief	Gets the number of observations. 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		 *
		 * \return	The number of observations. 
		*************************************************************************************************/
		int getNumberOfObservations() const
		{
			return (int)m_Points.size();
		}

		/**********************************************************************************************//**
		 * \fn	std::vector<ValueType>& getObservationValues()
		 *
		 * \brief	Gets the observation values. 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		 *
		 * \return	The values. 
		*************************************************************************************************/
		std::vector<ValueType>& getObservationValues() 
		{
			return m_Values;
		}

		/**********************************************************************************************//**
		 * \fn	void clearObservations()
		 *
		 * \brief	Clears the observations. 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		*************************************************************************************************/
		void clearObservations() 
		{
			m_Points.clear();
			m_Values.clear();
		}

		/**********************************************************************************************//**
		 * \fn	bool computeInterpolator()
		 *
		 * \brief	Calculates the interpolator. 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		 *
		 * \return	true if it succeeds, false if it fails. 
		*************************************************************************************************/
		bool computeInterpolator()
		{
	        if ( getNumberOfObservations() < SpaceDimension )
	             return false;

		    if ( !computeKrigeageMatrix() ) 
				return false;

			computeInterpolatorCoefficients();

		    return m_InterpolatorEstimated;
		}

		/**********************************************************************************************//**
		 * \fn	bool computeKrigeageMatrix(void)
		 *
		 * \brief	Calculates the krigeage matrix. 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		 *
		 * \return	true if it succeeds, false if it fails. 
		*************************************************************************************************/
		bool computeKrigeageMatrix(void)
		{
			int N = getNumberOfObservations();
			int i, j;

			m_KrigeageMatrix.set_size(N+SpaceDimension+1, N+SpaceDimension+1);

			//------------------------------------
			// Copy the observation points
			// Fill the last line and the last column with 1
			//------------------------------------
			for ( j=0 ; j<N ; j++ )
			{
				m_KrigeageMatrix[N][j] = 1.0;
				m_KrigeageMatrix[j][N] = 1.0;

				for ( unsigned int k=0; k<SpaceDimension;k++ )
				{
					m_KrigeageMatrix[N+1+k][j] = m_Points[j][k] ;
					m_KrigeageMatrix[j][N+1+k] = m_Points[j][k] ; 
				}
			}

			//------------------------------------
			// Sub matrix to zero
			//------------------------------------
			for ( unsigned i=N ; i<N+SpaceDimension+1 ; i++ )
				for ( unsigned j=N ; j<N+SpaceDimension+1 ; j++ )
					m_KrigeageMatrix[i][j] = 0.0;

			//------------------------------------
			// Now computes the distances for each points (Kii)
			//------------------------------------
			double val_g;
			for ( i=0 ; i<N ; i++ )
				for ( j=0 ; j<=i ; j++ )
				{
					val_g = computeHofDist(m_Points[i], m_Points[j]);
					
					m_KrigeageMatrix[i][j] = val_g;
					m_KrigeageMatrix[j][i] = val_g;
				}

			//-------------------------------
			// Compute the inverse
			//-------------------------------
			m_KrigeageMatrixInv = vnl_matrix_inverse<double>(m_KrigeageMatrix);
			return true;
		}

		/**********************************************************************************************//**
		 * \fn	virtual double computeHofDist(const SpacePointType& p1, const SpacePointType& p2) const
		 *
		 * \brief	Distance function. Calculates /f$ h( || pj - pi || ) /f$ for two points. 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		 *
		 * \param	p1	The first point. 
		 * \param	p2	The second point. 
		 *
		 * \return	The calculated distance. 
		*************************************************************************************************/
		virtual double computeHofDist(const SpacePointType& p1, const SpacePointType& p2) const
		{
			double dist=0;
			for ( unsigned int k=0;k<SpaceDimension; k++ )
				dist += (p1[k]-p2[k])*(p1[k]-p2[k]);

			dist=sqrt(dist);

			return (dist+3); // dist*dist*dist*log(dist);
		}

		/**********************************************************************************************//**
		 * \fn	bool computeInterpolatorCoefficients()
		 *
		 * \brief	Calculates the interpolator coefficients. 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		 *
		 * \return	true if it succeeds, false if it fails. 
		*************************************************************************************************/
		bool computeInterpolatorCoefficients()
		{
			int N = getNumberOfObservations();
			vnl_matrix<TInternalMatrixType> V(N+SpaceDimension+1, ValueDimension);

			for ( int i=0; i<N; i++ )
			{
				for ( unsigned int k=0; k<ValueDimension;k++ )
					V[i][k] = m_Values[i][k];
			}

			for ( unsigned int i=0; i<SpaceDimension+1; i++ )
				for ( unsigned int k=0; k<ValueDimension;k++ )
					V[N+i][k]=0;

			m_KrigeageCoeffs = m_KrigeageMatrixInv*V;

			m_InterpolatorEstimated=true;
			return true;
		}

		/**********************************************************************************************//**
		 * \fn	ValueType interpolate( const SpacePointType& pt ) const
		 *
		 * \brief	Interpolates. Computes /f$ \sum_{j=1}^N b_j h( || p - pj || ) + a_1 + < \left[ a_2, \
		 * 			ldots, a_{SpaceDimension+1} \right], p > /f$. 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		 *
		 * \param	pt	The point. 
		 *
		 * \return	. 
		*************************************************************************************************/
		ValueType interpolate( const SpacePointType& pt ) const
		{
			int N = getNumberOfObservations();

			ValueType val;
			val.fill(0);

			//------------------------------------
			// Compute the 'derive'
			//------------------------------------
			for ( unsigned int k=0; k<ValueDimension; k++ )
			{
				val[k]=m_KrigeageCoeffs[N][k];

				for ( unsigned int j=0;j<SpaceDimension; j++ )
					val[k] += m_KrigeageCoeffs[N+j+1][k]*pt[j];
			}

			//------------------------------------
			// Compute the 'random fluctuation'
			//------------------------------------
			for ( int i=0; i<N; i++ )
			{
				for ( unsigned int k=0; k<ValueDimension; k++ )
					val[k] += m_KrigeageCoeffs[i][k] * computeHofDist(m_Points[i], pt);
			}

			return val;
		}

		/**********************************************************************************************//**
		 * \fn	TValueType scalarInterpolate( const SpacePointType& pt ) const
		 *
		 * \brief	Provided for convenience. Interpolate and returns a scalar value (ValueDimension
		 * 			should be 1). 
		 *
		 * \author	Benoit Scherrer
		 * \date	February 2011
		 *
		 * \param	pt	The point. 
		 *
		 * \return	. 
		*************************************************************************************************/
		TValueType scalarInterpolate( const SpacePointType& pt ) const
		{
			ValueType val = interpolate(pt);
			return val[0];
		}

	protected:
		std::vector<SpacePointType> m_Points;
		std::vector<ValueType> m_Values;

		bool								m_InterpolatorEstimated;
		vnl_matrix<TInternalMatrixType>		m_KrigeageMatrix;
		vnl_matrix<TInternalMatrixType>		m_KrigeageMatrixInv;
		vnl_matrix<TInternalMatrixType>		m_KrigeageCoeffs;
	};
} // end of namespace crl

#endif

