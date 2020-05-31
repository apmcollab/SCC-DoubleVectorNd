/*
 * SCC_DoubleVectorNdSkin.h
 *
 * Classes : SCC::DoubleVector1dSkin
 *           SCC::DoubleVector2dSkin
 *           SCC::DoubleVector3dSkin
 *
 *
 * These classes provide a DoubleVectorNd wrapper (or skin) for existing data.
 *
 * It is assumed that the data associated with an M by N by K array
 * are stored in a double array of size MxNxK with data storage
 * by ROWS (C convention).
 *
 * There is no internal validation on the assumption about the size and storage
 * format associated with the data pointer used to initialize an instance
 * of the DoubleVectorNdSkin class.
 *
 * The methods are primarily provided as a utility for utilizing
 * methods whose input and output arguments are SCC::DoubleVectorNd
 * instances.
 *
 * Typical usage
 *
 * G = some type of double vector with data stored in a double* objectData
 *
 * (1) DoubleVectorNdSkin GS(G.objectData,M,N,K)
 *
 * (2) Call member function requiring DoubleVectorNd arguments and specify GS instead of G
 *
 * Transformations of the data associated with GS will will be performed on the
 * data of the associated grid object G.
 *
 * No data copying is required.
 *
 *
 *  Created on: Sep 30, 2016
 *      Author: anderson
 */

/*
#############################################################################
#
# Copyright 2015-16 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/

#include "SCC_DoubleVector1d.h"
#include "SCC_DoubleVector2d.h"
#include "SCC_DoubleVector3d.h"


#ifndef SCC_DOUBLE_VECTOR_ND_SKIN_
#define SCC_DOUBLE_VECTOR_ND_SKIN_

namespace SCC
{

class DoubleVector1dSkin : public DoubleVector1d
{
public:

	 DoubleVector1dSkin() : DoubleVector1d()
	 {
		 this->dataPtr = nullptr;
	 };


	 DoubleVector1dSkin(const DoubleVector1dSkin& G)
	 {
		 initialize(G.dataPtr, G.index1Size);
	 };


	 DoubleVector1dSkin(double* data1d, long M)
	 {
		 initialize(data1d, M);
	 }

	 void initialize(double* data1d, long M)
	 {
	    // Set internal data pointer to nullptr then call base class initialize

		this->dataPtr = nullptr;
		DoubleVector1d::initialize();

		// Set internal data pointer to supplied data pointer

		this->dataPtr = data1d;

		// Pack structure informaiton

		// DoubleVector1d structure information

		this->index1Size = M;

	 }


	 virtual ~DoubleVector1dSkin()
	 {
		// Don't delete the data associated with the instance

		 this->dataPtr = nullptr;
	 }
};


class DoubleVector2dSkin : public DoubleVector2d
{
public:

	 DoubleVector2dSkin() : DoubleVector2d()
	 {
		 this->dataPtr = nullptr;
	 };


	 DoubleVector2dSkin(const DoubleVector2dSkin& G)
	 {
		 initialize(G.dataPtr, G.index1Size, G.index2Size);
	 };


	 DoubleVector2dSkin(double* data2d, long M, long N)
	 {
		 initialize(data2d, M, N);
	 }


	 void initialize(double* data2d, long M, long N)
	 {
	    // Set internal data pointer to null ptr then call base class initialize

		this->dataPtr = nullptr;
		DoubleVector2d::initialize();

		// Set internal data pointer to supplied data pointer

		this->dataPtr = data2d;

		// Pack structure informaiton

		// DoubleVector3d structure information

		this->index1Size = M;
    	this->index2Size = N;;
	 }


	 virtual ~DoubleVector2dSkin()
	 {
		// Don't delete the data associated with the instance

		 this->dataPtr = nullptr;
	 }
};



class DoubleVector3dSkin : public DoubleVector3d
{
public:

	 DoubleVector3dSkin() : DoubleVector3d()
	 {
		 this->dataPtr = nullptr;
	 };


	 DoubleVector3dSkin(const DoubleVector3dSkin& G)
	 {
		 initialize(G.dataPtr, G.index1Size, G.index2Size, G.index3Size);
	 };


	 DoubleVector3dSkin(double* data3d, long M, long N, long K)
	 {
		 initialize(data3d,M,N,K);
	 }


	 void initialize(double* data3d, long M, long N, long K)
	 {
	    // Set internal data pointer to null ptr then call base class initialize

		this->dataPtr = nullptr;
		DoubleVector3d::initialize();

		// Set internal data pointer to supplied data pointer

		this->dataPtr = data3d;

		// Pack structure informaiton

		// DoubleVector3d structure information

		this->index1Size = M;
    	this->index2Size = N;
    	this->index3Size = K;
	 }


	 virtual ~DoubleVector3dSkin()
	 {
		// Don't delete the data associated with the instance

		 this->dataPtr = nullptr;
	 }
};
}




#endif /*SCC_DoubleVectorNdSkin_*/
