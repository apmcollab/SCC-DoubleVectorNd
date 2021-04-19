/*
 * SCC_DoubleVector3d.h
 *
 *  Created on: Jun 26, 2015
 *      Author: anderson
 *
 * Decisions:
 * Not making the class backward compatible with C++ versions
 *
 * Use of move semantics
 * Use of nullptr instead of 0 for null pointers (or NULL)
 * Use of std::copy to copy data instead of low level loop
 * Bounds checking completely turned off when _DEBUG is not set
 * Use of assert to facilitate index bounds error location
 *
 * Revised: Nov. 26, 2015
 *        : Jan. 3,  2016   If MS compiler, use std::memcpy instead of std:copy
 *                          to avoid MS compiler warnings about unchecked index ranges.
 *        : Jan. 13,  2016  Added iostream output to facilitate debugging
 *
*/
/*
#############################################################################
#
# Copyright 2015-17 Chris Anderson
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



#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <vector>

#ifdef  _DEBUG
#include <cstdio>
#else
#ifndef NDEBUG
#define NDEBUG
#endif
#endif
#include <cassert>

#undef VERBOSE_OPS

#ifndef SCC_DOUBLE_VECTOR_3D_
#define SCC_DOUBLE_VECTOR_3D_

namespace SCC
{
class DoubleVector3d
{
    public:

    DoubleVector3d()
    {
    dataPtr    = nullptr;
    index1Size = 0;
    index2Size = 0;
    index3Size = 0;
    }

    explicit DoubleVector3d(long m, long n, long p)
    {
    if(m*n*p > 0)
    {
        dataPtr    = new double[m*n*p];
        index1Size = m;
        index2Size = n;
        index3Size = p;
    }
    else {dataPtr    = nullptr; index1Size = 0; index2Size = 0; index3Size = 0;};
    }

   DoubleVector3d(const DoubleVector3d& V)
   {
      #ifdef VERBOSE_OPS
      std::cout << "Standard Copy " << std::endl;
      #endif

      if(V.dataPtr == nullptr) {dataPtr = nullptr; index1Size = 0; index2Size = 0; index3Size = 0; return;}

      dataPtr     = new double[V.index1Size*V.index2Size*V.index3Size];
      index1Size = V.index1Size;
      index2Size = V.index2Size;
      index3Size = V.index3Size;

#ifdef _MSC_VER
     std::memcpy(dataPtr, V.dataPtr, (sizeof(double))*index1Size*index2Size*index3Size);
#else
     std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size*index3Size, dataPtr);
#endif
   }

    DoubleVector3d(DoubleVector3d&& V)
    {
      #ifdef VERBOSE_OPS
      std::cout << "Move Copy " << std::endl;
      #endif

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      index3Size   = V.index3Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
      V.index3Size = 0;
    }

       virtual ~DoubleVector3d()
    {
    if(dataPtr != nullptr) delete [] dataPtr;
    }

    void initialize()
    {
    if(dataPtr != nullptr) delete [] dataPtr;
    dataPtr    = nullptr;
    index1Size = 0;
    index2Size = 0;
    index3Size = 0;
    }

    void initialize(long m, long n, long p)
    {
        if((index1Size != m)||(index2Size != n)||(index3Size != p))
        {
            if(dataPtr != nullptr) delete [] dataPtr;
            if(m*n*p > 0)
            {
                dataPtr    = new double[m*n*p];
                index1Size = m;
                index2Size = n;
                index3Size = p;
            }
            else {dataPtr    = nullptr; index1Size = 0; index2Size = 0; index3Size = 0;};
        }
    }

    void initialize(const DoubleVector3d& V)
    {
      if(V.dataPtr == nullptr)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr = nullptr;
      index1Size = 0;
      index2Size = 0;
      index3Size = 0;
      return;
      }

      if((index1Size != V.index1Size)||(index2Size != V.index2Size)||(index3Size != V.index3Size))
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr     = new double[V.index1Size*V.index2Size*V.index3Size];
      index1Size = V.index1Size;
      index2Size = V.index2Size;
      index3Size = V.index3Size;
      }

#ifdef _MSC_VER
     std::memcpy(dataPtr, V.dataPtr, (sizeof(double))*index1Size*index2Size*index3Size);
#else
     std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size*index3Size, dataPtr);
#endif

    }

    void initialize(DoubleVector3d&& V)
    {
      if(V.dataPtr == nullptr)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr = nullptr;
      index1Size = 0;
      index2Size = 0;
      index3Size = 0;
      return;
      }

      if(dataPtr != nullptr) delete [] dataPtr;

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      index3Size   = V.index3Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
      V.index3Size = 0;
    }

    // Assignment operators : Being careful with nullptr instances

    DoubleVector3d& operator=(const DoubleVector3d& V)
    {
      #ifdef VERBOSE_OPS
      std::cout << "Standard Assignment" << std::endl;
      #endif

      if (this != &V)
      {
         if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
         {
         index1Size  = V.index1Size;
         index2Size  = V.index2Size;
         index3Size  = V.index3Size;
         dataPtr     = new double[index1Size*index2Size*index3Size];

#ifdef _MSC_VER
         std::memcpy(dataPtr, V.dataPtr, (sizeof(double))*index1Size*index2Size*index3Size);
#else
         std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size*index3Size, dataPtr);
#endif
         }
         else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; index1Size = 0; index2Size = 0; index3Size = 0;return *this;}
         else
         {
         assert(sizeCheck(this->index1Size,V.index1Size,1));
         assert(sizeCheck(this->index2Size,V.index2Size,2));
         assert(sizeCheck(this->index3Size,V.index3Size,3));
#ifdef _MSC_VER
         std::memcpy(dataPtr, V.dataPtr, (sizeof(double))*index1Size*index2Size*index3Size);
#else
         std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size*index3Size, dataPtr);
#endif
         }
      }
      return *this;
    }

    DoubleVector3d& operator=(DoubleVector3d&& V)
    {
    #ifdef VERBOSE_OPS
    std::cout << "Move Assignment" << std::endl;
    #endif

    if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
    {
      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      index3Size   = V.index3Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
      V.index3Size = 0;
    }
    else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; index1Size = 0; index2Size = 0; index3Size = 0; return *this;}
    else
    {
      assert(sizeCheck(this->index1Size,V.index1Size,1));
      assert(sizeCheck(this->index2Size,V.index2Size,2));
      assert(sizeCheck(this->index3Size,V.index3Size,3));
      // Remove existing data

      delete [] dataPtr;

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      index3Size   = V.index3Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
      V.index3Size = 0;
    }
    return *this;
    }


    inline void transformValues(std::function<double(double)> F)
    {
    for(long k =0; k < index1Size*index2Size*index3Size; k++)
    {dataPtr[k] = F(dataPtr[k]);}
    }

    DoubleVector3d applyFunction(std::function<double(double)> F)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "F(*this)" << std::endl;
    #endif

    DoubleVector3d R(*this);
    R.transformValues(F);
    return std::move(R);
    }


    inline void operator+=(const  DoubleVector3d& D)
    {
      assert(sizeCheck(this->index1Size,D.index1Size,1));
      assert(sizeCheck(this->index2Size,D.index2Size,2));
      assert(sizeCheck(this->index3Size,D.index3Size,3));

      for(long i = 0; i < index1Size*index2Size*index3Size; i++)
      {
            dataPtr[i] += D.dataPtr[i];
      }
    }

    friend DoubleVector3d operator+(const DoubleVector3d& A, const DoubleVector3d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A + &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));
    DoubleVector3d R(A);
    R += B;
    return std::move(R);
    }

    friend DoubleVector3d operator+(const DoubleVector3d& A, DoubleVector3d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A + &&B " << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));

    B += A;
    return std::move(B);
    }

    friend DoubleVector3d operator+(DoubleVector3d&& A, const DoubleVector3d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A +  &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));

    A += B;
    return std::move(A);
    }

    friend DoubleVector3d operator+(DoubleVector3d&& A, DoubleVector3d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A +  &&B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));

    A += B;
    return std::move(A);
    }


    inline void operator-=(const  DoubleVector3d& D)
    {
    assert(sizeCheck(this->index1Size,D.index1Size,1));
    assert(sizeCheck(this->index2Size,D.index2Size,2));
    assert(sizeCheck(this->index3Size,D.index3Size,3));

    for(long i = 0; i < index1Size*index2Size*index3Size; i++)
    {
        dataPtr[i] -= D.dataPtr[i];
    }
    }

    friend DoubleVector3d operator-(const DoubleVector3d& A, const DoubleVector3d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A - &B" << std::endl;
    #endif
    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));

    DoubleVector3d R(A);
    R -= B;
    return std::move(R);
    }

    friend DoubleVector3d operator-(const DoubleVector3d& A, DoubleVector3d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A - &&B " << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));

    for(long i = 0; i < B.index1Size*B.index2Size*B.index3Size; i++)
    {
    B.dataPtr[i] = A.dataPtr[i] - B.dataPtr[i];
    }
    return std::move(B);
    }

    friend DoubleVector3d operator-(DoubleVector3d&& A, const DoubleVector3d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A -  &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));

    A -= B;
    return std::move(A);
    }

    friend DoubleVector3d operator-(DoubleVector3d&& A, DoubleVector3d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A -  &&B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));

    A -= B;
    return std::move(A);
    }

    friend DoubleVector3d operator-(const DoubleVector3d& A)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "-&A" << std::endl;
    #endif

    DoubleVector3d R(A);
    R *= -1.0;
    return std::move(R);
    }

    friend DoubleVector3d operator-(DoubleVector3d&& A)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "- &&A" << std::endl;
    #endif

    A *= -1.0;
    return std::move(A);
    }

    friend DoubleVector3d operator+(const DoubleVector3d& A)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "+&A" << std::endl;
    #endif

    DoubleVector3d R(A);
    return std::move(R);
    }

    friend DoubleVector3d operator+(DoubleVector3d&& A)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "+&&A" << std::endl;
    #endif

    return std::move(A);
    }


    inline void operator*=(const DoubleVector3d& B)
    {
        assert(sizeCheck(index1Size,B.index1Size,1));
        assert(sizeCheck(index2Size,B.index2Size,2));
        assert(sizeCheck(index3Size,B.index3Size,2));
        for(long i = 0; i < index1Size*index2Size*index3Size; i++)
        {
            dataPtr[i] *= B.dataPtr[i];
        }
    }

    friend DoubleVector3d operator*(const DoubleVector3d& A, const DoubleVector3d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A * &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));
    DoubleVector3d R(A);
    R *= B;
    return std::move(R);
    }

    friend DoubleVector3d operator*(const DoubleVector3d& A, DoubleVector3d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A * &&B " << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));

    B *= A;
    return std::move(B);
    }

    friend DoubleVector3d operator*(DoubleVector3d&& A, const DoubleVector3d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A *  &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));

    A *= B;
    return std::move(A);
    }

    friend DoubleVector3d operator*(DoubleVector3d&& A, DoubleVector3d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A *  &&B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));
    assert(A.sizeCheck(A.index3Size,B.index3Size,3));

    A *= B;
    return std::move(A);
    }



    inline void operator*=(const double alpha)
    {
    for(long i = 0; i < index1Size*index2Size*index3Size; i++)
    {
            dataPtr[i] *= alpha;
    }
    }




    friend DoubleVector3d operator*(const double alpha, const DoubleVector3d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "alpha* &B" << std::endl;
    #endif

    DoubleVector3d R(B);
    R *= alpha;
    return std::move(R);
    }

    friend DoubleVector3d operator*(const double alpha, DoubleVector3d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "alpha*+ &&B " << std::endl;
    #endif

    B *= alpha;
    return std::move(B);
    }

    DoubleVector3d operator*(const double alpha) const
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A*alpha" << std::endl;
    #endif

    DoubleVector3d R(*this);
    R *= alpha;
    return std::move(R);
    }

    friend DoubleVector3d operator*(DoubleVector3d&& A,const double alpha)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A * alpha" << std::endl;
    #endif

    A *= alpha;
    return std::move(A);
    }

    inline void operator/=(const double alpha)
    {
    for(long i = 0; i < index1Size*index2Size*index3Size; i++)
    {
            dataPtr[i] /= alpha;
    }
    }


    inline void operator/=(const DoubleVector3d& B)
    {
        assert(sizeCheck(index1Size,B.index1Size,1));
        assert(sizeCheck(index2Size,B.index2Size,2));
        assert(sizeCheck(index3Size,B.index3Size,2));
        for(long i = 0; i < index1Size*index2Size*index3Size; i++)
        {
            dataPtr[i] /= B.dataPtr[i];
        }
    }

    friend DoubleVector3d operator/(DoubleVector3d& A, const double alpha)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A/alpha" << std::endl;
    #endif

    DoubleVector3d R(A);
    R /= alpha;
    return std::move(R);
    }

    friend DoubleVector3d operator/(DoubleVector3d&& A,const double alpha)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A/alpha" << std::endl;
    #endif

    A /= alpha;
    return std::move(A);
    }
/*!  Sets all values of the vector to d. */
    void setToValue(double d)
    {
    for(long i = 0; i < index1Size*index2Size*index3Size; i++)
    {dataPtr[i] =  d;}
    }

    void addValue(double d)
    {
    for(long i = 0; i < index1Size*index2Size*index3Size; i++)
    {dataPtr[i] += d;}
    }


/*!  Standard vector dot product.  */

virtual double dot(const DoubleVector3d& v) const
{
    double dotVal = 0.0;
    for(long i = 0; i < index1Size*index2Size*index3Size; i++)
    {
    dotVal += (dataPtr[i]*v.dataPtr[i]);
    }
    return dotVal;
}

/*!  Maximal absolute value of the elements of the vector. */

double normInf() const
{
    double valMax = 0.0;
    for(long i = 0; i < index1Size*index2Size*index3Size; i++)
    {
    valMax = (valMax > std::abs(dataPtr[i])) ? valMax : std::abs(dataPtr[i]);
    }
    return valMax;
}

/*!  The Euclidean norm of the vector. */

virtual double norm2() const
{
    double val = 0.0;
    for(long i = 0; i < index1Size*index2Size*index3Size; i++)
    {
    val += (dataPtr[i]*dataPtr[i]);
    }
    return std::sqrt(std::abs(val));
}

// Selected BLAS interface

/*! BLAS Euclidean norm of the vector */

virtual double nrm2() const
{
   return norm2();
}

/*! BLAS scalar multiplication  */

void  scal(double alpha)
{
   (*this) *= alpha;
}

/*! BLAS copy : this  <-v   */

void copy(const DoubleVector3d& v)
{
    this->operator=(v);
}

void copy(DoubleVector3d&& v)
{
    this->operator=((DoubleVector3d&&)v);
}


/*! BLAS axpby : this  <- alpha*v + beta*this  */

void axpby(double alpha, const DoubleVector3d& v, double beta)
{
   assert(sizeCheck(this->index1Size,v.index1Size,1));
   assert(sizeCheck(this->index2Size,v.index2Size,2));
   assert(sizeCheck(this->index3Size,v.index3Size,3));

   for(long i = 0; i < index1Size*index2Size*index3Size; i++)
   {dataPtr[i]  = alpha*v.dataPtr[i] + beta*dataPtr[i];}
}

/*! BLAS axpy : this  <- alpha*v + this  */

void axpy(double alpha, const DoubleVector3d& v)
{
   assert(sizeCheck(this->index1Size,v.index1Size,1));
   assert(sizeCheck(this->index2Size,v.index2Size,2));
   assert(sizeCheck(this->index3Size,v.index3Size,3));

   for(long i = 0; i < index1Size*index2Size*index3Size; i++)
   {dataPtr[i]  += alpha*v.dataPtr[i];}
}


/*!  Returns the dimension of the vector */
virtual long getDimension() const
{
    return index1Size*index2Size*index3Size;
}

//
//###################################################################
//      Bounds checking enforced by defining  _DEBUG
//###################################################################
//

#ifdef _DEBUG
    double&  operator()(long i1, long i2, long i3)
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    assert(boundsCheck(i3, 0, index3Size-1,3));
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };

    const double&  operator()(long i1, long i2, long i3) const
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    assert(boundsCheck(i3, 0, index3Size-1,3));
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };
#else
    /*!
    Returns a reference to the element with index (i1,i2,i3) - indexing
    starting at (0,0,0).
    */
    inline double&  operator()(long i1, long i2, long i3)
    {
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };

    /*!
    Returns a reference to the element with index (i1,i2,i3) - indexing
    starting at (0,0,0).
    */
    inline const double&  operator()(long i1, long i2, long i3) const
    {
    return *(dataPtr + i3 + index3Size*(i2 + i1*index2Size));
    };
#endif


/*!  Outputs the vector values to a stream z-slices and using first quadrant indexing for (i,j);  (i,j) = (0,0) in lower left corner. */

    friend std::ostream& operator<<(std::ostream& outStream, const DoubleVector3d&V )
    {
            long i; long j; long k;

            for(k = 0; k < V.index3Size; k++)
            {
            for(j = V.index2Size-1; j >= 0; j--)
            {
            for(i = 0; i <  V.index1Size; i++)
            {
              outStream <<  std::scientific << std::setprecision(3) <<  std::right << std::setw(10) << V(i,j,k) << " ";
            }
            outStream << std::endl;
            }
            outStream << std::endl << std::endl << std::endl;
            }
            return outStream;
    }



    long getIndex1Size()  const {return index1Size;}
    long getIndex2Size()  const {return index2Size;}
    long getIndex3Size()  const {return index3Size;}

    double* getDataPointer(){return dataPtr;}

    const  double* getDataPointer() const {return dataPtr;};

    double*  dataPtr;
    long  index1Size;
    long  index2Size;
    long  index3Size;


//###################################################################
//                      Bounds Checking
//###################################################################
//
#ifdef _DEBUG
        bool boundsCheck(long i, long begin, long end, int coordinate) const
        {
        if((i < begin)||(i  > end))
        {
        std::cerr << "SCC::DoubleVector3d index " << coordinate << " out of bounds " << std::endl;
        std::cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "]" << std::endl;
        return false;
        }
        return true;
        }
#else
        bool boundsCheck(long, long, long, int) const {return true;}
#endif

#ifdef _DEBUG
    bool sizeCheck(long size1, long size2, int /* coordinate */)
    {
    if(size1 != size2)
    {
    std::cerr << "SCC::DoubleVector3d sizes are incompatible : " << size1 << " != " << size2  << std::endl;
    return false;
    }
    return true;
    }

    bool sizeCheck(long size1, long size2, int /* coordinate */) const
    {
    if(size1 != size2)
    {
    std::cerr << "SCC::DoubleVector3d sizes are incompatible : " << size1 << " != " << size2  << std::endl;
    return false;
    }
    return true;
    }
#else
    bool sizeCheck(long, long, int) {return true;}
    bool sizeCheck(long, long, int) const{return true;}
#endif

};
}

#endif /* SCC_DoubleVector3d_ */
