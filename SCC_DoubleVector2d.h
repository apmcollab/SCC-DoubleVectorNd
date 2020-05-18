/*
 * SCC_DoubleVector2d.h
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

#ifdef  _DEBUG
#include <cstdio>
#else
#ifndef NDEBUG
#define NDEBUG
#endif
#endif
#include <cassert>

#undef _VERBOSE_OPS_

#ifndef _SCC_DoubleVector2d_
#define _SCC_DoubleVector2d_

namespace SCC
{
class DoubleVector2d
{
    public:

    DoubleVector2d()
    {
    dataPtr    = nullptr;
    index1Size = 0;
    index2Size = 0;
    }

    explicit DoubleVector2d(long m, long n)
    {
    if(m*n > 0)
    {
    dataPtr    = new double[m*n];
    index1Size = m;
    index2Size = n;
    }
    else {dataPtr    = nullptr; index1Size = 0; index2Size = 0;}
    }

   DoubleVector2d(const DoubleVector2d& V)
   {
      #ifdef _VERBOSE_OPS_
      std::cout << "Standard Copy " << std::endl;
      #endif

      if(V.dataPtr == nullptr) {dataPtr = nullptr; index1Size = 0; index2Size = 0; return;}

      dataPtr     = new double[V.index1Size*V.index2Size];
      index1Size = V.index1Size;
      index2Size = V.index2Size;

#ifdef _MSC_VER
     std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*index1Size*index2Size);
#else
     std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size, dataPtr);
#endif

   }

    DoubleVector2d(DoubleVector2d&& V)
    {
      #ifdef _VERBOSE_OPS_
      std::cout << "Move Copy " << std::endl;
      #endif

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
    }

       virtual ~DoubleVector2d()
    {
    if(dataPtr != nullptr) delete [] dataPtr;
    }

    void initialize()
    {
    if(dataPtr != nullptr) delete [] dataPtr;
    dataPtr    = nullptr;
    index1Size = 0;
    index2Size = 0;
    }

    void initialize(long m, long n)
    {
        if((index1Size != m)||(index2Size != n))
        {
            if(dataPtr != nullptr) delete [] dataPtr;
            if(m*n > 0)
            {
                dataPtr    = new double[m*n];
                index1Size = m;
                index2Size = n;
            }
            else {dataPtr    = nullptr; index1Size = 0; index2Size = 0;}
        }
    }

    void initialize(const DoubleVector2d& V)
    {
      if(V.dataPtr == nullptr)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr = nullptr;
      index1Size = 0;
      index2Size = 0;
      return;
      }

      if((index1Size != V.index1Size)||(index2Size != V.index2Size))
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr     = new double[V.index1Size*V.index2Size];
      index1Size = V.index1Size;
      index2Size = V.index2Size;
      }

#ifdef _MSC_VER
     std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*index1Size*index2Size);
#else
     std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size, dataPtr);
#endif

    }

    void initialize(DoubleVector2d&& V)
    {
      if(V.dataPtr == nullptr)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr = nullptr;
      index1Size = 0;
      index2Size = 0;
      return;
      }

      if(dataPtr != nullptr) delete [] dataPtr;

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
    }

    // Assignment operators : Being careful with nullptr instances

    DoubleVector2d& operator=(const DoubleVector2d& V)
    {
      #ifdef _VERBOSE_OPS_
      std::cout << "Standard Assignment" << std::endl;
      #endif

      if (this != &V)
      {
         if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
         {
         index1Size  = V.index1Size;
         index2Size  = V.index2Size;
         dataPtr     = new double[index1Size*index2Size];

#ifdef _MSC_VER
         std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*index1Size*index2Size);
#else
         std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size, dataPtr);
#endif

         }
         else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; index1Size = 0; index2Size = 0; return *this;}
         else
         {
         assert(sizeCheck(this->index1Size,V.index1Size,1));
         assert(sizeCheck(this->index2Size,V.index2Size,2));
#ifdef _MSC_VER
         std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*index1Size*index2Size);
#else
         std::copy(V.dataPtr, V.dataPtr + index1Size*index2Size, dataPtr);
#endif
         }
      }
      return *this;
    }

    DoubleVector2d& operator=(DoubleVector2d&& V)
    {
    #ifdef _VERBOSE_OPS_
    std::cout << "Move Assignment" << std::endl;
    #endif

    if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
    {
      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
    }
    else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; index1Size = 0; index2Size = 0; return *this;}
    else
    {
      assert(sizeCheck(this->index1Size,V.index1Size,1));
      assert(sizeCheck(this->index2Size,V.index2Size,2));

      // Remove existing data

      delete [] dataPtr;

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      index2Size   = V.index2Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
      V.index2Size = 0;
    }
    return *this;
    }

    inline void transformValues(std::function<double(double)> F)
    {
    for(long k =0; k < index1Size*index2Size; k++)
    {dataPtr[k] = F(dataPtr[k]);}
    }

    DoubleVector2d applyFunction(std::function<double(double)> F)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "F(*this)" << std::endl;
    #endif

    DoubleVector2d R(*this);
    R.transformValues(F);
    return std::move(R);
    }


    inline void operator+=(const  DoubleVector2d& D)
    {
      assert(sizeCheck(this->index1Size,D.index1Size,1));
      assert(sizeCheck(this->index2Size,D.index2Size,2));
      for(long i = 0; i < index1Size*index2Size; i++)
      {
            dataPtr[i] += D.dataPtr[i];
      }
    }

    friend DoubleVector2d operator+(const DoubleVector2d& A, const DoubleVector2d& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&A + &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    DoubleVector2d R(A);
    R += B;
    return std::move(R);
    }

    friend DoubleVector2d operator+(const DoubleVector2d& A, DoubleVector2d&& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&A + &&B " << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    B += A;
    return std::move(B);
    }

    friend DoubleVector2d operator+(DoubleVector2d&& A, const DoubleVector2d& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&&A +  &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    A += B;
    return std::move(A);
    }

    friend DoubleVector2d operator+(DoubleVector2d&& A, DoubleVector2d&& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&&A +  &&B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    A += B;
    return std::move(A);
    }


    inline void operator-=(const  DoubleVector2d& D)
    {
    assert(sizeCheck(this->index1Size,D.index1Size,1));
    assert(sizeCheck(this->index2Size,D.index2Size,2));

    for(long i = 0; i < index1Size*index2Size; i++)
    {
        dataPtr[i] -= D.dataPtr[i];
    }
    }

    friend DoubleVector2d operator-(const DoubleVector2d& A, const DoubleVector2d& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&A - &B" << std::endl;
    #endif
    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    DoubleVector2d R(A);
    R -= B;
    return std::move(R);
    }

    friend DoubleVector2d operator-(const DoubleVector2d& A, DoubleVector2d&& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&A - &&B " << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    for(long i = 0; i < B.index1Size*B.index2Size; i++)
    {
    B.dataPtr[i] = A.dataPtr[i] - B.dataPtr[i];
    }
    return std::move(B);
    }

    friend DoubleVector2d operator-(DoubleVector2d&& A, const DoubleVector2d& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&&A -  &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    A -= B;
    return std::move(A);
    }

    friend DoubleVector2d operator-(DoubleVector2d&& A, DoubleVector2d&& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&&A -  &&B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    A -= B;
    return std::move(A);
    }

    friend DoubleVector2d operator-(const DoubleVector2d& A)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "-&A" << std::endl;
    #endif

    DoubleVector2d R(A);
    R *= -1.0;
    return std::move(R);
    }

    friend DoubleVector2d operator-(DoubleVector2d&& A)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "- &&A" << std::endl;
    #endif

    A *= -1.0;
    return std::move(A);
    }

    friend DoubleVector2d operator+(const DoubleVector2d& A)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "+&A" << std::endl;
    #endif

    DoubleVector2d R(A);
    return std::move(R);
    }

    friend DoubleVector2d operator+(DoubleVector2d&& A)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "+&&A" << std::endl;
    #endif

    return std::move(A);
    }


    inline void operator*=(const DoubleVector2d& B)
    {
        assert(sizeCheck(index1Size,B.index1Size,1));
        assert(sizeCheck(index2Size,B.index2Size,2));
        for(long i = 0; i < index1Size*index2Size; i++)
        {
            dataPtr[i] *= B.dataPtr[i];
        }
    }

    friend DoubleVector2d operator*(const DoubleVector2d& A, const DoubleVector2d& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&A * &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    DoubleVector2d R(A);
    R *= B;
    return std::move(R);
    }

    friend DoubleVector2d operator*(const DoubleVector2d& A, DoubleVector2d&& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&A * &&B " << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    B *= A;
    return std::move(B);
    }

    friend DoubleVector2d operator*(DoubleVector2d&& A, const DoubleVector2d& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&&A *  &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    A *= B;
    return std::move(A);
    }

    friend DoubleVector2d operator*(DoubleVector2d&& A, DoubleVector2d&& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&&A *  &&B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size,1));
    assert(A.sizeCheck(A.index2Size,B.index2Size,2));

    A *= B;
    return std::move(A);
    }



    inline void operator*=(const double alpha)
    {
    for(long i = 0; i < index1Size*index2Size; i++)
    {
            dataPtr[i] *= alpha;
    }
    }

    friend DoubleVector2d operator*(const double alpha, const DoubleVector2d& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "alpha* &B" << std::endl;
    #endif

    DoubleVector2d R(B);
    R *= alpha;
    return std::move(R);
    }

    friend DoubleVector2d operator*(const double alpha, DoubleVector2d&& B)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "alpha*+ &&B " << std::endl;
    #endif

    B *= alpha;
    return std::move(B);
    }

    DoubleVector2d operator*(const double alpha) const
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&A*alpha" << std::endl;
    #endif

    DoubleVector2d R(*this);
    R *= alpha;
    return std::move(R);
    }

    friend DoubleVector2d operator*(DoubleVector2d&& A,const double alpha)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&&A * alpha" << std::endl;
    #endif

    A *= alpha;
    return std::move(A);
    }


    inline void operator/=(const DoubleVector2d& B)
    {
        assert(sizeCheck(index1Size,B.index1Size,1));
        assert(sizeCheck(index2Size,B.index2Size,2));
        for(long i = 0; i < index1Size*index2Size; i++)
        {
            dataPtr[i] /= B.dataPtr[i];
        }
    }


    inline void operator/=(const double alpha)
    {
    for(long i = 0; i < index1Size*index2Size; i++)
    {
            dataPtr[i] /= alpha;
    }
    }

    friend DoubleVector2d operator/(DoubleVector2d& A, const double alpha)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&A/alpha" << std::endl;
    #endif

    DoubleVector2d R(A);
    R /= alpha;
    return std::move(R);
    }

    friend DoubleVector2d operator/(DoubleVector2d&& A,const double alpha)
    {
    #ifdef _VERBOSE_OPS_
    std::cout  << "&&A/alpha" << std::endl;
    #endif

    A /= alpha;
    return std::move(A);
    }
/*!  Sets all values of the vector to d. */
    void setToValue(double d)
    {
    for(long i = 0; i < index1Size*index2Size; i++)
    {dataPtr[i] =  d;}
    }

    void addValue(double d)
    {
    for(long i = 0; i < index1Size*index2Size; i++)
    {dataPtr[i] += d;}
    }


/*!  Standard vector dot product.  */

    double dot(const DoubleVector2d& v) const
    {
    double dotVal = 0.0;
    for(long i = 0; i < index1Size*index2Size; i++)
    {
    dotVal += (dataPtr[i]*v.dataPtr[i]);
    }
    return dotVal;
    }

/*!  Maximal absolute value of the elements of the vector. */

    double normInf() const
    {
    double valMax = 0.0;
    for(long i = 0; i < index1Size*index2Size; i++)
    {
    valMax = (valMax > std::abs(dataPtr[i])) ? valMax : std::abs(dataPtr[i]);
    }
    return valMax;
    }

/*!  The Euclidean norm of the vector. */

    double norm2() const
    {
    double val = 0.0;
    for(long i = 0; i < index1Size*index2Size; i++)
    {
    val += (dataPtr[i]*dataPtr[i]);
    }
    return std::sqrt(std::abs(val));
}

// Selected BLAS interface

/*! BLAS Euclidean norm of the vector */

double nrm2() const
{
   return norm2();
}

/*! BLAS scalar multiplication  */

void  scal(double alpha)
{
   (*this) *= alpha;
}

/*! BLAS copy : this  <-v   */

void copy(const DoubleVector2d& v)
{
    this->operator=(v);
}

void copy(DoubleVector2d&& v)
{
    this->operator=((DoubleVector2d&&)v);
}

/*! BLAS axpby : this  <- alpha*v + beta*this  */

void axpby(double alpha, const DoubleVector2d& v, double beta)
{
   assert(sizeCheck(this->index1Size,v.index1Size,1));
   assert(sizeCheck(this->index2Size,v.index2Size,2));

   for(long i = 0; i < index1Size*index2Size; i++)
   {dataPtr[i]  = alpha*v.dataPtr[i] + beta*dataPtr[i];}
}

/*! BLAS axpy : this  <- alpha*v + this  */

void axpy(double alpha, const DoubleVector2d& v)
{
   assert(sizeCheck(this->index1Size,v.index1Size,1));
   assert(sizeCheck(this->index2Size,v.index2Size,2));

   for(long i = 0; i < index1Size*index2Size; i++)
   {dataPtr[i]  += alpha*v.dataPtr[i];}
}


/*!  Returns the dimension of the vector */

virtual long getDimension()
{
    return index1Size*index2Size;
}


//
//###################################################################
//      Bounds checking enforced by defining  _DEBUG
//###################################################################
//

#ifdef _DEBUG
    double&  operator()(long i1, long i2)
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    return *(dataPtr +  i2 + i1*index2Size);
    };

    const double&  operator()(long i1, long i2) const
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    assert(boundsCheck(i2, 0, index2Size-1,2));
    return *(dataPtr +   i2  + i1*index2Size);
    };
#else
    /*!
    Returns a reference to the element with index (i1,i2) - indexing
    starting at (0,0).
    */
    inline double&  operator()(long i1, long i2)
    {
    return *(dataPtr +  i2 + i1*index2Size);
    };

    /*!
    Returns a reference to the element with index (i1,i2) - indexing
    starting at (0,0).
     */
    inline const double&  operator()(long i1, long i2) const
    {
    return *(dataPtr +  i2  + i1*index2Size);
    };
#endif


/*!  Outputs the vector values to a stream using first quadrant indexing: (i,j) = (0,0) in lower left corner. */

friend std::ostream& operator<<(std::ostream& outStream, const DoubleVector2d&V )
{
        long i; long j;

        for(j = V.index2Size-1; j >= 0; j--)
        {
        for(i = 0; i <  V.index1Size; i++)
        {
          outStream <<   std::scientific << std::setprecision(3) <<  std::right << std::setw(10) << V(i,j) << " ";
        }
        outStream << std::endl;
        }
        return outStream;
}



    long getIndex1Size()  const {return index1Size;}
    long getIndex2Size()  const {return index2Size;}


    double* getDataPointer(){return dataPtr;}
    const  double* getDataPointer() const {return dataPtr;};

    double*       dataPtr;
    long index1Size;
    long index2Size;


//###################################################################
//                      Bounds Checking
//###################################################################
//
#ifdef _DEBUG
        bool boundsCheck(long i, long begin, long end, int coordinate) const
        {
        if((i < begin)||(i  > end))
        {
        std::cerr << "SCC::DoubleVector2d index " << coordinate << " out of bounds " << std::endl;
        std::cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "]" << std::endl;
        return false;
        }
        return true;
        }
#else
        bool boundsCheck(long, long, long, int) const {return true;}
#endif

#ifdef _DEBUG
    bool sizeCheck(long size1, long size2, int coordinate)
    {
    if(size1 != size2)
    {
    std::cerr << "SCC::DoubleVector2d sizes are incompatible : " << size1 << " != " << size2;
    return false;
    }
    return true;
    }

    bool sizeCheck(long size1, long size2, int coordinate) const
    {
    if(size1 != size2)
    {
    std::cerr << "SCC::DoubleVector2d sizes are incompatible : " << size1 << " != " << size2;
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

#endif /* SCC_DoubleVector2d_ */
