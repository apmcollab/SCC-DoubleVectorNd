/*
 * SCC_DoubleVector1d.h
 *
 *  Created on: Jun 26, 2015
 *      Author: anderson
 *
 *
 *  A minimal 1d double vector class with move semantic implementation.
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

#ifndef SCC_DOUBLE_VECTOR_1D_
#define SCC_DOUBLE_VECTOR_1D_

namespace SCC
{
class DoubleVector1d
{
    public:

    DoubleVector1d()
    {
    dataPtr    = nullptr;
    index1Size = 0;
    }

    DoubleVector1d(long n)
    {
    if(n > 0)
    {
        dataPtr    = new double[n];
        index1Size = n;
    }
    else {dataPtr    = nullptr; index1Size = 0;}
    }

    DoubleVector1d(const DoubleVector1d& V)
    {
      #ifdef VERBOSE_OPS
      std::cout << "Standard Copy " << std::endl;
      #endif

      if(V.dataPtr == nullptr)
      {dataPtr = nullptr; index1Size = 0; return;}

      dataPtr     = new double[V.index1Size];
      index1Size = V.index1Size;

#ifdef _MSC_VER
      std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*index1Size);
#else
     std::copy(V.dataPtr, V.dataPtr + index1Size, dataPtr);
#endif
   }

    DoubleVector1d(DoubleVector1d&& V)
    {
      #ifdef VERBOSE_OPS
      std::cout << "Move Copy " << std::endl;
      #endif

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;;
    }

    // Conversion operators for std::vector<double>

    DoubleVector1d(const std::vector<double>& V)
    {
    	dataPtr    = nullptr;
    	index1Size = 0;
    	initialize(V);
    }

    operator std::vector<double>() const
	{
    std::vector<double> V;
    if(index1Size == 0){return V;}

    V.resize(index1Size);
#ifdef _MSC_VER
    std::memcpy(V.data(),dataPtr,(sizeof(double))*index1Size);
#else
     std::copy(dataPtr, dataPtr + index1Size, V.data());
#endif
    return V;
	}



    virtual ~DoubleVector1d()
    {
    if(dataPtr != nullptr) delete [] dataPtr;
    }

    void initialize()
    {
    if(dataPtr != nullptr) delete [] dataPtr;
    dataPtr    = nullptr;
    index1Size = 0;
    }

    void initialize(long n)
    {
      if(index1Size != n)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      if(n > 0)
      {
        dataPtr    = new double[n];
        index1Size = n;
      }
      else {dataPtr    = nullptr; index1Size = 0;}
      }
    }

    void initialize(const DoubleVector1d& V)
    {
      if(V.dataPtr == nullptr)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr = nullptr;
      index1Size = 0;
      return;
      }

      if(index1Size != V.index1Size)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr     = new double[V.index1Size];
      index1Size = V.index1Size;
      }

#ifdef _MSC_VER
      std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*index1Size);
#else
     std::copy(V.dataPtr, V.dataPtr + index1Size, dataPtr);
#endif

     }

    void initialize(DoubleVector1d&& V)
    {
      if(V.dataPtr == nullptr)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr = nullptr;
      index1Size = 0;
      return;
      }

      if(dataPtr != nullptr) delete [] dataPtr;

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      V.dataPtr    = 0;
      V.index1Size = 0;
    }

    void initialize(const std::vector<double>& V)
    {
      if(V.size() == 0)
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      dataPtr = nullptr;
      index1Size = 0;
      return;
      }

      if(index1Size != (long)V.size())
      {
      if(dataPtr != nullptr) delete [] dataPtr;
      index1Size  = (long)V.size();
      dataPtr     = new double[index1Size];
      }

#ifdef _MSC_VER
      std::memcpy(dataPtr,  V.data(), (sizeof(double))*index1Size);
#else
     std::copy(V.data(), V.data()+ index1Size, dataPtr);
#endif
    }

    bool isNull() const
    {
    if(dataPtr == nullptr){ return true;}
    return false;
    }


    // Assignment operators : Being careful with nullptr instances

    DoubleVector1d& operator=(const DoubleVector1d& V)
    {
      #ifdef VERBOSE_OPS
      std::cout << "Standard Assignment" << std::endl;
      #endif

      if (this != &V)
      {
         if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
         {
         index1Size  = V.index1Size;
         dataPtr     = new double[index1Size];
#ifdef _MSC_VER
         std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*index1Size);
#else
         std::copy(V.dataPtr, V.dataPtr + index1Size, dataPtr);
#endif
         }
         else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; index1Size = 0; return *this;}
         else
         {
         assert(sizeCheck(this->index1Size,V.index1Size));
#ifdef _MSC_VER
         std::memcpy(dataPtr,  V.dataPtr, (sizeof(double))*index1Size);
#else
         std::copy(V.dataPtr, V.dataPtr + index1Size, dataPtr);
#endif
         }
      }
      return *this;
    }

    DoubleVector1d& operator=(DoubleVector1d&& V)
    {
    #ifdef VERBOSE_OPS
    std::cout << "Move Assignment" << std::endl;
    #endif

    if((dataPtr == nullptr)&&(V.dataPtr != nullptr))
    {
      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
    }
    else if((dataPtr == nullptr)&&(V.dataPtr == nullptr)){dataPtr = nullptr; index1Size = 0; return *this;}
    else
    {
      assert(sizeCheck(this->index1Size,V.index1Size));

      // Remove existing data

      delete [] dataPtr;

      dataPtr      = V.dataPtr;
      index1Size   = V.index1Size;
      V.dataPtr    = nullptr;
      V.index1Size = 0;
    }
    return *this;
    }

    inline void transformValues(std::function<double(double)> F)
    {
    for(long k =0; k < index1Size; k++)
    {dataPtr[k] = F(dataPtr[k]);}
    }

    DoubleVector1d applyFunction(std::function<double(double)> F)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "F(*this)" << std::endl;
    #endif

    DoubleVector1d R(*this);
    R.transformValues(F);
    return R;
    }


    inline void operator+=(const  DoubleVector1d& D)
    {
    assert(sizeCheck(this->index1Size,D.index1Size));
    for(long i = 0; i < this->index1Size; i++)
    {
        dataPtr[i] += D.dataPtr[i];
    }
    }


    friend DoubleVector1d operator+(const DoubleVector1d& A, const DoubleVector1d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A + &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size));
    DoubleVector1d R(A);
    R += B;
    return R;
    }

    friend DoubleVector1d operator+(const DoubleVector1d& A, DoubleVector1d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A + &&B " << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size));
    B += A;
    return std::move(B);
    }

    friend DoubleVector1d operator+(DoubleVector1d&& A, const DoubleVector1d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A +  &B" << std::endl;
    #endif

    assert(B.sizeCheck(A.index1Size,B.index1Size));
    A += B;
    return std::move(A);
    }

    friend DoubleVector1d operator+(DoubleVector1d&& A, DoubleVector1d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A +  &&B" << std::endl;
    #endif

    assert(B.sizeCheck(A.index1Size,B.index1Size));
    A += B;
    return std::move(A);
    }

    inline void operator-=(const  DoubleVector1d& D)
    {
        assert(sizeCheck(this->index1Size,D.index1Size));
        for(long i = 0; i < this->index1Size; i++)
        {
            dataPtr[i] -= D.dataPtr[i];
        }
    }

    friend DoubleVector1d operator-(const DoubleVector1d& A, const DoubleVector1d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A - &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size));
    DoubleVector1d R(A);
    R -= B;
    return R;
    }

    friend DoubleVector1d operator-(const DoubleVector1d& A, DoubleVector1d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A - &&B " << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size));
    for(long i = 0; i < B.index1Size; i++)
    {
    B.dataPtr[i] = A.dataPtr[i] - B.dataPtr[i];
    }
    return std::move(B);
    }

    friend DoubleVector1d operator-(DoubleVector1d&& A, const DoubleVector1d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A -  &B" << std::endl;
    #endif

    assert(B.sizeCheck(A.index1Size,B.index1Size));
    A -= B;
    return std::move(A);
    }

    friend DoubleVector1d operator-(DoubleVector1d&& A, DoubleVector1d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A -  &&B" << std::endl;
    #endif

    assert(B.sizeCheck(A.index1Size,B.index1Size));
    A -= B;
    return std::move(A);
    }

    friend DoubleVector1d operator-(const DoubleVector1d& A)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "-&A" << std::endl;
    #endif

    DoubleVector1d R(A);
    R *= -1.0;
    return R;
    }

    friend DoubleVector1d operator-(DoubleVector1d&& A)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "- &&A" << std::endl;
    #endif

    A *= -1.0;
    return std::move(A);
    }

    friend DoubleVector1d operator+(const DoubleVector1d& A)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "+&A" << std::endl;
    #endif

    DoubleVector1d R(A);
    return R;
    }

    friend DoubleVector1d operator+(DoubleVector1d&& A)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "+&&A" << std::endl;
    #endif

    return std::move(A);
    }

    inline void operator*=(const DoubleVector1d& B)
    {
        assert(sizeCheck(index1Size,B.index1Size));
        for(long i = 0; i < this->index1Size; i++)
        {
            dataPtr[i] *= B.dataPtr[i];
        }
    }

    friend DoubleVector1d operator*(const DoubleVector1d& A, const DoubleVector1d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A * &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size));
    DoubleVector1d R(A);
    R *= B;
    return R;
    }

    friend DoubleVector1d operator*(const DoubleVector1d& A, DoubleVector1d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A * &&B " << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size));
    B *= A;
    return std::move(B);
    }

    friend DoubleVector1d operator*(DoubleVector1d&& A, const DoubleVector1d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A *  &B" << std::endl;
    #endif

    assert(B.sizeCheck(A.index1Size,B.index1Size));
    A *= B;
    return std::move(A);
    }

    friend DoubleVector1d operator*(DoubleVector1d&& A, DoubleVector1d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A *  &&B" << std::endl;
    #endif

    assert(B.sizeCheck(A.index1Size,B.index1Size));
    A *= B;
    return std::move(A);
    }

    inline void operator*=(const double alpha)
    {
        for(long i = 0; i < this->index1Size; i++)
        {
            dataPtr[i] *= alpha;
        }
    }

    friend DoubleVector1d operator*(const double alpha, const DoubleVector1d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "alpha* &B" << std::endl;
    #endif

    DoubleVector1d R(B);
    R *= alpha;
    return R;
    }

    friend DoubleVector1d operator*(const double alpha, DoubleVector1d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "alpha*+ &&B " << std::endl;
    #endif

    B *= alpha;
    return std::move(B);
    }

    DoubleVector1d operator*(const double alpha) const
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A*alpha" << std::endl;
    #endif

    DoubleVector1d R(*this);
    R *= alpha;
    return R;
    }

    friend DoubleVector1d operator*(DoubleVector1d&& A,const double alpha)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A * alpha" << std::endl;
    #endif

    A *= alpha;
    return std::move(A);
    }


    inline void operator/=(const DoubleVector1d& B)
    {
        assert(sizeCheck(index1Size,B.index1Size));
        for(long i = 0; i < this->index1Size; i++)
        {
            dataPtr[i] /= B.dataPtr[i];
        }
    }

    inline void operator/=(const double alpha)
    {
    for(long i = 0; i < this->index1Size; i++)
    {
        dataPtr[i] /= alpha;
    }
    }

    friend DoubleVector1d operator/(DoubleVector1d& A, const double alpha)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A/alpha" << std::endl;
    #endif

    DoubleVector1d R(A);
    R /= alpha;
    return R;
    }

    friend DoubleVector1d operator/(DoubleVector1d&& A,const double alpha)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A/alpha" << std::endl;
    #endif

    A /= alpha;
    return std::move(A);
    }

    friend DoubleVector1d operator/(const DoubleVector1d& A, const DoubleVector1d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A * &B" << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size));
    DoubleVector1d R(A);
    R /= B;
    return R;
    }

    friend DoubleVector1d operator/(const DoubleVector1d& A, DoubleVector1d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&A * &&B " << std::endl;
    #endif

    assert(A.sizeCheck(A.index1Size,B.index1Size));
    DoubleVector1d R(A);
    R /= B;
    return R;
    }

    friend DoubleVector1d operator/(DoubleVector1d&& A, const DoubleVector1d& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A *  &B" << std::endl;
    #endif

    assert(B.sizeCheck(A.index1Size,B.index1Size));
    A /= B;
    return std::move(A);
    }

    friend DoubleVector1d operator/(DoubleVector1d&& A, DoubleVector1d&& B)
    {
    #ifdef VERBOSE_OPS
    std::cout  << "&&A *  &&B" << std::endl;
    #endif

    assert(B.sizeCheck(A.index1Size,B.index1Size));
    A /= B;
    return std::move(A);
    }
/*!  Sets all values of the vector to d. */

    void setToValue(double d)
    {
    for(long i = 0; i < index1Size; i++)
    {dataPtr[i] =  d;}
    }

    void addValue(double d)
    {
    for(long i = 0; i < index1Size; i++)
    {dataPtr[i] += d;}
    }


/*!  Standard vector dot product.  */

    virtual double dot(const DoubleVector1d& v) const
    {
    double dotVal = 0.0;
    for(long i = 0; i < index1Size; i++)
    {
    dotVal += (dataPtr[i]*v.dataPtr[i]);
    }
    return dotVal;
    }

/*!  Maximal absolute value of the elements of the vector. */

    virtual double normInf() const
    {
    double valMax = 0.0;
    for(long i = 0; i < index1Size; i++)
    {
    valMax = (valMax > std::abs(dataPtr[i])) ? valMax : std::abs(dataPtr[i]);
    }
    return valMax;
    }

/*!  The Euclidean norm of the vector. */

    virtual double norm2() const
    {
    double val = 0.0;
    for(long i = 0; i < index1Size; i++)
    {
    val += (dataPtr[i]*dataPtr[i]);
    }
    return std::sqrt(std::abs(val));
    }

// Selected Level 1 BLAS interface

/*! BLAS Euclidean norm of the vector */

   virtual double nrm2() const
   {
   return norm2();
   }

/*! BLAS scalar multiplication  */

   virtual void scal(double alpha)
   {
   (*this) *= alpha;
   }

/*! BLAS copy : this  <-v   */

    void copy(const DoubleVector1d& v)
    {
    this->operator=(v);
    }

    void copy(DoubleVector1d&& v)
    {
    this->operator=((DoubleVector1d&&)v);
    }


/*! BLAS axpby : this  <- alpha*v + beta*this  */

   virtual void axpby(double alpha, const DoubleVector1d& v, double beta)
   {
   assert(sizeCheck(this->index1Size,v.index1Size));
   for(long i = 0; i < index1Size; i++)
   {dataPtr[i]  = alpha*v.dataPtr[i] + beta*dataPtr[i];}
   }

/*! BLAS axpy : this  <- alpha*v + this  */

   virtual void axpy(double alpha, const DoubleVector1d& v)
   {
   assert(sizeCheck(this->index1Size,v.index1Size));
   for(long i = 0; i < index1Size; i++)
   {dataPtr[i]  += alpha*v.dataPtr[i];}
   }

/*!  Returns the dimension of the vector */

    virtual long getDimension() const
    {
    return index1Size;
    }

#ifndef NDEBUG
    double&  operator()(long i1)
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    return *(dataPtr +  i1);
    };

    const double&  operator()(long i1) const
    {
    assert(boundsCheck(i1, 0, index1Size-1,1));
    return *(dataPtr +  i1);
    };

#else
    inline double&  operator()(long i1)
    {
    return *(dataPtr + i1);
    };

    inline const double&  operator()(long i1) const
    {
    return *(dataPtr + i1);
    };
#endif


 /*!  Outputs vector values to a stream. */

friend std::ostream& operator<<(std::ostream& outStream, const DoubleVector1d& V)
{

        long i;
        for(i = 0; i <  V.index1Size; i++)
        {
          outStream << std::setprecision(3) <<  std::right << std::setw(10) << V(i) << " ";
          outStream << std::endl;
        }
        return outStream;
}


    long getSize()  const {return index1Size;}

    long getIndex1Size()  const {return index1Size;}

    double* getDataPointer(){return dataPtr;}
    const  double* getDataPointer()  const  {return dataPtr;}

    double*       dataPtr;
    long index1Size;




//###################################################################
//                      Bounds Checking
//###################################################################
//
#ifdef _DEBUG
        bool boundsCheck(long i, long begin, long end, int coordinate) const
        {
        if((i < begin)||(i  > end))
        {
        std::cerr << "SCC::DoubleVector1d index " << coordinate << " out of bounds " << std::endl;
        std::cerr << "Offending index value : " << i << " Acceptable Range [" << begin << "," << end << "] " << std::endl;
        return false;
        }
        return true;
        }
#else
        bool boundsCheck(long, long, long, int) const {return true;}
#endif

#ifdef _DEBUG
    bool sizeCheck(long size1, long size2)
    {
    if(size1 != size2)
    {
    std::cerr << "SCC::DoubleVector1d sizes are incompatible : " << size1 << " != " << size2 << "  \n";
    return false;
    }
    return true;
    }

    bool sizeCheck(long size1, long size2) const
    {
    if(size1 != size2)
    {
    std::cerr << "SCC::DoubleVector1d sizes are incompatible : " << size1 << " != " << size2 << "  \n";
    return false;
    }
    return true;
    }
#else
    bool sizeCheck(long, long) {return true;}
    bool sizeCheck(long, long) const{return true;}
#endif

};
}

#endif /* SCC_DoubleVector1d_ */
