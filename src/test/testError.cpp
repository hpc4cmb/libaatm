/*******************************************************************************
 * ALMA - Atacama Large Millimeter Array
 * (c) Instituto de Estructura de la Materia, 2011
 * (in the framework of the ALMA collaboration).
 * All rights reserved.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *******************************************************************************/

#include <iostream>
#include <vector>
#include <sstream>

#include "ATMError.h"

using namespace std;
using namespace atm;


const int DefaultSize = 10;

// define the exception class
class xBoundary
{
public:
  xBoundary() {}
  ~xBoundary() {}
private:
};



class Array
{
public:
  // constructors
  Array(int itsSize = DefaultSize);
  Array(const Array &rhs);
  ~Array() { delete [] pType;}

  // operators
  Array& operator=(const Array&);
  int& operator[](int offSet);
  const int& operator[](int offSet) const;

  // accessors
  int GetitsSize() const { return itsSize; }

  // friend function
  friend ostream& operator<< (ostream&, const Array&);

private:
  int *pType;
  int  itsSize;
};


Array::Array(int size):
  itsSize(size)
{
  pType = new int[size];
  for (int i = 0; i<size; i++)
    pType[i] = 0;
}


Array& Array::operator=(const Array &rhs)
{
  if (this == &rhs)
    return *this;
  delete [] pType;
  itsSize = rhs.GetitsSize();
  pType = new int[itsSize];
  for (int i = 0; i<itsSize; i++)
    pType[i] = rhs[i];
  return *this;
}

Array::Array(const Array &rhs)
{
  itsSize = rhs.GetitsSize();
  pType = new int[itsSize];
  for (int i = 0; i<itsSize; i++)
    pType[i] = rhs[i];
}


int& Array::operator[](int offSet)
{
  // int size = GetitsSize();
  if (offSet >= 0 && offSet < GetitsSize())
    return pType[offSet];
  //throw xBoundary();
  //Error(SERIOUS,"Out of range");
  ostringstream msg;
  msg << "Out of range; index must be in the range [" << 0 << ", " << GetitsSize() - 1 << "]";
  Error(SERIOUS, msg.str());
  return pType[offSet]; // to appease MSC!

}


const int& Array::operator[](int offSet) const
{
  //  int mysize = GetitsSize();
  if (offSet >= 0 && offSet < GetitsSize())
    return pType[offSet];
  throw xBoundary();
  return pType[offSet]; // to appease MSC!
}

ostream& operator<< (ostream& output, const Array& theArray)
{
  for (int i = 0; i<theArray.GetitsSize(); i++)
    output << "[" << i << "] " << theArray[i] << endl;
  return output;
}

// testbed

int main()
{
  cout<<"The default acceptable error level is "<<Error::getAcceptableLevelToString()<<endl;
  Array intArray(20);
  try
    {
      for (int j = 0; j< 100; j++)
	{
	  intArray[j] = j;
	  cout << "intArray[" << j << "] okay..." << endl;
	}
    }
  catch (Error){
      cout << "This level of the error, "<<Error::getLevelToString()
	   << ", is not acceptable!" << endl;
      cout << "Some action should be taken here!" << endl;}
  catch (xBoundary){
      cout << "Unable to process your input!";}

  if(Error::getLevel()<Error::getAcceptableLevel())
    cout << "The level of the error, "<<Error::getLevelToString()
	 << ", is acceptable; we can proceed" << endl;

  cout<<"We now give less tolerance on the acceptable level of an error."<<endl;
  Error::setAcceptableLevel(FATAL);
  cout<<"The acceptable error level has been set to "<<Error::getAcceptableLevelToString()<<endl;
  cout<<"Test: the same error is now generated a second time"<<endl;
  try
    {
      for (int j = 0; j< 100; j++)
	{
	  intArray[j] = j;
	  cout << "intArray[" << j << "] okay..." << endl;
	}
    }
  catch (Error)
    {
      cout << "An error has been thrown because its level ("
	   << Error::getLevelToString()<<") is severe enough relative"<< endl;
      cout <<"to the current tolerence ("<<Error::getAcceptableLevelToString()<<")"<<endl;
    }

  cout<<"Test completed successfuly"<<endl;
  return 0;

}
