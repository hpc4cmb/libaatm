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
#include <string>
using namespace std;

#include <limits>
#include <math.h>


#include "ATMSpectralGrid.h"



#include <iostream>
using namespace atm;
  /** \brief A C++ main code to test the <a href="classatm_1_1SpectralGrid.html">SpectralGrid</a> Class
   *
   *   Test 1 is structured as follows:
   *         - Creates a pointer called "sgPtr1" of objects belonging to the class <a href="classatm_1_1SpectralGrid.html">SpectralGrid</a>
   *         - Initializes that pointer with a first object created with the <a href="classatm_1_1SpectralGrid.html#a1">this</a> constructor.
   *         - The pointer has then one element . On this element, several things are checked using the following operators of the class:
   *           <a href="classatm_1_1SpectralGrid.html#z14_6">getRefFreq(unsigned int spwId)</a>,
   *           <a href="classatm_1_1SpectralGrid.html#z14_8">getChanSep(unsigned int spwId)</a>,
   *           <a href="classatm_1_1SpectralGrid.html#z14_0">getNumSpectralWindow()</a>, and
   *           <a href="classatm_1_1SpectralGrid.html#z14_2">getNumChan(unsigned int spwId)</a>.
   *         - An rrror message is expected when trying getNumChan(1) because there is no spectral window number 1.
   *         - Finally, the operators <a href="classatm_1_1SpectralGrid.html#z14_22">isRegular(unsigned int spwId)</a>,
   *           and <a href="classatm_1_1SpectralGrid.html#z14_25">getAssocSpwId(unsigned int spwId)</a> are tested.
   *
   * The ouput of this test should be as follows:
   *
   * <b>
   * SpectralGridTest: Test 1: <br>
   * SpectralGridTest: Create a pointer with first element construted with SpectralGrid(usigned int numChan, unsigned int refChan, Frequency refFreq, Frequency chanSep):<br>
   * SpectralGridTest: Number of channels retrieved:  64 (Value entered to constructor:64)<br>
   * SpectralGridTest: Reference frequency retrieved: 9e+10 Hz;  SpectralGridTest:  Input:90 GHz<br>
   * SpectralGridTest: Reference frequency retrieved: 90GHz  SpectralGridTest: Input:90 GHz<br>
   * SpectralGridTest: Reference channel retrieved:   32 (Value entered to constructor:32)<br>
   * SpectralGridTest: Channel separation retrieved:  1e+07 Hz;  SpectralGridTest: Input:0.01 GHz<br>
   * SpectralGridTest: Channel separation retrieved:  10000 kHz  SpectralGridTest: Input:0.01 GHz<br>
   * SpectralGridTest: Number of spectral windows: 1<br>
   * SpectralGridTest: Number of channels for spectral window identifier 0: 64<br>
   * SpectralGrid: ERROR: 1 is a wrong spectral window identifier<br>
   * SpectralGridTest: Number of channels for spectral window identifier 1: 0<br>
   * SpectralGridTest: the first spectral window is regularily sampled<br>
   * SpectralGridTest: First spectral window has no associated spectral windows<br>
   * SpectralGridTest: End of Test 1<br>
   * </b>
   *
   * Test 2 starts from the pointer created in Test 1:
   *       - A new spectral window is added with 128 channels, channel number 32 as the reference channel, reference frequency at 215 GHz, and regular channel separation of 0.02 GHz.
   *       This new spectral window is added using <a href="classatm_1_1SpectralGrid.html#a11">add(unsigned int numChan, unsigned int refChan, Frequency refFreq, Frequency chanSep)</a>
   *       - The next step is to verify that nothing has change for spectral window number 0 from Test 1.
   *       - Finally, spectral window number 1 is tested similarly to spectral window #0.
   *       - It is also checked that spectral windows with numbers > 1 do not exist.
   *
   * The ouput of Test 2 should be as follows:
   *
   * <b>
   * SpectralGridTest: Test 2<br>
   * SpectralGridTest: New spectral window using add(unsigned int numChan, unsigned int refChan, Frequency refFreq, Frequency chanSep):<br>
   * SpectralGridTest: A new spectral window has been appended and got the identifier number:1<br>
   * SpectralGridTest: Number of spectral windows: 2<br>
   * SpectralGridTest: Number of channels retrieved for spwId 0: 64<br>
   * SpectralGridTest: Reference frequency retrieved: 90 GHz<br>
   * SpectralGridTest: Channel separation retrieved:  0.01 GHz GHz<br>
   * SpectralGridTest: Number of channels retrieved for spwId 1: 128<br>
   * SpectralGridTest: Reference frequency retrieved: 215 GHz<br>
   * SpectralGridTest: Channel separation retrieved:  0.02 GHz GHz<br>
   * SpectralGridTest: the spectral window with id 1 is regularily sampled<br>
   * SpectralGridTest: Number of spectral windows: 2<br>
   * SpectralGridTest: As expected this spectral window with spwid=1 has no sideband specification<br>
   * SpectralGrid: ERROR: 10 is a wrong spectral window identifier<br>
   * SpectralGridTest: Spectral window with id=10 does not exist!<br>
   * SpectralGridTest: End of Test 2<br>
   * </b>
   */

int main()
{


  unsigned int     numChan         = 64;
  unsigned int     refChan         = 32;

  Frequency myRefFreq(90.0,"GHz");
  Frequency myChanSep(0.01,"GHz");

  SpectralGrid* sgPtr1;

  cout << " SpectralGridTest: Test 1:" <<endl;
  cout << " SpectralGridTest: Create a pointer with first element construted with SpectralGrid(usigned int numChan, unsigned int refChan, Frequency refFreq, Frequency chanSep):" << endl;
  sgPtr1 = new SpectralGrid(numChan, refChan, myRefFreq, myChanSep);
  cout << " SpectralGridTest: Number of channels retrieved:  " << sgPtr1->getNumChan() << " (Value entered to constructor:" << numChan << ")" << endl;
  cout << " SpectralGridTest: Reference frequency retrieved: "
       << sgPtr1->getRefFreq().get() << " Hz; "
       << " SpectralGridTest:  Input:" << myRefFreq.get("GHz") << " GHz" << endl;
  cout << " SpectralGridTest: Reference frequency retrieved: "
       << sgPtr1->getRefFreq().get("GHz")<< "GHz "
       <<" SpectralGridTest: Input:" << myRefFreq.get("GHz") << " GHz" << endl;
  cout << " SpectralGridTest: Reference channel retrieved:   " << sgPtr1->getRefChan() << " (Value entered to constructor:" << refChan << ")" << endl;
  cout << " SpectralGridTest: Channel separation retrieved:  "
       << sgPtr1->getChanSep().get() << " Hz; "
       << " SpectralGridTest: Input:" << myChanSep.get("GHz") << " GHz" << endl;
  cout << " SpectralGridTest: Channel separation retrieved:  "
       << sgPtr1->getChanSep().get("kHz") << " kHz "
       << " SpectralGridTest: Input:" << myChanSep.get("GHz") << " GHz" << endl;
  cout << " SpectralGridTest: Number of spectral windows: "
       << sgPtr1->getNumSpectralWindow() << endl;
  cout << " SpectralGridTest: Number of channels for spectral window identifier 0: "
       << sgPtr1->getNumChan(0) << endl;
  cout << " SpectralGridTest: Number of channels for spectral window identifier 1: "
       << sgPtr1->getNumChan(1) << endl;
  if(sgPtr1->isRegular(0)){
    cout << " SpectralGridTest: the first spectral window is regularily sampled" << endl;
  }else{
    cout << " SpectralGridTest: the first spectral window is not regularily sampled" << endl;
  }
  if(sgPtr1->getAssocSpwId(0).size()==0){
    cout<<" SpectralGridTest: First spectral window has no associated spectral windows" << endl;
  }else{
    cout<<" SpectralGridTest: First spectral window has associated spectral windows "<<endl;
  }
  cout << " SpectralGridTest: End of Test 1" <<endl;
  cout << "    " << endl;




  cout << " SpectralGridTest: Test 2" << endl;
   cout << " SpectralGridTest: New spectral window using add(unsigned int numChan, unsigned int refChan, Frequency refFreq, Frequency chanSep):" << endl;
  unsigned int numChan1         = 128;
  unsigned int refChan1         = 32;
  Frequency myNewRefFreq(215.0,"GHz");
  Frequency myNewChanSep(0.02,"GHz");
  unsigned int spwId = sgPtr1->add( numChan1, refChan1, myNewRefFreq, myNewChanSep);
  cout << " SpectralGridTest: A new spectral window has been appended and got the identifier number:" << spwId << endl;
  cout << " SpectralGridTest: Number of spectral windows: "
       << sgPtr1->getNumSpectralWindow() << endl;

  spwId=0;
  cout << " SpectralGridTest: Number of channels retrieved for spwId "<<spwId<<": "
       << sgPtr1->getNumChan(spwId) << endl;
  cout << " SpectralGridTest: Reference frequency retrieved: "
       << sgPtr1->getRefFreq(spwId).get("GHz") << " GHz" << endl;
  cout << " SpectralGridTest: Channel separation retrieved:  "
       << sgPtr1->getChanSep(spwId).get("GHz") << " GHz " << "GHz" << endl;

  spwId=1;
  cout << " SpectralGridTest: Number of channels retrieved for spwId "<<spwId<<": "
       << sgPtr1->getNumChan(spwId) << endl;
  cout << " SpectralGridTest: Reference frequency retrieved: "
       << sgPtr1->getRefFreq(spwId).get("GHz") << " GHz" << endl;
  cout << " SpectralGridTest: Channel separation retrieved:  "
       << sgPtr1->getChanSep(spwId).get("GHz") << " GHz " << "GHz" << endl;

  if(sgPtr1->isRegular(spwId)){
    cout << " SpectralGridTest: the spectral window with id "<<spwId<<" is regularily sampled" << endl;
  }else{
    cout << " SpectralGridTest: the spectral window with id "<<spwId<<" is not regularily sampled" << endl;
  }

  cout << " SpectralGridTest: Number of spectral windows: "
       << sgPtr1->getNumSpectralWindow() << endl;

  if(sgPtr1->getSideband(spwId).size()==0)
    cout << " SpectralGridTest: As expected this spectral window with spwid="<<spwId
	 << " has no sideband specification" << endl;
  unsigned int id=10;
  if(sgPtr1->getSideband(id).size()==0)
    cout << " SpectralGridTest: Spectral window with id="<<id
	 << " does not exist!" << endl;
  cout << " SpectralGridTest: End of Test 2" <<endl;







  cout << " SpectralGridTest: Channel frequency and number for the first spectral window: " << endl;
  double chFreq[sgPtr1->getNumChan()];           // one dynamic alloc
  vector<double> chanFreq;
  chanFreq.reserve(sgPtr1->getNumChan());        // a more versatil dynamic alloc (allowing eg resizing)

  for(int i=0; i<(int)sgPtr1->getNumChan(); i++){
    chanFreq[i] = sgPtr1->getChanFreq(i).get();
    chFreq[i] = chanFreq[i];
    cout << "SpectralGridTest: " << i << " channel: " << i-(int)refChan+1 << " freq: " << chanFreq[i] << endl;
  }
  cout << endl;

  delete sgPtr1; sgPtr1=0;



  double  refFreq         = 90.0E9;

  SpectralGrid* sgPtr2;
  cout << " SpectralGridTest: Test 2:" <<endl;
  cout << " SpectralGridTest: Build using SpectralGrid( unsigned int numChan, unsigned int refChan, double* chFreq, string units):" << endl;
  sgPtr2 = new SpectralGrid( numChan, refChan, chFreq, "Hz");
  cout << " SpectralGridTest: Number of channels retrieved: " << sgPtr2->getNumChan()      << "    Input:  " << numChan << endl;
  cout << " SpectralGridTest: Reference frequency retrieved:" << sgPtr2->getRefFreq().get()      << "Hz  Initial: none" << endl;
  cout << " SpectralGridTest: Reference frequency retrieved:" << sgPtr2->getRefFreq().get("MHz") << "MHz Initial: none" << endl;
  cout << " SpectralGridTest: Reference channel retrieved:  " << sgPtr2->getRefChan()      << "    Input:  " << refChan << endl;
  cout << " SpectralGridTest: Channel separation retrieved: " << sgPtr2->getChanSep().get()      << "Hz  Initial: none" << endl;
  cout << " SpectralGridTest: Channel separation retrieved: " << sgPtr2->getChanSep().get("MHz") << "MHz Initial: none" << endl;
  if(sgPtr2->isRegular()){
    cout << " SpectralGridTest: the first spectral window with id 0 is regularily sampled" << endl;
  }else{
    cout << " SpectralGridTest: the first spectral window with id 0 is not regularily sampled" << endl;
  }
  chFreq[sgPtr2->getNumChan()/4]=chFreq[sgPtr2->getNumChan()/4]+1.;
  cout << " SpectralGridTest: Add a second irregular spectral window using add( unsigned int numChan, unsigned int refChan, double* chFreq, string units):" << endl;
  sgPtr2->add( numChan, refChan, chFreq, "Hz");
  if(sgPtr2->isRegular()){
    cout << " SpectralGridTest: the first spectral window with id 0 is regularily sampled as expected" << endl;
  }else{
    cout << " SpectralGridTest: the first spectral window with id 0 is not regularily sampled ==> ERROR in the code" << endl;
  }
  if(sgPtr2->isRegular(spwId)){
    cout << " SpectralGridTest: the spectral window with id "<<spwId<<" is regularily sampled ==> ERROR in the code" << endl;
  }else{
    cout << " SpectralGridTest: the spectral window with id "<<spwId<<" is not regularily sampled as expected" << endl;
  }
  delete sgPtr2; sgPtr2=0;

  cout << endl;
  cout << endl;

  SpectralGrid* sgPtr3;

  cout << " SpectralGridTest: Test 3:" << endl;
  cout << " SpectralGridTest: Build using SpectralGrid( unsigned int numChan, double refFreq, double* chFreq, string freqUnits):" << endl;
  sgPtr3 = new SpectralGrid( numChan, refFreq, chFreq, "Hz");
  cout << " SpectralGridTest: Number of channels retrieved: " << sgPtr3->getNumChan() << " Input: " << numChan << endl;
  cout << " SpectralGridTest: Reference frequency retrieved:" << sgPtr3->getRefFreq().get()      << "Hz  Initial:" << refFreq << "Hz" << endl;
  cout << " SpectralGridTest: Reference frequency retrieved:" << sgPtr3->getRefFreq().get("MHz") << "MHz Initial:" << refFreq << "Hz" << endl;
  cout << " SpectralGridTest: Reference channel retrieved:  " << sgPtr3->getRefChan() << " Initial: " << refChan << endl;
  cout << " SpectralGridTest: Channel separation retrieved: " << sgPtr3->getChanSep().get() << " Initial: none"<< endl;

  { double chan=16.123456;
    cout << " SpectralGridTest: Position (GU) retrieved: "
         << sgPtr3->getChanNum(refFreq+(sgPtr3->getChanSep().get())*chan)
         << " SpectralGridTest:  Exact: " << chan << endl;
  }
  cout << " SpectralGridTest: Total bandwidth retrieved: " << sgPtr3->getBandwidth().get()
       << " SpectralGridTest: Initial: "<< (sgPtr3->getChanSep().get())*(numChan-1) << endl;

  cout << " SpectralGridTest: Frequency range: from "<< sgPtr3->getMinFreq().get() <<" to "<< sgPtr3->getMaxFreq().get() <<"Hz"<< endl;
  cout << " SpectralGridTest: Frequency range: from "<< sgPtr3->getMinFreq().get("GHz") <<" to "<< sgPtr3->getMaxFreq().get("GHz") <<"GHz"<< endl;

  delete sgPtr3;

  cout << endl;
  cout << endl;

  numChan         = 128;
  refChan         = 64;
  Frequency refFreq2(215.0,"GHz");
  Frequency chanSep2(0.02,"GHz");
  Frequency intermediateFreq(2.0,"GHz");
  Frequency bandWidth(1.0,"GHz");
  //  SidebandSide sbSide=LSB;
  //  SidebandType sbType=SSB;

  sgPtr1 = new SpectralGrid(numChan, refChan, refFreq2, chanSep2,
                            intermediateFreq, LSB, SSB);

  cout   << " SpectralGridTest: Number of spectral windows:            " << sgPtr1->getNumSpectralWindow() << " Expected: 2" << endl;

  for(unsigned int spwId=0; spwId<sgPtr1->getNumSpectralWindow(); spwId++){
    cout << " SpectralGridTest: Sideband:                              " << sgPtr1->getSideband(spwId) << endl;
    cout << " SpectralGridTest: LO frequency:                          " << sgPtr1->getLoFrequency(spwId) << "Hz " <<  endl;
    cout << " SpectralGridTest: Number of channels retrieved:          " << sgPtr1->getNumChan(spwId) << " for spwId " <<  spwId<<": "
         << " Input:" << numChan << endl;
    cout << " SpectralGridTest: Reference frequency retrieved:         " << sgPtr1->getRefFreq(spwId).get()
         << "  Input:" << refFreq2.get("GHz") << "GHz" << endl;
    cout << " SpectralGridTest: Reference frequency retrieved:         " << sgPtr1->getRefFreq(spwId).get("GHz")<< "GHz "
         << " Input:" << refFreq2.get("GHz") << "GHz" << endl;
    cout << " SpectralGridTest: Reference channel retrieved:           " << sgPtr1->getRefChan()
         << " Input:" << refChan << endl;
    cout << " SpectralGridTest: Channel separation retrieved:          " << sgPtr1->getChanSep(spwId).get() << "Hz "
         << " Input: |" << chanSep2.get("GHz") << "| GHz" << endl;
    cout << " SpectralGridTest: Channel separation retrieved:          " << sgPtr1->getChanSep(spwId).get("kHz") << "kHz "
         << " Input: |" << chanSep2.get("GHz") << "| GHz" << endl;
    cout << " SpectralGridTest: minFreq:                               " << sgPtr1->getMinFreq(spwId).get("GHz") << " GHz" << endl;
    cout << " SpectralGridTest: maxFreq:                               " << sgPtr1->getMaxFreq(spwId).get("GHz") << " GHz" << endl;
    cout << " SpectralGridTest: Channel (grid units) for the min:      " << sgPtr1->getChanNum(spwId,sgPtr1->getMinFreq(spwId).get()) << endl;
    cout << " SpectralGridTest: Channel (grid units) for the max:      " << sgPtr1->getChanNum(spwId,sgPtr1->getMaxFreq(spwId).get()) << endl;

    if(sgPtr1->isRegular(spwId)){
      cout << " SpectralGridTest: the spectral window with id "<<spwId<<" is regularily sampled" << endl;
    }else{
      cout << " SpectralGridTest: the spectral window with id "<<spwId<<" is not regularily sampled" << endl;
    }

    if(sgPtr1->getAssocSpwId(spwId).size()==0){
      cout << " SpectralGridTest: the spectral window with id "<< spwId <<" has no associated spectral window" << endl;
    }else{
      for(unsigned int n=0; n<sgPtr1->getAssocSpwId(spwId).size(); n++){
	unsigned int assocSpwId = sgPtr1->getAssocSpwId(spwId)[n];
	cout << " SpectralGridTest: the spectral window with id "<< spwId
	     << " has the associated spec. win. with id " <<  assocSpwId
	     << " (" <<  sgPtr1->getAssocNature(spwId)[n] << ")" << endl;

	for(unsigned int i=0; i<sgPtr1->getNumChan(spwId); i++){
	  cout << " SpectralGridTest: chan index:" << i << " "
	       <<  sgPtr1->getSideband(spwId) <<" "<<sgPtr1->getChanFreq(spwId,i).get("GHz")<<"GHz  "
	       <<  sgPtr1->getAssocNature(spwId)[n] <<" "<<sgPtr1->getChanFreq(assocSpwId,i).get("GHz")<<"GHz"<<endl;
	}
      }
    }
    cout << endl;



  }

  cout << " SpectralGridTest: TESTBED done" << endl;
  return 0;
}
