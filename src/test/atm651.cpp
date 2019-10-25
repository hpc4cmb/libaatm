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

// C++ main code corresponding to atm651.py (cf CSV-526)
//
// atm651.py (freq=651Ghz) is a short script created by T. Hunter to generate atmospheric
// ozone line profiles for comparison to observation

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

using namespace std;

#include "ATMPercent.h"
#include "ATMPressure.h"
#include "ATMNumberDensity.h"
#include "ATMMassDensity.h"
#include "ATMTemperature.h"
#include "ATMLength.h"
#include "ATMInverseLength.h"
#include "ATMOpacity.h"
#include "ATMAngle.h"
#include "ATMHumidity.h"
#include "ATMFrequency.h"
#include "ATMWaterVaporRadiometer.h"
#include "ATMWVRMeasurement.h"
#include "ATMProfile.h"
#include "ATMSpectralGrid.h"
#include "ATMRefractiveIndex.h"
#include "ATMRefractiveIndexProfile.h"
#include "ATMSkyStatus.h"

int main()
{
  // Initialize the Atmospheric model

  size_t      atmType = 3;       // 1=tropical (ALMA site), 2=midlatSummer, 3=midlatWinter
  atm::Temperature  T(269. ,"K" );     // Ground temperature
  atm::Pressure     P(55450. ,"Pa");   // Ground Pressure
  atm::Humidity     H(7.,"%" );       // Ground Relative Humidity (indication)
  atm::Length       Alt(5059,"m" );    // Altitude of the site
  atm::Length       WVL(   2.0,"km");  // Water vapor scale height
  double            TLR=  -5.6      ;  // Tropospheric lapse rate (must be in K/km)
  atm::Length       topAtm(  48.0,"km");  // Upper atm. boundary for calculations
  atm::Pressure     Pstep(  10.0,"mb"); // Primary pressure step (10.0 mb)
  double   PstepFact=         1.2;     // Pressure step ratio between two consecutive layers

  atm::AtmProfile myProfile( Alt, P, T, TLR, H, WVL, Pstep, PstepFact,  topAtm, atmType );

  cout<<"BASIC ATMOSPHERIC PARAMETERS TO GENERATE REFERENCE ATMOSPHERIC PROFILE"<<endl;
  printf("Ground temperature T:          %.1f K \n",T.get());
  printf("Ground pressure P:             %.0f Pa \n",P.get("Pa"));
  printf("Relative humidity rh:          %.0f percent \n",H.get("%"));
  printf("Scale height h0:               %.0f m \n",WVL.get("m"));
  printf("Pressure step dp:              %.0f mbar \n",Pstep.get("mbar"));
  printf("Altitude alti:                 %.0f m \n",Alt.get());
  printf("Altitude top atm profile atmh: %.0f m \n",topAtm.get("m"));
  printf("Pressure step factordp1:       %.2f \n",PstepFact);
  printf("Tropospherique lapse rate:     %.1f K/km \n" , TLR);

  cout<<"Number of layers :"<<myProfile.getNumLayer()<<endl;
  cout<<"First guess precipitable water vapor content: " << myProfile.getGroundWH2O().get("mm") << " mm" << endl;

  // -------------------------------------------------------------------
  // Create a spectralGrid containing 1 astronomical band
  // This SpectralGrid is a member of a RefrativeIndexProfile object
  // -------------------------------------------------------------------

  // Astronomical  spectral window
  //
  double refFreqAstro   = 650.75*1.e9;
  //double bandWidthAstro = 2.*1.e9;
  double chanWidthAstro = 0.002*1e9; // 2Mhz

  size_t numChanAstro=1500;
  size_t refChanAstro=numChanAstro/2+1;


  atm::SpectralGrid * spGrid_p = new atm::SpectralGrid(numChanAstro,   // size_t numChan
                                                                       refChanAstro,   // size_t refChan
                                                                       atm::Frequency(refFreqAstro), // Frequency refFreq
                                                                       atm::Frequency(chanWidthAstro));



    // Print description of spectral windows
    //
  cout<<"-----------------------------------------------"<<endl;
  for(size_t j=0; j<spGrid_p->getNumSpectralWindow(); j++){
    cout << "Spectral Window " << j
         << " Reference Frequency: " <<  spGrid_p->getRefFreq(j).get("GHz") << " GHz, "
         << " Reference Channel: " <<  spGrid_p->getRefChan(j)
         << " Freq. Resolution: " <<  spGrid_p->getChanSep(j).get("MHz") << " MHz, "
         << " LO frequency: " <<  spGrid_p->getLoFrequency(j)/1.e9 << " GHz, "
         << " Num. of channels: " << spGrid_p->getNumChan(j)<<endl
         << " "<<spGrid_p->getChanFreq(j,0).get()
         << " "<<spGrid_p->getChanFreq(j,numChanAstro-1).get()
         <<endl;

    //         for (size_t i=0;i<spGrid_p->getNumChan(j);i++)
    //           cout<<spGrid_p->getChanFreq(j,i).get("GHz")<<" ";
    //         cout<<endl;


  }
  cout<<"-----------------------------------------------"<<endl;



  cout << "Create SkyStatus object"<<endl; // and associates a WaterVaporRadiometer to it " << endl;

  // Define absorption profile with SpectralGridand Profile
  atm::RefractiveIndexProfile absProfile(*spGrid_p, myProfile);

  // -------------------------------------------------------------------
  // - Create a SkyStatus object
  // -------------------------------------------------------------------

  // Define skystatus instance
  atm::SkyStatus skyStatus(absProfile);
  skyStatus.setAirMass(1.);
  skyStatus.setUserWH2O(myProfile.getGroundWH2O());

  // -------------------------------------------------------------------
  // - Compute the transmission at 651Ghz
  // -------------------------------------------------------------------

  string filename="atm651.dat";
  ofstream resultFile(filename.c_str());
  int spwid=0;
  for (size_t i=0; i<numChanAstro; i++) {
    double dryOpa = skyStatus.getDryOpacity(spwid,i).get();
    double wetOpa = skyStatus.getWetOpacity(spwid,i).get();
    resultFile<<spGrid_p->getChanFreq(spwid,i).get()<<" "<<exp(-wetOpa-dryOpa)<<" "<<exp(-wetOpa)<<" "<<exp(-dryOpa)<<endl;
  }
  resultFile.close();

  cout<<endl;
  cout<<"Transmission from "<<spGrid_p->getChanFreq(spwid,0).get()<<"Ghz to "<<
    spGrid_p->getChanFreq(spwid,numChanAstro-1).get()<<"Ghz stored in file "<<filename<<endl;
  return 0;
}
