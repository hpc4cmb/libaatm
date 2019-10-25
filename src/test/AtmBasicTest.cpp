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

  size_t      atmType = 1;       // 1=tropical (ALMA site), 2=midlatSummer, 3=midlatWinter
  atm::Temperature  T(273. ,"K" );     // Ground temperature
  atm::Pressure     P(55000. ,"Pa");   // Ground Pressure
  atm::Humidity     H(8.,"%" );       // Ground Relative Humidity (indication)
  atm::Length       Alt(5000,"m" );    // Altitude of the site
  atm::Length       WVL(   3.0,"km");  // Water vapor scale height
  double            TLR=  -6.5      ;  // Tropospheric lapse rate (must be in K/km)
  atm::Length       topAtm(  48.0,"km");  // Upper atm. boundary for calculations
  atm::Pressure     Pstep(  10.0,"mb"); // Primary pressure step (10.0 mb)
  double   PstepFact=         1.2;     // Pressure step ratio between two consecutive layers

  atm::AtmProfile myProfile( Alt, P, T, TLR, H, WVL, Pstep, PstepFact,  topAtm, atmType );

  cout<<"BASIC ATMOSPHERIC PARAMETERS TO GENERATE REFERENCE ATMOSPHERIC PROFILE"<<endl;
  printf("Ground temperature T:          %.1f K \n",T.get());
  printf("Ground pressure P:             %.0f Pa \n",P.get("Pa"));
  printf("Relative humidity rh:          %.0f percent \n",H.get("%"));
  printf("Scale height h0:               %.0f m \n",WVL.get("m"));
  printf("Pressure step dp:              %.0f Pa \n",Pstep.get("Pa"));
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
  double refFreqAstro   = 86.3*1.e9;
  double bandWidthAstro = 2.*1.e9;
  double chanWidthAstro = 0.1*1e9; // 100Mhz

  size_t numChanAstro=int(bandWidthAstro/chanWidthAstro);
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
      // << " Central Frequency: " <<  spGrid_p->getRefFreq(j).get("GHz") << " GHz, "
         << " Central Frequency: " <<  spGrid_p->getChanFreq(j,(spGrid_p->getNumChan(j)/2)).get("GHz") << " GHz, "
         << " Freq. Resolution: " <<  spGrid_p->getChanSep(j).get("MHz") << " MHz, "
         << " LO frequency: " <<  spGrid_p->getLoFrequency(j)/1.e9 << " GHz, "
         << " Num. of channels: " << spGrid_p->getNumChan(j)<<endl;

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

  cout<<"end"<<endl;

  // -------------------------------------------------------------------
  // - Retrieve path length
  // -------------------------------------------------------------------

  double deltaPressure = 100.;

  double startPressure = 49000.;
  double endPressure   = 61000.;
  int numPressure = int((endPressure-startPressure)/deltaPressure);

  int spwid = spGrid_p->getNumSpectralWindow()-1; // last spw = astronomical band

  string filename="path.dat";
  ofstream resultFile(filename.c_str());
  resultFile <<"! Pressure TotalPath H2OPathLength DispersiveDryPathLength NonDispersiveDryPathLength"<<endl;
  for (size_t i=0; i<numPressure; i++) {
    skyStatus.setBasicAtmosphericParameters(atm::Pressure(startPressure+i*deltaPressure,"Pa"));

    cout<<"Number of layers :"<<skyStatus.getNumLayer()<<endl;
    cout<<"First guess precipitable water vapor content: " << skyStatus.getGroundWH2O().get("mm") << " mm" << endl;
    cout<<"Ground Pressure: " << skyStatus.getGroundPressure().get("Pa") << " Pa" << endl;
    cout<<"Ground Temperature: " << skyStatus.getGroundTemperature().get("K") << " K" << endl;
    for (size_t ii=0; ii<skyStatus.getNumLayer(); ii++) {
      cout << "Pressure in Layer " << ii << " = " <<  skyStatus.getLayerPressure(ii).get("Pa") << " Pa " <<
	"Temperature: " <<  skyStatus.getLayerTemperature(ii).get("K") << " K " << endl;
    }


    atm::Length path = skyStatus.getAverageH2OPathLength(spwid)
      + skyStatus.getAverageDispersiveDryPathLength(spwid)
      + skyStatus.getAverageNonDispersiveDryPathLength(spwid);

    resultFile<<startPressure+i*deltaPressure<<" "
              <<path.get()<<" "
              <<skyStatus.getAverageH2OPathLength(spwid).get()<<" "
              <<skyStatus.getAverageDispersiveDryPathLength(spwid).get()<<" "
              <<skyStatus.getAverageNonDispersiveDryPathLength(spwid).get()<<endl;
  }
  resultFile.close();

  return 0;

}
