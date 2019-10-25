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
using namespace atm;
  /** \brief A C++ main code to test the <a href="classatm_1_1SkyStatus.html">SkyStatus</a> Class
   *
   *   The test 1 is structured as follows:
   *         - Creates an object called "myProfile" belonging to the <a href="classatm_1_1AtmProfile.html">AtmProfile</a> Class. This object
   *           is the same created in <a href="AtmProfileTest_8cpp.html">AtmProfileTest.cpp</a>
   *         - A single frequency is defined (850 GHz).
   *         - A refractive index profile, called "myRefractiveIndexProfile", for 850 GHz and "myProfile"
   *        is constructed using  a constructor in the form
   *        <a href="classatm_1_1RefractiveIndexProfile.html#z8_0">RefractiveIndexProfile (Frequency frequency, AtmProfile atmProfile)</a>.
   *         - An object of the <a href="classatm_1_1SkyStatus.html">SkyStatus</a> class is created using
   *         <a href="classatm_1_1SkyStatus.html#z12_0">the basic constructor</a> of this class. This object is called "mySky_siglefreq"
   *         - Radiative transfer results are obtained for "mySky_siglefreq" forcing some parameters (such as the water vapor column) to different
   *           values external to the object.
   *         - Several parameters are changed internally in the object "mySky_siglefreq" using setters, e.g.
   *           <a href="classatm_1_1SkyStatus.html#z13_1">setUserWH2O</a> and
   *           <a href="classatm_1_1SkyStatus.html#z13_3">setAirMass</a>, and new radiative transfer results are obtained.
   *         - In the final part of this test, a new object, called "mySky_1band_8channels", of the
   *           <a href="classatm_1_1SkyStatus.html">SkyStatus</a> class is
   *           created with again the <a href="classatm_1_1SkyStatus.html#z12_0">the basic constructor</a> but with an object
   *         of the <a href="classatm_1_1RefractiveIndexProfile.html">RefractiveIndexProfile</a> class that has a
   *         <a href="classatm_1_1SpectralGrid.html">SpectralGrid</a> with one band and eight frequecy channels.
   *         - Radiative transfer results are obtained for the eight frequency channels under different conditions specified by setters.
   *
   *  The output of this test is the following:
   *
   * <b>
   * SkyStatusTest: BASIC ATMOSPHERIC PARAMETERS TO GENERATE REFERENCE ATMOSPHERIC PROFILE
   * SkyStatusTest:   <br>
   * SkyStatusTest: Ground temperature T:         270 K<br>
   * SkyStatusTest: Ground pressure P:            560 mb<br>
   * SkyStatusTest: Relative humidity rh:         20 %<br>
   * SkyStatusTest: Scale height h0:              2 km<br>
   * SkyStatusTest: Pressure step dp:             5 mb<br>
   * SkyStatusTest: Altitude alti:                5000 m<br>
   * SkyStatusTest: Attitude top atm profile atmh:48 km<br>
   * SkyStatusTest: Pressure step factordp1:      1.1 <br>
   * SkyStatusTest: Tropospherique lapse rate:    -5.6 K/km<br>
   * SkyStatusTest: Atmospheric type:             tropical<br>
   * SkyStatusTest:   <br>
   * SkyStatusTest:   <br>
   * SkyStatusTest: Object myProfile built with the AtmProfile CONSTRUCTOR and the above entries<br>
   * SkyStatusTest:   <br>
   * SkyStatusTest: Number of layers returned:  22<br>
   * SkyStatusTest: Layer parameters:  <br>
   * SkyStatusTest: size of Temperature Profile: 22<br>
   * SkyStatusTest:  P: 554.977 mb T: 269.598 K Thickness: 143.444 m WaterVapor: 0.000762199 kg m-3 WaterVapor: 2.5504e+22 m-3 CO: 1.93876e+18 m-3 O3: 5.62901e+17 m-3 N2O: 4.76212e+18 m-3<br>
   * SkyStatusTest:  P: 543.967 mb T: 268.706 K Thickness: 175.101 m WaterVapor: 0.000703854 kg m-3 WaterVapor: 2.35517e+22 m-3 CO: 1.90629e+18 m-3 O3: 5.58954e+17 m-3 N2O: 4.68575e+18 m-3<br>
   * SkyStatusTest:  P: 530.751 mb T: 267.615 K Thickness: 214.579 m WaterVapor: 0.000638519 kg m-3 WaterVapor: 2.13656e+22 m-3 CO: 1.86583e+18 m-3 O3: 5.54106e+17 m-3 N2O: 4.59355e+18 m-3<br>
   * SkyStatusTest:  P: 514.888 mb T: 266.275 K Thickness: 264.252 m WaterVapor: 0.000566481 kg m-3 WaterVapor: 1.89551e+22 m-3 CO: 1.81512e+18 m-3 O3: 5.48123e+17 m-3 N2O: 4.48208e+18 m-3<br>
   * SkyStatusTest:  P: 495.844 mb T: 264.618 K Thickness: 327.475 m WaterVapor: 0.000488584 kg m-3 WaterVapor: 1.63486e+22 m-3 CO: 1.75116e+18 m-3 O3: 5.40691e+17 m-3 N2O: 4.34711e+18 m-3<br>
   * SkyStatusTest:  P: 472.979 mb T: 262.555 K Thickness: 409.147 m WaterVapor: 0.000406408 kg m-3 WaterVapor: 1.35988e+22 m-3 CO: 1.66905e+18 m-3 O3: 5.31317e+17 m-3 N2O: 4.18305e+18 m-3<br>
   * SkyStatusTest:  P: 445.521 mb T: 259.963 K Thickness: 516.736 m WaterVapor: 0.00032243 kg m-3 WaterVapor: 1.07889e+22 m-3 CO: 1.56532e+18 m-3 O3: 5.19589e+17 m-3 N2O: 3.98414e+18 m-3<br>
   * SkyStatusTest:  P: 412.536 mb T: 256.662 K Thickness: 662.272 m WaterVapor: 0.000240119 kg m-3 WaterVapor: 8.03464e+21 m-3 CO: 1.43556e+18 m-3 O3: 5.01384e+17 m-3 N2O: 3.74538e+18 m-3<br>
   * SkyStatusTest:  P: 372.891 mb T: 252.381 K Thickness: 866.58 m WaterVapor: 0.000163845 kg m-3 WaterVapor: 5.48243e+21 m-3 CO: 1.26268e+18 m-3 O3: 4.89541e+17 m-3 N2O: 3.44722e+18 m-3<br>
   * SkyStatusTest:  P: 325.19 mb T: 246.68 K Thickness: 1169.46 m WaterVapor: 9.84855e-05 kg m-3 WaterVapor: 3.29543e+21 m-3 CO: 1.03631e+18 m-3 O3: 4.87238e+17 m-3 N2O: 3.07023e+18 m-3<br>
   * SkyStatusTest:  P: 267.67 mb T: 238.76 K Thickness: 1658.92 m WaterVapor: 4.85606e-05 kg m-3 WaterVapor: 1.62489e+21 m-3 CO: 7.6885e+17 m-3 O3: 5.03793e+17 m-3 N2O: 2.58146e+18 m-3<br>
   * SkyStatusTest:  P: 256.935 mb T: 230.765 K Thickness: 1000 m WaterVapor: 2.49802e-05 kg m-3 WaterVapor: 8.35865e+20 m-3 CO: 5.53847e+17 m-3 O3: 5.37164e+17 m-3 N2O: 2.16403e+18 m-3<br>
   * SkyStatusTest:  P: 257.26 mb T: 224.015 K Thickness: 1000 m WaterVapor: 1.51513e-05 kg m-3 WaterVapor: 5.06978e+20 m-3 CO: 3.99903e+17 m-3 O3: 5.63587e+17 m-3 N2O: 1.87717e+18 m-3<br>
   * SkyStatusTest:  P: 211.112 mb T: 215.315 K Thickness: 1625 m WaterVapor: 7.86037e-06 kg m-3 WaterVapor: 2.63017e+20 m-3 CO: 2.48942e+17 m-3 O3: 5.7276e+17 m-3 N2O: 1.56008e+18 m-3<br>
   * SkyStatusTest:  P: 158.866 mb T: 203.315 K Thickness: 2000 m WaterVapor: 3.17588e-06 kg m-3 WaterVapor: 1.06268e+20 m-3 CO: 1.23924e+17 m-3 O3: 5.91703e+17 m-3 N2O: 1.1682e+18 m-3<br>
   * SkyStatusTest:  P: 114.088 mb T: 191.865 K Thickness: 2000 m WaterVapor: 1.16834e-06 kg m-3 WaterVapor: 3.90939e+19 m-3 CO: 5.57264e+16 m-3 O3: 1.47198e+18 m-3 N2O: 7.62395e+17 m-3<br>
   * SkyStatusTest:  P: 81.0467 mb T: 190.615 K Thickness: 2000 m WaterVapor: 4.29808e-07 kg m-3 WaterVapor: 1.43819e+19 m-3 CO: 2.61142e+16 m-3 O3: 2.77652e+18 m-3 N2O: 4.64509e+17 m-3<br>
   * SkyStatusTest:  P: 50.4577 mb T: 201.065 K Thickness: 3750 m WaterVapor: 1.02088e-07 kg m-3 WaterVapor: 3.41599e+18 m-3 CO: 1.54239e+16 m-3 O3: 3.93168e+18 m-3 N2O: 2.3445e+17 m-3<br>
   * SkyStatusTest:  P: 27.4549 mb T: 212.515 K Thickness: 4000 m WaterVapor: 1.47072e-08 kg m-3 WaterVapor: 4.9212e+17 m-3 CO: 1.04872e+16 m-3 O3: 4.55929e+18 m-3 N2O: 1.03582e+17 m-3<br>
   * SkyStatusTest:  P: 15.0592 mb T: 221.315 K Thickness: 4000 m WaterVapor: 1.99041e-09 kg m-3 WaterVapor: 6.66012e+16 m-3 CO: 7.05928e+15 m-3 O3: 3.23216e+18 m-3 N2O: 4.5226e+16 m-3<br>
   * SkyStatusTest:  P: 8.2977 mb T: 230.215 K Thickness: 4250 m WaterVapor: 2.53052e-10 kg m-3 WaterVapor: 8.46739e+15 m-3 CO: 4.4967e+15 m-3 O3: 1.76311e+18 m-3 N2O: 1.69912e+16 m-3<br>
   * SkyStatusTest:  P: 4.09633 mb T: 241.315 K Thickness: 6000 m WaterVapor: 1.95133e-11 kg m-3 WaterVapor: 6.52935e+14 m-3 CO: 2.56843e+15 m-3 O3: 6.47884e+17 m-3 N2O: 3.88669e+15 m-3<br>
   * SkyStatusTest: First guess precipitable water vapor content: 1.57182mm<br>
   * SkyStatusTest: (This value is estimated from the relative humidity at ground level and the water vapor scale height)<br>
   * SkyStatusTest:   <br>
   * SkyStatusTest:   <br>
   * SkyStatusTest:    TEST FOR JUST 1 FREQUENCY  <br>
   * SkyStatusTest:    =========================  <br>
   * SkyStatusTest:  <br>
   * SkyStatusTest: Frequency      : 850 GHz<br>
   * SkyStatusTest: Wet opacity    : 1.09385 nepers, for 1 mm H2O<br>
   * SkyStatusTest: Dry opacity    : 0.12149 nepers <br>
   * SkyStatusTest: Sky brightness : 189.623 K<br>
   * SkyStatusTest: (EXTERNAL CHANGE) water vapor column:           0.4 mm<br>
   * SkyStatusTest: (NEW OUTPUT) T_EBB =                         121.052 K <br>
   * SkyStatusTest: Current water vapor column in Radiance Class:1 mm<br>
   * SkyStatusTest:  <br>
   * SkyStatusTest: (INTERNAL CHANGE) Air mass:                     2   <br>
   * SkyStatusTest: (NEW OUTPUT) T_EBB =                         242.61 K <br>
   *<br>
   * SkyStatusTest: (INTERNAL CHANGE) water vapor column:           0.8 mm<br>
   * SkyStatusTest: (NEW OUTPUT) T_EBB =                         229.859 K <br>
   *<br>
   * SkyStatusTest:    TEST FOR ONE SPECTRAL WINDOW WITH SEVERAL CHANNELS  <br>
   * SkyStatusTest:    =====================================================  <br>
   * SkyStatusTest:  <br>
   * SkyStatusTest: water vapor column: 1 mm<br>
   * SkyStatusTest: Air mass          : 1<br>
   * <br>
   * SkyStatusTest: Freq: 847 GHz  /  T_EBB=190.099 K <br>
   * SkyStatusTest: Freq: 848 GHz  /  T_EBB=189.865 K <br>
   * SkyStatusTest: Freq: 849 GHz  /  T_EBB=196.705 K <br>
   * SkyStatusTest: Freq: 850 GHz  /  T_EBB=189.623 K <br>
   * SkyStatusTest: Freq: 851 GHz  /  T_EBB=189.39 K <br>
   * SkyStatusTest: Freq: 852 GHz  /  T_EBB=189.422 K <br>
   * SkyStatusTest: Freq: 853 GHz  /  T_EBB=189.933 K <br>
   * SkyStatusTest: Freq: 854 GHz  /  T_EBB=189.59 K <br>
   * SkyStatusTest:  <br>
   * SkyStatusTest: water vapor column: 1 mm<br>
   * SkyStatusTest: Air mass          : 1.2<br>
   * <br>
   * SkyStatusTest: Freq: 847 GHz  /  T_EBB=206.063 K <br>
   * SkyStatusTest: Freq: 848 GHz  /  T_EBB=205.841 K <br>
   * SkyStatusTest: Freq: 849 GHz  /  T_EBB=212.186 K <br>
   * SkyStatusTest: Freq: 850 GHz  /  T_EBB=205.61 K <br>
   * SkyStatusTest: Freq: 851 GHz  /  T_EBB=205.389 K <br>
   * SkyStatusTest: Freq: 852 GHz  /  T_EBB=205.419 K <br>
   * SkyStatusTest: Freq: 853 GHz  /  T_EBB=205.9 K <br>
   * SkyStatusTest: Freq: 854 GHz  /  T_EBB=205.578 K <br>
   * </b>
   */

int main()
{
  //  double h_div_k=0.04799274551;        // plank/boltz in un of K/GHz

  //  Atmospheretype   atmType = tropical; // Atmospheric type (to reproduce behavior above the tropopause)
  unsigned int atmType = 1; // TROPICAL
  Temperature      T( 270.0,"K" );     // Ground temperature
  Pressure         P( 560.0,"mb");     // Ground Pressure
  Humidity         H(  20.0,"%" );     // Ground Relative Humidity (indication)
  Length         Alt(  5000,"m" );     // Altitude of the site
  Length         WVL(   2.0,"km");     // Water vapor scale height
  double         TLR=  -5.6      ;     // Tropospheric lapse rate (must be in K/km)
  Length      topAtm(  48.0,"km");     // Upper atm. boundary for calculations
  Pressure     Pstep(   5.0,"mb");     // Primary pressure step (10.0 mb)
  double   PstepFact=         1.1;     // Pressure step ratio between two consecutive layers


  cout<<" SkyStatusTest: BASIC ATMOSPHERIC PARAMETERS TO GENERATE REFERENCE ATMOSPHERIC PROFILE"<<endl;
  cout<<" SkyStatusTest:   "<<endl;
  cout<<" SkyStatusTest: Ground temperature T:         " << T.get()         << " K"    <<endl;
  cout<<" SkyStatusTest: Ground pressure P:            " << P.get("mb")     << " mb"   <<endl;
  cout<<" SkyStatusTest: Relative humidity rh:         " << H.get("%")      << " %"    <<endl;
  cout<<" SkyStatusTest: Scale height h0:              " << WVL.get("km")   << " km"   <<endl;
  cout<<" SkyStatusTest: Pressure step dp:             " << Pstep.get("mb") << " mb"   <<endl;
  cout<<" SkyStatusTest: Altitude alti:                " << Alt.get()       << " m"    <<endl;
  cout<<" SkyStatusTest: Attitude top atm profile atmh:" << topAtm.get("km")<< " km"   <<endl;
  cout<<" SkyStatusTest: Pressure step factordp1:      " << PstepFact          << " "    <<endl;
  cout<<" SkyStatusTest: Tropospherique lapse rate:    " << TLR                << " K/km" <<endl;

  AtmProfile myProfile( Alt, P, T, TLR, H, WVL, Pstep, PstepFact,  topAtm, atmType );

  cout<<" SkyStatusTest: Atmospheric type:             " << myProfile.getAtmosphereType() <<endl;
  cout<<" SkyStatusTest:   "<<endl;
  cout<<" SkyStatusTest:   "<<endl;

  cout<<" SkyStatusTest: Object myProfile built with the AtmProfile CONSTRUCTOR and the above entries"<<endl;
  cout<<" SkyStatusTest:   "<<endl;
  cout<<" SkyStatusTest: Number of layers returned:  " << myProfile.getNumLayer() <<endl;
  cout<<" SkyStatusTest: Layer parameters:  " <<endl;


  cout << " SkyStatusTest: size of Temperature Profile: " << myProfile.getTemperatureProfile().size() << endl;
  for(unsigned int i=0; i<myProfile.getNumLayer(); i++){
    cout << " SkyStatusTest:  P: "          << myProfile.getLayerPressure(i).get("mb")    << " mb"
	 << " T: "          << myProfile.getLayerTemperature(i).get("K")   << " K"
	 << " Thickness: "  << myProfile.getLayerThickness(i).get("m")   << " m"
	 << " WaterVapor: " << myProfile.getLayerWaterVaporMassDensity(i).get("kgm**-3")  << " kg m-3"
	 << " WaterVapor: " << myProfile.getLayerWaterVaporNumberDensity(i).get("m**-3")  << " m-3"
	 << " CO: "         << myProfile.getLayerCO(i).get("m**-3")      << " m-3"
         << " O3: "         << myProfile.getLayerO3(i).get("m**-3")          << " m-3"
         << " N2O: "        << myProfile.getLayerN2O(i).get("m**-3")         << " m-3" << endl;
  }

  cout << " SkyStatusTest: First guess precipitable water vapor content: " << myProfile.getGroundWH2O().get("mm") << "mm" << endl;
  cout << " SkyStatusTest: (This value is estimated from the relative humidity at ground level and the water vapor scale height)" << endl;
  cout<<" SkyStatusTest:   "<<endl;
  cout<<" SkyStatusTest:   "<<endl;

  cout << " SkyStatusTest:    TEST FOR JUST 1 FREQUENCY  " << endl;
  cout << " SkyStatusTest:    =========================  " << endl;
  cout << " SkyStatusTest:  " << endl;

  Frequency  mySingleFreq(850,"GHz");

  RefractiveIndexProfile myRefractiveIndexProfile(mySingleFreq, myProfile);

  SkyStatus mySky_siglefreq(myRefractiveIndexProfile);

  Length   new_wh2o0(1.0,"mm");
  mySky_siglefreq.setUserWH2O(new_wh2o0);


  for(unsigned int i=0; i<mySky_siglefreq.getNumSpectralWindow(); i++){

    for(unsigned int j=0; j<mySky_siglefreq.getNumChan(i); j++){

      cout << " SkyStatusTest: Frequency      : " << mySky_siglefreq.getChanFreq(j).get("GHz") << " GHz" << endl;

      cout << " SkyStatusTest: Wet opacity    : " << mySky_siglefreq.getWetOpacity(j).get() << " nepers, for " << mySky_siglefreq.getUserWH2O().get("mm") << " mm H2O" << endl;

      cout << " SkyStatusTest: Dry opacity    : " << mySky_siglefreq.getDryOpacity(j).get() << " nepers " << endl;

      cout << " SkyStatusTest: Sky brightness : " << mySky_siglefreq.getAverageTebbSky(j).get("K") << " K" << endl;

    }

  }

  Length   new_wh2o(0.45,"mm");
  cout << " SkyStatusTest: (EXTERNAL CHANGE) water vapor column:           " << new_wh2o.get("mm")                       << " mm" << endl;
  cout << " SkyStatusTest: (NEW OUTPUT) T_EBB =                         " << mySky_siglefreq.getAverageTebbSky(new_wh2o).get("K") << " K " << endl;
  cout << " SkyStatusTest: Current water vapor column in Radiance Class:" << mySky_siglefreq.getUserWH2O().get("mm")           << " mm" << endl;
  cout << " SkyStatusTest:  " <<endl;


  mySky_siglefreq.setAirMass(2.0);

  cout << " SkyStatusTest: (INTERNAL CHANGE) Air mass:                     " <<  mySky_siglefreq.getAirMass()              << "   " <<endl;
  cout << " SkyStatusTest: (NEW OUTPUT) T_EBB =                         " <<  mySky_siglefreq.getAverageTebbSky().get("K")     << " K " <<endl;
  cout << endl;

  Length   new_wh2o1(0.8,"mm" );

  mySky_siglefreq.setUserWH2O(new_wh2o1);

  cout<< " SkyStatusTest: (INTERNAL CHANGE) water vapor column:           " << mySky_siglefreq.getUserWH2O().get("mm")        << " mm" <<endl;
  cout<< " SkyStatusTest: (NEW OUTPUT) T_EBB =                         " << mySky_siglefreq.getAverageTebbSky().get("K")      << " K " <<endl;
  cout<<endl;

  cout << " SkyStatusTest:    TEST FOR ONE SPECTRAL WINDOW WITH SEVERAL CHANNELS  " << endl;
  cout << " SkyStatusTest:    =====================================================  " << endl;
  cout << " SkyStatusTest:  " << endl;

  //  unsigned int numchan=4;     /* 400; */      /* 500; */   /* 64; */
  //  unsigned int refchan=2;     /* 200;  */    /* 250; */    /* 32; */

  //  Frequency reffreq(899.55,"GHz");
  // Frequency chansep(  0.5,"GHz");


  unsigned int numchan=5;     /* 400; */      /* 500; */   /* 64; */
  unsigned int refchan=3;      /* 200;  */    /* 250; */    /* 32; */

  Frequency reffreq(616.50,"GHz");
  Frequency chansep(  0.01,"GHz");


  SpectralGrid band(numchan, refchan, reffreq, chansep);
  RefractiveIndexProfile abs_band(band, myProfile);
  SkyStatus mySky_1band_8channels(abs_band);


  cout << " SkyStatusTest: water vapor column: " << mySky_1band_8channels.getUserWH2O().get("mm") << " mm" << endl;
  cout << " SkyStatusTest: Air mass          : " << mySky_1band_8channels.getAirMass() << endl;
  cout << " " << endl;
  for(unsigned int i=0; i<mySky_1band_8channels.getNumChan(0); i++){

    cout << " SkyStatusTest: Freq: " <<  mySky_1band_8channels.getChanFreq(i).get("GHz") << " GHz  /  T_EBB=" << mySky_1band_8channels.getTebbSky(i).get("K")  << " K " <<endl;

  }
  cout << " SkyStatusTest:  " << endl;

  mySky_1band_8channels.setAirMass(1.0);
  mySky_1band_8channels.setUserWH2O(0.45,"mm");
  mySky_siglefreq.setUserWH2O(0.45,"mm");

  cout << " SkyStatusTest: water vapor column: " << mySky_1band_8channels.getUserWH2O().get("mm") << " mm" << endl;
  cout << " SkyStatusTest: Air mass          : " << mySky_1band_8channels.getAirMass() << endl;
  cout << " " << endl;

  for(unsigned int i=0; i<mySky_1band_8channels.getNumChan(0); i++){

    cout << " SkyStatusTest: Freq: " <<  mySky_1band_8channels.getChanFreq(i).get("GHz") << " GHz  /  T_EBB=" << mySky_1band_8channels.getTebbSky(i).get("K")  << " Dry opacity: " << mySky_1band_8channels.getDryOpacity(i).get ("np") << " np "
	 << " Wet opacity: " << mySky_1band_8channels.getWetOpacity(i).get ("np") << " np "
	 << " Total opacity: " << mySky_1band_8channels.getTotalOpacity(i).get ("np") << " np "
<<endl;

  }


  for(unsigned int i=0; i<mySky_1band_8channels.getNumChan(0); i++){

    //    cout << " SkyStatusTest: Freq: " <<  mySky_1band_8channels.getChanFreq(i).get("GHz") << " GHz  /  T_EBB=" << mySky_1band_8channels.getTebbSky(i).get("K")  << " K " <<endl;
    /*  cout << mySky_1band_8channels.getChanFreq(i).get("GHz") << " " << mySky_1band_8channels.getTebbSky(i).get("K")  << " " <<
      mySky_1band_8channels.getH2OLinesOpacity(i).get() << " " <<
      mySky_1band_8channels.getH2OContOpacity(i).get() << " " <<
      mySky_1band_8channels.getO2LinesOpacity(i).get() << " " <<
      mySky_1band_8channels.getDryContOpacity(i).get() << " " <<
      mySky_1band_8channels.getO3LinesOpacity(i).get() << " " <<
      endl; */

  }

  return 0;

}
