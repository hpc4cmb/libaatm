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
#include <math.h>
using namespace std;


#include "ATMPercent.h"
#include "ATMPressure.h"
#include "ATMAngle.h"
#include "ATMNumberDensity.h"
#include "ATMMassDensity.h"
#include "ATMTemperature.h"
#include "ATMLength.h"
#include "ATMInverseLength.h"
#include "ATMOpacity.h"
#include "ATMHumidity.h"
#include "ATMFrequency.h"
#include "ATMProfile.h"
#include "ATMSpectralGrid.h"
#include "ATMRefractiveIndex.h"
#include "ATMRefractiveIndexProfile.h"

using namespace atm;
  /** \brief A C++ main code to test the <a href="classatm_1_1RefractiveIndexProfile.html">RefractiveIndexProfile</a> Class
   *
   *   The test is structured as follows:
   *         - Creates an object called "myProfile" belonging to the <a href="classatm_1_1AtmProfile.html">AtmProfile</a> Class. This object
   *           is the same created in <a href="AtmProfileTest_8cpp.html">AtmProfileTest.cpp</a>
   *         - A single frequency is defined (850 GHz).
   *         - A refractive index profile for 850 GHz and "myProfile" is constructed using  a constructor in the form
   *        <a href="classatm_1_1RefractiveIndexProfile.html#z8_0">RefractiveIndexProfile (Frequency frequency, AtmProfile atmProfile)</a>.
   *         - Several test are then preformed on this created RefractiveIndexProfile object (called "myRefractiveIndexProfile")
   *           are then performed and written on the screen.
   *
   *  The output of this test is the following:
   *
   *  <b>
   * RefractiveIndexProfileTest: BASIC ATMOSPHERIC PARAMETERS TO GENERATE REFERENCE ATMOSPHERIC PROFILE<br>
   *  <br>
   * RefractiveIndexProfileTest: Ground temperature T:         270 K<br>
   * RefractiveIndexProfileTest: Ground pressure P:            560 mb<br>
   * RefractiveIndexProfileTest: Relative humidity rh:         20 %<br>
   * RefractiveIndexProfileTest: Scale height h0:              2 km<br>
   * RefractiveIndexProfileTest: Pressure step dp:             5 mb<br>
   * RefractiveIndexProfileTest: Altitude alti:                5000 m<br>
   * RefractiveIndexProfileTest: Attitude top atm profile atmh:48 km<br>
   * RefractiveIndexProfileTest: Pressure step factordp1:      1.1 <br>
   * RefractiveIndexProfileTest: Tropospherique lapse rate:    -5.6 K/km<br>
   * RefractiveIndexProfileTest: Atmospheric type:             tropical<br>
   * RefractiveIndexProfileTest:   <br>
   * RefractiveIndexProfileTest:   <br>
   * RefractiveIndexProfileTest: Object myProfile built with the AtmProfile CONSTRUCTOR and the above entries<br>
   * RefractiveIndexProfileTest:   <br>
   * RefractiveIndexProfileTest: Number of layers returned:  22<br>
   * RefractiveIndexProfileTest: Layer parameters:  <br>
   * RefractiveIndexProfileTest:  P: 554.977 mb T: 269.598 K Thickness: 143.444 m WaterVapor: 0.000762199 kg m-3 WaterVapor: 2.5504e+22 m-3 CO: 1.93876e+18 m-3 O3: 5.62901e+17 m-3 N2O: 4.76212e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 543.967 mb T: 268.706 K Thickness: 175.101 m WaterVapor: 0.000703854 kg m-3 WaterVapor: 2.35517e+22 m-3 CO: 1.90629e+18 m-3 O3: 5.58954e+17 m-3 N2O: 4.68575e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 530.751 mb T: 267.615 K Thickness: 214.579 m WaterVapor: 0.000638519 kg m-3 WaterVapor: 2.13656e+22 m-3 CO: 1.86583e+18 m-3 O3: 5.54106e+17 m-3 N2O: 4.59355e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 514.888 mb T: 266.275 K Thickness: 264.252 m WaterVapor: 0.000566481 kg m-3 WaterVapor: 1.89551e+22 m-3 CO: 1.81512e+18 m-3 O3: 5.48123e+17 m-3 N2O: 4.48208e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 495.844 mb T: 264.618 K Thickness: 327.475 m WaterVapor: 0.000488584 kg m-3 WaterVapor: 1.63486e+22 m-3 CO: 1.75116e+18 m-3 O3: 5.40691e+17 m-3 N2O: 4.34711e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 472.979 mb T: 262.555 K Thickness: 409.147 m WaterVapor: 0.000406408 kg m-3 WaterVapor: 1.35988e+22 m-3 CO: 1.66905e+18 m-3 O3: 5.31317e+17 m-3 N2O: 4.18305e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 445.521 mb T: 259.963 K Thickness: 516.736 m WaterVapor: 0.00032243 kg m-3 WaterVapor: 1.07889e+22 m-3 CO: 1.56532e+18 m-3 O3: 5.19589e+17 m-3 N2O: 3.98414e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 412.536 mb T: 256.662 K Thickness: 662.272 m WaterVapor: 0.000240119 kg m-3 WaterVapor: 8.03464e+21 m-3 CO: 1.43556e+18 m-3 O3: 5.01384e+17 m-3 N2O: 3.74538e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 372.891 mb T: 252.381 K Thickness: 866.58 m WaterVapor: 0.000163845 kg m-3 WaterVapor: 5.48243e+21 m-3 CO: 1.26268e+18 m-3 O3: 4.89541e+17 m-3 N2O: 3.44722e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 325.19 mb T: 246.68 K Thickness: 1169.46 m WaterVapor: 9.84855e-05 kg m-3 WaterVapor: 3.29543e+21 m-3 CO: 1.03631e+18 m-3 O3: 4.87238e+17 m-3 N2O: 3.07023e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 267.67 mb T: 238.76 K Thickness: 1658.92 m WaterVapor: 4.85606e-05 kg m-3 WaterVapor: 1.62489e+21 m-3 CO: 7.6885e+17 m-3 O3: 5.03793e+17 m-3 N2O: 2.58146e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 256.935 mb T: 230.765 K Thickness: 1000 m WaterVapor: 2.49802e-05 kg m-3 WaterVapor: 8.35865e+20 m-3 CO: 5.53847e+17 m-3 O3: 5.37164e+17 m-3 N2O: 2.16403e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 257.26 mb T: 224.015 K Thickness: 1000 m WaterVapor: 1.51513e-05 kg m-3 WaterVapor: 5.06978e+20 m-3 CO: 3.99903e+17 m-3 O3: 5.63587e+17 m-3 N2O: 1.87717e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 211.112 mb T: 215.315 K Thickness: 1625 m WaterVapor: 7.86037e-06 kg m-3 WaterVapor: 2.63017e+20 m-3 CO: 2.48942e+17 m-3 O3: 5.7276e+17 m-3 N2O: 1.56008e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 158.866 mb T: 203.315 K Thickness: 2000 m WaterVapor: 3.17588e-06 kg m-3 WaterVapor: 1.06268e+20 m-3 CO: 1.23924e+17 m-3 O3: 5.91703e+17 m-3 N2O: 1.1682e+18 m-3<br>
   * RefractiveIndexProfileTest:  P: 114.088 mb T: 191.865 K Thickness: 2000 m WaterVapor: 1.16834e-06 kg m-3 WaterVapor: 3.90939e+19 m-3 CO: 5.57264e+16 m-3 O3: 1.47198e+18 m-3 N2O: 7.62395e+17 m-3<br>
   * RefractiveIndexProfileTest:  P: 81.0467 mb T: 190.615 K Thickness: 2000 m WaterVapor: 4.29808e-07 kg m-3 WaterVapor: 1.43819e+19 m-3 CO: 2.61142e+16 m-3 O3: 2.77652e+18 m-3 N2O: 4.64509e+17 m-3<br>
   * RefractiveIndexProfileTest:  P: 50.4577 mb T: 201.065 K Thickness: 3750 m WaterVapor: 1.02088e-07 kg m-3 WaterVapor: 3.41599e+18 m-3 CO: 1.54239e+16 m-3 O3: 3.93168e+18 m-3 N2O: 2.3445e+17 m-3<br>
   * RefractiveIndexProfileTest:  P: 27.4549 mb T: 212.515 K Thickness: 4000 m WaterVapor: 1.47072e-08 kg m-3 WaterVapor: 4.9212e+17 m-3 CO: 1.04872e+16 m-3 O3: 4.55929e+18 m-3 N2O: 1.03582e+17 m-3<br>
   * RefractiveIndexProfileTest:  P: 15.0592 mb T: 221.315 K Thickness: 4000 m WaterVapor: 1.99041e-09 kg m-3 WaterVapor: 6.66012e+16 m-3 CO: 7.05928e+15 m-3 O3: 3.23216e+18 m-3 N2O: 4.5226e+16 m-3<br>
   * RefractiveIndexProfileTest:  P: 8.2977 mb T: 230.215 K Thickness: 4250 m WaterVapor: 2.53052e-10 kg m-3 WaterVapor: 8.46739e+15 m-3 CO: 4.4967e+15 m-3 O3: 1.76311e+18 m-3 N2O: 1.69912e+16 m-3<br>
   * RefractiveIndexProfileTest:  P: 4.09633 mb T: 241.315 K Thickness: 6000 m WaterVapor: 1.95133e-11 kg m-3 WaterVapor: 6.52935e+14 m-3 CO: 2.56843e+15 m-3 O3: 6.47884e+17 m-3 N2O: 3.88669e+15 m-3<br>
   * RefractiveIndexProfileTest: First guess precipitable water vapor content: 1.57182mm<br>
   * RefractiveIndexProfileTest: (This value is estimated from the relative humidity at ground level and the water vapor scale height)<br>
   * RefractiveIndexProfileTest:   <br>
   * RefractiveIndexProfileTest:   <br>
   * RefractiveIndexProfileTest: Example 1: Absorption profile for a single frequency: 850 GHz<br>
   *  <br>
   * RefractiveIndexProfileTest: Absorption Profile built from RefractiveIndexProfile CONSTRUCTOR. Summary of results:<br>
   *<br>
   * RefractiveIndexProfileTest: Total Dry Opacity at      850 GHz for 1.0 air mass: 0.12149<br>
   *<br>
   * RefractiveIndexProfileTest: Total Dry Cont Opacity at      850 GHz for 1.0 air mass: 0.102063<br>
   * RefractiveIndexProfileTest: Total O2 lines Opacity at      850 GHz for 1.0 air mass: 0.0136942<br>
   * RefractiveIndexProfileTest: Total O3 lines Opacity at      850 GHz for 1.0 air mass: 0.00506383<br>
   * RefractiveIndexProfileTest: Total CO lines Opacity at      850 GHz for 1.0 air mass: 2.21902e-06<br>
   * RefractiveIndexProfileTest: Total N2O lines Opacity at      850 GHz for 1.0 air mass: 0.000666054<br>
   *<br>
   * RefractiveIndexProfileTest: Total Wet Opacity at      850 GHz for 1.0 air mass: 1.71934<br>
   *<br>
   * RefractiveIndexProfileTest: Total H2O lines Opacity at      850 GHz for 1.0 air mass: 1.1558<br>
   * RefractiveIndexProfileTest: Total H2O Cont Opacity at      850 GHz for 1.0 air mass: 0.563542<br>
   *<br>
   *<br>
   * RefractiveIndexProfileTest: Total Dispersive Delay at 850 GHz for 1.0 air mass: 0.0018242 meters <br>
   * RefractiveIndexProfileTest: Total Non-Dispersive Delay at 850 GHz for 1.0 air mass: 0.0110324 meters <br>
   * RefractiveIndexProfileTest: Total Dry Delay at 850 GHz for 1.0 air mass: 1.47498 meters <br>
   * RefractiveIndexProfileTest: Total O2 lines Delay at 850 GHz for 1.0 air mass: -3.91084e-05 meters <br>
   * RefractiveIndexProfileTest: Total O3 lines Delay at 850 GHz for 1.0 air mass: 2.49173e-07 meters <br>
   * RefractiveIndexProfileTest: Total CO lines Delay at 850 GHz for 1.0 air mass: 8.40603e-09 meters <br>
   * RefractiveIndexProfileTest: Total N2O lines Delay at 850 GHz for 1.0 air mass: 1.4191e-07 meters <br>
   *<br>
   *<br>
   * RefractiveIndexProfileTest: (your actual water vapor column is 1.57182 mm; 1.57182 mm<br>
   *  </b>
  */



int main()
{
  //  double h_div_k=0.04799274551;        // plank/boltz in units of K/GHz

  // Atmospheretype   atmType = tropical; // Atmospheric type (to reproduce behavior above the tropopause)
  size_t atmType = 1;  // TROPICAL
  Temperature      T( 270.0,"K" );     // Ground temperature
  Pressure         P( 560.0,"mb");     // Ground Pressure
  Humidity         H(  20.0,"%" );     // Ground Relative Humidity (indication)
  Length         Alt(  5000,"m" );     // Altitude of the site
  Length         WVL(   2.0,"km");     // Water vapor scale height
  double         TLR=  -5.6      ;     // Tropospheric lapse rate (must be in K/km)
  Length      topAtm(  48.0,"km");     // Upper atm. boundary for calculations
  Pressure     Pstep(   5.0,"mb");     // Primary pressure step
  double   PstepFact=         1.1;     // Pressure step ratio between two consecutive layers


  cout<<" RefractiveIndexProfileTest: BASIC ATMOSPHERIC PARAMETERS TO GENERATE REFERENCE ATMOSPHERIC PROFILE"<<endl;
  cout<<"  "<<endl;
  cout<<" RefractiveIndexProfileTest: Ground temperature T:         " << T.get()         << " K"    <<endl;
  cout<<" RefractiveIndexProfileTest: Ground pressure P:            " << P.get("mb")     << " mb"   <<endl;
  cout<<" RefractiveIndexProfileTest: Relative humidity rh:         " << H.get("%")      << " %"    <<endl;
  cout<<" RefractiveIndexProfileTest: Scale height h0:              " << WVL.get("km")   << " km"   <<endl;
  cout<<" RefractiveIndexProfileTest: Pressure step dp:             " << Pstep.get("mb") << " mb"   <<endl;
  cout<<" RefractiveIndexProfileTest: Altitude alti:                " << Alt.get()       << " m"    <<endl;
  cout<<" RefractiveIndexProfileTest: Attitude top atm profile atmh:" << topAtm.get("km")<< " km"   <<endl;
  cout<<" RefractiveIndexProfileTest: Pressure step factordp1:      " << PstepFact          << " "    <<endl;
  cout<<" RefractiveIndexProfileTest: Tropospherique lapse rate:    " << TLR                << " K/km" <<endl;

  AtmProfile myProfile( Alt, P, T, TLR, H, WVL, Pstep, PstepFact,  topAtm, atmType );

  cout<<" RefractiveIndexProfileTest: Atmospheric type:             " << myProfile.getAtmosphereType() <<endl;
  cout<<" RefractiveIndexProfileTest:   "<<endl;
  cout<<" RefractiveIndexProfileTest:   "<<endl;

  cout<<" RefractiveIndexProfileTest: Object myProfile built with the AtmProfile CONSTRUCTOR and the above entries"<<endl;
  cout<<" RefractiveIndexProfileTest:   "<<endl;
  cout<<" RefractiveIndexProfileTest: Number of layers returned:  " << myProfile.getNumLayer() <<endl;
  cout<<" RefractiveIndexProfileTest: Layer parameters:  " <<endl;


  for(size_t i=0; i<myProfile.getNumLayer(); i++){
    cout << " RefractiveIndexProfileTest:  P: "          << myProfile.getLayerPressure(i).get("mb")    << " mb"
	 << " T: "          << myProfile.getLayerTemperature(i).get("K")   << " K"
	 << " Thickness: "  << myProfile.getLayerThickness(i).get("m")   << " m"
	 << " WaterVapor: " << myProfile.getLayerWaterVaporMassDensity(i).get("kgm**-3")  << " kg m-3"
	 << " WaterVapor: " << myProfile.getLayerWaterVaporNumberDensity(i).get("m**-3")  << " m-3"
	 << " CO: "         << myProfile.getLayerCO(i).get("m**-3")          << " m-3"
         << " O3: "         << myProfile.getLayerO3(i).get("m**-3")          << " m-3"
         << " N2O: "        << myProfile.getLayerN2O(i).get("m**-3")         << " m-3" << endl;
  }


  //  int mylayers=myProfile.getNumLayer();


  cout << " RefractiveIndexProfileTest: First guess precipitable water vapor content: " << myProfile.getGroundWH2O().get("mm") << "mm" << endl;
  cout << " RefractiveIndexProfileTest: (This value is estimated from the relative humidity at ground level and the water vapor scale height)" << endl;
  cout<<" RefractiveIndexProfileTest:   "<<endl;
  cout<<" RefractiveIndexProfileTest:   "<<endl;

  Frequency          mySingleFreq(850,"GHz");
  cout << " RefractiveIndexProfileTest: Example 1: Absorption profile for a single frequency: " << mySingleFreq.get("GHz") << " GHz" << endl;
  cout<<"  "<<endl;


  RefractiveIndexProfile myRefractiveIndexProfile(mySingleFreq, myProfile);


  cout<<" RefractiveIndexProfileTest: Absorption Profile built from RefractiveIndexProfile CONSTRUCTOR. Summary of results:"<<endl;
  cout<<endl;
  cout<<" RefractiveIndexProfileTest: Total Dry Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " << myRefractiveIndexProfile.getDryOpacity().get() <<endl;
  cout<<endl;
  cout<<" RefractiveIndexProfileTest: Total Dry Cont Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " << myRefractiveIndexProfile.getDryContOpacity().get() <<endl;
  cout<<" RefractiveIndexProfileTest: Total O2 lines Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " << myRefractiveIndexProfile.getO2LinesOpacity().get() <<endl;
  cout<<" RefractiveIndexProfileTest: Total O3 lines Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " << myRefractiveIndexProfile.getO3LinesOpacity().get() <<endl;
  cout<<" RefractiveIndexProfileTest: Total CO lines Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " << myRefractiveIndexProfile.getCOLinesOpacity().get() <<endl;
  cout<<" RefractiveIndexProfileTest: Total N2O lines Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " << myRefractiveIndexProfile.getN2OLinesOpacity().get() <<endl;
  cout<<endl;
  cout<<" RefractiveIndexProfileTest: Total Wet Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<  myRefractiveIndexProfile.getWetOpacity().get() << endl;
  cout<<endl;
  cout<<" RefractiveIndexProfileTest: Total H2O lines Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<  myRefractiveIndexProfile.getH2OLinesOpacity().get() << endl;
  cout<<" RefractiveIndexProfileTest: Total H2O Cont Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<  myRefractiveIndexProfile.getH2OContOpacity().get() << endl;
  cout<<endl;
  cout<<endl;
   cout<<" RefractiveIndexProfileTest: Total Dispersive PathLength at "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<
     myRefractiveIndexProfile.getDispersiveH2OPathLength().get() << " meters " << endl;

   cout<<" RefractiveIndexProfileTest: Total Non-Dispersive PathLength at "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<
     myRefractiveIndexProfile.getNonDispersiveH2OPathLength().get() << " meters " << endl;

   cout<<" RefractiveIndexProfileTest: Ratio Dispersive/Non-Dispersive PathLength at "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<
     myRefractiveIndexProfile.getDispersiveH2OPathLength().get()/myRefractiveIndexProfile.getNonDispersiveH2OPathLength().get() << endl;

   cout<<" RefractiveIndexProfileTest: Total Dry PathLength at "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<
     myRefractiveIndexProfile.getNonDispersiveDryPathLength().get() << " meters " << endl;

   cout<<" RefractiveIndexProfileTest: Total O2 lines PathLength at "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<
     myRefractiveIndexProfile.getO2LinesPathLength().get() << " meters " << endl;
   cout<<" RefractiveIndexProfileTest: Total O3 lines PathLength at "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<
     myRefractiveIndexProfile.getO3LinesPathLength().get() << " meters " << endl;
   cout<<" RefractiveIndexProfileTest: Total CO lines PathLength at "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<
     myRefractiveIndexProfile.getCOLinesPathLength().get() << " meters " << endl;
   cout<<" RefractiveIndexProfileTest: Total N2O lines PathLength at "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<
     myRefractiveIndexProfile.getN2OLinesPathLength().get() << " meters " << endl;


  cout<<endl;
  cout<<endl;
  cout << " RefractiveIndexProfileTest: (your actual water vapor column is " << (myProfile.getGroundWH2O()).get("mm") << " mm; " << (myRefractiveIndexProfile.getGroundWH2O()).get("mm") << " mm" <<endl;
  cout<<endl;



  /*
  cout << "change the basic parameter"<<endl;
  cout << "=========================="<<endl;
  cout << "Old ground temperature: "<< T.get() << " K"    <<endl;
  Temperature newT(275.0,"K");
  cout << "New ground temperature:" << newT.get() << " K"    <<endl;

  myRefractiveIndexProfile.setBasicAtmosphericParameters(Alt, P, newT, TLR, H, WVL);

  cout << "(your actual water vapor column is " << (myRefractiveIndexProfile.getGroundWH2O()).get("mm") << " mm" <<endl;


  cout<<"Absorption Profile with this new temperature. Summary of results:"<<endl;
  cout<<endl;
  cout<<"Total Dry Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " << myRefractiveIndexProfile.getDryOpacity().get() <<endl;
  cout<<"Total Wet Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<
    myRefractiveIndexProfile.getWetOpacity().get()/(myProfile.getGroundWH2O()).get("mm") << " per mm " << endl;


  cout<<"Total Dispersive Delay at "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " <<
    (myRefractiveIndexProfile.getDispersivePathLength().get())/((myProfile.getGroundWH2O()).get("mm")) << " meters per mm of water vapor (" <<
    (100*(myRefractiveIndexProfile.getDispersivePathLength().get())/(myProfile.getGroundWH2O().get("mm")))/
    ((myRefractiveIndexProfile.getNonDispersivePathLength().get())/((myProfile.getGroundWH2O().get("mm")))) <<" % of the Non-dispersive one)" << endl;



  cout << "(your actual water vapor column is " << (myProfile.getGroundWH2O()).get("mm") << " mm)" <<endl;
  cout<<endl;

  cout << "Add a spectral window"<<endl;
  cout << "====================="<<endl;

  int numChan = 4;
  int refChan = 2;
  Frequency refFreq(284.97346,"GHz");   // 350
  Frequency chanSep(2.0,"MHz");

  myRefractiveIndexProfile.add(numChan, refChan, refFreq, chanSep);

  int numSpw = myRefractiveIndexProfile.getNumSpectralWindow();

  cout << "There are now " << numSpw << " spectral windows" << endl;


  cout<<"Absorption profiles including this new spectral window. Summary of results:"<<endl;
  cout<<endl;
  cout<<"Total Dry Opacity at      "<<mySingleFreq.get("GHz") << " GHz for 1.0 air mass: " << myRefractiveIndexProfile.getDryOpacity().get() <<endl;
  double freq;
  int    spwid=0;
  int    numCh;
  int    k=0;
  for(spwid=0; spwid<numSpw; spwid++){

    numCh = myRefractiveIndexProfile.getNumChan(spwid); cout <<"Spectral window "<<spwid<<" has "<<numCh<<" frequency channels"<<endl;
    for(int n=0; n<numCh; n++){

      freq = myRefractiveIndexProfile.getChanFreq(spwid,n).get("GHz");
      cout<<"Total Wet Opacity at      "<< freq << " GHz for 1.0 air mass: " <<
	myRefractiveIndexProfile.getWetOpacity(k).get()/(myRefractiveIndexProfile.getGroundWH2O()).get("mm") << " per mm " << endl;

      cout<<"Total Non-Dispersive Delay at "<< freq << " GHz for 1.0 air mass: " <<
	(myRefractiveIndexProfile.getNonDispersivePathLength(k).get())/(myProfile.getGroundWH2O().get("mm")) << " meters per mm of water vapor " << endl;

      cout<<"Total Dispersive Delay at "<< freq << " GHz for 1.0 air mass: " <<
	(myRefractiveIndexProfile.getDispersivePathLength(k).get())/(myProfile.getGroundWH2O().get("mm")) << " meters per mm of water vapor (" <<
	(100*(myRefractiveIndexProfile.getDispersivePathLength(k).get())/(myRefractiveIndexProfile.getGroundWH2O().get("mm")))/
	((myRefractiveIndexProfile.getNonDispersivePathLength(k).get())/(myRefractiveIndexProfile.getGroundWH2O().get("mm"))) <<" % of the Non-dispersive one)" << endl;

      cout<<"Total Dry Delay at "<< freq << " GHz for 1.0 air mass: " <<
	(myRefractiveIndexProfile.getDryPathLength(k).get("micron")) << " microns " << endl;


      cout << "(your actual water vapor column is " << (myRefractiveIndexProfile.getGroundWH2O()).get("mm") << " mm)." <<endl;
      cout<<endl;
      k++;
    }
  }
  */


}
