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
using namespace std;


#include "ATMPercent.h"
#include "ATMPressure.h"
#include "ATMNumberDensity.h"
#include "ATMMassDensity.h"
#include "ATMTemperature.h"
#include "ATMLength.h"
#include "ATMHumidity.h"
#include "ATMProfile.h"
using namespace atm;

  /** \brief A C++ main code to test the <a href="classatm_1_1AtmProfile.html">AtmProfile</a> Class
   *
   *   The test is structured as follows:
   *         - Creates the object "myProfile" using the <a href="classatm_1_1AtmProfile.html#z7_1">longest constructor </a>
   *           of the Class <a href="classatm_1_1AtmProfile.html">AtmProfile</a>.
   *         - Gets the number of layers of "myProfile" using the operator
   *           <a href="classatm_1_1AtmProfile.html#z7_21">.getNumLayer()</a> of the Class <a href="classatm_1_1AtmProfile.html">AtmProfile</a>.
   *         - Lists different physical and chemical parameters layer by layer, using operators such as
   *           <a href="classatm_1_1AtmProfile.html#z7_29">.getLayerWaterVaporMassDensity(size_t i)</a> of the Class <a href="classatm_1_1AtmProfile.html">AtmProfile</a>.
   *
   * The ouput of this test should be as follows:
   *
   * <b>
   * AtmProfileTest: BASIC ATMOSPHERIC PARAMETERS TO GENERATE REFERENCE ATMOSPHERIC PROFILE  <br>
   * AtmProfileTest:    <br>
   * AtmProfileTest: Ground temperature T:         270 K <br>
   * AtmProfileTest: Ground pressure P:            560 mb <br>
   * AtmProfileTest: Relative humidity rh:         20 % <br>
   * AtmProfileTest: Scale height h0:              2 km <br>
   * AtmProfileTest: Pressure step dp:             5 mb <br>
   * AtmProfileTest: Altitude alti:                5000 m <br>
   * AtmProfileTest: Attitude top atm profile atmh:48 km <br>
   * AtmProfileTest: Pressure step factordp1:      1.1  <br>
   * AtmProfileTest: Tropospheric lapse rate:    -5.6 K/km <br>
   * AtmProfileTest: Atmospheric type:             tropical <br>
   * AtmProfileTest:    <br>
   * AtmProfileTest:    <br>
   * AtmProfileTest: Object myProfile built with the AtmProfile CONSTRUCTOR and the above entries <br>
   *   <br>
   * AtmProfileTest: Number of layers returned:  22 <br>
   * AtmProfileTest: Layer parameters:   <br>
   * AtmProfileTest:  P: 554.977 mb  T: 269.598 K  Thickness: 143.444 m  WaterVapor: 0.000762199 kg m-3  WaterVapor: 2.5504e+22 m-3  CO: 1.93876e+18 m-3  O3: 5.62901e+17 m-3  N2O: 4.76212e+18 m-3  <br>
   * AtmProfileTest:  P: 543.967 mb  T: 268.706 K  Thickness: 175.101 m  WaterVapor: 0.000703854 kg m-3  WaterVapor: 2.35517e+22 m-3  CO: 1.90629e+18 m-3  O3: 5.58954e+17 m-3  N2O: 4.68575e+18 m-3 <br>
   * AtmProfileTest:  P: 530.751 mb  T: 267.615 K  Thickness: 214.579 m  WaterVapor: 0.000638519 kg m-3  WaterVapor: 2.13656e+22 m-3  CO: 1.86583e+18 m-3  O3: 5.54106e+17 m-3  N2O: 4.59355e+18 m-3 <br>
   * AtmProfileTest:  P: 514.888 mb  T: 266.275 K  Thickness: 264.252 m  WaterVapor: 0.000566481 kg m-3  WaterVapor: 1.89551e+22 m-3  CO: 1.81512e+18 m-3  O3: 5.48123e+17 m-3  N2O: 4.48208e+18 m-3 <br>
   * AtmProfileTest:  P: 495.844 mb  T: 264.618 K  Thickness: 327.475 m  WaterVapor: 0.000488584 kg m-3  WaterVapor: 1.63486e+22 m-3  CO: 1.75116e+18 m-3  O3: 5.40691e+17 m-3  N2O: 4.34711e+18 m-3 <br>
   * AtmProfileTest:  P: 472.979 mb  T: 262.555 K  Thickness: 409.147 m  WaterVapor: 0.000406408 kg m-3  WaterVapor: 1.35988e+22 m-3  CO: 1.66905e+18 m-3  O3: 5.31317e+17 m-3  N2O: 4.18305e+18 m-3 <br>
   * AtmProfileTest:  P: 445.521 mb  T: 259.963 K  Thickness: 516.736 m  WaterVapor: 0.00032243 kg m-3  WaterVapor: 1.07889e+22 m-3  CO: 1.56532e+18 m-3  O3: 5.19589e+17 m-3  N2O: 3.98414e+18 m-3 <br>
   * AtmProfileTest:  P: 412.536 mb  T: 256.662 K  Thickness: 662.272 m  WaterVapor: 0.000240119 kg m-3  WaterVapor: 8.03464e+21 m-3  CO: 1.43556e+18 m-3  O3: 5.01384e+17 m-3  N2O: 3.74538e+18 m-3 <br>
   * AtmProfileTest:  P: 372.891 mb  T: 252.381 K  Thickness: 866.58 m  WaterVapor: 0.000163845 kg m-3  WaterVapor: 5.48243e+21 m-3  CO: 1.26268e+18 m-3  O3: 4.89541e+17 m-3  N2O: 3.44722e+18 m-3 <br>
   * AtmProfileTest:  P: 325.19 mb  T: 246.68 K  Thickness: 1169.46 m  WaterVapor: 9.84855e-05 kg m-3  WaterVapor: 3.29543e+21 m-3  CO: 1.03631e+18 m-3  O3: 4.87238e+17 m-3  N2O: 3.07023e+18 m-3 <br>
   * AtmProfileTest:  P: 267.67 mb  T: 238.76 K  Thickness: 1658.92 m  WaterVapor: 4.85606e-05 kg m-3  WaterVapor: 1.62489e+21 m-3  CO: 7.6885e+17 m-3  O3: 5.03793e+17 m-3  N2O: 2.58146e+18 m-3 <br>
   * AtmProfileTest:  P: 256.935 mb  T: 230.765 K  Thickness: 1000 m  WaterVapor: 2.49802e-05 kg m-3  WaterVapor: 8.35865e+20 m-3  CO: 5.53847e+17 m-3  O3: 5.37164e+17 m-3  N2O: 2.16403e+18 m-3 <br>
   * AtmProfileTest:  P: 257.26 mb  T: 224.015 K  Thickness: 1000 m  WaterVapor: 1.51513e-05 kg m-3  WaterVapor: 5.06978e+20 m-3  CO: 3.99903e+17 m-3  O3: 5.63587e+17 m-3  N2O: 1.87717e+18 m-3 <br>
   * AtmProfileTest:  P: 211.112 mb  T: 215.315 K  Thickness: 1625 m  WaterVapor: 7.86037e-06 kg m-3  WaterVapor: 2.63017e+20 m-3  CO: 2.48942e+17 m-3  O3: 5.7276e+17 m-3  N2O: 1.56008e+18 m-3 <br>
   * AtmProfileTest:  P: 158.866 mb  T: 203.315 K  Thickness: 2000 m  WaterVapor: 3.17588e-06 kg m-3  WaterVapor: 1.06268e+20 m-3  CO: 1.23924e+17 m-3  O3: 5.91703e+17 m-3  N2O: 1.1682e+18 m-3 <br>
   * AtmProfileTest:  P: 114.088 mb  T: 191.865 K  Thickness: 2000 m  WaterVapor: 1.16834e-06 kg m-3  WaterVapor: 3.90939e+19 m-3  CO: 5.57264e+16 m-3  O3: 1.47198e+18 m-3  N2O: 7.62395e+17 m-3 <br>
   * AtmProfileTest:  P: 81.0467 mb  T: 190.615 K  Thickness: 2000 m  WaterVapor: 4.29808e-07 kg m-3  WaterVapor: 1.43819e+19 m-3  CO: 2.61142e+16 m-3  O3: 2.77652e+18 m-3  N2O: 4.64509e+17 m-3 <br>
   * AtmProfileTest:  P: 50.4577 mb  T: 201.065 K  Thickness: 3750 m  WaterVapor: 1.02088e-07 kg m-3  WaterVapor: 3.41599e+18 m-3  CO: 1.54239e+16 m-3  O3: 3.93168e+18 m-3  N2O: 2.3445e+17 m-3 <br>
   * AtmProfileTest:  P: 27.4549 mb  T: 212.515 K  Thickness: 4000 m  WaterVapor: 1.47072e-08 kg m-3  WaterVapor: 4.9212e+17 m-3  CO: 1.04872e+16 m-3  O3: 4.55929e+18 m-3  N2O: 1.03582e+17 m-3 <br>
   * AtmProfileTest:  P: 15.0592 mb  T: 221.315 K  Thickness: 4000 m  WaterVapor: 1.99041e-09 kg m-3  WaterVapor: 6.66012e+16 m-3  CO: 7.05928e+15 m-3  O3: 3.23216e+18 m-3  N2O: 4.5226e+16 m-3 <br>
   * AtmProfileTest:  P: 8.2977 mb  T: 230.215 K  Thickness: 4250 m  WaterVapor: 2.53052e-10 kg m-3  WaterVapor: 8.46739e+15 m-3  CO: 4.4967e+15 m-3  O3: 1.76311e+18 m-3  N2O: 1.69912e+16 m-3 <br>
   * AtmProfileTest:  P: 4.09633 mb  T: 241.315 K  Thickness: 6000 m  WaterVapor: 1.95133e-11 kg m-3  WaterVapor: 6.52935e+14 m-3  CO: 2.56843e+15 m-3  O3: 6.47884e+17 m-3  N2O: 3.88669e+15 m-3  <br>
   */
int main()
{
  // double h_div_k=0.04799274551;        // plank/boltz in units of K/GHz

  // Atmospheretype   atmType = tropical; // Atmospheric type (to reproduce behavior above the tropopause)
  size_t atmType = 1; // TROPICAL
  Temperature      T( 270.0,"K" );     // Ground temperature
  Pressure         P( 560.0,"mb");     // Ground Pressure
  Humidity         H(  20.0,"%" );     // Ground Relative Humidity (indication)
  Length         Alt(  5000,"m" );     // Altitude of the site
  Length         WVL(   2.0,"km");     // Water vapor scale height
  double         TLR=  -5.6      ;     // Tropospheric lapse rate (must be in K/km)
  Length      topAtm(  48.0,"km");     // Upper atm. boundary for calculations
  Pressure     Pstep(  5.0,"mb");     // Primary pressure step
  double   PstepFact=         1.1;     // Pressure step ratio between two consecutive layers


  cout<<" AtmProfileTest: BASIC ATMOSPHERIC PARAMETERS TO GENERATE REFERENCE ATMOSPHERIC PROFILE"<<endl;
  cout<<" AtmProfileTest:   "<<endl;
  cout<<" AtmProfileTest: Ground temperature T:         " << T.get()         << " K"    <<endl;
  cout<<" AtmProfileTest: Ground pressure P:            " << P.get("mb")     << " mb"   <<endl;
  cout<<" AtmProfileTest: Relative humidity rh:         " << H.get("%")      << " %"    <<endl;
  cout<<" AtmProfileTest: Scale height h0:              " << WVL.get("km")   << " km"   <<endl;
  cout<<" AtmProfileTest: Pressure step dp:             " << Pstep.get("mb") << " mb"   <<endl;
  cout<<" AtmProfileTest: Altitude alti:                " << Alt.get()       << " m"    <<endl;
  cout<<" AtmProfileTest: Attitude top atm profile atmh:" << topAtm.get("km")<< " km"   <<endl;
  cout<<" AtmProfileTest: Pressure step factordp1:      " << PstepFact          << " "    <<endl;
  cout<<" AtmProfileTest: Tropospheric lapse rate:    " << TLR                << " K/km" <<endl;

  AtmProfile myProfile( Alt, P, T, TLR, H, WVL, Pstep, PstepFact,  topAtm, atmType );

  cout<<" AtmProfileTest: Atmospheric type:             " << myProfile.getAtmosphereType() <<endl;
  cout<<" AtmProfileTest:   "<<endl;
  cout<<" AtmProfileTest:   "<<endl;

  cout<<" AtmProfileTest: Object myProfile built with the AtmProfile CONSTRUCTOR and the above entries"<<endl;
  cout<<"  "<<endl;
  cout<<" AtmProfileTest: Number of layers returned:  " << myProfile.getNumLayer() <<endl;
  cout<<" AtmProfileTest: Layer parameters:  " <<endl;


  for(size_t i=0; i<myProfile.getNumLayer(); i++){
    cout << " AtmProfileTest:  P: "          << myProfile.getLayerPressure(i).get("mb")    << " mb"
	 << "  T: "          << myProfile.getLayerTemperature(i).get("K")   << " K"
	 << "  Thickness: "  << myProfile.getLayerThickness(i).get("m")   << " m"
	 << "  WaterVapor: " << myProfile.getLayerWaterVaporMassDensity(i).get("kgm**-3")  << " kg m-3"
	 << "  WaterVapor: " << myProfile.getLayerWaterVaporNumberDensity(i).get("m**-3")  << " m-3"
	 << "  CO: "         << myProfile.getLayerCO(i).get("m**-3")          << " m-3"
         << "  O3: "         << myProfile.getLayerO3(i).get("m**-3")          << " m-3"
         << "  N2O: "        << myProfile.getLayerN2O(i).get("m**-3")         << " m-3" << endl;
  }

}
