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
#include <math.h>
#include <string.h>
#include <stdlib.h>

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
#define STRLEN  40        // Max length of a row in a tpoint file
using namespace atm;
  /** \brief This is a C++ main code to test THE AB-INITIO SCHEME proposed by Robert Lucas:        <br>
   *  1. Get the ATM model parameters that reproduce:                                              <br>
   *     - The observed WVR channel temperatures                                                   <br>
   *     - Additional data from observed astronomical band                                         <br>
   *     - Additional weather data (P,T, sounding profile?)                                        <br>
   *  2. Compute from ATM the increased WVR channel temperatures for a water vapor increment w     <br>
   *  3. Get best linear stimator of w from these as a function of observed changes                <br>
   *     - Partial derivatives in (2) as coefficients                                              <br>
   *     -Weight approximately taking noise into account                                           <br>
   *  4. Compute form ATM model the path length prediction in observed astronomical band for an increment of water content w. <br>
   *  5. Compute output coefficients (combine 4 and 5)                                             <br>
   *                                                                                               <br>
   * The test is structured in xx steps and the output should be:                                  <br>
   *
   * <b>
   * AbInitioTest: <br>
   * </b>
   */
int main()
{
  #ifdef HAVE_WINDOWS
  fprintf(stdout, "Skipping test under Windows.\n");
  return 0;
  #endif
  cout << " AbInitioTest: STEP 1: CREATES REFERENCE ATMOSPHERIC PROFILE CORRESTONDING TO ATMOSPHERIC CONDITIONS AT THE GROUND:" << endl;
  cout << " AbInitioTest: Atmosphere Type: TROPICAL" << endl;
  cout << " AbInitioTest: Site Altitude: 4100 m above sea level" << endl;
  cout << " AbInitioTest: Ground Temperature: 268.15 K" << endl;
  cout << " AbInitioTest: Ground Pressure: 623.0 mb" << endl;
  cout << " AbInitioTest: Relative Humidity at ground level: 11.3 % " << endl;
  cout << " AbInitioTest: Water vapor scale height: 2.2 km" << endl;
  cout << " AbInitioTest: Tropospheric Temperature lapse rate: -5.6 K/km" << endl;
  cout << " AbInitioTest: Top of the atmosphere for calculations: 48 km" << endl;
  cout << " AbInitioTest: Primary pressure step: 10 mb (please check the above reference)" << endl;
  cout << " AbInitioTest: Pressure step factor: 1.2 (please check the above reference)" << endl;

  // Atmospheretype   atmType = tropical; // Atmospheric type (to reproduce behavior above the tropopause)
  size_t atmType = 1; // TROPICAL
  Temperature      T( 268.15,"K" );    // Ground temperature
  Pressure         P( 623.0,"mb");     // Ground Pressure
  Humidity         H(  11.30,"%" );    // Ground Relative Humidity (indication)
  Length         Alt(  4100,"m" );     // Altitude of the site
  Length         WVL(   2.2,"km");     // Water vapor scale height
  double         TLR=  -5.6      ;     // Tropospheric lapse rate (must be in K/km)
  Length      topAtm(  48.0,"km");     // Upper atm. boundary for calculations
  Pressure     Pstep(   5.0,"mb");     // Primary pressure step (5.0 mb)
  double   PstepFact=         1.1;     // Pressure step ratio between two consecutive layers

  AtmProfile myProfile( Alt, P, T, TLR, H, WVL, Pstep, PstepFact,  topAtm, atmType );

  cout << " AbInitioTest: First guess precipitable water vapor content: " << myProfile.getGroundWH2O().get("mm") << " mm" << endl;
  cout << " AbInitioTest:  " << endl;

  cout << " AbInitioTest: STEP 2: CREATES SpectralGrid and RefractiveIndexProfile OBJECTS CONTAINING" << endl;
  cout << " AbInitioTest:       3 WATER VAPOR RADIOMETER CHANNELS + 1 ASTRONOMICAL BAND       " << endl;
  cout << " AbInitioTest:         (Negative frequency resolutions are used for image side bands)    " << endl;

  vector<size_t> WVR_signalId;    // IDs of the Signal Side Band of each Water Vapor Radiometer Channel in the SpectralGrid object.

  size_t numchan1=11;  size_t refchan1=6;
  Frequency reffreq1(182.11,"GHz"); Frequency chansep1(  0.04,"GHz"); Frequency intfreq1(  1.20,"GHz");
  SidebandSide sidebandside1=LSB; SidebandType sidebandtype1=DSB;
  SpectralGrid alma_SpectralGrid(numchan1, refchan1, reffreq1, chansep1, intfreq1, sidebandside1, sidebandtype1);
  WVR_signalId.push_back(0); // this is the Id of the 1st spectral window of the 1st pair

  RefractiveIndexProfile alma_RefractiveIndexProfile(alma_SpectralGrid, myProfile);

  size_t numchan2=11;  size_t refchan2=6;
  Frequency reffreq2(179.11,"GHz"); Frequency chansep2(  0.09,"GHz"); Frequency intfreq2(  4.20,"GHz");
  SidebandSide sidebandside2=LSB; SidebandType sidebandtype2=DSB;
  WVR_signalId.push_back(alma_RefractiveIndexProfile.getNumSpectralWindow());  // this will be the Id of the 1st spw of the 2nd pair
  alma_RefractiveIndexProfile.addNewSpectralWindow(numchan2, refchan2, reffreq2, chansep2, intfreq2, sidebandside2, sidebandtype2);

  size_t numchan3=11; size_t refchan3=6;
  Frequency reffreq3(175.51,"GHz"); Frequency chansep3(  0.10,"GHz"); Frequency intfreq3(  7.80,"GHz");
  SidebandSide sidebandside3=LSB; SidebandType sidebandtype3=DSB;
  WVR_signalId.push_back(alma_RefractiveIndexProfile.getNumSpectralWindow());
  alma_RefractiveIndexProfile.addNewSpectralWindow(numchan3, refchan3, reffreq3, chansep3, intfreq3, sidebandside3, sidebandtype3);

  Frequency  mySingleFreq_astro2(467.75,"GHz");  Frequency  chanSep_astro2(0.050,"GHz");
  size_t numChan_astro2=31;  size_t refChan_astro2=16;


  alma_RefractiveIndexProfile.addNewSpectralWindow(numChan_astro2, refChan_astro2, mySingleFreq_astro2, chanSep_astro2);

  size_t astro_band=alma_RefractiveIndexProfile.getNumSpectralWindow()-1;

  for(size_t j=0; j<alma_RefractiveIndexProfile.getNumSpectralWindow(); j++){
    cout << " AbInitioTest: Spectral Window " << j
	 << " Central Frequency: " <<  alma_RefractiveIndexProfile.getRefFreq(j).get("GHz") << " GHz, "
	 << " Freq. Resolution: " <<  alma_RefractiveIndexProfile.getChanSep(j).get("MHz") << " MHz, "
	 << " Num. of channels: " << alma_RefractiveIndexProfile.getNumChan(j)    << endl;
  }

  cout << " AbInitioTest: Spectral windows associations:    " << endl;
  for(size_t j=0; j<alma_RefractiveIndexProfile.getNumSpectralWindow(); j++){
    if(alma_RefractiveIndexProfile.getAssocSpwId(j).size()==1){
      cout << " AbInitioTest: Spectral Window " << j << " associated to spectral window: "
	   <<  alma_RefractiveIndexProfile.getAssocSpwId(j)[0] << " (double band)" << endl;
    }else{
      cout << " AbInitioTest: Spectral Window " << j << " associated to spectral window: "
	   <<  j << " (single band)" << endl;
    }
  }
  cout << " AbInitioTest:  " << endl;




  cout << " AbInitioTest: STEP 4: CREATES SkyStatus object and associates a WaterVaporRadiometer to it " << endl;

  SkyStatus skyAntenna1(alma_RefractiveIndexProfile);
  vector<double> skycoupling183; // Sky couplings the WVR Channels (set to 0.7 below)
  vector<Percent> signalgain183; // Signal Side Band Gain of the WVR Channels (set to 50% below)
  double skyCoupling_1stGuess = 0.8;
  Temperature tspill(285.15,"K"); // Spillover Temperature see by the WVR
  for(size_t i=0; i<WVR_signalId.size(); i++){
    skycoupling183.push_back(skyCoupling_1stGuess);
    signalgain183.push_back(Percent(50.0,"%"));
  }
  WaterVaporRadiometer wvr183ghz(WVR_signalId,skycoupling183,signalgain183,tspill);
  skyAntenna1.setWaterVaporRadiometer(wvr183ghz);

  cout << " AbInitioTest:  " << endl;

  cout << " AbInitioTest: WaterVaporRadiometer characteristics: " << endl;
  size_t Ids;
  for(size_t i=0; i<skyAntenna1.getWaterVaporRadiometer().getIdChannels().size(); i++){
    Ids=skyAntenna1.getWaterVaporRadiometer().getIdChannels()[i];
    cout << " AbInitioTest: WVR Channel " << i << " SpectralGrid Id of Signal sideband: " <<
      Ids << " / SpectralGrid Id of Image sideband: " << skyAntenna1.getAssocSpwId(Ids)[0] << endl;
    cout << " AbInitioTest: " << " Sky Coupling: " << skyAntenna1.getWaterVaporRadiometerSkyCoupling(i) <<
      " / Gain of Signal sideband: " << skyAntenna1.getWaterVaporRadiometerSignalGain(i).get("%") << " %" << endl;
  }

  cout << " AbInitioTest:  " << endl;
  cout << " AbInitioTest: STEP 5: Reads a block of real WVRMeasurements taken at Mauna Kea on March 3, 2002 (local code)" << endl;

  vector<WVRMeasurement> RadiometerData;
  WVRMeasurement singleRadiometerData;
  vector<double> time;
  Angle aaa;
  FILE*  fp;
  #ifdef HAVE_WINDOWS
  fp = fopen("WVR_MAUNA_KEA\\radiometer_data.dat", "r");
  #else
  fp = fopen("WVR_MAUNA_KEA/radiometer_data.dat", "r");
  #endif
  if (fp != 0) {
    char  aRow[STRLEN+1];
    char* token;
    vector<Temperature> v_tsky;
    size_t numWVRChannels = 0;
    char * fgrow = fgets( aRow, STRLEN, fp );
    size_t inilen=strlen(aRow);
    size_t lacum=0;
    token = strtok(aRow,","); time.push_back(atof(token));
    lacum=lacum+strlen(token)+1;
    token = 0; token = strtok(token,","); aaa=Angle(atof(token),"deg");
    lacum=lacum+strlen(token)+1;
    while (lacum<=inilen){
      numWVRChannels++;
      token = 0; token = strtok(token,","); Temperature tt(atof(token),"K"); v_tsky.push_back(tt);
      lacum=lacum+strlen(token)+1;
    }
    singleRadiometerData=WVRMeasurement(aaa,v_tsky);
    RadiometerData.push_back(singleRadiometerData);
    fgrow = fgets( aRow, STRLEN, fp );
    while (feof(fp)==0){
      if(strncmp(aRow," ",1)==0){
	token = strtok(aRow,","); time.push_back(atof(token));
	token = 0; token = strtok(token,","); aaa=Angle(atof(token),"deg");
	for(size_t j=0; j<numWVRChannels-1; j++){
	  token = 0; token = strtok(token,","); v_tsky[j]=Temperature(atof(token),"K");
	}
	token = 0; token = strtok(token,"\n"); v_tsky[numWVRChannels-1]=Temperature(atof(token),"K");
      }
      singleRadiometerData=WVRMeasurement(aaa,v_tsky);
      RadiometerData.push_back(singleRadiometerData);
      fgrow = fgets( aRow, STRLEN, fp );
    }
  }
  fclose( fp );

  cout << " AbInitioTest: Total number of WVR data: " << RadiometerData.size() << endl;
  cout << " AbInitioTest:  " << endl;

  size_t FirstMeasurementAnalyzed=9300;
  size_t NumberofMeasurementsAnalyzed=1; // 7
  size_t NumberofMeasurementsforSkyCoupligRetrieval=5;


  cout << " AbInitioTest: STEP 6: Performing Water Vapor Retrieval over " << NumberofMeasurementsAnalyzed <<
    " WVR measurements starting at " << FirstMeasurementAnalyzed/3600 << " UT on March/3/2002 " << endl;

  skyAntenna1.WaterVaporRetrieval_fromWVR(RadiometerData,FirstMeasurementAnalyzed,FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed);
  cout << " AbInitioTest: The average Sigma of this ensemble of fits is: " <<
    skyAntenna1.getWVRAverageSigmaTskyFit(RadiometerData,FirstMeasurementAnalyzed,FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed).get("K")
       << " K" << endl;
  cout << " AbInitioTest: User Water Vapor Column: " << skyAntenna1.getUserWH2O().get("mm") << " mm" << endl;

  cout << " AbInitioTest:  " << endl;


  cout << " AbInitioTest: STEP 7: RETRIEVES THE BEST SKY COUPLING OF THE WVR CHANNELS USING DATA " <<  FirstMeasurementAnalyzed <<
    " TO " << FirstMeasurementAnalyzed+NumberofMeasurementsforSkyCoupligRetrieval << endl;

  skyAntenna1.updateSkyCoupling_fromWVR(RadiometerData,FirstMeasurementAnalyzed,FirstMeasurementAnalyzed+NumberofMeasurementsforSkyCoupligRetrieval);

  cout << " AbInitioTest: The best sky coupling is: " << endl;
  for(size_t i=0; i<skyAntenna1.getWaterVaporRadiometer().getIdChannels().size(); i++){
    cout << " AbInitioTest: WVR Channel " << i << ": " << skyAntenna1.getWaterVaporRadiometerSkyCoupling(i) << ": " << endl;
  }
  cout << " AbInitioTest: For which the AverageSigmaTskyFit over those measurements is: " <<
    skyAntenna1.getWVRAverageSigmaTskyFit(RadiometerData,FirstMeasurementAnalyzed,FirstMeasurementAnalyzed+NumberofMeasurementsforSkyCoupligRetrieval).get("K")
       << " K" << endl;
  cout << " AbInitioTest: User Water Vapor Column: " << skyAntenna1.getUserWH2O().get("mm") << " mm" << endl;
  cout << " AbInitioTest:" << endl;


  cout << " AbInitioTest: STEP 8: DATA ANALYSIS USING THE NEW SKY COUPLING " << endl;

  skyAntenna1.WaterVaporRetrieval_fromWVR(RadiometerData,FirstMeasurementAnalyzed,FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed);

  for(size_t i=FirstMeasurementAnalyzed; i<FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed; i++){
    cout << " AbInitioTest: Data point " << i << " (UT time: " << time[i]/3600 << " hours on 2002, March 3)" << endl;
    cout << " AbInitioTest: Measured and fitted Sky Tebb's (in K): " << endl;
    for(size_t j=0; j<RadiometerData[i].getmeasuredSkyBrightness().size(); j++){
      cout << " AbInitioTest:   Channel " << j << ": " << RadiometerData[i].getmeasuredSkyBrightness()[j].get("K") <<
	"  "<< RadiometerData[i].getfittedSkyBrightness()[j].get("K") << endl;
    }
    cout << " AbInitioTest:   Retrieved Zenith Water Vapor Column: " << RadiometerData[i].getretrievedWaterVaporColumn().get("mm") << " mm /"
	 << " Sigma Fit: " << RadiometerData[i].getSigmaFit().get("K") << " K" << endl;


    skyAntenna1.setUserWH2O(Length(2.0,"mm"));

    cout << " AbInitioTest: Alternative method: Retrieved Zenith Water Vapor Column: " <<
      skyAntenna1.WaterVaporRetrieval_fromTEBB(WVR_signalId,Percent(50.0,"%"),RadiometerData[i].getmeasuredSkyBrightness(),
					       RadiometerData[i].getAirMass(),skyAntenna1.getWaterVaporRadiometerSkyCoupling(0),
					       Temperature(285.15,"K")).get("mm") << " mm" << endl;

  }

  cout << " AbInitioTest: " << endl;
  cout << " AbInitioTest: STEP 9: ATMOSPHERE OBSERVABLES IN ASTRONOMICAL BAND (ID=" << astro_band << ")" << " CENTRAL FREQUENCY=" << skyAntenna1.getRefFreq(astro_band).get("GHz") << " GHz" << endl;

  cout << " AbInitioTest: " << endl;
  cout << " AbInitioTest: Zenith water vapor column above Antenna 1: " << skyAntenna1.getUserWH2O().get("mm") << " mm" << endl;
  cout << " AbInitioTest: Air Mass for Antenna 1: " << skyAntenna1.getAirMass() << endl;
  cout << " AbInitioTest: AVERAGE TEBB=" << skyAntenna1.getAverageTebbSky(astro_band).get("K") << " K" << endl;
  cout << " AbInitioTest: TEBB DERIVATIVE=" <<
    skyAntenna1.getAverageTebbSky(astro_band,Length(skyAntenna1.getUserWH2O().get("mm")+0.001,"mm")).get("K")-
    skyAntenna1.getAverageTebbSky(astro_band).get("K") << " K/(micron H2O)" << endl;
  cout << " AbInitioTest: AVERAGE WET OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageWetOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE H2O Lines OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageH2OLinesOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE H2O Cont. OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageH2OContOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest: AVERAGE DRY OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageDryOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE O2  Lines OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageO2LinesOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE Dry Cont. OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageDryContOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE O3  Lines OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageO3LinesOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE N2O Lines OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageN2OLinesOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE CO  Lines OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageCOLinesOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest: AVERAGE H2O DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageDispersiveH2OPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest: AVERAGE H2O NON-DISPER PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageNonDispersiveH2OPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest: AVERAGE H2O PATHLENGTH DERIVATIVE (ZENITH VALUE)=" << skyAntenna1.getAverageH2OPathLengthDerivative(astro_band) << " microns/micron_H2O" << endl;
  cout << " AbInitioTest: AVERAGE DRY DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageDispersiveDryPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest:        AVERAGE O2  DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageO2LinesPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest:        AVERAGE O3  DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageO3LinesPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest:        AVERAGE N2O DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageN2OLinesPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest:        AVERAGE CO  DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageCOLinesPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest: AVERAGE DRY DISPERSIVE PATHLENGTH Pr_DERIVATIVE  (ZENITH VALUE microns/mb)=" << skyAntenna1.getAverageDispersiveDryPathLength_GroundPressureDerivative(astro_band) << " microns/mb" << endl;
  cout << " AbInitioTest: AVERAGE DRY DISPERSIVE PATHLENGTH Temp_DERIVATIVE (ZENITH VALUE microns/K)=" << skyAntenna1.getAverageDispersiveDryPathLength_GroundTemperatureDerivative(astro_band) << " microns/K" << endl;
  cout << " AbInitioTest: AVERAGE DRY NON-DISPER PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageNonDispersiveDryPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest: AVERAGE DRY NON-DISPER PATHLENGTH Pr_DERIVATIVE  (ZENITH VALUE microns/mb)=" << skyAntenna1.getAverageNonDispersiveDryPathLength_GroundPressureDerivative(astro_band) << " microns/mb" << endl;
  cout << " AbInitioTest: AVERAGE DRY NON-DISPER PATHLENGTH Temp_DERIVATIVE (ZENITH VALUE microns/K)=" << skyAntenna1.getAverageNonDispersiveDryPathLength_GroundTemperatureDerivative(astro_band) << " microns/K" << endl;




  skyAntenna1.setUserWH2O(0.5,"mm");
  skyAntenna1.setAirMass(2.0);
  cout << " AbInitioTest: " << endl;
  cout << " AbInitioTest: Zenith water vapor column above Antenna 1: " << skyAntenna1.getUserWH2O().get("mm") << " mm" << endl;
  cout << " AbInitioTest: Air Mass for Antenna 1: " << skyAntenna1.getAirMass() << endl;
  cout << " AbInitioTest: AVERAGE TEBB=" << skyAntenna1.getAverageTebbSky(astro_band).get("K") << " K" << endl;
  cout << " AbInitioTest: TEBB DERIVATIVE=" <<
    skyAntenna1.getAverageTebbSky(astro_band,Length(skyAntenna1.getUserWH2O().get("mm")+0.001,"mm")).get("K")-
    skyAntenna1.getAverageTebbSky(astro_band).get("K") << " K/(micron H2O)" << endl;
  cout << " AbInitioTest: AVERAGE WET OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageWetOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE H2O Lines OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageH2OLinesOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE H2O Cont. OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageH2OContOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest: AVERAGE DRY OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageDryOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE O2  Lines OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageO2LinesOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE Dry Cont. OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageDryContOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE O3  Lines OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageO3LinesOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE N2O Lines OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageN2OLinesOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest:        AVERAGE CO  Lines OPACITY ALONG LINE OF SIGHT=" << skyAntenna1.getAverageCOLinesOpacity(astro_band).get("np")*skyAntenna1.getAirMass() << " np" << endl;
  cout << " AbInitioTest: AVERAGE H2O DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageDispersiveH2OPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest: AVERAGE H2O NON-DISPER PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageNonDispersiveH2OPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest: AVERAGE H2O PATHLENGTH DERIVATIVE (ZENITH VALUE)=" << skyAntenna1.getAverageH2OPathLengthDerivative(astro_band) << " microns/micron_H2O" << endl;
  cout << " AbInitioTest: AVERAGE DRY DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageDispersiveDryPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest:        AVERAGE O2  DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageO2LinesPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest:        AVERAGE O3  DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageO3LinesPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest:        AVERAGE N2O DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageN2OLinesPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest:        AVERAGE CO  DISPERSIVE PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageCOLinesPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest: AVERAGE DRY DISPERSIVE PATHLENGTH Pr_DERIVATIVE  (ZENITH VALUE microns/mb)=" << skyAntenna1.getAverageDispersiveDryPathLength_GroundPressureDerivative(astro_band) << " microns/mb" << endl;
  cout << " AbInitioTest: AVERAGE DRY DISPERSIVE PATHLENGTH Temp_DERIVATIVE (ZENITH VALUE microns/K)=" << skyAntenna1.getAverageDispersiveDryPathLength_GroundTemperatureDerivative(astro_band) << " microns/K" << endl;
  cout << " AbInitioTest: AVERAGE DRY NON-DISPER PATHLENGTH ALONG LINE OF SIGHT=" << skyAntenna1.getAverageNonDispersiveDryPathLength(astro_band).get("microns")*skyAntenna1.getAirMass() << " microns" << endl;
  cout << " AbInitioTest: AVERAGE DRY NON-DISPER PATHLENGTH Pr_DERIVATIVE  (ZENITH VALUE microns/mb)=" << skyAntenna1.getAverageNonDispersiveDryPathLength_GroundPressureDerivative(astro_band) << " microns/mb" << endl;
  cout << " AbInitioTest: AVERAGE DRY NON-DISPER PATHLENGTH Temp_DERIVATIVE (ZENITH VALUE microns/K)=" << skyAntenna1.getAverageNonDispersiveDryPathLength_GroundTemperatureDerivative(astro_band) << " microns/K" << endl;





  return 0;

}
