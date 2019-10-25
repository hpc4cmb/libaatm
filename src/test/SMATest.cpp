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
#define STRLEN  200        // Max length of a row in a tpoint file
using namespace atm;

int main()
{

  Frequency rff;
  Frequency if1;
  Frequency chansep;
  cout << " SMATest:" << endl;
  cout << " SMATest: THIS PROPOSED TEST OF THE ATM INTERFACE SOFTWARE IS BASED ON SMA WVR DATA (17/July/2006)" << endl;
  cout << " " << endl;

  cout << " SmaTest: STEP 1: CREATES REFERENCE ATMOSPHERIC PROFILE CORRESTONDING TO THE FOLLOWING BASIC PARAMETERS:" << endl;

  //  Atmospheretype   atmType = tropical; // Atmospheric type (to reproduce behavior above the tropopause)
  unsigned int atmType = 1; // TROPICAL
  Temperature      T( 268.15,"K" );    // Ground temperature
  Pressure         P( 623.0,"mb");     // Ground Pressure
  Humidity         H(  11.30,"%" );    // Ground Relative Humidity (indication)
  Length         Alt(  4100,"m" );     // Altitude of the site
  Length         WVL(   2.2,"km");     // Water vapor scale height
  double         TLR=  -5.6      ;     // Tropospheric lapse rate (must be in K/km)
  Length      topAtm(  48.0,"km");     // Upper atm. boundary for calculations
  Pressure     Pstep(   5.0,"mb");     // Primary pressure step (10.0 mb)
  double   PstepFact=         1.1;     // Pressure step ratio between two consecutive layers

  Length      wvr_retrieved;
  Length      wvr_retrieved_astro;
  Temperature sigma_wvr_retrieved_astro;


  AtmProfile sma_AtmProfile( Alt, P, T, TLR, H, WVL, Pstep, PstepFact,  topAtm, atmType );

  cout << " SMATest: First guess precipitable water vapor content: " << sma_AtmProfile.getGroundWH2O().get("mm") << "mm" << endl;
  cout << " SMATest:  " << endl;


  cout << " SMATest: STEP 2: CREATES SpectralGrid and RefractiveIndexProfile objects" << endl;

  vector<unsigned int> WVR_signalId;

  unsigned int numchan0=25;  unsigned int refchan0=13;
   Frequency reffreq0(183.310-5.225,"GHz"); Frequency chansep0(2.650/25.0,"GHz"); Frequency intfreq0(0,"GHz");
   //  SidebandSide sidebandside0=LSB; SidebandType sidebandtype0=SSB;
  SpectralGrid sma_SpectralGrid(numchan0, refchan0, reffreq0, chansep0); //, intfreq0, sidebandside0, sidebandtype0);
  WVR_signalId.push_back(0); // this is the Id of the 1st spectral window of the 1st pair (WVR channel 1)
  // Creates RefractiveIndexProfile for the current atm profile (sma_AtmProfile) and spectral grid
  // (sma_SpectralGrid, so far with only the first WVR channel). Later on we will add new spectral
  // windows directly at the level of this RefractiveIndexProfile object.
  RefractiveIndexProfile sma_RefractiveIndexProfile(sma_SpectralGrid, sma_AtmProfile);

  unsigned int numchan1=25;  unsigned int refchan1=13;
  Frequency reffreq1(183.310-3.18,"GHz"); Frequency chansep1(1.400/25.0,"GHz"); Frequency intfreq1(0,"GHz");
  //  SidebandSide sidebandside1=LSB; SidebandType sidebandtype1=SSB;
  WVR_signalId.push_back(sma_RefractiveIndexProfile.getNumSpectralWindow()); // this is the Id of the 1st spectral window of the 2nd pair (WVR channel 2)
  sma_RefractiveIndexProfile.addNewSpectralWindow(numchan1, refchan1, reffreq1, chansep1); //, intfreq1, sidebandside1, sidebandtype1);

  unsigned int numchan2=25; unsigned int refchan2=13;
  Frequency reffreq2(183.310-1.9475,"GHz"); Frequency chansep2(0.845/25.0,"GHz"); Frequency intfreq2(0,"GHz");
  //  SidebandSide sidebandside2=LSB; SidebandType sidebandtype2=SSB;
  WVR_signalId.push_back(sma_RefractiveIndexProfile.getNumSpectralWindow()); // this is the Id of the 1st spectral window of the 3rd pair (WVR channel 3)
  sma_RefractiveIndexProfile.addNewSpectralWindow(numchan2, refchan2, reffreq2, chansep2);  //, intfreq2, sidebandside2, sidebandtype2);

  unsigned int numchan3=25; unsigned int refchan3=13;
  Frequency reffreq3(183.310-0.882,"GHz"); Frequency chansep3(0.206/25.0,"GHz"); Frequency intfreq3(0,"GHz");
  //  SidebandSide sidebandside3=LSB; SidebandType sidebandtype3=SSB;
  WVR_signalId.push_back(sma_RefractiveIndexProfile.getNumSpectralWindow()); // this is the Id of the 1st spectral window of the 3rd pair (WVR channel 3)
  sma_RefractiveIndexProfile.addNewSpectralWindow(numchan3, refchan3, reffreq3, chansep3);  //, intfreq3, sidebandside3, sidebandtype3);

  unsigned int numchan4=25; unsigned int refchan4=13;
  Frequency reffreq4(183.310+0.882,"GHz"); Frequency chansep4(0.206/25.0,"GHz"); Frequency intfreq4(0,"GHz");
  //  SidebandSide sidebandside3=LSB; SidebandType sidebandtype3=SSB;
  WVR_signalId.push_back(sma_RefractiveIndexProfile.getNumSpectralWindow()); // this is the Id of the 1st spectral window of the 3rd pair (WVR channel 3)
  sma_RefractiveIndexProfile.addNewSpectralWindow(numchan4, refchan4, reffreq4, chansep4);  //, intfreq3, sidebandside3, sidebandtype3);

  unsigned int numchan5=25; unsigned int refchan5=13;
  Frequency reffreq5(183.310+1.9475,"GHz"); Frequency chansep5(0.845/25.0,"GHz"); Frequency intfreq5(0,"GHz");
  //  SidebandSide sidebandside2=LSB; SidebandType sidebandtype2=SSB;
  WVR_signalId.push_back(sma_RefractiveIndexProfile.getNumSpectralWindow()); // this is the Id of the 1st spectral window of the 3rd pair (WVR channel 3)
  sma_RefractiveIndexProfile.addNewSpectralWindow(numchan5, refchan5, reffreq5, chansep5);  //, intfreq2, sidebandside2, sidebandtype2);

  unsigned int numchan6=25;  unsigned int refchan6=13;
  Frequency reffreq6(183.310+3.18,"GHz"); Frequency chansep6(1.400/25.0,"GHz"); Frequency intfreq6(0,"GHz");
  //  SidebandSide sidebandside1=LSB; SidebandType sidebandtype1=SSB;
  WVR_signalId.push_back(sma_RefractiveIndexProfile.getNumSpectralWindow()); // this is the Id of the 1st spectral window of the 2nd pair (WVR channel 2)
  sma_RefractiveIndexProfile.addNewSpectralWindow(numchan6, refchan6, reffreq6, chansep6); //, intfreq1, sidebandside1, sidebandtype1);

  unsigned int numchan7=25;  unsigned int refchan7=13;
  Frequency reffreq7(183.310+5.225,"GHz"); Frequency chansep7(2.650/25.0,"GHz"); Frequency intfreq7(0,"GHz");
  //  SidebandSide sidebandside1=LSB; SidebandType sidebandtype1=SSB;
  WVR_signalId.push_back(sma_RefractiveIndexProfile.getNumSpectralWindow()); // this is the Id of the 1st spectral window of the 2nd pair (WVR channel 2)
  sma_RefractiveIndexProfile.addNewSpectralWindow(numchan7, refchan7, reffreq7, chansep7); //, intfreq1, sidebandside1, sidebandtype1);


  for(unsigned int j=0; j<sma_RefractiveIndexProfile.getNumSpectralWindow(); j++){
    cout << " SMATest: Spectral Window " << j
	 << " Central Frequency: " <<  sma_RefractiveIndexProfile.getRefFreq(j).get("GHz") << " GHz, "
	 << " Freq. Resolution: " <<  sma_RefractiveIndexProfile.getChanSep(j).get("MHz") << " MHz, "
	 << " Num. of channels: " << sma_RefractiveIndexProfile.getNumChan(j)    << endl;
  }
  cout << " " << endl;

  cout << " SMATest: Spectral windows associations:    " << endl;
  for(unsigned int j=0; j<sma_RefractiveIndexProfile.getNumSpectralWindow(); j++){
    if(sma_RefractiveIndexProfile.getAssocSpwId(j).size()==1){
      cout << " SMATest: Spectral Window " << j << " associated to spectral window: "
	   <<  sma_RefractiveIndexProfile.getAssocSpwId(j)[0] << " (double band)" << endl;
    }else{
      cout << " SMATest: Spectral Window " << j << " associated to spectral window: "
	   <<  j << " (single band)" << endl;
    }
  }
  cout << " SMATest:  " << endl;


  cout << " SMATest:  STEP 3: CREATES SkyStatus object " << endl;

  SkyStatus sma_SkyStatus(sma_RefractiveIndexProfile);

  vector<double> skycouplingastro; // Sky couplings the astro channel taken as WWR
  vector<Percent> signalgainastro; // Signal Side Band Gain of the astro channel taken as WWR
  for(unsigned int i=0; i<1; i++){
    skycouplingastro.push_back(1.0);
    signalgainastro.push_back(Percent(50.0,"%"));
  }


  vector<double> skycoupling183; // Sky couplings the WVR Channels (set to 0.7 below)
  vector<Percent> signalgain183; // Signal Side Band Gain of the WVR Channels (set to 50% below)
  double skyCoupling_1stGuess = 0.8;    //0.9752;

  Temperature tspill(268.15,"K"); // Spillover Temperature seen by the WVR



  for(unsigned int i=0; i<WVR_signalId.size(); i++){
    skycoupling183.push_back(skyCoupling_1stGuess);
    signalgain183.push_back(Percent(100.0,"%"));
  }
  WaterVaporRadiometer wvr183ghz(WVR_signalId,skycoupling183,signalgain183,tspill);



  /*
  int numchan=25;  int refchan=13; Frequency chansep(  0.04,"GHz");
  Frequency if1(  6.00,"GHz");
  SidebandSide sidebandside_1=LSB; SidebandType sidebandtype_1=DSB;
  Frequency rff(342.505469,"GHz");
  sma_SkyStatus.addNewSpectralWindow(numchan, refchan, rff, chansep, if1, sidebandside_1, sidebandtype_1);
  */

  cout << " SMATest: skycoupling183.size()=" << skycoupling183.size() << endl;
  cout << " SMATest: wvr183ghz.getSkyCoupling().size()=" << wvr183ghz.getSkyCoupling().size() << endl;
  cout << " SMATest: wvr183ghz.getIdChannels().size()=" << wvr183ghz.getIdChannels().size() << endl;

  cout << " SMATest:  " << endl;

  cout << " SMATest: WaterVaporRadiometer characteristics: " << endl;
  unsigned int Ids;


  for(unsigned int i=0; i<sma_SkyStatus.getWaterVaporRadiometer().getIdChannels().size(); i++){
    cout << " SMATest: i=" << i << endl;
    Ids=sma_SkyStatus.getWaterVaporRadiometer().getIdChannels()[i];
    cout << " SMATest: WVR Channel " << i << " SpectralGrid Id of Signal sideband: " <<
      Ids << " / SpectralGrid Id of Image sideband: " << sma_SkyStatus.getAssocSpwId(Ids)[0] << endl;
    cout << "              " << " Sky Coupling: " << sma_SkyStatus.getWaterVaporRadiometerSkyCoupling(i) <<
      " / Gain of Signal sideband: " << sma_SkyStatus.getWaterVaporRadiometerSignalGain(i).get("%") << " %" << endl;
  }


  cout << " SMATest:  " << endl;
  cout << " SMATest: STEP 5: Reads a block of real WVRMeasurements taken with SMA on May 18, 2006" << endl;

  vector<WVRMeasurement> RadiometerData;
  WVRMeasurement singleRadiometerData;
  vector<double> time_mjd;
  vector<Angle> elevation;
  vector<Frequency> astrosignalfreq;
  vector<Frequency> astroimagefreq;
  vector<Frequency> astrobandwidth;
  vector<Temperature> groundtemp;
  vector<Pressure> groundpressure;
  vector<Percent> groundhumidity;
  vector<Temperature> astrotebb;
  vector<Length> astrowater;
  vector<Length> wvrwater;

  Angle aaa;
  FILE*  fp;
  fp = fopen("SMA/SMA_17JUL2006_skydip1.dat", "r");
  if (fp != 0) {
    cout << " SMATest: file open" << endl;
    char  aRow[STRLEN+1];
    char* token;
    vector<Temperature> v_tsky;
    unsigned int numWVRChannels = 0;
    char * fgrow = fgets( aRow, STRLEN, fp );
    unsigned int inilen=strlen(aRow);
    unsigned int lacum=0;

    token = strtok(aRow,","); lacum=lacum+strlen(token)+1; time_mjd.push_back(atof(token));
    token = 0; token = strtok(token,","); lacum=lacum+strlen(token)+1; aaa=Angle(atof(token),"rad"); elevation.push_back(aaa);
    token = 0; token = strtok(token,","); lacum=lacum+strlen(token)+1;

    while (lacum<=inilen){
      numWVRChannels++;
      token = 0; token = strtok(token,","); Temperature tt(atof(token)+0.0,"K"); v_tsky.push_back(tt);
      lacum=lacum+strlen(token)+1;
    }
    singleRadiometerData=WVRMeasurement(aaa,v_tsky);
    RadiometerData.push_back(singleRadiometerData);

    cout << " SMATest: WVR Channels with valid data: " << numWVRChannels <<endl;


    fgrow = fgets( aRow, STRLEN, fp );

    while (feof(fp)==0){

      if(strncmp(aRow," ",1)==0){              // Data file has a blank at the beginning of each row

	token = strtok(aRow,","); time_mjd.push_back(atof(token));
	token = 0; token = strtok(token,","); aaa=Angle(atof(token),"rad");     elevation.push_back(aaa);
	token = 0; token = strtok(token,",");

	for(unsigned int j=0; j<numWVRChannels-1; j++){
	  token = 0; token = strtok(token,","); v_tsky[j]=Temperature(atof(token)+0.0,"K");
	  // cout << "v_tsky["<<j<<"]=" << v_tsky[j].get("K") <<endl;
	}
	token = 0; token = strtok(token,"\n"); v_tsky[numWVRChannels-1]=Temperature(atof(token),"K");
	// cout << "v_tsky[numWVRChannels-1]=" << v_tsky[numWVRChannels-1].get("K") <<endl;
	singleRadiometerData=WVRMeasurement(aaa,v_tsky);
	RadiometerData.push_back(singleRadiometerData);
      }
      fgrow = fgets( aRow, STRLEN, fp );
    }
  }
  fclose( fp );

  cout << " SMATest: Total number of WVR data: " << RadiometerData.size() << endl;
  cout << " SMATest: Elevation of last measurement: " << RadiometerData[RadiometerData.size()-1].getElevation().get("deg") << " deg" << endl;
  cout << " " << endl;

  unsigned int FirstMeasurementAnalyzed=0;
  unsigned int NumberofMeasurementsAnalyzed=5;

  //  sma_SkyStatus.setUserWH2O(Length(1.18686,"mm"));
  // sma_SkyStatus.setAirMass(1.213137262);   ELV 55.5185
  // sma_SkyStatus.setAirMass(1.063907694);   // ELV 70.04
  // sma_SkyStatus.setAirMass(1.2569);   // ELV 52.71

    cout << "STEP 6: Performing Water Vapor Retrieval over " << NumberofMeasurementsAnalyzed <<
    " WVR measurements starting " << FirstMeasurementAnalyzed/3600 << " on May/18/2002 " << endl;

  sma_SkyStatus.setWaterVaporRadiometer(wvr183ghz);

  cout << " SMATest: sma_SkyStatus.getWaterVaporRadiometer().getIdChannels().size()=" << sma_SkyStatus.getWaterVaporRadiometer().getIdChannels().size() << endl;

  cout << "Sky Coupling=" << sma_SkyStatus.getWaterVaporRadiometerSkyCoupling(0) << endl;

  sma_SkyStatus.WaterVaporRetrieval_fromWVR(RadiometerData,FirstMeasurementAnalyzed,FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed);

  cout << "The average Sigma of this ensemble of fits is: " <<
    sma_SkyStatus.getWVRAverageSigmaTskyFit(RadiometerData,FirstMeasurementAnalyzed,FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed).get("K")
       << " K" << endl;
  for(unsigned int i=FirstMeasurementAnalyzed; i<FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed; i++){
    cout << "Data point analyzed: " << i << "/ Measured (fitted) Sky Tebb's (in K): " <<
      RadiometerData[i].getmeasuredSkyBrightness()[0].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[0].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[1].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[1].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[2].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[2].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[3].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[3].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[4].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[4].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[5].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[5].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[6].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[6].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[7].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[7].get("K") << " K) "  << endl;

    cout << " Sigma Fit: " << RadiometerData[i].getSigmaFit().get("K") << " K / Retrieved Water Vapor Column: " <<
      RadiometerData[i].getretrievedWaterVaporColumn().get("mm") << " mm " << endl;



  }



  /*  for(unsigned int i=FirstMeasurementAnalyzed; i<FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed; i++){
    cout << "Measurement analyzed " << i << ": " << RadiometerData[i].getretrievedWaterVaporColumn().get("mm") << " mm" << endl;
} */

  cout << " " << endl;

  FirstMeasurementAnalyzed=100;
  NumberofMeasurementsAnalyzed=5;

  cout << " SMATest: STEP 7: RETRIEVES THE BEST SKY COUPLING AND WATER VAPOR FOR EACH MEASUREMENT INDIVIDUALLY " <<  FirstMeasurementAnalyzed <<
    " TO " << FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed << endl;


  sma_SkyStatus.updateSkyCoupling_fromWVR(RadiometerData,FirstMeasurementAnalyzed,FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed);

  cout << " SMATest: sma_SkyStatus.getWaterVaporRadiometer().getIdChannels().size()=" << sma_SkyStatus.getWaterVaporRadiometer().getIdChannels().size() << endl;


  cout << "Sky Coupling=" << sma_SkyStatus.getWaterVaporRadiometerSkyCoupling(0) << endl;

  FirstMeasurementAnalyzed=501;
  NumberofMeasurementsAnalyzed=1;
  sma_SkyStatus.WaterVaporRetrieval_fromWVR(RadiometerData,FirstMeasurementAnalyzed,FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed);

  cout << "The average Sigma of this ensemble of fits is: " <<
    sma_SkyStatus.getWVRAverageSigmaTskyFit(RadiometerData,FirstMeasurementAnalyzed,FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed).get("K")
       << " K" << endl;
  for(unsigned int i=FirstMeasurementAnalyzed; i<FirstMeasurementAnalyzed+NumberofMeasurementsAnalyzed; i++){
    /* cout << "Data point analyzed: " << i << "/ Measured (fitted) Sky Tebb's (in K): " <<
      RadiometerData[i].getmeasuredSkyBrightness()[0].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[0].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[1].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[1].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[2].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[2].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[3].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[3].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[4].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[4].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[5].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[5].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[6].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[6].get("K") << " K) "  <<
      RadiometerData[i].getmeasuredSkyBrightness()[7].get("K") <<"  K ("<< RadiometerData[i].getfittedSkyBrightness()[7].get("K") << " K) "  << endl;
    */
    /*    cout << i << " Time : " <<   time_mjd[i]
         << " elevation: " << elevation[i].get("deg")
         << " Sigma Fit: " << RadiometerData[i].getSigmaFit().get("K") << " K /"
         << " Retrieved Water Vapor Column: " << RadiometerData[i].getretrievedWaterVaporColumn().get("mm") << " mm " << endl;
    */

    cout << i << " " << time_mjd[i] << " " << elevation[i].get("deg") << " "
	 << RadiometerData[i].getfittedSkyBrightness()[0].get("K") << " "
	 << RadiometerData[i].getfittedSkyBrightness()[1].get("K")  << " "
	 << RadiometerData[i].getfittedSkyBrightness()[2].get("K")  << " "
	 << RadiometerData[i].getfittedSkyBrightness()[3].get("K")  << " "
	 << RadiometerData[i].getfittedSkyBrightness()[4].get("K")  << " "
	 << RadiometerData[i].getfittedSkyBrightness()[5].get("K")  << " "
	 << RadiometerData[i].getfittedSkyBrightness()[6].get("K")  << " "
	 << RadiometerData[i].getfittedSkyBrightness()[7].get("K")  << " "
	 << RadiometerData[i].getretrievedWaterVaporColumn().get("mm") << " " << RadiometerData[i].getSigmaFit().get("K") << endl;

  }




  return 0;


}
