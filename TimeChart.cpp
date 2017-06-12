// BY STEVEN DORSHER, FALL 2012
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TimeChart.h"
#include "AstronomicalTime.h"
#include <string>
#include <math.h>
using namespace std;

int main(void) 
{
  TimeChart TC = TimeChart();
  double longitude, lattitude;
  string startDate;
  int numDaysInTable;
  int hourInterval;
  int starType;
  cout << "This program will print a table of GMUT, Julian Date," << endl;
  cout << "and local sidereal time to the file timeChart.txt" << endl<< endl;
  cout << "Please enter a longitude in hours: ";
  cin >> longitude;
  cout << "Please enter a lattitude in degrees: ";
  cin >> lattitude;
  cout << "Please enter a start date for the table in mm/dd/yyyy form: ";
  cin >> startDate;
  cout << "Please enter a number of days you wish the table to show: ";
  cin >> numDaysInTable;
  cout << "Please choose a star:\n";
  cout << "1. V471 Tauri\n";
  cout << "2. FF Aquarii\n";
  cout << "3. Algol\n";
  cout << "Enter a number: ";
  cin >> starType;

  // Call method printChart to generate the chart, with the arguments which
  // have been read in. 18 and 6 are the start and end hours in local time.

  TC.printChart(numDaysInTable, 17, 7, startDate, longitude, lattitude, starType);
   
  return 0;

}  

void TimeChart::printChart(int dayRequested, int firstHourShown, int lastHourShown, string startDateString, double longitude, double lattitude, int starType)
{

  double RAepoch[3];
  double* decEpock = new decEpoch[3];
  double* ephemeris = new double[2];
  string starString;
  double* phaseEclipse = new double[4]; // phase enter, reach min, end min, leave

  switch (starType)
    {
    case 1: // V471 Tauri
      starString = "V471 Tauri";
      RAepoch[0] = 3;
      RAepoch[1] = 50;
      RAepoch[2] = 24.968;
      decEpoch[0] = 17;
      decEpoch[1] = 14;
      decEpoch[2] = 47.42;
      // from SIMBAD, J2000

      ephemeris[0] = 2454028.452551; //551; //HJD0
      ephemeris[1] = 0.521183439; // period in days
      //L.Hric, E. Kundra, P. Dubovsky
      //"The fast brightness decline eclipses study-- the case of V471 Tau"
      //Contrib. Astron. Obs. Skalnate Pleso 41, 39-53 (2011)

      //This enters the eclipse and exits it about 0.04 in phase to either side of the minimum.
      //V. Miranda, T. Vaccaro, T.D. Oswalt
      // "V471 Tauri light curves and spot modelling"
      // Journal of the Southeastern Association for Research in Astronomy, 1, 17-20 (2007)
      break;
    case 2: //FF Aqr
      starString = "FF Aquarii";
      RAepoch[0] = 22;
      RAepoch[1] = 00;
      RAepoch[2] = 36.418;
      decEpoch[0] = -2;
      decEpoch[1] = 44;
      decEpoch[2] = 26.86;
      // from SIMBAD, J2000

      ephemeris[0] = 2452844.8186;
      ephemeris[1] = 9.207763;
      // E. Sipahi, S. Evern, G. Tas, C. Ibanoglu
      // "Photoelectric photometry of the unusual eclipsing binary system FF Aquarii"
      // Mem. SA. It. Vol 76, 627 (2005)

      // This also enters the eclipse and exits it about 0.04 in phase to either side of the minimum.
      // It was difficult to estimate from this light curve and this may be a slight overestimate.
      // E. Sipahi, S. Evren, G. Tas, and C. Ibanoglu
      // "Photoelectric photometry of the unusual eclipsing binary system FF Aquarii"
      // Mem. S.A. It. Vol. 76, 627 (2005)
      break;
    case 3:
      starString= "Algol";
      RAepoch[0] = 3;
      RAepoch[1] = 8;
      RAepoch[2] = 10.132;
      decEpoch[0] = 40;
      decEpoch[1] = 57;
      decEpoch[2] = 20.33;
      ephemeris[0] = 2452500.156;
      ephemeris[1] = 2.867357;
    }


  // Instantiate the AstronomicalTime object. See AstronomicalTime.cpp 
  // for methods called by AT.method(argument)
  AstronomicalTime AT = AstronomicalTime();
  
  // The date is read in as 09/29/2012, for example. This retrieves it and
  // stores it as {9, 29, 2012} in a 3 element array. Anything declared with
  // a star and a new is an array (technically a pointer to an array) 
  // of the dimensions given after the data type on the 
  // right hand side.
  
  int* startDate = new int[3];
  startDate = AT.processDateString(startDateString);
  
  int* yearDayLeap = new int[3];
  int* UTdate = new int[3];
  UTdate = startDate;

  double JD, ST, ST2,  UTtime;

  int* UTclockTime = new int[3];
  UTclockTime[1] = 0;
  UTclockTime[2] = 0;
  
  int* localDateTime = new int[6];
  int* localTime = new int[3];
  int* localDate = new int[3];

  double HJD;

  int timeZone;
  if(longitude>0)
    {
      timeZone = floor(longitude);
    }
  else
    {
      timeZone = ceil(longitude);
    }
  
    double RAdecEpoch[2];

  // convert RA and dec to decimal form
  double pi = 3.14159265;
  RAdecEpoch[0] = (RAepoch[0] + RAepoch[1]/60.0 + RAepoch[2]/60.0/60.0);
  RAdecEpoch[1]= (decEpoch[0] + decEpoch[1]/60.0 + decEpoch[2]/60.0/60.0);
    
  double* RAdecPrecessed = new double[2];
  double* RAtriplet = new double[3];
  double* decTriplet = new double[3];  
  
  double RAprecessed, decPrecessed;
  string RAstring, decString;
  string RAepochString, decEpochString;
  
  RAepochString = AT.makeTripletString(RAepoch, " ", 3);
  decEpochString = AT.makeTripletString(decEpoch, " ", 2);

  // Open the file stream output to timeChart.txt and send it some 
  // strings written below. 
  ofstream fs;
  fs.open("timeChart.txt");
  fs << "Steven Dorsher\n\n";
  fs << starString << endl;
  fs << "For longitude = " << longitude << " h and lattitude ";
  fs << lattitude << " degrees" << endl;
  fs << "For RA = " << RAepochString << ", dec = " << decEpochString << endl;
  fs << "For ephemeris: Min = HDJ" << fixed << setprecision(10); 
  fs <<ephemeris[0] << " + " << ephemeris[1] << " E" << endl;
  fs << endl;
  
  fs << "Local time            Universal Time        Julian Date";
  fs << "   HJD             ST            Prec. RA  Prec. Dec   Az.       Alt.\t\tAir M\tH Ang\tPhase\n";
  
  string localDateString, localTimeString, STstring, ST2string;
  string UTdateString, UTtimeString, RAdecString;
  int localHour;
  int hoursPerDay = 24;
  bool visible, tempBool1, tempBool2, tempBool3;
  double phase;

  double hourAngle;
  string altString, azString;

  double* azTriplet = new double[3];
  double* altTriplet = new double[3];
  double* altAz = new double[2];
  double airmass;
  string inEclipseString;

  // Loop over all the days within the range requested
  for (int day = 0; day < dayRequested; day++)
    {
      // For every hour in the day, going by the UT hour. 
      // Later we will check if it appears on the chart in local time.

      for(int UThour =0; UThour < hoursPerDay; UThour++)
	
	{
	  
	  /// THIS IS THE CODE WHICH DOES THE BULK OF THE CALCULATION
	  // see AstronomicalTime.cpp for the subroutines written here

	  yearDayLeap = AT.calendarDateToYearDayLeap(UTdate);
	  UTclockTime[0] = UThour;
	  UTtime = AT.clockTimeToDecimalHours(UTclockTime);
	  JD = AT.julianDate(yearDayLeap, UTtime);
	  ST = AT.siderealTime(JD, UTtime, longitude);
      	  //ST2 = AT.siderealTime2(JD, UTtime, longitude);
	  RAdecPrecessed = AT.precessCoordinates(RAdecEpoch, JD);
	  HJD = AT.HJD(JD, RAdecPrecessed);
	  localDateTime = AT.UTtoLocalTime(UTclockTime,UTdate,timeZone);
	  phase = AT.phase(HJD, ephemeris);
	  hourAngle = AT.hourAngle(ST, RAdecPrecessed);
	  altAz = AT.altitudeAzimuth(hourAngle, RAdecPrecessed, lattitude);
	  airmass = AT.airMass(hourAngle, RAdecPrecessed, lattitude);

	  ///////////////////////
	  //do some string processessing to get the results in printable form
	  altTriplet = AT.decimalHourToTriplet(altAz[0]);
	  azTriplet = AT.decimalHourToTriplet(altAz[1]);
	  altString = AT.makeTripletString(altTriplet, " ", 0);
	  azString = AT.makeTripletString(azTriplet, " ", 0);
	  
	  RAtriplet = AT.decimalHourToTriplet(RAdecPrecessed[0]);
	  RAstring = AT.makeTripletString(RAtriplet, " ", 0);
	  decTriplet = AT.decimalHourToTriplet(RAdecPrecessed[1]);
	  // the hour routine also works for degrees since they are 
	  // subdivided in the same manor
	  decString = AT.makeTripletString(decTriplet, " ", 0);
  	  
	  localDate[0]=localDateTime[0];
	  localDate[1]=localDateTime[1];
	  localDate[2]=localDateTime[2];
	  localTime[0]=localDateTime[3];
	  localTime[1]=localDateTime[4];
	  localTime[2]=localDateTime[5];
	  
	  UTdateString = AT.makeDateOrTimeString(UTdate, "date");
	  UTtimeString = AT.makeDateOrTimeString(UTclockTime, "time");
	  localDateString = AT.makeDateOrTimeString(localDate, "date");
	  localTimeString = AT.makeDateOrTimeString(localTime, "time");
	  STstring = AT.decimalHourToString(ST);
	  ST2string = AT.decimalHourToString(ST2);
	  ///////////////////////
	  
	  // logic for what appears on the table. 
	  // Note that the first night may be truncated based on what 
	  // local time corresponds to UT 00:00:00
	  
	  visible = (lastHourShown>firstHourShown);
	  visible = visible&&(localTime[0]<=lastHourShown);
	  visible = visible&&(localTime[0]>=firstHourShown);
	  tempBool1 = (localTime[0]<hoursPerDay);
	  tempBool1 = tempBool1&&(localTime[0]>=firstHourShown);
	  tempBool2 = (localTime[0]>=0);
	  tempBool2 = tempBool2&&(localTime[0]<=lastHourShown);
	  tempBool3 = tempBool1||tempBool2;
	  tempBool3 = tempBool3&&(lastHourShown<firstHourShown);
	  visible = visible||tempBool3;

	  //////////////
	  
	  // output to the file for each row
	  
	  if(visible)
	    {
	      
	      fs << localDateString << " " << localTimeString << "   " ;
	      fs << UTdateString << " " << UTtimeString << "   ";
	      fs << fixed << setprecision(5) <<JD << " ";
	      fs << HJD << "   ";
	      fs << STstring  << "    " ;
	      fs << RAstring << "  " << decString;
	      fs << "    "; 
	      fs << azString << "  " << altString; 
	      fs << setprecision(2); 
	      fs << "   " << airmass << " " <<hourAngle << "  ";
	      fs << setprecision(3) <<phase << endl;
	    }
	  //////////////////
	}  

      // at the end of the day, increment the day then reconcile the 
      // date triplet to make sure that it is still a consistent calendar day

      UTdate[1]++;
      UTdate = AT.reconcileCalendarDate(UTdate);

    } // return to start of loop
  
  fs.close();

}
	 
