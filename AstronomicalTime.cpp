//BY STEVEN DORSHER, FALL 2012
#include "AstronomicalTime.h"
#include <string>
#include <iostream>
#include <iomanip>
#include<sstream>
#include<math.h>
using namespace std;

// Note that JDepoch is different for each method, since each one uses a 
// different reference Julian Date in its calculation


double AstronomicalTime::julianDate(int* yearDayLeap, double UTtime)
{
  int years = yearDayLeap[0];
  int daysIntoYear = yearDayLeap[1];
  int numberLeapYears = yearDayLeap[2];
  double epochOffset = JDepoch -0.5; // midnight of Jan 1 2000


  double JD = ((double) intDaysPerYear) * (years - Yepoch) + numberLeapYears;
  JD += epochOffset + (double) (daysIntoYear + UTtime/((double) hoursPerDay));

  return JD;
}



double AstronomicalTime::siderealTime(double JDdate, double UTtime, double longitude)
{

  double SToffset = 6.6460556; // ST at midnight at 0 longitude at the epoch
  double drift = 2400.0512617; // sidereal hours drifted per year
  double hoursSTperUT = 1.0027379;
  double epoch = 2415020;

  double ST = drift/daysPerCentury * (JDdate - epoch);
  ST += hoursSTperUT * UTtime +SToffset +longitude;
  ST = fmod(ST, (double) hoursPerDay);
  return ST;

}

double AstronomicalTime::HJD(double JD, double* RAdecPrecessed)
{
  double lightTimeToEarth = 8.316; //minutes
 
  double delta = RAdecPrecessed[1]*radiansPerDegree;
  double alpha = RAdecPrecessed[0]*radiansPerHour;
   
  double T = (JD - JDepoch)/daysPerCentury;
  double L = 280.460 + 36000.772*T;
  double M = 357.528+ 35999.050*T;
  double rate1 = 1.915 + 0.0048 *T;
  double rate2 = 0.02;

  //convert to radians
  L *= radiansPerDegree;
  M *= radiansPerDegree;
  rate1 *= radiansPerDegree;
  rate2 *= radiansPerDegree;

  double sunLong =  L + rate1*sin(M)+rate2*sin(2.0*M);

  double cosu = sin(delta)*sin(sunLong)*sin(obliqEcliptic);
  cosu += cos(delta)*cos(alpha)*cos(sunLong);
  cosu += cos(delta)*sin(alpha)*sin(sunLong)*cos(obliqEcliptic);
  
  double deltaT = lightTimeToEarth * cosu;
  double HJD = JD - deltaT/minutesPerDay;
  return HJD;
}



double* AstronomicalTime::precessCoordinates(double* RAdec, double JD)
{// calcualtion in radians, returns an answer in degrees and hours
   double RAepoch = RAdec[0]*radiansPerHour;
  double decEpoch = RAdec[1]*radiansPerDegree;
  
  double RAprecessed, decPrecessed;
  double precessConst = 50.4/(secondsPerMinute*minutesPerDegree); // 50.4 arcsec/year
  double time = (JD - JDepoch)/daysPerYear;
  precessConst *=radiansPerDegree;

  decPrecessed = decEpoch + precessConst*sin(obliqEcliptic)*cos(RAepoch)*time;
  RAprecessed=cos(obliqEcliptic)+sin(obliqEcliptic)*sin(RAepoch)*tan(decEpoch);
  RAprecessed *= time*precessConst;
  RAprecessed += RAepoch;

  double* RAdecPrecessed = new double[2];
  RAdecPrecessed[0] = RAprecessed/radiansPerHour;
  RAdecPrecessed[1] = decPrecessed/radiansPerDegree;
  return RAdecPrecessed;

}

double AstronomicalTime::hourAngle(double localSiderealTime, double* RAdecPrecessed)
{
  double hourAngle = localSiderealTime - RAdecPrecessed[0];
  if(hourAngle>hoursPerDay/2)
    {
      hourAngle = hourAngle-hoursPerDay;
    }

  return hourAngle;
}

double* AstronomicalTime::altitudeAzimuth(double hourAngle, double* RAdecPrecessed, double lattitudeDeg)
{
  double alpha = radiansPerHour*RAdecPrecessed[0]; // RA
  double delta = radiansPerDegree*RAdecPrecessed[1]; // declination
  double phi = lattitudeDeg * radiansPerDegree; // lattitude in radians
  double H = radiansPerHour*hourAngle; // 
  
  double sinAlt = sin(phi)*sin(delta)+cos(phi)*cos(delta)*cos(H);
  
  double altitude = asin(sinAlt);

  double cosAz = (sin(delta)-sin(phi)*sinAlt)/(cos(phi)*cos(altitude));
  
  double azimuth = acos(cosAz);
  if(H>0)
    {
      azimuth = 2.0*pi-azimuth;
    }

  azimuth *= degreesPerRadian;
  altitude *= degreesPerRadian;

  double* altAz = new double[2];
  altAz[0] = altitude;
  altAz[1] = azimuth;
  return altAz;

}

double AstronomicalTime::airMass(double hourAngle, double* RAdecPrecessed, double lattitudeDeg)
{
  double phi = radiansPerDegree*lattitudeDeg;
  double delta = radiansPerDegree*RAdecPrecessed[1]; //declination
  double H = radiansPerHour * hourAngle;
  double secz = 1.0/(sin(phi)*sin(delta) + cos(phi) *cos(delta)*cos(H));
  
  double airmass;
  
  //simple implementation, valid for hour angle less than 4. 
  airmass = secz;

  double approxLimitH = 4.0; //this approximation is only valid for low angle
  double q = airmass -1.0;
  
  double c1 = -0.0018167;
  double c2 = -0.002875;
  double c3 = -0.0008083;

  if(fabs(H)>approxLimitH)
    {
      airmass = airmass + c1*q + c2*q*q + c3*q*q*q;
    }
  
  return airmass;
}

double AstronomicalTime::phase(double HJD, double* ephemeris)
{
  double HJD0 = ephemeris[0];
  double period = ephemeris[1]; //days
  double phase = (HJD - HJD0)/period;
  phase = phase - floor(phase);
  //phase = fmod(phase, 1.0);
  return phase;
}



/////////////////////////////
// Time and date utilities



int* AstronomicalTime::calendarDateToYearDayLeap(int* calendarDate)
{// tested for 01/01/2000, 01/01/2001, 03/01/2012, 10/01/2012, 01/01/2099, 01/01/2101, 04/12/2013

  int month = calendarDate[0]; // mm
  int day = calendarDate[1]; // day within the month, dd
  int year = calendarDate[2]; // yyyy
  int startOfMonthDay; // 0 to 365. Will be added to day to get 
  // the day of the year.
  int dayOfTheYear;

  startOfMonthDay = 0;
  for(int m = 1; m<month; m++)
    {
      startOfMonthDay += daysInMonth(m, year);
    }

   

  bool currentYearIsLeap, leapYear;
  int numberLeapYears =0; // number of leap years since 2001

  currentYearIsLeap = false;
  for (int y=Yepoch+1; y<year; y++)
    {
      leapYear = checkLeapYear(y);
      if(leapYear)
	{
	  numberLeapYears ++;
	  if (y==year) currentYearIsLeap = true;
	}
    }

  dayOfTheYear = startOfMonthDay + day;
  
  int* dateTriplet = new int[3];
  dateTriplet[0] = year;
  dateTriplet[1]= dayOfTheYear;
  dateTriplet[2] = numberLeapYears;
  return dateTriplet;
}

int AstronomicalTime::daysInMonth(int month, int year)
{

  int monthEnd;
  switch(month)
    {
    case 1:
      monthEnd = 31;
      break;
    case 2:
      if(checkLeapYear(year))
	{
	  monthEnd = 29;
	}
      else
	{
	  monthEnd = 28;
	}
      break;
    case 3:
      monthEnd = 31;
      break;
    case 4:
      monthEnd = 30;
      break;
    case 5:
      monthEnd = 31;
      break;
    case 6:
      monthEnd = 30;
      break;
    case 7:
      monthEnd = 31;
      break;
    case 8:
      monthEnd = 31;
      break;
    case 9:
      monthEnd = 30;
      break;
    case 10:
      monthEnd = 31;
      break;
    case 11:
      monthEnd =30;
      break;
    case 12:
      monthEnd = 31;
      break;
    }

  return monthEnd;
}


bool AstronomicalTime::checkLeapYear(int year)
{//tested
  int leapYearInterval = 4;
  int leapYearSkipped =100;
  int leapYearNotSkipped = 400;
  bool leapYearRegular= ((year%leapYearInterval)==0);
  bool exception = ((year%leapYearSkipped)==0)&&((year%leapYearNotSkipped)!=0);
  bool leapYear = leapYearRegular&&(!exception);
  return leapYear;
}

bool AstronomicalTime::checkDaylightSavings(int* clockTime, int* calendarDate, int timeZone)
{// only works for 2012 currently and ignores the repeated hour. clocktime in ut

  int month = calendarDate[0];
  int date = calendarDate[1];
  int hour = clockTime[0]+timeZone;
  
  int* localDate = new int[3];
  localDate[0]=month;
  localDate[2] = calendarDate[2];
  
  if(hour>hoursPerDay)
    {
      hour-=hoursPerDay;
      date++;
    }
  else if (hour < 0)
    {
      hour +=hoursPerDay;
      date--; 
    }

  localDate[1]=date;
  localDate = reconcileCalendarDate(localDate);
  month = localDate[0];
  date = localDate[1];
  int year = localDate[2];
  
  bool daylightTime;

  int startMonth, startDate, endMonth, endDate, startHour, endHour;

  switch (year)
    {
    case 2012:

      startMonth = 3;
      startDate = 11;
      endMonth = 11;
      endDate = 4;
      startHour = 2;
      endHour = 2;
      break;
    case 2005:
      startMonth = 4;
      startDate = 3;
      startHour = 2;
      endMonth = 10;
      endDate = 30;
      endHour = 2;
      break;
    }

  daylightTime =false;
  if ((startMonth<month)&&(month<endMonth))
    {
      daylightTime =true;
    }
  else if (month==startMonth)
    {
      if(date>startDate)
	{
	  daylightTime=true;
	}
      else if(date==startDate)
	{
	  if(hour>=startHour)
	    {
	      daylightTime = true;
	    }
	}
    }
  else if (month == endMonth)
    {
      if(date<endDate)
	{
	  daylightTime = true;
	}
      else if (date == endDate)
	{
	  if(hour <endHour)
	    {
	      daylightTime = true;
	    }
	}
    }

  return daylightTime;
}



int* AstronomicalTime::UTtoLocalTime(int* UTclockTime, int* calendarDate, int timeZone)
{// only works for October 2012 currently
  int UTmonth = calendarDate[0];
  int UTdate = calendarDate[1];
  int UTyear = calendarDate[2];
  
  int UThour = UTclockTime[0];

  int localHour = UThour +timeZone;
  int localDate = UTdate;

  bool daylightTime = checkDaylightSavings(UTclockTime, calendarDate, timeZone);
  if (daylightTime)
    {
      localHour++;
    }


  if (localHour > hoursPerDay)
    {
      localHour -=hoursPerDay;
      localDate ++;
    }
  else if(localHour<0)
    {
      localHour +=hoursPerDay;
      localDate--;
    }
  
  int* localCalendarDate = new int[3];
  localCalendarDate[0] = UTmonth;
  localCalendarDate[1] = localDate;
  localCalendarDate[2] = UTyear;
  localCalendarDate = reconcileCalendarDate(localCalendarDate);

  int* localDateTime = new int[6];
  localDateTime[0] = localCalendarDate[0];
  localDateTime[1] = localCalendarDate[1];
  localDateTime[2] = localCalendarDate[2];
  localDateTime[3] = localHour;
  localDateTime[4] = UTclockTime[1];
  localDateTime[5] = UTclockTime[2];
  
  return localDateTime;
}



int* AstronomicalTime::reconcileCalendarDate(int* calendarDate)
{//only works for incrementing or decrementing date
  //tested and works
  int month = calendarDate[0];
  int date = calendarDate[1];
  int year = calendarDate[2];
  

  int* newCalendarDate = new int[3];
  int monthEnd;

  monthEnd = daysInMonth(month, year);

  if(date>monthEnd) 
    {
      month++;
      date =1;
    }
   if(month >monthsPerYear)
     {
       year++;
       month =1;
     }
 
  
  if(date<=0)
    {
      month--;
      if(month<=0)
	{
	  year --;
	  month=monthsPerYear;
	}
      date = daysInMonth(month,year);
    }

  newCalendarDate[0] = month;
  newCalendarDate[1] = date;
  newCalendarDate[2] = year;
  return newCalendarDate;

}

///////////////////////
// Simple arithmetic on dates or times, or string processing below here   


string AstronomicalTime::makeDateOrTimeString(int* dateOrTime, string type)
{// tested
  
  string separator;
  if (type == "date")
    {
      separator = "/";
    }
  else 
    {
      separator = ":";
    }

  int x = dateOrTime[0];
  int y = dateOrTime[1];
  int z = dateOrTime[2];

  string dateOrTimeString;
  ostringstream convert;
  
  if (x<10)
    {
      convert << "0" << x;
    }
  else 
    {
      convert << x;
    }
  convert << separator;
  if (y < 10)
    {
      convert << "0" <<y;
    }
  else
    {
      convert << y;
    }
  convert << separator;
  if (z<10) // only true for times
    {
      convert << "0" << z;
    }
  else
    {
      convert << z;
    }

  
  dateOrTimeString = convert.str();

  return dateOrTimeString;
}
  

string AstronomicalTime::makeTripletString(double* triplet, string separator, int precision)
{// tested
  
  /* date /, time :, RA and dec space
     date and time precision = 0, RA and dec precision = 1 */

  int x = (int) triplet[0];
  int y = (int) triplet[1];
  double zd = triplet[2];

  string tripletString;
  ostringstream convert;
  
  if (x < 0)
    {
      convert << "-" ;
    }
  if (abs(x)<10)
    {
      convert << "0";
    }
 
  convert << abs(x);
 
  convert << separator;
  if (y < 10)
    {
      convert << "0" <<y;
    }
  else
    {
      convert << y;
    }
  convert << separator;

  if (zd<10) // only true for times
    {
      convert << "0";
    }

  convert << fixed << setprecision(precision) << zd;
  
  tripletString = convert.str();

  return tripletString;
}
  


int* AstronomicalTime::processDateString(string dateString)
{//tested
  string mm = dateString.substr(0,2); // 0 and 1
  string dd = dateString.substr(3,2); //  3 and 4
  string yyyy = dateString.substr(6,4); // 6 through 9
  
  istringstream convertmm(mm); 
  istringstream convertdd(dd);
  istringstream convertyyyy(yyyy);
  
  int month, day, year;

  convertmm >> month;
  convertdd >> day;
  convertyyyy >> year;
  
  int* calendarDate = new int[3];

  calendarDate[0] = month;
  calendarDate[1] = day;
  calendarDate[2] = year;
  
  return calendarDate;
}




string AstronomicalTime::decimalHourToString(double time)
{// REMOVE
  int hour = floor(time);
  double remainder = time - hour;
  int minutes = floor(60*remainder);
  remainder = 60.0*remainder - minutes;
  double seconds = remainder*60.0;

  ostringstream convert;
  if(hour<10)
    {
      convert << "0" << hour;
    }
  else
    {
      convert << hour;
    }
  convert << ":";
  if(minutes<10)
    {
      convert << "0" << minutes;
    }
  else 
    {
      convert << minutes;
    }
  convert << ":";
  convert << fixed << setprecision(1);
  if(seconds<10)
    {
      convert << "0"<< seconds;
    }
  else
    {
      convert << seconds;
    }
  
  string timeString = convert.str();
  return timeString;
}


double AstronomicalTime::clockTimeToDecimalHours(int* clockTime)
{//tested
  int hours = clockTime[0];
  int minutes = clockTime[1];
  int seconds = clockTime[2];
  
  double decimalSeconds = seconds + secondsPerMinute*minutes;
  decimalSeconds += secondsPerHour*hours;
  double decimalHours = (decimalSeconds/secondsPerDay)*((double) hoursPerDay);
  return decimalHours;
}


double* AstronomicalTime::decimalHourToTriplet(double time)
{// also works for degrees
  double sgn = time/fabs(time);
  time = fabs(time);
  int hour = floor(time);
  double remainder = time - hour;
  int minutes = floor(minutesPerHour*remainder);
  remainder = minutesPerHour*remainder-minutes;
  double seconds = remainder*secondsPerMinute;

  double* timeTriplet = new double[3];
  timeTriplet[0] = sgn*hour;
  timeTriplet[1] = minutes;
  timeTriplet[2] = seconds;
  return timeTriplet;
}
