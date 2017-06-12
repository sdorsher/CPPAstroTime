#include <string>
using namespace std;

class AstronomicalTime {
 private:
  static const double pi = 3.14159265;
  static const double radiansPerHour = 3.14159265/12.0;
  static const double radiansPerDegree = 3.14159265/180.0;
  static const double degreesPerRadian = 180.0/3.14159265;
  static const double hoursPerRadian = 12.0/3.14159265;
  static const double secondsPerMinute = 60.0;
  static const double minutesPerHour = 60.0;
  static const double secondsPerHour = 60.0*60.0;
  static const int hoursPerDay = 24;
  static const double secondsPerDay = 24.0*60.0*60.0;
  static const double minutesPerDay = 60.0*24.0;
  static const int monthsPerYear = 12;
  static const double daysPerYear = 365.25;
  static const int intDaysPerYear = 365;
  static const double minutesPerDegree = 60.0;
  static const double Yepoch = 2000.0; // epoch 2000. default in program
  static const double JDepoch = 2451545.0; //epoch 2000, Jan 1 at noon 2000
  static const double daysPerCentury = 36525.0;
  static const double obliqEcliptic = 23.439*3.14159265/180.0; // from usno website
    //(23 + 27.0/60.0+8.0/60.0/60.0)*3.14159265/180.0;// 23d27m8s

 public:
  ///////////
  // Astronomical methods
  double julianDate(int* yearDayLeap, double UTtime);
  double siderealTime(double JDate, double UTtime, double longitude);
  double HJD(double JD, double* RAdecEpoch);
  double* precessCoordinates(double* RAdec, double JD);
  double hourAngle(double siderealTime, double* RAdecPrecessed);
  double* altitudeAzimuth(double hourAngle, double* RAdecPrecessed, double lattitude);
  double airMass(double hourAngle, double* RAdecPrecessed, double lattitude);
  double phase(double JD, double* ephemeris);
  //////////////////
  /////// Time and date computation methods
  int* calendarDateToYearDayLeap(int* calendarDate);
  int daysInMonth(int month, int year);
  bool checkLeapYear(int year);
  bool checkDaylightSavings(int* clockTime, int* calendarDate, int timeZone);
  int* UTtoLocalTime(int* UTclockTime, int* UTdate, int timeZone);
  int* reconcileCalendarDate(int* calendarDate);
  /////////////////
  /// String and simple time or date arithmetic methods
  string makeDateOrTimeString(int* dateOrTime, string type);
  string makeTripletString(double* triplet, string separator, int precision);
  int* processDateString(string dateString);
  string decimalHourToString(double time);
  double clockTimeToDecimalHours(int* clockTime);
  double* decimalHourToTriplet(double time);
};
  

  
  
