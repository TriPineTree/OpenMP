#ifndef _TRIP_H_
#define _TRIP_H_

#define CHROMOSOMES    50000 // 50000 different trips                              (DO NOT CHANGE)
#define CITIES         36    // 36 cities = ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789   (DO NOT CHANGE)
#define MAX_GENERATION 150   //                                                    (DO NOT CHANGE)
#define TOP_X          25000 // top 50%                                            (DO NOT CHANGE)
#define MUTATE_RATE    40    // Most optimal rate                                 
                                                                    
#define DEBUG          false // for debugging   

// Each trip (or chromosome) informatioin
class Trip {
public:
  char itinerary[CITIES + 1];  // a route through all 36 cities from (0, 0) 
  float fitness;               // the distance of this entire route
};

#endif
