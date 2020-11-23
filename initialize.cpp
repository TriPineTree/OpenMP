#include <iostream>  // cout
#include <fstream>   // ofstream
#include <string.h>  // strncmp
#include <stdlib.h>  // rand

#define CHROMOSOMES    50000 // 50000
#define CITIES         36    // Cities = ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789
#define DEBUG          false // for debugging

using namespace std;

char getCity( );
void initialize( char trip[][CITIES + 1], int coordinates[CITIES][2], 
		 int nChromosomes );

int main( int argc, char* argv[] ) {
  // default values
  int nChromosomes = CHROMOSOMES;

  // argument verification
  if ( argc == 2 ) {
    nChromosomes = atoi( argv[1] );
  }
  else {
    cout << "usage: initialize nChromosomes" << endl;
    if ( argc != 1 )
      exit( -1 );
  }
  cout << "# chromosomes = " << nChromosomes 
       << ", # cities = " << CITIES << endl;

  // declare chromosomes and cities
  char trip[nChromosomes][CITIES + 1];
  int coordinates[CITIES][2]; 

  // open two files to store chromosomes and cities
  ofstream chromosome_file( "chromosome.txt" );
  ofstream cities_file( "cities.txt" );

  // initialize chormosomes and cities
  initialize( trip, coordinates, nChromosomes );
  
  // write data to the two files
  // chromosome.txt:
  //   T8JHFKM7BO5XWYSQ29IP04DL6NU3ERVA1CZG
  //   FWLXU2DRSAQEVYOBCPNI608194ZHJM73GK5T
  //   HU93YL0MWAQFIZGNJCRV12TO75BPE84S6KXD
  for ( int i = 0; i < nChromosomes; i++ )
    chromosome_file << trip[i] << endl;

  // cities.txt:
  // name    x       y
  // A       83      99
  // B       77      35
  // C       14      64
  for ( int i = 0; i < CITIES; i++ ) {
    char city_name = ( i < 26 ) ? i + 'A' : i - 26 + '0';
    cities_file << city_name 
		<< "\t" << coordinates[i][0] << "\t" << coordinates[i][1] 
		<< endl;
  }

  // close the files.
  chromosome_file.close( );
  cities_file.close( );
    
  return 0;
}

/*
 * Returns a randome city name, one of: ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789
 */
char getCity( ) {
  char city = rand( ) % CITIES;
  if ( city < 26 )
    city += 'A';
  else
    city = city - 26 + '0';
  return city;
}

/*
 * Initialize nChoromosomes number of trips and nCities of coordinates.
 *
 * @param trip:         all chromosomes to be initialized
 * @param coordinates:  all cities to be initialized for x and y
 * @param nChromosomes: # of chromosomes
 */
void initialize( char trip[][CITIES + 1], int coordinates[CITIES][2], 
		 int nChromosomes ) {
  // initialize chromosomes
  for ( int chrom = 0; chrom < nChromosomes; chrom++ ) {

    // show the progress
    if ( chrom % 10000 == 0 )
      cout << chrom << " initialized" << endl;

    // add the very last delimiter
    trip[chrom][CITIES] = 0;
    
    while ( true ) {
      // initialize each trip
      for ( int cities = 0; cities < CITIES; cities++ ) {
	while ( true ) {
	  char candidate = getCity( );
	  
	  // check if there is the same city
	  if ( cities > 0 ) {
	    for ( int prev = 0; prev < cities; prev++ )
	      if ( trip[chrom][prev] == candidate ) {
		// found the same city on the route
		candidate = 0; // make it void
		break;
	      }
	  }

	  if ( candidate != 0 ) {
	    trip[chrom][cities] = candidate;
	    break;          // go to the next city
	  }
	  // else get another city
	}
      }
      
      // check if there is the same trip
      bool found = false;
      if ( chrom > 0 ) {
	for ( int prev = 0; prev < chrom; prev++ )
	  if ( found = ( strncmp( trip[prev], trip[chrom], CITIES ) == 0 ) ) {
	    break;
	  }
      }
      if ( found )
	continue;        // get another trip
      else
	break;           // go to the next trip
    }

    if ( DEBUG ) 
      cout << "chrom[" << chrom << "] = " << trip[chrom] << endl;
  }

  // initialize each city's x and y coordinates
  for ( int cities = 0; cities < CITIES; cities++ ) {
    while ( true ) {
      coordinates[cities][0] = rand( ) % 100;
      coordinates[cities][1] = rand( ) % 100;

      // check if there is the same coordiante
      bool found = false;
      if ( cities > 0 ) {
	for ( int prev = 0; prev < cities; prev++ )
	  if ( found = ( coordinates[prev][0] == coordinates[cities][0] && 
			 coordinates[prev][1] == coordinates[cities][1] ) )
	    // found the same coordinate
	    break;
      }
      if ( found )
	continue; // get another coordinate
      else
	break;    // go to the next city
    }
    if ( DEBUG )
      cout << "coordinates[" 
	   << ( char )( ( cities < 26 ) ? cities + 'A' : cities - 26 + '0' ) 
	   << "] = ("
	   << coordinates[cities][0] << ", " << coordinates[cities][1] << ")" 
	   << endl;
  }
}

