#include <iostream>  // cout
#include <fstream>   // ifstream
#include <string.h>  // strncpy
#include <stdlib.h>  // rand
#include <math.h>    // sqrt, pow
#include <omp.h>     // OpenMP
#include <algorithm>
#include "Timer.h"
#include "Trip.h"


using namespace std;

// Already implemented. see the actual implementations below
void initialize(Trip trip[CHROMOSOMES], int coordinates[CITIES][2]);
void select(Trip trip[CHROMOSOMES], Trip parents[TOP_X]);
void populate(Trip trip[CHROMOSOMES], Trip offsprings[TOP_X]);


// need to implement for your program 1
extern void evaluate(Trip trip[CHROMOSOMES], int coordinates[CITIES][2]);
extern void crossover(Trip parents[TOP_X], Trip offsprings[TOP_X], int coordinates[CITIES][2]);
extern void mutate(Trip offsprings[TOP_X]);


// Helper functions
bool visited(Trip trip, char city, int maxIndex);
int compare(const void* a, const void* b) {
	Trip* tripA = (Trip*)a;
	Trip* tripB = (Trip*)b;
	return (tripA->fitness - tripB->fitness);
}
int getIndex(char city);
char getChar(int index);
float distance(int coordinates[CITIES][2], char cityA, char cityB);

/*
 * MAIN: usage: Tsp #threads
 */
int main(int argc, char* argv[]) {
	Trip trip[CHROMOSOMES];       // all 50000 different trips (or chromosomes)
	Trip shortest;                // the shortest path so far
	int coordinates[CITIES][2];   // (x, y) coordinates of all 36 cities:
	int nThreads = 1;

	// verify the arguments
	if (argc == 2)
		nThreads = atoi(argv[1]);
	else {
		cout << "usage: Tsp #threads" << endl;
		if (argc != 1)
			return -1; // wrong arguments
	}
	cout << "# threads = " << nThreads << endl;

	// shortest path not yet initialized
	shortest.itinerary[CITIES] = 0;  // null path
	shortest.fitness = -1.0;         // invalid distance

	// initialize 5000 trips and 36 cities' coordinates
	initialize(trip, coordinates);

	// start a timer 
	Timer timer;
	timer.start();

	// change # of threads
	omp_set_num_threads(nThreads);

	// find the shortest path in each generation
	for (int generation = 0; generation < MAX_GENERATION; generation++) {

		// evaluate the distance of all 50000 trips
		evaluate(trip, coordinates);
		// just print out the progress
		if (generation % 20 == 0)
			cout << "generation: " << generation << endl;

		// whenever a shorter path was found, update the shortest path
		if (shortest.fitness < 0 || shortest.fitness > trip[0].fitness) {

			strncpy(shortest.itinerary, trip[0].itinerary, CITIES);
			shortest.fitness = trip[0].fitness;

			cout << "generation: " << generation
				<< " shortest distance = " << shortest.fitness
				<< "\t itinerary = " << shortest.itinerary << endl;
		}
		// define TOP_X parents and offsprings.
		Trip parents[TOP_X], offsprings[TOP_X];

		// choose TOP_X parents from trip
		select(trip, parents);

		// generates TOP_X offsprings from TOP_X parenets
		crossover(parents, offsprings, coordinates);

		// mutate offsprings
		mutate(offsprings);

		// populate the next generation.
		populate(trip, offsprings);
	}

	// stop a timer
	cout << "elapsed time = " << timer.lap() << endl;
	return 0;
}

/*
 * Initializes trip[CHROMOSOMES] with chromosome.txt and coordiantes[CITIES][2] with cities.txt
 *
 * @param trip[CHROMOSOMES]:      50000 different trips
 * @param coordinates[CITIES][2]: (x, y) coordinates of 36 different cities: ABCDEFGHIJKLMNOPQRSTUVWXYZ
 */
void initialize(Trip trip[CHROMOSOMES], int coordinates[CITIES][2]) {
	// open two files to read chromosomes (i.e., trips)  and cities
	ifstream chromosome_file("chromosome.txt");
	ifstream cities_file("cities.txt");

	// read data from the files
	// chromosome.txt:                                                                                           
	//   T8JHFKM7BO5XWYSQ29IP04DL6NU3ERVA1CZG                                                                    
	//   FWLXU2DRSAQEVYOBCPNI608194ZHJM73GK5T                                                                    
	//   HU93YL0MWAQFIZGNJCRV12TO75BPE84S6KXD
	for (int i = 0; i < CHROMOSOMES; i++) {
		chromosome_file >> trip[i].itinerary;
		trip[i].fitness = 0.0;
	}

	// cities.txt:                                                                                               
	// name    x       y                                                                                         
	// A       83      99                                                                                        
	// B       77      35                                                                                        
	// C       14      64                                                                                        
	for (int i = 0; i < CITIES; i++) {
		char city;
		cities_file >> city;
		int index = (city >= 'A') ? city - 'A' : city - '0' + 26;
		cities_file >> coordinates[index][0] >> coordinates[index][1];
	}

	// close the files.
	chromosome_file.close();
	cities_file.close();

	// just for debugging
	if (DEBUG) {
		for (int i = 0; i < CHROMOSOMES; i++)
			cout << trip[i].itinerary << endl;
		for (int i = 0; i < CITIES; i++)
			cout << coordinates[i][0] << "\t" << coordinates[i][1] << endl;
	}
}
/*
 * Select the first TOP_X parents from trip[CHROMOSOMES]
 *
 * @param trip[CHROMOSOMES]: all trips
 * @param parents[TOP_X]:    the firt TOP_X parents
 */
void select(Trip trip[CHROMOSOMES], Trip parents[TOP_X]) {
	// just copy TOP_X trips to parents
#pragma omp parallel for
	for (int i = 0; i < TOP_X; i++)
		strncpy(parents[i].itinerary, trip[i].itinerary, CITIES + 1);
}

/*
 * Replace the bottom TOP_X trips with the TOP_X offsprings
 */
void populate(Trip trip[CHROMOSOMES], Trip offsprings[TOP_X]) {
	// just copy TOP_X offsprings to the bottom TOP_X trips.
#pragma omp parallel for
	for (int i = 0; i < TOP_X; i++)
		strncpy(trip[CHROMOSOMES - TOP_X + i].itinerary, offsprings[i].itinerary, CITIES + 1);
}

/*
 * Evaluate the distance of all 50000 trips
 *
 * @param trip[CHROMOSOMES]: all trips
 * @param int coordinates[CITIES][2]: coordinates of all cities
 */
void evaluate(Trip trip[CHROMOSOMES], int coordinates[CITIES][2]) {
#pragma omp parallel for
	for (int i = 0; i < CHROMOSOMES; i++) {
		// Reset the trip fitness to 0 to reevaluate
		trip[i].fitness = 0;
		// Add the distance from the origin to the first city
		trip[i].fitness += sqrt(pow(coordinates[getIndex(trip[i].itinerary[0])][0], 2)
			+ pow(coordinates[getIndex(trip[i].itinerary[0])][1], 2));
		// Adding distance to the total (fitness)
		for (int j = 1; j < CITIES; j++) {
			trip[i].fitness += distance(coordinates, trip[i].itinerary[j - 1], trip[i].itinerary[j]);
		}
	}
	qsort(trip, CHROMOSOMES, sizeof(*trip), compare);
}

/*
 * generates TOP_X offsprings from TOP_X parenets
 *
 * @param parents[TOP_X]: the first TOP_X parents
 * @param offsprings[TOP_X]: the first TOP_X childs
 * @param int coordinates[CITIES][2]: coordinates of all cities
 */
void crossover(Trip parents[TOP_X], Trip offsprings[TOP_X], int coordinates[CITIES][2]) {
#pragma omp parallel for
	for (int i = 0; i < TOP_X; i += 2) {
		// assign parent's first city as child[i]'s first city
		offsprings[i].itinerary[0] = parents[i].itinerary[0];
		for (int j = 1; j < CITIES; j++) {
			// If 2 cities are not visited, choose the shorter distance one
			if (!visited(offsprings[i], parents[i].itinerary[j], j) &&
				!visited(offsprings[i], parents[i + 1].itinerary[j], j)) {
				// compare distance
				if (distance(coordinates, offsprings[i].itinerary[j - 1], parents[i].itinerary[j]) <
					distance(coordinates, offsprings[i].itinerary[j - 1], parents[i + 1].itinerary[j])) {
					// assign next city to child[i] from parent[i]
					offsprings[i].itinerary[j] = parents[i].itinerary[j];
				}
				else {
					// assign next city to child[i] from parent[i+1]
					offsprings[i].itinerary[j] = parents[i + 1].itinerary[j];
				}
			}
			// If parent[i]'s next is visited and parent[i+1] not
			else if (visited(offsprings[i], parents[i].itinerary[j], j) &&
				!visited(offsprings[i], parents[i + 1].itinerary[j], j)) {
				// child[i] 
				offsprings[i].itinerary[j] = parents[i + 1].itinerary[j];
			}
			// If parent[i+1]'s next is visited and parent[i] not
			else if (visited(offsprings[i], parents[i + 1].itinerary[j], j) &&
				!visited(offsprings[i], parents[i].itinerary[j], j)) {
				// child[i]
				offsprings[i].itinerary[j] = parents[i].itinerary[j];
			}
			else {
				// Choose a random city from parents[i]
				int index = rand() % CITIES;
				for (int z = index; z < CITIES; z++) {
					if (!visited(offsprings[i], parents[i].itinerary[z], j)) {
						// child[i]
						offsprings[i].itinerary[j] = parents[i].itinerary[z];
						break;
					}
					// Loop to the beginning to check for the skipped cities
					if (z == CITIES - 1) {
						z = -1;
					}
					if (z == index - 1) {
						break;
					}

				}
			}
		}
		// Create child[i+1] as compliment to child[i]
		for (int c = 0; c < CITIES; c++) {
			offsprings[i + 1].itinerary[c] = getChar(CITIES - getIndex(offsprings[i].itinerary[c]) - 1);
		}
	}
}

/*
 * mutate offsprings
 *
 * @param offsprings[TOP_X]: the first TOP_X childs
 */
void mutate(Trip offsprings[TOP_X]) {
	for (int i = 0; i < TOP_X; i++) {
		int percentage = rand() % 100;
		if (percentage <= MUTATE_RATE) {
			int rand1 = rand() % (CITIES - 1);
			int rand2 = rand() % (CITIES - 1);
			char temp = offsprings[i].itinerary[rand1];
			offsprings[i].itinerary[rand1] = offsprings[i].itinerary[rand2];
			offsprings[i].itinerary[rand2] = temp;
		}
	}
}

/*
 * Iterate through the current trip given and check if the city is visited
 * 
 * @param trip: the current trip
 * @param city: the city to check
 * @param maxIndex: the current length of the child
 */
bool visited(Trip trip, char city, int maxIndex) {
	for (int i = 0; i < maxIndex; i++) {
		if (trip.itinerary[i] == city) {
			return true;
		}
	}
	return false;
}
/*
 * Return the index of the given city in the cities list
 * 
 * @param city: the city to get index of
 */
int getIndex(char city) {
	return (int)((city >= 'A') ? city - 'A' : city - '0' + 26);
}

/*
 * Return the city name of the given index in the cities list
 * 
 * @param: the index to get city of
 */
char getChar(int index) {
	return (char)((index < 26) ? index + 'A' : index - 26 + '0');
}
/*
 * Calculate the distance between 2 points using pythagorean theorem
 * 
 * @param coordinates[CITIES][2]: cities coordinates
 * @param cityA: the current city
 * @param cityB: the next city
 */
float distance(int coordinates[CITIES][2], char cityA, char cityB) {
	int a = getIndex(cityA);
	int b = getIndex(cityB);

	return  sqrt(pow(coordinates[b][0] - coordinates[a][0], 2)
		+ pow(coordinates[b][1] - coordinates[a][1], 2));
}
