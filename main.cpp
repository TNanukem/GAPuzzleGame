#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <numeric>
#include <ctime>
#include <cassert>

#define BOARD_DIM 3

#define t_fitness float

using std::cout;
using std::cin;
using std::cerr;
using std::array;
using std::vector;
using std::ifstream;

typedef array<array<int, BOARD_DIM>, BOARD_DIM> board_array;

typedef vector<bool> chromosome_vector;

typedef vector<chromosome_vector> pop_vector;

typedef vector<t_fitness> fitness_vector;

typedef vector<float> probability_vector;


// Generate and return an initial random population of chromosomes
pop_vector generatePopulation(size_t initial_pop_size){

	pop_vector population;
	population.resize(initial_pop_size);

	for(auto& chromosome : population)
	{
		chromosome.resize(4);
		for(auto&& cell : chromosome)
		{
			if(rand() % 100 > 50)
				cell = 1;
			else
				cell = 0;
		}
	}

	return population;
}

// Check if any chromosome is unsolvable and, if so, removes it from the population
void check_solvability(pop_vector& population){

	// quem fizer essa função vai usar o método "erase" do vector

	// por exemplo, o código abaixo remove o cromossomo de posição igual a 5
	population.erase(population.begin() + 5);

	// outro exemplo, o código abaixo itera sobre a população e remove os cromossomos cujo tamanho é maior que 2
	for(auto elem = population.begin(); elem != population.end(); elem++){
		if((*elem).size() > 2)
			population.erase(elem);
	}

}

// Calculate and return the fitness of each chromosome
fitness_vector fitnessCalculation(pop_vector& population){

}

// Calculate the probability of selection for each chromosome based on its fitness value
probability_vector calculateProbabilities(const fitness_vector& fitness){

	t_fitness sum_of_fitness = std::accumulate(fitness.begin(), fitness.end(), 0);
	probability_vector probabilities;
	probabilities.resize(fitness.size());
	float previous_probability = 0;

	// Calculate the probability of each chromosome being selected,
	// i.e., the distance between his value and the previous one
	for(size_t i = 0; i < fitness.size(); i++){
		previous_probability += (fitness[i] / sum_of_fitness);
		probabilities[i] = previous_probability;
	}

	return probabilities;
}

// Select pairs of reproducers with probability proportional to its fitness values
pop_vector selectReproducers(const pop_vector& population, const fitness_vector& fitness){

	pop_vector parents;
	parents.resize(population.size());
	probability_vector probabilities = calculateProbabilities(fitness);

	// The population size must be even for the next "for" to work properly
	assert(population.size() % 2 == 0);

	// Select the parents based on previous probability calculation
	for(size_t i = 0; i < population.size(); i += 2){

		// Selecting the first parent
		float random = (float) rand() / RAND_MAX; // Random number between 0 and 1
		float previous = 0;
		for(size_t candidate = 0; candidate < probabilities.size(); candidate++){
			if(random >= previous && random < probabilities[candidate]){
				parents[i] = population[candidate];
				break;
			}
			if(random == probabilities[candidate])
				parents[i] = population[candidate];
		}

		// Selecting the second parent
		random = (float) rand() / RAND_MAX;
		previous = 0;
		for(size_t candidate = 0; candidate < probabilities.size(); candidate++){
			if(random >= previous && random < probabilities[candidate]){
				parents[i+1] = population[candidate];
				break;
			}
			if(random == probabilities[candidate])
				parents[i+1] = population[candidate];
		}

		// The parents need to be different
		while(parents[i] == parents[i+1]){
			random = (float) rand() / RAND_MAX;
			previous = 0;
			for(size_t candidate = 0; candidate < probabilities.size(); candidate++){
				if(random >= previous && random < probabilities[candidate]){
					parents[i+1] = population[candidate];
					break;
				}
				if(random == probabilities[candidate])
					parents[i+1] = population[candidate];
			}
		}
	}

	return parents;
}

// Cross every pair of consecutive parents, which will generate new pairs of chromosomes (a new population)
// The "generation" variable is used to control the number of board moves made by each chromosome
pop_vector reproducePopulation(const pop_vector& parents, float crossover_probability, unsigned generation){

	// para concatenar bits aos cromossomos, utilizem o método "push_back" do vector

	// por exemplo, o código abaixo cria uma população igual aos pais e concatena um 0 no final do 5º cromossomo:
	pop_vector population = parents;
	population[4].push_back(0);

}

// Flip one random bit of each chromosome according to a fixed probability rate
void mutatePopulation(pop_vector& population, float probability){

}

// Verifies if any chromosome on population have achivied a solution.
// If a solution have been found, fill the "solution" vector with the chosen solution
bool isSolved(const pop_vector& population, chromosome_vector& solution){

}

// Print the solution found
void printSolution(const chromosome_vector& solution){

}

int main (int argc, char* argv[]){

	srand(time(NULL));

	// Requires the puzzle game
	if(argc < 2){
		cerr << "Please, insert the path to the puzzle" << "\n";
		return 1;
	}

	// Reads the data and puts it on the board array
	ifstream file(argv[1]);
	if(!file){
		cerr << "There was an error loading the puzzle file" << "\n";
		return 1;
	}

	board_array board;

	for(auto& line : board){
		for(auto& col : line){
			file >> col;
			cout << col;
		}
		cout << "\n";
	}
	cout << "\n";

	size_t initial_pop_size = 100;

	// Generate a random population
	pop_vector population = generatePopulation(initial_pop_size);

	bool solved = false;
	chromosome_vector solution;

	const float crossover_probability = 0.04;
	const float mutation_probability = 0.004;

	unsigned generation = 1;

	while(solved == false){

		// Check if any chromosome is unsolvable and, if so, removes it from the population
		check_solvability(population);

		// Function to calculate the fitness of each candidate
		fitness_vector fitness = fitnessCalculation(population);

		// Selects the candidates to reproduce
		pop_vector parents = selectReproducers(population, fitness);

		// Candidates reproduction
		population = reproducePopulation(parents, crossover_probability, generation);

		// Verifies if some chromosome has achieved the solution
		// If so, the solution will be returned by reference
		solved = isSolved(population, solution);

		// Mutate the population only if a solution hasn't been found
		if(!solved)
			mutatePopulation(population, mutation_probability);

		generation++;
	}

	printSolution(solution);

	return 0;
}
