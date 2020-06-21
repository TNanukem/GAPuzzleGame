#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <numeric>
#include <ctime>
#include <cassert>
#include <algorithm>

#define BOARD_DIM 3
#define DOWN 0
#define UP 1
#define RIGHT 2
#define LEFT 3

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

typedef vector<array<int, 2>> board_position;


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

board_array generateFinalBoard(vector<int> moves, board_array board){

	// Passar isso como parâmetro depois
	//board_array board = {{5,4,2},{8,6,7},{1,3,0}};
	int old_value = -1;
	array<int, 2> blank_position;

	// Finding the blank
	for(int i = 0; i < BOARD_DIM; i++){
		for(int j = 0; j < BOARD_DIM; j++){
			if(board[i][j] == 0){
				blank_position = {i, j};
			}
		}
	}

	for(int i = 0; i < moves.size(); i++){
		switch (moves[i]){
		
		case UP:
			old_value = board[blank_position[0]-1][blank_position[1]];
			board[blank_position[0]-1][blank_position[1]] = 0;
			board[blank_position[0]][blank_position[1]] = old_value;
			blank_position[0] -= 1; 
			break;
		
		case DOWN:
			old_value = board[blank_position[0]+1][blank_position[1]];
			board[blank_position[0]+1][blank_position[1]] = 0;
			board[blank_position[0]][blank_position[1]] = old_value;
			blank_position[0] += 1; 
			break;

		case LEFT:
			old_value = board[blank_position[0]][blank_position[1]-1];
			board[blank_position[0]][blank_position[1]-1] = 0;
			board[blank_position[0]][blank_position[1]] = old_value;
			blank_position[1] -= 1; 
			break;

		case RIGHT:
			old_value = board[blank_position[0]][blank_position[1]+1];
			board[blank_position[0]][blank_position[1]+1] = 0;
			board[blank_position[0]][blank_position[1]] = old_value;
			blank_position[1] += 1;
			break;
		
		}
	}
	return board;
}

// Given a chromosome, goes down on the tree to grab the final board for him
board_array downOnTree(chromosome_vector& chromosome, board_array board){

	// Corner -> Only one possibility ((0,2),(2,0),(0,0),(2,2))
	// Edge -> Two possibilities ((0,1), (1,0), (2,1), (1,2))
	// Middle -> Three possibilities ((1,1))

	board_position corner = {{0,2}, {2,0}, {0,0}, {2,2}};
	board_position edge = {{0,1}, {1,0}, {2,1}, {1,2}};
	board_position middle = {{1,1}};

	array<int, 2> blank_position = {2, 2};
	int previous_move = -1;

	vector<int> moves = {};

	for(int i = 1; i < chromosome.size(); i += 2){
		
		vector<int> possible_moves = {UP, DOWN, LEFT, RIGHT};

		int node = -1;
		int cells = 2 * chromosome[i] + chromosome[i-1];

		// Verify where the blank is so we know which number to module for
		// And if this is the first move, the module sums 1
		if (find(corner.begin(), corner.end(), blank_position) != corner.end())
			
			if(previous_move != -1)
				node = cells % 1;
			else node = cells % 2;
			
		else if (find(edge.begin(), edge.end(), blank_position) != edge.end())
			
			if(previous_move != -1)
				node = cells % 2;
			else node = cells % 3;
		
		else 
			
			if(previous_move != -1)
				node = cells % 3;
			else node = cells % 4;

		// Filters the possible moves
		if(previous_move != -1){
			// If not the first move, we do not allow the opposite of the last move to happen
			if(previous_move == UP)
				possible_moves.erase(possible_moves.begin()+1);
			else if (previous_move == DOWN)
				possible_moves.erase(possible_moves.begin());
			else if (previous_move == LEFT)
				possible_moves.erase(possible_moves.begin()+3);
			else if (previous_move == RIGHT)
				possible_moves.erase(possible_moves.begin()+2);
		}

		// Out of the possible moves, see which can be made given blank position
		if(blank_position[0] == 0){
			// We cannot go up
			cout << "Cannot go up\n";
			possible_moves.erase(std::remove(possible_moves.begin(), possible_moves.end(), UP), possible_moves.end());
		}
		if(blank_position[0] == BOARD_DIM-1){
			// We cannot go down
			cout << "Cannot go down\n";
			possible_moves.erase(std::remove(possible_moves.begin(), possible_moves.end(), DOWN), possible_moves.end());
		}
		if(blank_position[1] == 0){
			// We cannot go left
			cout << "Cannot go left\n";
			possible_moves.erase(std::remove(possible_moves.begin(), possible_moves.end(), LEFT), possible_moves.end());
		}
		if(blank_position[1] == BOARD_DIM-1){
			// We cannot go right
			cout << "Cannot go right\n";
			possible_moves.erase(std::remove(possible_moves.begin(), possible_moves.end(), RIGHT), possible_moves.end());
		}

		// The move becomes the possible move of index node

		previous_move = possible_moves[node];
		moves.push_back(previous_move);

		if(previous_move == UP)
			blank_position[0] -= 1; 
		if(previous_move == DOWN)
			blank_position[0] += 1;
		if(previous_move == LEFT)
			blank_position[1] -= 1;
		if(previous_move == RIGHT)
			blank_position[1] += 1;
	}

	cout << "Chromosome sequence of events: ";
	for (int i = 0; i < moves.size(); i++){
		cout << moves[i] << " ";
	}
	cout << "\n";

	board_array final_board = generateFinalBoard(moves, board);

	for(int i = 0; i < BOARD_DIM; i++){
		for(int j = 0; j < BOARD_DIM; j++){
			cout << final_board[i][j] << " ";
		}
		cout << "\n";
	}

	return final_board;
}

// Calculate and return the fitness of each chromosome
fitness_vector fitnessCalculation(pop_vector& population, board_array board){

	// The fitness of the chromosome is the inverse of the number of wrong tiles its final board
	// has in relation to the target.

	fitness_vector fitness_;
	int target[BOARD_DIM][BOARD_DIM] = {{1,2,3},{4,5,6},{7,8,0}};

	for (auto &chromosome : population){
		board_array chromosome_final_board = downOnTree(chromosome, board);

		float error = 0.0;
		for(int i = 0; i < BOARD_DIM; i++){
			for(int j = 0; j < BOARD_DIM; j++){
				if(chromosome_final_board[i][j] != target[i][j])
					error += 1;
			}
		}
		float fitness = 1/error;
		//cout << fitness << "\n";
		fitness_.push_back(fitness);
		cout << "Fitness: " << fitness << "\n\n";
	}

	return fitness_;
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
	
	for (auto &chromosome : population){
		for (auto &&cell : chromosome)
		{
			if (rand() % 100 > (probability*100)){
				// Bit flipping
				cell = 1 - cell;
				// If flips, then no other change is made on the chromosome
				break;
			}
		}
	}
}

// Verifies if any chromosome on population have achieved a solution.
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
	cout << "Teste\n";
	// Generate a random population
	pop_vector population = generatePopulation(initial_pop_size);
	
	bool solved = false;
	chromosome_vector solution;

	const float crossover_probability = 0.04;
	const float mutation_probability = 0.004;

	unsigned generation = 1;


	while(solved == false){

		cout << "Running generation " << generation << "\n";

	// 	// Check if any chromosome is unsolvable and, if so, removes it from the population
	// 	check_solvability(population);

	// 	// Function to calculate the fitness of each candidate
		fitness_vector fitness = fitnessCalculation(population, board);

	// 	// Selects the candidates to reproduce
	// 	pop_vector parents = selectReproducers(population, fitness);

	// 	// Candidates reproduction
	// 	population = reproducePopulation(parents, crossover_probability, generation);

	// 	// Verifies if some chromosome has achieved the solution
	// 	// If so, the solution will be returned by reference
	// 	solved = isSolved(population, solution);

	// 	// Mutate the population only if a solution hasn't been found
	// 	if(!solved)
	 		mutatePopulation(population, mutation_probability);

		generation++;
		break;
	}

	// printSolution(solution);

	return 0;
}
