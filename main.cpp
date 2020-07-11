#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <numeric>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include <random>

// Board dimensions
#define BOARD_DIM 3

// Possible moves of the blank tile
#define DOWN 0
#define UP 1
#define RIGHT 2
#define LEFT 3

// Assert with messages
#define assertm(exp, msg) assert(((void)msg, exp))

using std::cout;
using std::cin;
using std::cerr;
using std::array;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::string;

typedef double t_fitness;

typedef double t_probability;

typedef array<array<int, BOARD_DIM>, BOARD_DIM> board_array;

typedef vector<bool> chromosome_vector;

typedef vector<chromosome_vector> pop_vector;

typedef vector<t_fitness> fitness_vector;

typedef vector<t_probability> probability_vector;

typedef vector<array<int, 2>> board_position;

bool solved = false;
vector<int> final_moves;
board_array target = {{{1,2,3},{4,5,6},{7,8,0}}};

// To generate random numbers
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

// Generate and return an initial random population of chromosomes
pop_vector generatePopulation(const size_t initial_pop_size, const size_t initial_chromosome_size){

	// The initial chromosome size must be even because we use 2 bits for each gene/move
	assertm(initial_chromosome_size % 2 == 0, "The initial chromosome size must be even");

	pop_vector population;
	population.resize(initial_pop_size);

	// Generate random int numbers between 1 and 100
	std::uniform_int_distribution<> distrib(1, 100);

	for(auto& chromosome : population)
	{
		chromosome.resize(initial_chromosome_size);
		for(auto&& cell : chromosome)
		{
			float random = distrib(gen);
			if(random > 50)
				cell = 1;
			else
				cell = 0;
		}
	}
	return population;
}

// Given an initial state and a list of movements, returns the final board obtained
board_array generateFinalBoard(const vector<int>& moves, board_array board){

	int old_value = -1;
	array<int, 2> blank_position;
	vector<int> aux;

	// Finding the blank position on the board
	for(int i = 0; i < BOARD_DIM; i++){
		for(int j = 0; j < BOARD_DIM; j++){
			if(board[i][j] == 0){
				blank_position = {i, j};
			}
		}
	}

	// For each move, we move the blank position to the right spot
	for(size_t i = 0; i < moves.size(); i++){
		aux.push_back(moves[i]);
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

		// Verifies if the obtained board is a solution by comparing it to the target board
		bool aux_sol = true;
		bool outer_loop = true;
		for(int j = 0; j < BOARD_DIM && outer_loop; j++){
			for(int k = 0; k < BOARD_DIM; k++){
				if( board[j][k] != target[j][k] ){
					aux_sol = false;
					outer_loop = false;
					break;
				}
			}
		}

		// If the solution is found, the execution stops and the board is returned as well as the list of movements
		if(aux_sol == true){

			for(size_t k = 0; k < aux.size(); k++){
				final_moves.push_back(aux[k]);
			}

			solved = true;
			return board;
		}
	}
	return board;
}

// Given a chromosome, goes down on the tree keeping each movement made
board_array downOnTree(const chromosome_vector& chromosome, const board_array& board, vector<int>* mov){

	// If the solution was found, we return
	if(solved == true){
		return board;
	}

	// Corner -> Only one possibility ((0,2),(2,0),(0,0),(2,2))
	// Edge -> Two possibilities ((0,1), (1,0), (2,1), (1,2))
	// Middle -> Three possibilities ((1,1))

	board_position corner = {{0,2}, {2,0}, {0,0}, {2,2}};
	board_position edge = {{0,1}, {1,0}, {2,1}, {1,2}};
	board_position middle = {{1,1}};

	array<int, 2> blank_position;

	// Finding the blank position
	for (int i = 0; i < BOARD_DIM; i++)
	{
		for (int j = 0; j < BOARD_DIM; j++)
		{
			if (board[i][j] == 0)
			{
				blank_position = {i, j};
			}
		}
	}

	int previous_move = -1;

	vector<int> moves = {};

	for(size_t i = 1; i < chromosome.size(); i += 2){

		vector<int> possible_moves = {UP, DOWN, LEFT, RIGHT};

		int node = -1;
		int cells = 2 * chromosome[i] + chromosome[i-1];

		// Verify where the blank is so we know which number to module for
		// And if this is the first move, the module needs to be summed by 1 because there is no
		// previous movement
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
			possible_moves.erase(std::remove(possible_moves.begin(), possible_moves.end(), UP), possible_moves.end());
		}
		if(blank_position[0] == BOARD_DIM-1){
			// We cannot go down
			possible_moves.erase(std::remove(possible_moves.begin(), possible_moves.end(), DOWN), possible_moves.end());
		}
		if(blank_position[1] == 0){
			// We cannot go left
			possible_moves.erase(std::remove(possible_moves.begin(), possible_moves.end(), LEFT), possible_moves.end());
		}
		if(blank_position[1] == BOARD_DIM-1){
			// We cannot go right
			possible_moves.erase(std::remove(possible_moves.begin(), possible_moves.end(), RIGHT), possible_moves.end());
		}

		// The move becomes the possible move of index node

		previous_move = possible_moves[node];
		moves.push_back(previous_move);

		// Updates the blank position
		if(previous_move == UP)
			blank_position[0] -= 1;
		if(previous_move == DOWN)
			blank_position[0] += 1;
		if(previous_move == LEFT)
			blank_position[1] -= 1;
		if(previous_move == RIGHT)
			blank_position[1] += 1;
	}


	// Retrieves the final board after all the movements
	board_array final_board = generateFinalBoard(moves, board);

	if(solved == true){
		*mov = final_moves;
	}

	return final_board;
}

// Calculate and return the fitness of each chromosome in a Manhattan fashion
fitness_vector fitnessCalculationManhattan(const pop_vector& population, const board_array& board, vector<int>* mov){

	// The fitness of the chromosome is the sum of manhattan distances of every instance of the final board
	// to the target position.

	fitness_vector fitness_;

	int pos_x, pos_y;

	for (const auto &chromosome : population){
		board_array chromosome_final_board = downOnTree(chromosome, board, mov);

		float error = 0.0;
		for(int i = 0; i < BOARD_DIM; i++){
			for(int j = 0; j < BOARD_DIM; j++){
				pos_x = pos_y = -1;

				// Finds the position of the target board
				for(int k = 0; k < BOARD_DIM; k++){
					for(int w = 0; w < BOARD_DIM; w++){
						if(chromosome_final_board[i][j] == target[k][w]){
							pos_x = k;
							pos_y = w;
						}
					}
				}

				error += abs((i - pos_x)) + abs((j - pos_y));
			}
		}
		if(error == 0)
			error = 1;
		float fitness = 1/error;
		fitness_.push_back(fitness);
	}

	return fitness_;
}

// Calculate and return the fitness of each chromosome in a simple fashion
fitness_vector fitnessCalculationSimple(const pop_vector& population, const board_array& board, vector<int>* mov){

	// The fitness of the chromosome is the inverse of the number of wrong tiles its final board
	// has in relation to the target.

	fitness_vector fitness_;

	for (auto &chromosome : population){
		board_array chromosome_final_board = downOnTree(chromosome, board, mov);

		float error = 0.0;
		for(int i = 0; i < BOARD_DIM; i++){
			for(int j = 0; j < BOARD_DIM; j++){
				if(chromosome_final_board[i][j] != target[i][j])
					error += 1;
			}
		}
		if(error == 0)
			error = 1;
		float fitness = 1/error;
		fitness_.push_back(fitness);
	}

	return fitness_;
}

// Calculate the probability of selection for each chromosome based on its fitness value
probability_vector calculateProbabilities(const fitness_vector& fitness){

	t_fitness sum_of_fitness = std::accumulate(fitness.begin(), fitness.end(), 0.0);
	probability_vector probabilities;
	probabilities.resize(fitness.size());
	t_probability previous = 0;

	// Calculate the probability of each chromosome being selected,
	// i.e., the distance between his value and the previous one
	for(size_t i = 0; i < fitness.size(); i++){
		probabilities[i] = previous + (fitness[i] / sum_of_fitness);
		previous = probabilities[i];
	}

	return probabilities;
}

// Select pairs of reproducers with probability proportional to its fitness values
pop_vector selectReproducers(const pop_vector& population, const fitness_vector& fitness){

	pop_vector parents;
	parents.resize(population.size());
	probability_vector probabilities = calculateProbabilities(fitness);

	// The population size must be even for the next "for" to work properly
	assertm(population.size() % 2 == 0, "The population size must be even");

	// Generate random real numbers between 0 and 1
	static std::uniform_real_distribution<> distrib(0.0, 1.0);

	// Select the parents based on previous probability calculation
	for(size_t i = 0; i < population.size(); i += 2){

		// Selecting the first parent
		t_probability random = distrib(gen);
		t_probability previous = 0;
		size_t candidate;
		for(candidate = 0; candidate < probabilities.size(); candidate++){
			if(random >= previous && random < probabilities[candidate]){
				parents[i] = population[candidate];
				break;
			}
			previous = probabilities[candidate];
		}
		if(random == probabilities[candidate])
			parents[i] = population[candidate];

		// Selecting the second parent
		random = distrib(gen);
		previous = 0;
		for(candidate = 0; candidate < probabilities.size(); candidate++){
			if(random >= previous && random < probabilities[candidate]){
				parents[i+1] = population[candidate];
				break;
			}
			previous = probabilities[candidate];
		}
		if(random == probabilities[candidate])
			parents[i+1] = population[candidate];

		// The parents need to be different
		while(parents[i] == parents[i+1]){
			random = distrib(gen);
			previous = 0;
			for(candidate = 0; candidate < probabilities.size(); candidate++){
				if(random >= previous && random < probabilities[candidate]){
					parents[i+1] = population[candidate];
					break;
				}
				previous = probabilities[candidate];
			}
			if(random == probabilities[candidate])
				parents[i+1] = population[candidate];
		}
	}

	return parents;
}

//Function that adds random binary numbers at the end of each child
//Duplicates the size of the chromossome for each child
chromosome_vector generateTail(chromosome_vector child)
{
	// bit number that will be appended to tail
	bool bin;

	size_t tail_size = 2;

	// Generate random bit numbers (0 or 1)
	static std::uniform_int_distribution<> distrib(0, 1);

	for (size_t i = 0; i < tail_size; i++)
	{
		bin = distrib(gen);
		child.push_back(bin);
	}
	return child;
}

//Generates the crossover for the next generation
pop_vector crossOver(const chromosome_vector& chromossomeX, const chromosome_vector& chromossomeY, size_t crossover_point)
{
	chromosome_vector childX;
	chromosome_vector childY;
	pop_vector child;
	size_t tam = chromossomeX.size();

	//Create copies of the chromossome until the crossover point
	for(size_t i = 0; i < crossover_point; i++)
	{
		childX.push_back(chromossomeX[i]);
		childY.push_back(chromossomeY[i]);
	}
	//Complete the children adding binaries from the parent after the crossover point
	for (size_t i = crossover_point; i < tam; i++)
	{
		childX.push_back(chromossomeY[i]);
		childY.push_back(chromossomeX[i]);
	}

	child.push_back(childX);
	child.push_back(childY);

	return child;
}

// Cross every pair of consecutive parents, which will generate new pairs of chromosomes (a new population)
// The crossover is done by a single point of division
pop_vector reproducePopulation(const pop_vector& parents, float crossover_probability, bool grow_chromosome){

	float random;
	size_t tam = parents[0].size();
	size_t crossover_point;

	//Array that will store the new population
	pop_vector children;

	//Auxiliary array to receive children from crossover function
	pop_vector aux;

	// Generate random real numbers between 0 and 1
	static std::uniform_real_distribution<> distrib_float(0.0, 1.0);

	// Generate random int numbers between 0 and tam-1
	std::uniform_int_distribution<> distrib_int(0, tam-1);

	for(size_t i = 0; i < parents.size(); i += 2)
	{
		// Random number between 0 and 1
		random = distrib_float(gen);

		//Check if the crossover will happen
		if(random < crossover_probability)
		{
			crossover_point = distrib_int(gen);
			aux = crossOver(parents[i], parents[i+1], crossover_point);
			chromosome_vector child1 = grow_chromosome ? generateTail(aux[0]) : aux[0];
			chromosome_vector child2 = grow_chromosome ? generateTail(aux[1]) : aux[1];
			children.push_back(child1);
			children.push_back(child2);
		}
		//The child is the parent
		else
		{
			chromosome_vector child1 = grow_chromosome ? generateTail(parents[0]) : parents[0];
			chromosome_vector child2 = grow_chromosome ? generateTail(parents[1]) : parents[1];
			children.push_back(child1);
			children.push_back(child2);
		}
	}

	return children;
}

// Replace a random children with the best parent
void elitism(const pop_vector& population, pop_vector& children, size_t best_id)
{
	size_t tam_population = population.size();

	// Generate random int numbers between 0 and tam_population-1
	std::uniform_int_distribution<> distrib_int_p(0, tam_population-1);

	size_t to_replace = distrib_int_p(gen);
	children[to_replace] = population[best_id];
}

// Flip one random bit of each chromosome according to a fixed probability rate
pop_vector mutatePopulation(pop_vector population, float probability){

	// Generate random real numbers between 0 and 1
	static std::uniform_real_distribution<> distrib_float(0.0, 1.0);

	// The cell that will be mutated
	size_t cell;
	for (auto &chromosome : population){
		// Generate random int numbers between 0 and the size of the chromosome - 1
		std::uniform_int_distribution<> distrib_int(0, chromosome.size()-1);
		float random = distrib_float(gen);

		if (random < probability){
			cell = distrib_int(gen);
			// Flip the bit at the "cell" position
			chromosome[cell] = 1 - chromosome[cell];
		}
	}
	return population;
}

// Print the solution found
void printSolution(const vector<int>& moves){

	cout << "\nSolution was found with " << moves.size() << " movements:\n";
	for(size_t i = 0; i < moves.size(); i++){

		switch (moves[i]){
			case UP:
				cout << "[UP] ";
				break;
			case DOWN:
				cout << "[DOWN] ";
				break;
			case LEFT:
				cout << "[LEFT] ";
				break;
			case RIGHT:
				cout << "[RIGHT] ";
				break;
		}
		if( i % 10 == 0 && i != 0){
			cout << "\n";
		}
	}
	cout << "\n";
}
//Function that executes each movement from the final solution
//Used to test if the movements are correct
void printBoardMoves(const vector<int> moves,board_array board)
{
	//Gets the row and collumn for blanck space
	int linha,coluna,aux;
	bool outer_loop = true;
	for (size_t i = 0; i < BOARD_DIM && outer_loop; i++)
	{
		for (size_t j = 0; j < BOARD_DIM; j++)
		{
			if(board[i][j]==0)
			{
				linha = i;
				coluna = j;
				outer_loop = false;
				break;
			}
		}
	}
	//Executes the movements
	for(size_t i = 0; i < moves.size(); i++){

		if(moves[i] == UP)
		{
			cout <<"UP\n";
			aux = board[linha-1][coluna];
			board[linha-1][coluna] = board[linha][coluna];
			board[linha][coluna] = aux;
			linha -=1;
		}
		if(moves[i] == DOWN)
		{
			cout <<"DOWN\n";
			aux = board[linha+1][coluna];
			board[linha+1][coluna] = board[linha][coluna];
			board[linha][coluna] = aux;
			linha+=1;
		}
		if(moves[i] == LEFT)
		{
			cout <<"LEFT\n";
			aux = board[linha][coluna-1];
			board[linha][coluna-1] = board[linha][coluna];
			board[linha][coluna] = aux;
			coluna-=1;
		}
		if(moves[i] == RIGHT)
		{
			cout <<"RIGHT\n";
			aux = board[linha][coluna+1];
			board[linha][coluna+1] = board[linha][coluna];
			board[linha][coluna] = aux;
			coluna+=1;
		}
		for (size_t i = 0; i < BOARD_DIM; i++)
		{
			for (size_t j = 0; j < BOARD_DIM; j++)
			{
				cout << board[i][j];
			}
			cout<<"\n";
		}
		cout<<"\n";
	}

}
// Calculate the sum inversion for each position of the board
unsigned calculateInversion(const board_array& board, size_t x, size_t y)
{
	unsigned sum = 0;
	size_t aux = y;
	for (size_t i = x; i < BOARD_DIM; i++)
	{
		for (size_t j = aux; j < BOARD_DIM; j++)
		{
			// If the tile is not blank (zero)
			if(board[i][j] != 0)
			{
				if(board[i][j] < board[x][y])
				{
					sum += 1;
				}
			}
		}
		aux = 0;
	}
	return sum;
}

//Calculates the total sum of the board
//Returns if the board is solvable or not
bool isSolvable(const board_array& board)
{
	unsigned sum = 0;
	for(size_t i = 0; i < BOARD_DIM; i++)
	{
		for (size_t j = 0; j < BOARD_DIM; j++)
		{
			// If the tile is not blank (zero)
			if(board[i][j] != 0)
				sum += calculateInversion(board, i, j);
		}
	}

	//If the total sum is even the board is solvable
	//Otherwise, it's unsolvable
	return (sum % 2 == 0);
}

// Read the initial board inside file with the given filename
// and store it inside the "board" array.
// Returns true if the read was successful or false if it wasn't
bool readBoard(const string& filename, board_array& board)
{
	ifstream file(filename);
	if(!file){
		return false;
	}

	for(auto& line : board) {
		for(auto& col : line) {
			file >> col;
		}
	}

	return true;
}

// Check if the board is well formed, i.e. have all the numbers between 0 and 8 and they all appear just once
bool checkBoard(const board_array& board)
{
	array<bool, BOARD_DIM*BOARD_DIM> contain;
	contain.fill(false);

	for(const auto& line : board)
	{
		for(const auto& col : line)
		{
			contain[col] = true;
		}
	}

	for(const auto& elem : contain)
	{
		if(elem == false)
			return false;
	}

	return true;
}

int main (int argc, char* argv[]){

	// Print every float number with a precision of 3 digits
	std::cout << std::setprecision(3) << std::fixed;

	// Requires the puzzle game
	if(argc != 6){
		cerr << "Wrong number of arguments. The arguments must be given in the following order:\n"
			"(1) path to the file containing the initial board\n"
			"(2) initial population size\n"
			"(3) initial chromosome size\n"
			"(4) crossover probability\n"
			"(5) mutation probability\n";
		return 1;
	}

	// Getting the parameters of the GA
	const size_t initial_pop_size = atol(argv[2]);
	const size_t initial_chromosome_size = atol(argv[3]);
	const float crossover_probability = atof(argv[4]);
	const float mutation_probability = atof(argv[5]);

	// Asserting that all probabilities are inside the required range
	assertm(0. <= crossover_probability && crossover_probability <= 1., "The crossover probability value must be between 0 and 1!");
	assertm(0. <= mutation_probability && mutation_probability <= 1., "The mutation probability value must be between 0 and 1!");

	// Reading the initial board configuration
	board_array board;
	const string filename = argv[1];
	if(!readBoard(filename, board))
	{
		cerr << "There was an error loading the puzzle file" << "\n";
		return 1;
	}

	cout << "Initial Board is: \n";
	for(const auto& line : board){
		for(const auto& col : line){
			cout << col << " ";
		}
		cout << "\n";
	}
	cout << "\n";

	if(!checkBoard(board))
	{
		cerr << "This initial board is invalid. Either one value appears twice or does not appear in the board\n";
		return 1;
	}

	if(!isSolvable(board))
	{
		cout << "This initial board cannot be solved.\n";
		return 0;
	}

	// File in which the results of each generation will be writen
	string out_filename = "results.txt";
	ofstream out_file(out_filename);
	if(!out_file){
		cerr << "There was a error creating the results file\n";
		return 1;
	}
	out_file.precision(3);

	out_file << "Generation | Avg fitness | Max fitness\n\n";

	// Generate a random population
	pop_vector population = generatePopulation(initial_pop_size, initial_chromosome_size);

	vector<int> mov;

	t_fitness best_fitness = 0;
	bool improved = false;
	unsigned n_generations_no_improve = 0;
	bool grow_chromosome = false;

	chromosome_vector the_best;

	unsigned generation = 1;

	while(!solved){
		cout << "________________________________________________________________________________________________________________________________________\n";
		cout << "|\t\t\t\t\t\t"
			<< "Running generation " << generation << "\t\t\t\t\t\t\t\t\t|"
			<< "\n";

		out_file << generation << " ";

		// Function to calculate the fitness of each candidate
		fitness_vector fitness = fitnessCalculationManhattan(population, board, &mov);
		cout << "|\t" << "Size of chromosome" << "\t|";
		cout << "\tAverage fitness" << "\t\t|";
		cout << "\tMax fitness" << "\t|";
		cout << "\tBest Board" << "\t|";
		cout << "\tBest so far" << "\t|";
		cout << "\n";

		double total = 0;
		double max = 0;
		size_t max_id = -1;

		// Retrieving the average and the maximum fitness of each generation
		for(size_t i = 0; i < fitness.size(); i++){
			total += fitness[i];
			if(fitness[i] > max){
				max = fitness[i];
				max_id = i;
			}
		}
		float average_fitness = total/fitness.size();

		cout << "|\t" << population[0].size() << "\t\t\t|";
		cout << "\t\t" << average_fitness << "\t\t|";
		cout << "\t" << max << "\t\t|";
		cout << "   ";

		out_file << average_fitness << " ";
		out_file << max << "\n";

		// The fitness have improved from the previous best one
		if(max > best_fitness)
		{
			best_fitness = max;
			the_best = population[max_id];
			improved = true;
			n_generations_no_improve = 0;
		}
		else
		{
			improved = false;
			n_generations_no_improve++;
		}

		// Retrieves the final board of the maximum fitness of this generation
		if(!solved){
			board_array best_board_of_generation = downOnTree(population[max_id], board, nullptr);
			for (int j = 0; j < BOARD_DIM; j++){
				for (int k = 0; k < BOARD_DIM; k++){
					cout << best_board_of_generation[j][k] << " ";
				}
			}

			cout << "  |   ";

			board_array best_board_so_far = downOnTree(the_best, board, nullptr);
			for (int j = 0; j < BOARD_DIM; j++){
				for (int k = 0; k < BOARD_DIM; k++){
					cout << best_board_so_far[j][k] << " ";
				}
			}

			cout << "  |";

			cout << "\n";

			// Selects the candidates to reproduce
			pop_vector parents = selectReproducers(population, fitness);

			if((!improved) && (n_generations_no_improve > 4))
			{
				grow_chromosome = true;
				n_generations_no_improve = 0;
			}

			// Candidates reproduction
			pop_vector children = reproducePopulation(parents, crossover_probability, grow_chromosome);

			// Keeps the best chromosome in the new population if it didn't grow
			if(!grow_chromosome)
				elitism(population, children, max_id);

			// The children becomes the new population
			population = children;

			grow_chromosome = false;

			// Mutate the population only if a solution hasn't been found
			population = mutatePopulation(population, mutation_probability);
		}
		else {
			cout << "1 2 3 4 5 6 7 8 0   |\n";
			cout << "________________________________________________________________________________________________________________________________________\n";
		}

		generation++;

	}
	printSolution(final_moves);
	printBoardMoves(final_moves,board);

	out_file << "Moves: " << final_moves.size();

	return 0;
}

