#include <stdio.h>
#include <stdlib.h>

int* generatePopulation(){

}

int* fitnessCalculation(int *population){

}

int* reproducePopulation(int *population, int *fitness_vector){

}

int isSolved(int *population){

}

int main (int argc, char* argv[]){

    int board[9], i;
    int solved = 0;
    int *population, *fitness_vector;

    // Requires the puzzle game
    if(argc < 2){
        printf("Please, insert the path to the puzzle");
        return 1;
    }

    // Reads the data and puts it on the board array
    FILE *p = fopen(argv[1], "r");
    if(p == NULL){
        printf("There was an error loading the puzzle file");
        return 1;
    }

    for(i = 0; i < 9; i++){
        fscanf (p, "%d", &board[i]);
        printf("%d ", board[i]);
    }
    
    // Generate a random population
    population = generatePopulation();

    while(solved == 0){

        // Function to calculate the fitness of each candidate
        fitness_vector = fitnessCalculation(population);

        // Candidates reproduction
        population = reproducePopulation(population, fitness_vector);
            // crossOver()
            // mutation()

        // Verifies if some chromosome has achieved the solution
        solved = isSolved(population);
        solved = 1;
    }

    return 0;
}