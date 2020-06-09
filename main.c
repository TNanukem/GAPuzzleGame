#include <stdio.h>
#include <stdlib.h>
#include<time.h> 

int* generatePopulation(){

}

int* fitnessCalculation(int *population){

}

int* select_reproducers(int *population, int *fitness_vector, int population_size, int new_pop_size){
    int i, candidate_1, candidate_2;
    int *new_population;
    int fitter, non_fitter;
    float r, k = 0.75;

    new_population = (int*) malloc(sizeof(int) * new_pop_size);
    
    for(i = 0; i < new_pop_size; i++){
        candidate_1 = rand() % population_size;
        candidate_2 = rand() % population_size;
        r = (float) rand() / population_size;

        // They need to be different
        while(candidate_2 == candidate_1){
            candidate_2 = rand() % population_size;
        }

        // Verifying which is fitter
        if(fitness_vector[candidate_1] > fitness_vector[candidate_2]){
            fitter = candidate_1;
            non_fitter = candidate_2;
        }
        else{
            fitter = candidate_2;
            non_fitter = candidate_1;
        }

        if(r < k){
            // Select fitter
            new_population[i] = fitter;
        }
        else {
            // Select non-fitter
            new_population[i] = non_fitter;
        }
    }
    return new_population;
}

int* reproducePopulation(int *population, int *fitness_vector, int *selected){

}

int isSolved(int *population){

    int goalState[9] = {1, 2, 3, 4, 5, 6, 7, 8, 0};
    int i, solved = 1;

    // This is not right, it depends on how the population is set
    for(i = 0; i < 9; i++){
        if(goalState[i] != population[i]){
            solved = 0;
            break;
        }
    }
    printf("%d\n", solved);
    return solved;
}

int main (int argc, char* argv[]){

    int board[9], i;
    int solved = 0;
    int *population, *fitness_vector, *selected;
    time_t t;

    int new_pop_size, population_size = 10; // Isso aqui o certo é mudar dentro do loop, mas ainda n sei como

    srand((unsigned) time(&t));

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
        
        // Também não tá certo, precisa ver o jeito certo
        new_pop_size = population_size - 1;

        // Function to calculate the fitness of each candidate
        fitness_vector = fitnessCalculation(population);

        // Selects the candidates to reproduce
        selected = select_reproducers(population, fitness_vector, population_size, new_pop_size);

        // Candidates reproduction
        population = reproducePopulation(population, fitness_vector, selected);
            // crossOver()
            // mutation()

        population_size--;

        // Verifies if some chromosome has achieved the solution
        solved = isSolved(population);
        solved = 1;
    }

    return 0;
}