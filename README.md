# GAPuzzleGame
An implementation of an genetic algorithm to solve the 8-Puzzle Game.

## Description
The 8-Puzzle Game consists of a 3x3 board in which are placed, in random order, all the numbers between 1 and 8, represented by tiles. In the end there will be a blank tile which is used to slide the tiles that are next to it. The objective of the game is to slide the tiles in consecutive moves until until the board reaches a target configuration, usually defined by the following configuration:

1 2 3<br>
4 5 6<br>
7 8 0

where the zero represents the blank tile.

## Puzzles
Inside the "puzzles" folder there are 20 differents board configurations which can be used as an initial board. There are, however, some configurations that can't be solved by any number of moves, which is the case of the board inside the "puzzle3.txt" file; in these cases, our program will print that this board cannot be solved and terminate itself.<br>
Aditionally, inside this folder there's a "goalstate.txt" file which defines the target configuration of the board.<br>
Both the initial and target boards can be customized. The "goalstate.txt" can be modified and other initial boards files can be created, as long as both have all the numbers between 0 and 8 in any order and without repetition.

## Population and chromosomes
The population of this genetic algorithm is consisted by N chromosomes - N is called "population size" -, each chromosome is consisted by n bits - n is called "chromosome size".<br>
Each chromosome represents a sequence of moves of the blank tile starting from the initial board. Each move is represented by a pair of non-overlapping bits of the chromosome. For example, if the chromosome starts with the pair of bits 10, it moves the blank tile to the left; if, instead, it's the pair 01, it moves the blank tile down, and so on.

## Crossover
The crossover is done by selecting one random point of the chromosome of the parents such that the first children takes the left size of the first parent and the right size of the second parent, and the second children takes the right size of the first parent and the left size of the second parent.<br>
This is done according to a probability rate defined by the user.

## Mutation
The mutation is done by flipping a random bit of the chromosome and also follows a probability rate to decide if it will happen or not for each chromosome.

## Compiling
Inside the project's root, execute:
```
make all
```

This will create an executable called "8puzzle.exe" in the same folder.

## Running
The executable takes 5 arguments from the command line: the path to the initial board (puzzle), the size of the population of chromosomes, the initial size of the chromosome, the probability that a crossover will happen at each reproduction (a pair of parents generating a pair of childrens) and the probability that a mutation will happen for each chromosome.

To run, inside the project's root execute:
```
./8puzzle.exe <path to the puzzle> <population size> <initial chromosome size> <crossover probability> <mutation probability>
```
