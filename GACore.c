#include "GACore.h"
#include <stdio.h>
#include "RKF78.h"
#include "fitness.h"

unsigned char random_dosage() {
    unsigned char lower = 0;
    unsigned char upper = 15;
    unsigned char value =((rand() % (upper - lower + 1)) + lower);
    unsigned char result = value;
    //printf("%d|", value);
    return result;
}

unsigned char random_bit() {
    unsigned char lower = 0;
    unsigned char upper = 3;
    unsigned char value =((rand() % (upper - lower + 1)) + lower);
    unsigned char result = value;
    //printf("%d|", value);
    return result;
}

unsigned char * createGenotype(){
    unsigned char * g = malloc(sizeof(unsigned char) * n_treatments * n_drugs);
    for(int i = 0; i < n_treatments * n_drugs; i++){
        *(g + i) = random_dosage();
    }
    return g;
}

unsigned char * noDosageGenotype(){
    unsigned char * g = malloc(sizeof(unsigned char) * n_treatments * n_drugs);
    for(int i = 0; i < n_treatments * n_drugs; i++){
        *(g + i) = 0;
    }
    return g;
}


Individual createIndividual(){
    Individual individual;
    individual.genotype = createGenotype();
    individual.fitness = Curative_Fitness(individual.genotype);
    return  individual;
}

int fitnessCompare (const void * a, const void * b) {
    Individual *individualA = (Individual *)a;
    Individual *individualB = (Individual *)b;
    return individualB->fitness > individualA->fitness ? 0 : 1;
}

Individual * createPopulation(){
    Individual * population = malloc(sizeof(Individual) * population_size);
    for(int i = 0; i < population_size; i++){
        *(population + i ) = createIndividual();
    }
    return population;
}

/**
 * Uniform Crossover.
 * A 0 0 1 1 0 0 1 1
 * B 0 1 0 1 0 1 0 1
 * C 0 0 0 0 1 1 1 1
 * R 0 0 1 1 0 1 0 1
 *
 * a && (( ! b && ! c ) || ( b && ! c )) || (b && c)
 *(a && not b && not c) || (not a && b && c) || (a && b)
 *
 * */

unsigned char crossOverGene(unsigned char a, unsigned char b, unsigned char c){
    return (a & ~b & ~c) | (~a & b & c) | (a & b);
}

unsigned char * crossOver(unsigned char * genoTypeA, unsigned char * genoTypeB){
    // Uniform crossover
    unsigned char * childGenes = malloc(sizeof(unsigned char) * n_treatments * n_drugs);
    for (int i = 0; i < n_treatments * n_drugs; i++){
        unsigned char randomMask = random_dosage();
        unsigned char genA = genoTypeA[i];
        unsigned char genB = genoTypeB[i];
        *(childGenes + i) = crossOverGene(genA, genB, randomMask);
    }
    return childGenes;
}

void mutation(unsigned char * genoType){
    for (int i = 0; i < n_treatments * n_drugs; i++){
        *(genoType + i)  ^= 1UL << random_bit();
    }
}

void printGenotype(unsigned char * genotype){
    for(int i =0; i < n_drugs * n_treatments; i++){
        printf("%02d|", genotype[i]);
    }
    printf("\n");
}

Individual offSpring(Individual parentA, Individual parentB) {
    Individual child;
    child.genotype = crossOver(parentA.genotype, parentB.genotype);
    mutation(child.genotype);
    child.fitness = Curative_Fitness(child.genotype);
    return child;
}