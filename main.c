#include <stdio.h>
#include <time.h>
#include "RKF78.h"
#include "fitness.h"
#include "GACore.h"

int main() {
    srand((unsigned) time(NULL));
    printf("Start Running!\n");
    Individual * population = createPopulation();
    Individual individual = *(population);

    printf("Checking first individual after population creation \n");
    printf("Fitness: %f \n", individual.fitness);
    printf("Genotype: ");
    for (int g = 0; g < n_treatments * n_drugs; g++ ){
        printf("%d|", individual.genotype[g]);
        //printf("%s|", randomBoolean() ? "true" : "false");
    }
    printf("\n");

    qsort(population, population_size, sizeof(Individual), fitnessCompare);

    for(int i = 0; i < population_size; i++){
        printf("%f \n", population[i].fitness);
    }

    Individual child = offSpring(population[0], population[1]);

    for(int i =0; i < n_drugs * n_treatments; i++){
        printf("%02d|", population[0].genotype[i]);
    }
    printf("\n");
    for(int i =0; i < n_drugs * n_treatments; i++){
        printf("%02d|", population[1].genotype[i]);
    }
    printf("\n");
    for(int i =0; i < n_drugs * n_treatments; i++){
        printf("%02d|", child.genotype[i]);
    }

    printf("\n");
    printf("%f \n", population[0].fitness);
    printf("%f \n", population[1].fitness);
    printf("%f \n", child.fitness);

    // Create a Non-treatment individual.
    Individual nonTreatment;
    nonTreatment.genotype = noDosageGenotype();
    nonTreatment.fitness = Curative_Fitness(nonTreatment.genotype);
    printf("Max Grown Cells: %f \n",nonTreatment.fitness);

    /**
     * Algorithm from here.
     * */
    int generationCount = 0;

    do {
        printf("Generation count: %d ", generationCount);
        printf("Best: %f ", population[0].fitness);
        printf("Second Best: %f ", population[1].fitness);
        printf("Worst: %f \n", population[population_size-1].fitness);
        Individual a = population[0];
        Individual b = population[1];
        Individual c = population[3];
        Individual d = population[4];
        Individual worstFitnessA = population[population_size - 1];
        Individual worstFitnessB = population[population_size - 2];

        free(worstFitnessA.genotype);
        free(worstFitnessB.genotype);

        Individual childA = offSpring(a, b);
        Individual childB = offSpring(c, d);

        population[population_size -1] = childA;
        population[population_size -2] = childB;

        qsort(population, population_size, sizeof(Individual), fitnessCompare);

        generationCount++;
    } while (generationCount < 5000 && population[0].fitness > 1000.0);

    printf("Best Individual fount fitness: %f \n", population[0].fitness);
    printGenotype(population[0].genotype);

    return 0;
}

