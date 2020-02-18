typedef struct {
    double fitness;
    unsigned char * genotype;
} Individual;

#define n_drugs 10
#define n_treatments 10
#define population_size 500

Individual * createPopulation();
int fitnessCompare (const void * a, const void * b);
Individual offSpring(Individual parentA, Individual parentB);
unsigned char * noDosageGenotype();
void printGenotype(unsigned char * genotype);