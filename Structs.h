/* 
 * File:   Structs.h
 * Author: James
 *
 * Created on January 21, 2014, 10:47 PM
 */
#ifndef STRUCTS_H
#define	STRUCTS_H
using namespace std;
#include <string>

    struct Training_Data{
        double data =0;
        double time_information =0; 
    };
    struct Genetic_Parameters{
        int popsize =0;                    //initial genetic population size
        int gen_cycle =0;                   //number of iterations of the genetics function
        int error_ok =0;                   //if our error is lower than this value, break.
    };
    
    struct NeuralNetwork_Parameters{     // always only 1 hidden layer, assume layers ==1 and for multiple layers we will be using multiple networks, anything else is dumb.
        int I_num =0;                     //number of total inputs
        int O_num =0;                      //number of total outputs
        int training_data_length =0;       //how many individual training data steps there are
        int hidden_neurons =0;                   //number of hidden neurons
        int contextual_neurons =0;         //number of contextual neurons (IE neurons that remember previous states)
        int wtqt =0;                       //number of weights, created by NeuralNetworkFab)
        int training_verify_percentage =0.5;    // percentage of total training data used for training, with the rest for verification, 50% default
    };
    
    struct Probability_Parameters{
        double prob_genesplit_avg =0;       //probability of the childs DNA being averaged from the parents
        double prob_genesplit_add =0;      //probability of the child's DNA being created from the parents genes added together
        double prob_genesplit_sub =0;      //probability of the parent 1 subtracting from parent 2 to create child
        double prob_genesplit_parallel =0; //probability of the parents being added and divided by the multiplication
        double prob_mutation =0;           //probability of mutation per reproduction
        double RNG_SD =0;                  //standard deviaviation of the random values
        double lamda =0;                   //rate of decay of overall population, lower = slower pop decrease
    };
    struct Weights{
        double sfitness =0;          //scaled fitness score for the particular individual (comes from RMS)
        double weight =0;
        bool eligible_rep =0;        //whether this individual is eligible for reproduction 
    };
    
    struct Parent{
        double gene =0;
        bool parentnum =0;         // false = parent 1, true = parent 2
        
    };
    
#endif	/* STRUCTS_H */

