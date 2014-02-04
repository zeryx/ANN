/* 
 * File:   Structs.h
 * Author: James
 *
 * Created on January 21, 2014, 10:47 PM
 */
#ifndef STRUCTS_H
#define	STRUCTS_H
using namespace std;
#include <iostream>
#include <string>

    struct Input{
        double data;
    };
    struct ini{
        int innum;              //number of total inputs
        int hid_neur;           //number of hidden neurons
        int layers;             //number of hidden neuron layers
        int wtqt;
        int popsize;            //initial genetic population size
        int gen_cycle;           //number of iterations of the genetics function
        int error_ok;             //if our error is lower than this value, break.
        int data_size;
    };
    struct probability{
        double prob_genesplit_avg;  //probability of the childs DNA being averaged from the parents
        double prob_genesplit_add;  //probability of the child's DNA being created from the parents genes added together
        double prob_genesplit_sub; //probability of the parent 1 subtracting from parent 2 to create child
        double prob_genesplit_parallel; //probability of the parents being added and divided by the multiplication
        double prob_mutation;     //probability of mutation per reproduction
        double RNG_SD;          //standard deviaviation of the random values
        double lamda;           //rate of decay of overall population, lower = slower pop decrease
    }prob;
    struct individual{
        double sfitness;          //scaled fitness score for the particular individual (comes from RMS)
        double weight;
        string genename;               //name of the gene (what it effects)
        bool eligible_rep;        //whether this individual is eligible for reproduction 
        bool eligible_sort;       //whether this individual is eligible for sorting
        int order;              //order 
    };
    
    struct parent{
        double gene;
        bool parentnum;         // false = parent 1, true = parent 2
        
    };
    
#endif	/* STRUCTS_H */

