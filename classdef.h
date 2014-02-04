/* 
 * File:   Neural.h
 * Author: James
 *
 * Created on January 21, 2014, 12:50 PM
 */

#ifndef CLASSDEF_H
#define	CLASSDEF_H
#include "structs.h"
#include <string>
static ini ini1;
    class Neural{
    public:
        typedef boost::multi_array<Input, 2> in_array;
        typedef boost::multi_array<Input, 1> out_array;
        in_array inputs;         //input data storage
        out_array real_output;
        Neural(int , int , int ,int , int , int ,int);
        Neural();
        virtual ~Neural();
        void callgen();         //calls the genetics class to iterate until weights have been set
        void increment();       //will increment the values of the inputs as we're now done with the current data
        void prepare();           //gets the inputs in the right format

    private:
    };

    class genetic{
    public:
        typedef boost::multi_array<individual, 2> current_array;
        current_array current;
        individual best_pop;
        genetic();
        bool Parent_Chooser(int order, int population_size, double lamda);
        
        /* main genetic algo loop. creates initial pop, creates chroma from randomized weights
                      * then throws the individual into the RNNtion and the RMS thats returned
                      * is turned into a fitness, which is then scaled.
                      * The best individuals probably breed and their chroma is switched around based on probability
                      * and they could be mutated (low percentage of it though)
                      * keeps iterating till either gen_cycle is reached, or RMS of current best individual is sufficiently low*/
        double RNG_normal(double standard_deviation);
        bool RANDOM(double probability);
        virtual ~genetic();
    private:
        
        
        
    };
    class RNN : public Neural, public genetic
    {
    public:
        double fitness;
        RNN(int person_num, int layer_num, int num_hidden, int innum, int data_size, int wtqty);
        virtual ~RNN();
        void createfitness();
        double outputfit();
    private:
        double Activation_sig(double X);
        
    };

#endif	/* CLASSDEF_H */

