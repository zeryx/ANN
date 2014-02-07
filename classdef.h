/* 
 * File:   Neural.h
 * Author: James
 *
 * Created on January 21, 2014, 12:50 PM
 */

#ifndef CLASSDEF_H
#define	CLASSDEF_H
#include "structs.h"
class Dash_Board{
private:
    static void TrainingDataInitializer(double data, int IO_number, bool IO_indicator, double time);      // adds IO data to the static arrays
    static void TrainingRegimenDesign();
    static void CreateStructures(NeuralNetwork_Parameters tNN, Probability_Parameters tProb, Genetic_Parameters tGene);
public:
    typedef boost::multi_array<Training_Data, 2> IO_array;
    static IO_array Itraining;
    static IO_array Otraining;
    static IO_array Iverify;
    static IO_array Overify;
    static NeuralNetwork_Parameters NNP;
    static Probability_Parameters PP;
    static Genetic_Parameters GP;
    virtual ~Dash_Board()=0;
    
};
    class Neural_Network: public Dash_Board{
    private:
        static void CreateWeights();
        static void WeightOrganization();
        static double Neural_Function(IO_array Input, IO_array Outputs, int training_itr, int individual_itr, weight_array weights);
    public:
        typedef boost::multi_array<Weights, 2> weight_array;
        static weight_array current;
        virtual ~Dash_Board()=0;
    };

    class Genetics : public Neural_Network{
    private:
        typedef boost::multi_array<Weights, 1> weight_array;
        static void MainLoop();
        static void FitnessNormalization();
        static void Sort();
        static void ParentSelection();
        static void Reproduction();
        static void Mutation();
    public:
        static weight_array best_pop;
        virtual ~Genetics()=0;
    };
    class Random
    {
    public:
       static double RNGNormalDistribution(double standard_deviation);
       static bool BooleanRandom(double probability);
    };

#endif	/* CLASSDEF_H */

