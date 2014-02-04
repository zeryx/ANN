/* 
 * File:   MEMBERS.h
 * Author: James
 *
 * Created on January 21, 2014, 4:53 PM
 */

#ifndef MEMBERS_H
#define	MEMBERS_H
#include "classdef.h"
#include "structs.h"
#include <string>
#include <random>
    


   
    
    Neural::Neural(int input_quantity, int hidden_neurons, int layers,
            int initial_population_size, int probability_of_mutation, int number_of_generations, int error_ok)
    {
        ini1.layers = layers;
        ini1.gen_cycle = number_of_generations;
        ini1.hid_neur = hidden_neurons;
        ini1.popsize = initial_population_size;
        ini1.innum=input_quantity;
        ini1.error_ok=error_ok;
        prob.RNG_SD=15;
        prob.lamda=20;
        prob.prob_genesplit_add=0.5;
        prob.prob_genesplit_avg=0.7;
        prob.prob_genesplit_parallel=0.5;
        prob.prob_genesplit_sub=0.5;
        if(ini1.layers==1)
            {
                if(ini1.hid_neur>0)
                {
                    ini1.wtqt=(ini1.innum+1)*ini1.hid_neur+(ini1.hid_neur+1);
                }
                else
                {
                    ini1.wtqt=ini1.innum+1;
                }
            }
            if(ini1.layers==2)
            {
                ini1.wtqt=ini1.innum*ini1.hid_neur+2+ini1.layers*ini1.hid_neur;
            }
    }
    

void Neural::prepare()
    {
    double values[] ={0.0,0.0,0.0,1.0,1.0,0.0,1.0,1.0};
        double outputs[] = {0.0,1.0,1.0,0.0};
        inputs.resize(boost::extents[2][4]);
    real_output.resize(boost::extents[4]);
    int j=0;
    for(int i=0;i<4;i++)
    {
        for(int k=0;k<2;k++)
        {
            inputs[k][i].data=values[j];
            j++;
        }
        real_output[i].data = outputs[i];
    }
    
    }

// now genetics stuff
void Neural::callgen(){
    genetic prime;
}

genetic::genetic()          // meat and potatoes function
{
    int iter;
    for( iter=0;iter<ini1.gen_cycle;iter++)
    {
        current.resize(boost::extents[ini1.popsize][ini1.wtqt]);
        if(iter==0)                // if this is the first iteration, populate with randoms
        {
            for(int j=0;j<ini1.popsize;j++)        // create individuals
            {
                for(int k=0;k<ini1.wtqt;k++)             // randomize the "genes" (weights in this case) for the entire chroma
               {
                    current[j][k].weight=RNG_normal( prob.RNG_SD );
                }
            }
        }
        else{    // ok, since we already did one iteration lets start THE CYCLE


                /* plug our genes(weights) into the NN and get our error results as a fitness (1/RMS)*/
                /*fitness normalization begins here TEST THIS*/
                double sfit_sum;
                int numerize=0;
                double best_fit;
                for(int j=0;j<ini1.popsize;j++)  //for every current 
                {
                    RNN calculate(j, ini1.layers, ini1.hid_neur, ini1.innum, ini1.data_size, ini1.wtqt);
                    current[j][0].sfitness=calculate.outputfit();        // invert the RMS values of the output to create fitness
                    if(current[j][0].sfitness<= ini1.error_ok){
                        best_pop = current[j][0];
                        break;
                    }
                    sfit_sum += current[j][0].sfitness;              //compile a sum of fitnesses for normalization
                }

                sfit_sum /= ini1.popsize;                               // get the average fitness

                for(int j=0;j<ini1.popsize;j++)
                { 
                    current[j][0].sfitness = current[j][0].sfitness/sfit_sum; //normalize the fitness results of each individual
                    std::cout<<"normalized fitness: "<<current[j][0].sfitness<<endl;
                    current[j][0].eligible_sort = true;                      //prepare all individuals for sorting
                }
                /*normalized fitness section ends here*/


                /*sorting algorithm starts here TESTED, WORKS AS INTENDED*/
                    for( int n=0;n<ini1.popsize;n++)
                    {
                         best_fit=0;
                        for( int x=0;x<ini1.popsize;x++)
                        {
                            if(current[x][0].eligible_sort==true)
                                if(current[x][0].sfitness>best_fit)
                                    best_fit=current[x][0].sfitness;
                        }
                        for ( int x=0;x<ini1.popsize;x++)        // go back through the list and find the best
                        {
                            if(current[x][0].sfitness==best_fit)
                            {
                                current[x][0].order = numerize;
                                current[x][0].eligible_sort = false; //set illegible for future sorting
                                numerize++;
                            }
                            if(ini1.popsize<=numerize){
                                std::cout<<"error: two or more genes are identical."<<endl;
                                        break;
                            }
                            
                        }
                    }
                /*sorting algorithim stops here*/
                
                
                
                /*parent selection starts here*/
                typedef boost::multi_array<parent, 2> parent_array;
                typedef parent_array::index parent_index;
                parent_array parent;
                int breed_qt=0;
                for ( int x=0;x<ini1.popsize;x++)/* lets cycle through the individuals, and give the best ones*/                          
                {
                        if(Parent_Chooser(current[x][0].order, ini1.popsize, prob.lamda)){ // *FIGURE THIS OUT*
                            current[x][0].eligible_rep=true;             //if it passes, let it breed
                        breed_qt++;
                        }
                }
                /*parent selection stops here*/
                
                
                /*reproduction starts here*/
                parent.resize(boost::extents[breed_qt][ini1.wtqt]);
                int parent_indexer=0;
                int current_indexer, j;//parent "couple" selection
                for( current_indexer=0;current_indexer<ini1.popsize;current_indexer++)
                {
                    if(current[current_indexer][0].eligible_rep)                //CREATE RANDOM
                    {

                        for(j=0;j<ini1.wtqt;j++)
                        {
                            parent[parent_indexer][j].gene=current[current_indexer][j].weight;
                        }
                        parent_indexer++;
                                                       //we've made another parent, increment
                    }
                }
                /*parent now coupled, time to make babies*/            // empty the current array of people
                ini1.popsize=breed_qt*breed_qt;                       //reset popsize, since we're about to make new people
                current.resize(boost::extents[ini1.popsize][ini1.wtqt]);
                volatile double rn=rand()/((double)RAND_MAX+1);
                int selection;
                int parent1_index,parent2_index;
                srand(time(NULL));
                for(int curr_in=0; curr_in<breed_qt*breed_qt; curr_in++){    
                                              /* every parent can reproduce with each other to make 1 child, so that the total
                                              * number of combinations are breed_qt^2, and the weights are averaged
                                             * for each child, just to make things super easy */
                    selection =rn*breed_qt;
                    parent1_index=selection;            //randomly select the parents for this child
                    selection = rn*breed_qt;
                    parent2_index=selection;
                    if(RANDOM(prob.prob_genesplit_avg)) // if succeed, create children based on avg of parents genes
                    {
                        for(int gena=0;gena<ini1.wtqt;gena++){

                         current[curr_in][gena].weight=(parent[parent1_index][gena].gene+parent[parent2_index][gena].gene)/2;
                         gena++;
                        }

                     }
                    
                    if(RANDOM(prob.prob_genesplit_sub)) //if succeed, create chlidren based on subtraction of parent1 from parent2
                    {
                        for(int gena=0;gena<ini1.wtqt;gena++){

                         current[curr_in][gena].weight=(parent[parent1_index][gena].gene-parent[parent2_index][gena].gene);
                        }
                    }
                    
                    else if(RANDOM(prob.prob_genesplit_parallel)) //if succeed, create children made from parallel equivalent of parents
                    {
                        for(int gena=0;gena<ini1.wtqt;gena++){

                         current[curr_in][gena].weight=(parent[parent1_index][gena].gene+parent[parent2_index][gena].gene)/(parent[parent1_index][gena].gene*parent[parent2_index][gena].gene);
                        }
                    }
                    
                    else                //if all else fails, add the parents together to get the genes right
                    {
                        for(int gena=0;gena<ini1.wtqt;gena++){

                         current[curr_in][gena].weight=parent[parent1_index][gena].gene+parent[parent2_index][gena].gene;
                        }
                        
                     }
                    
                    
                 }
        }
    }
}
                    /*whew! thats the end of one full iteration, now its time to
                     test our new guys and see if their better than their parent!*/



 genetic::~genetic()
 {
     
 }
 Neural::Neural(){                      // placeholder constructor
 }
Neural::~Neural()
{
    cout<<"deleted everything"<<endl;
}




bool genetic::Parent_Chooser(int order, int population_size, double lamda){          
    double map = (double)(order)/(double)population_size;
    double y = exp(-lamda*map);
    double p = RAND_MAX*y - rand();
    if(p>0)
        return(true);
    else
        return(false);
}

bool genetic::RANDOM(double probability){                                                                      
        srand(time(NULL));
    double FITTED_MAX = probability*(double)(RAND_MAX);
        double p_scaled = FITTED_MAX - rand();
        if (p_scaled >=1)
            return(true);
        else
            return(false);
}

double genetic::RNG_normal( double standard_deviation){
    std::default_random_engine generator;
    std::normal_distribution<double> d1(0.0,standard_deviation);
    return(d1(generator));
}

 RNN::RNN(int person_num, int num_layers, int num_hidden, int innum, int data_size, int wtqty){        //assumes we only have 1 output
        double MS=0, RMS;
        double output[data_size];
        double hidden_val[num_hidden];
        for(int x=0;x<data_size;x++)
        {
            int wt_ctr=0;
            for(int i=0;i<num_hidden;i++)
            {
                for(int j=0;j<innum;j++)
                {
                    hidden_val[i]+= inputs[j][x].data*current[person_num][wt_ctr++].weight;
                }
             hidden_val[i]+= 1*current[person_num][wt_ctr++].weight;
             hidden_val[i] = Activation_sig(hidden_val[i]);
             assert(wt_ctr!=wtqty);
            }
            for(int i=0;i<num_hidden;i++)
            {
                output[x] += hidden_val[i]*current[person_num][wt_ctr++].weight;
                assert(wt_ctr!=wtqty);
            }
            output[x] += 1*current[person_num][wt_ctr++].weight;
            cout<<"output value: "<<output[x]<<endl;
            output[x] = Activation_sig(output[x]);
        }
    for(int x=0;x<data_size;x++)
    {
        MS += pow((output[x]-real_output[x].data), 2);
    }
    RMS = sqrt(MS/(double)data_size);
    fitness = exp(-20*RMS);
 }
 
double RNN::Activation_sig(double X)
{
    return(1/(1+exp(-X)));
}



RNN::~RNN(){            // simple deconstructor
}



double RNN::outputfit(){
    return fitness;
}
#endif	/* MEMBERS_H */

