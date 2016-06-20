#include <iostream>
#include <sstream>
#include <string>
#include <vector>

void readme();
void split_string(std::string , int* , int* );

void readme(){
      printf("##################################\n");
      printf("*******     This code is used to compute Transitive Closure with OPENMP \n");
      printf("*******     Author: Jiankai Sun\n");
      printf("*******     Date: 2016.06.15\n");
      printf("*******     Compile: make\n");
      printf("*******     Submit job on OSC: qsub -I -l walltime==00:59:00 -l nodes=1:ppn=12 -l mem=32G \n");
      printf("##################################\n");
      printf("*******     Parameters are like: \n");
      printf("*******     ./testgraph number_of_cores graph.txt 0/1 0/1 edges_to_be_removed_file \n");
      printf("*******     number_of_cores: number of threads to be used \n");
      printf("*******     graph.txt: graph's edges list file\n");
      printf("*******     First 0/1: divide blocks with(1)/without(0) locality \n");
      printf("*******     Second 0/1: reverse the graph(1) or not(0)\n");
      printf("*******     edges_to_be_removed_file: edges_to_be_removed file, Otherwise use: None\n");
      printf("##################################\n");
}

long long accumulate_vector(std::vector<long> a_vector){
  long a_size = a_vector.size();
  long long sum = 0;
  for (long i = 0; i < a_size; ++i){
    sum += a_vector[i];
  }
  return sum;
}


std::vector<int> split_string(std::string s)
{   
    std::vector<int> edge;
    std::istringstream iss(s);
    std::string sub;
    iss >> sub;
    edge.push_back(atoi(sub.c_str()));
    iss >> sub;
    edge.push_back(atoi(sub.c_str()));
    return edge;

}

void split_string(std::string s,int* a, int* b)
{
    std::istringstream iss(s);
    std::string sub;
    iss >> sub;
    *a = atoi(sub.c_str());
    iss >> sub;
    *b = atoi(sub.c_str());
}