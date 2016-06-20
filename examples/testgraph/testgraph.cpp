#include "stdafx.h"
#undef max
#undef min 
#include <vector>
#include <stack>
#include <map>
#include <set>
#include <iostream>
#include <omp.h>
#include <cstdlib>
#include <sys/time.h>
//
#include <sstream>
#include <string>
#include <fstream>

#include "utility.h"

using namespace std;


void print_sep(){
  printf("******************\n");
}

void dump_graph(PNGraph G){
    // dump the graph
  printf("Graph (%d, %d)\n", G->GetNodes(), G->GetEdges());
  for (PNGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
    printf("  %d: ", NI.GetId());
    for (int e = 0; e < NI.GetDeg(); e++) {
      printf(" %d", NI.GetNbrNId(e)); }
    printf("\n");
  }
}

double rtclock(void);

double rtclock(void)
{
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday (&Tp, &Tzp);
  if (stat != 0) printf("Error return from gettimeofday: %d",stat);
  return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

long bfs_visit(PNGraph G, int node_id){

	   //PNGraph::TObj::TNodeI NI
	  long sum_of_elems = 0;
      std::vector<int> queue;
      std::set<int> traversed;
    //int node_id = NI.GetId();
     queue.push_back(node_id);
    
    while (!queue.empty()) {
    	
      	int current_node_id = queue.back();
      	//printf("visiting node: %d",current_node_id);
      	queue.pop_back();
      //current_node_id = current_node_NI.GetId()
      	if (traversed.find(current_node_id) == traversed.end()) {
      		traversed.insert(current_node_id);
      		PNGraph::TObj::TNodeI current_node_NI = G->GetNI(current_node_id);
          	int out_degree = current_node_NI.GetOutDeg();
      		for(int e = 0; e < out_degree; e++){
     	 		int next_node_id = current_node_NI.GetOutNId(e);
     	 		if(traversed.find(next_node_id) == traversed.end()){
     	 			queue.push_back(next_node_id);
     	 		}
      		}
     	}
	}	
	//printf("max size of  set is: %ld \n",traversed.max_size());
	sum_of_elems = traversed.size() - 1;
	//printf("node: %d 's reach ability: %ld \n",node_id,sum_of_elems);
	return sum_of_elems;
}


/*
void parallel_wo_locality_bsf(PNGraph G,int n){
  //int n = 4;
  //n: number of cores
  int ppn = 12;
  int num_threads = n * ppn;
  printf("number of threads is: %d\n",num_threads);
  std::vector<long> reachability;//= new std::vector<long>();
  std::vector<int> nodes_id;

  printf("add all nodes' ids to vector...\n");
  for (PNGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
  //PGraph::TObj::TNodeI NI = G->BegNI();
    int node_id = NI.GetId();
    //if(node_id < 10000){
    nodes_id.push_back(node_id);
    //printf("node: %d is added to vec\n",node_id);
    //}
  }
  long number_of_nodes_in_graph = nodes_id.size();
  printf("%ld nodes in Graph \n",number_of_nodes_in_graph);

  #pragma omp parallel num_threads(num_threads)
  {
    int id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    //printf("thread_id: %d, and number of threads: %d\n",id,nthreads);
    long number_of_children = 0;
    //long block_size = number_of_nodes_in_graph / nthreads;
    //long extra_block_size = number_of_nodes_in_graph % nthreads;

    //long start_index = id*block_size;
    //long end_index = (id + 1) * block_size;
    //if(id == nthreads - 1){
    //  end_index = number_of_nodes_in_graph;
    //}

    for (int i = id; i < number_of_nodes_in_graph; ++nthreads){

      number_of_children = bfs_visit(G,nodes_id[i]);
      
      #pragma omp critical
      {
        reachability.push_back(number_of_children);
        long computed = reachability.size();
        if(computed % 20000 ==0){
          printf("%ld nodes has been computed, total edges: %lld...\n",computed,accumulate_vector(reachability));
        }
      }
    }

  }

  long number_of_nodes_computed = reachability.size();
  printf("number of nodes have been computed: %ld \n",number_of_nodes_computed);
  long long total_edges = 0;
  for(int i = 0; i < number_of_nodes_computed; ++i){
    //printf("edges: %ld\n",reachability[i]);
    total_edges += reachability[i];
  }
  printf("total edges in transitive closure: %lld , %lld \n",total_edges,accumulate_vector(reachability));
}
*/

void parallel_wo_locality_bsf(PNGraph G,int n){
  //int n = 4;
  //n: number of cores
  //int ppn = 8;
  int num_threads = n;
  printf("number of threads is: %d\n",num_threads);
  //std::vector<long> reachability;//= new std::vector<long>();
  long long reachability = 0;
  long number_of_nodes_have_been_computed = 0; 
  std::vector<int> nodes_id;

  printf("add all nodes' ids to vector...\n");
  for (PNGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
  
    if(NI.GetOutDeg() != 0){
          int node_id = NI.GetId();
          nodes_id.push_back(node_id);
    }
  }

  long number_of_nodes_in_graph = nodes_id.size();
  printf("%ld nodes added to vector...\n",number_of_nodes_in_graph);

  #pragma omp parallel num_threads(num_threads)
  {
    int id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    //printf("thread_id: %d, and number of threads: %d\n",id,nthreads);
    long number_of_children = 0;

    for (int i = id; i < number_of_nodes_in_graph; i += nthreads){

      number_of_children = bfs_visit(G,nodes_id[i]);
      
      #pragma omp critical
      {
        reachability += number_of_children;
        number_of_nodes_have_been_computed += 1;
        /*
        if(number_of_nodes_have_been_computed % 10000 == 0){
          printf("Thread: %d, %ld nodes has been computed, total edges: %lld...\n",id,number_of_nodes_have_been_computed,reachability);
        }
        */
      }
    }
  }

  printf("number of nodes have been computed: %ld \n",number_of_nodes_have_been_computed);
  printf("total edges in transitive closure: %lld \n",reachability);
  printf("*****************************************\n");
}

void parallel_bsf_wo_vector(PNGraph G,int n){
  //int n = 4;
  //n: number of cores
  //int ppn = 12;
  int num_threads = n;
  printf("number of threads is: %d\n",num_threads);
  long long reachability = 0;//= new std::vector<long>();
  long number_of_nodes_have_been_computed = 0;
  std::vector<int> nodes_id;

  printf("add all nodes' ids to vector...\n");
  for (PNGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
  //PGraph::TObj::TNodeI NI = G->BegNI();
    if(NI.GetOutDeg() != 0){
          int node_id = NI.GetId();
          nodes_id.push_back(node_id);
    }

  }

  long number_of_nodes_in_graph = nodes_id.size();
  printf("%ld nodes added to vector \n",number_of_nodes_in_graph);

  #pragma omp parallel num_threads(num_threads)
  {
    int id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    long number_of_children = 0;
    long block_size = number_of_nodes_in_graph / nthreads;

    long start_index = id*block_size;
    long end_index = (id + 1) * block_size;

    if(id == nthreads - 1){
      end_index = number_of_nodes_in_graph;
    }

    for (int i = start_index; i < end_index; ++i){
      number_of_children = bfs_visit(G,nodes_id[i]);

      #pragma omp critical
      {
        reachability += number_of_children;
        number_of_nodes_have_been_computed += 1;
        /*
        if(number_of_nodes_have_been_computed % 10000 == 0){
          printf("ThreadID: %d, %ld nodes has been computed, total edges: %lld...\n",id,number_of_nodes_have_been_computed,reachability);
        }
        */
      }
    }

  }

  printf("number of nodes have been computed: %ld \n",number_of_nodes_have_been_computed);
  printf("total edges in transitive closure: %lld \n",reachability);
  printf("*****************************************\n");
}
/*
void parallel_bsf(PNGraph G,int n){
	//int n = 4;
  //n: number of cores
	int ppn = 12;
	int num_threads = n * ppn;
  printf("number of threads is: %d\n",num_threads);
	std::vector<long> reachability;//= new std::vector<long>();
	std::vector<int> nodes_id;

	printf("add all nodes' ids to vector...\n");
	for (PNGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
	//PGraph::TObj::TNodeI NI = G->BegNI();
		int node_id = NI.GetId();
    //if(node_id < 10000){
		nodes_id.push_back(node_id);
    //printf("node: %d is added to vec\n",node_id);
    //}
	}
	long number_of_nodes_in_graph = nodes_id.size();
	printf("%ld nodes in Graph \n",number_of_nodes_in_graph);

	#pragma omp parallel num_threads(num_threads)
	{
		int id = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
		//printf("thread_id: %d, and number of threads: %d\n",id,nthreads);
		long number_of_children = 0;
		long block_size = number_of_nodes_in_graph / nthreads;
		//long extra_block_size = number_of_nodes_in_graph % nthreads;

		long start_index = id*block_size;
		long end_index = (id + 1) * block_size;
		if(id == nthreads - 1){
			end_index = number_of_nodes_in_graph;
		}

		for (int i = start_index; i < end_index; ++i){

			number_of_children = bfs_visit(G,nodes_id[i]);
			#pragma omp critical
			{
				reachability.push_back(number_of_children);
        long computed = reachability.size();
        if(computed % 10000 ==0){
          printf("%ld nodes has been computed, total edges: %lld...\n",computed,accumulate_vector(reachability));
        }
			}
		}

	}

	long number_of_nodes_computed = reachability.size();
	printf("number of nodes have been computed: %ld \n",number_of_nodes_computed);
	long long total_edges = 0;
	for(int i = 0; i < number_of_nodes_computed; ++i){
		//printf("edges: %ld\n",reachability[i]);
		total_edges += reachability[i];
	}
	printf("total edges in transitive closure: %lld , %lld \n",total_edges,accumulate_vector(reachability));
}
*/


void remove_edges(PNGraph G,std::string edges_list_file, int reverse){

    ifstream myin;
    myin.open(edges_list_file.c_str(), std::ifstream::in);
    int *a = new int();
    int *b = new int();
    int edges_to_be_removed = 0;

    if(myin.is_open()){
        std::string line;
        //stringstream ss;
        while(getline(myin,line)){
          split_string(line,a,b);
          edges_to_be_removed += 1;
          if (reverse == 1){
            G->DelEdge(*b, *a);
            //printf("edges to be removed: (%d,%d) \n",*b,*a);
          }
          else{
            G->DelEdge(*a, *b);
            //printf("edges to be removed: (%d,%d) \n",*a,*b);
          }
        }
        printf("number of edges to be removed: %d\n",edges_to_be_removed);
        myin.close();
    }
    else{
        printf("**NO EDGES TO BE DELETED**\n");
    }
}


void remove_edges_wo_pointer(PNGraph G,std::string edges_list_file, int reverse){

    ifstream myin;
    myin.open(edges_list_file.c_str(), std::ifstream::in);

    int edges_to_be_removed = 0;

    if(myin.is_open()){
        std::string line;
        //stringstream ss;
        while(getline(myin,line)){
          std::vector<int> edge = split_string(line);
          edges_to_be_removed += 1;
          if (reverse == 1){
            G->DelEdge(edge[1], edge[0]);
            //printf("edges to be removed: (%d,%d) \n",*b,*a);
          }
          else{
            G->DelEdge(edge[0], edge[1]);
            //printf("edges to be removed: (%d,%d) \n",*a,*b);
          }
        }
        printf("number of edges to be removed: %d\n",edges_to_be_removed);
        myin.close();
    }
    else{
        printf("**NO EDGES TO BE DELETED**\n");
    }
}


int main(int argc, char* argv[]) {
 

   if(argc != 6){
    readme();
    return 0;
   }

   int number_of_cores = atoi(argv[1]);
   printf("Number of cores: %d\n",number_of_cores);
   int locality = atoi(argv[3]);
   printf("Locality or not: %d\n",locality);
   int reverse = atoi(argv[4]);
   printf("Reverse the graph: %d\n",reverse);

   printf("Load graph file from: %s \n",argv[2]);
   double clkbegin_1 = rtclock();
   int source_index = 0;
   int destination_index = 1;

   if(reverse == 1){
   	source_index = 1;
   	destination_index = 0;
   }

   PNGraph G = TSnap::LoadEdgeList<PNGraph>(argv[2], source_index, destination_index);
   printf("Graph (%d, %d)\n", G->GetNodes(), G->GetEdges());

   printf("Edges to be deleted, load from file: %s\n",argv[5]);
   //remove_edges(G,argv[5],reverse);
   remove_edges_wo_pointer(G,argv[5],reverse);
   
   printf("Graph (%d, %d)\n", G->GetNodes(), G->GetEdges());
   
   double clkend_1 = rtclock();
   double t_1 = clkend_1 - clkbegin_1;
   printf("Prepare Graph Done: %0.3f sec\n",t_1);


  printf("Start to compute transitive closure\n");
  double clkbegin_2 = rtclock();
  
  if(locality == 0){
    parallel_wo_locality_bsf(G,number_of_cores);
  }
  else{
    parallel_bsf_wo_vector(G,number_of_cores);
  }
  
  double clkend_2 = rtclock();
  double t_2 = clkend_2 - clkbegin_2;
  printf("compute transitive closure used: %0.3f sec\n",t_2);
  /*
      All graph and network datatypes define node and edge iterators. 
      In general graph/network data types use the following functions
       to return various iterators:

      BegNI(): iterator to first node
      EndNI(): iterator to one past last node
      GetNI(u): iterator to node with id u
      BegEI(): iterator to first edge
      EndEI(): iterator to one past last edge
      GetEI(u,v): iterator to edge (u,v)
      GetEI(e): iterator to edge with id e (only for multigraphs)

      GetRndNId(): get a random node id
  */
  
  //PGraph::TObj::TNodeI NI = G->GetNI(0);

  /*
      In general node iterators provide the following functionality:

      GetId(): return node id
      GetOutDeg(): return out-degree of a node
      GetInDeg(): return in-degree of a node
      GetOutNId(e): return node id of the endpoint of e-th out-edge
      GetInNId(e): return node id of the endpoint of e-th in-edge
      GetNbrNId(e): return node e's neighbors
      IsOutNId(int NId): do we point to node id n
      IsInNId(n): does node id n point to us
      IsNbhNId(n): is node n our neighbor
  */

 // printf("Delete edge %d -- %d\n", NI.GetId(), NI.GetOutNId(0));
 // G->DelEdge(NI.GetId(), NI.GetOutNId(0));

 // int RndNId = G->GetRndNId();
  //printf("Trying to Delete node %d\n", RndNId);
  //RndNId = G->GetRndNId();
  //printf("Delete node %d\n", RndNId);

  //G->DelNode(RndNId);
  //G->Dump();
  //IAssert(G->IsOk());


  /**
  TIntV NIdV;
  for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
    if (NIdV.Len() < G->GetNodes()/2) { NIdV.Add(NI.GetId()); }
  }
  PGraph SubG = TSnap::GetSubGraph(G, NIdV);
  //SubG->Dump();
  // get UNGraph
  { PUNGraph UNG = TSnap::ConvertGraph<PUNGraph>(SubG);
  UNG->Dump();
  IAssert(UNG->IsOk());
  TSnap::ConvertSubGraph<PNGraph>(G, NIdV)->Dump(); }
  // get NGraph
  { PNGraph NG = TSnap::ConvertGraph<PNGraph>(SubG);
  NG->Dump();
  IAssert(NG->IsOk());
  TSnap::ConvertSubGraph<PNGraph>(G, NIdV)->Dump(); }
  // get NEGraph
  { PNEGraph NEG = TSnap::ConvertGraph<PNEGraph>(SubG);
  NEG->Dump();
  IAssert(NEG->IsOk());
  TSnap::ConvertSubGraph<PNGraph>(G, NIdV)->Dump(); }

  TSnap::TestAnf<PUNGraph>();
  **/

  //test_omp();
  return 0;
}
