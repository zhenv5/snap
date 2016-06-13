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

//using namespace std;

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




long bfs_visit(PNGraph G, int node_id){

	//PNGraph::TObj::TNodeI NI
	long sum_of_elems = 0;
    std::vector<int> queue;
    std::set<int> traversed;
    //int node_id = NI.GetId();
    queue.push_back(node_id);
    
    while (queue.size() > 0) {
    	
      	int current_node_id = queue.back();
      	//printf("visiting node: %d",current_node_id);
      	queue.pop_back();
      //current_node_id = current_node_NI.GetId()
      	if (traversed.find(current_node_id) == traversed.end()) {
      		traversed.insert(current_node_id);
      		PNGraph::TObj::TNodeI current_node_NI = G->GetNI(current_node_id);
      		for(int e = 0; e < current_node_NI.GetOutDeg(); e++){
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
    if (node_id > 50000){
      break;
    }
		nodes_id.push_back(node_id);
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
	printf("total edges in transitive closure: %lld \n",total_edges);
}

int main(int argc, char* argv[]) {
 
  printf("Creating graph:\n");
  
  /**
  PGraph TSnap::LoadEdgeList ( const TStr &  InFNm,
                              const int &   SrcColId,
                            const int &   DstColId 
  ) 
  Whitespace separated file of several columns:
  <source node="" id>=""> ... <destination node="" id>=""> ... 
  SrcColId and DstColId are column indexes of source/destination (integer!) node ids.
   This means there is one edge per line and node IDs are assumed to be integers.
  **/

   if(argc != 3){
      printf("need more parameters...\n");
      return 0;

   }
   int number_of_cores = atoi(argv[1]);
   printf("number of cores: %d \n",number_of_cores);
   printf("load file from: %s \n",argv[2]);

   PNGraph G = TSnap::LoadEdgeList<PNGraph>(argv[2], 0, 1);

  //printf("Dump Graph:");
  //G->Dump();
  //printf("---");
  //dump_graph(G);

  printf("start to compute transitive closure\n");
  parallel_bsf(G,number_of_cores);

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
