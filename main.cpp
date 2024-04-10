/************************************************************************/
/* $Id: Goal-Obfuscation 2019-17-07 created by Santiago Franco$                */
/* This is an implementation of "An Optimization Approach to Robust Goal Obfuscation",KR 2020 by:*/
/* Prof. Sara Bernardini, Prof. Fabio Fagnani, Dr. Santiago Franco */
/*  $Id: Original yen's algorithm implementation  10-09-08 06:48:36 yan.qi.asu$*/
/************************************************************************/

#include <limits>
#include <set>
#include <map>
#include <queue>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <memory>
#include "GraphElements.h"
#include "Graph.h"
#include "DijkstraShortestPathAlg.h"
//#include "Astar.h"
#include "YenTopKShortestPathsAlg.h"
#include <boost/numeric/ublas/matrix.hpp>
#include "colors.h"
#include <tuple>
#include <boost/algorithm/string.hpp>
#include <ctime>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/graph/graphviz.hpp>
#include <sys/time.h>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <math.h>    // for sqrt

//#include<graphics.h>
//#include <curses.h>
//#include <boost/timer/timer.hpp>
using namespace std;

bool debug = false;
float input_lambda = 0;
int N = 3;
int C = -1;
int K = -1;
bool min_cal = false;
bool random_destinations = false;
int random_positions = 0;
int random_seed = 1;
float optimal_limit = 1.20;
map<int, int> optimal_distances;
int start = 0;
vector<int> Destinations;
vector<bool> DestAdjacentVec;
float YenTime = 0;
unordered_map<int, pair<int, int>> node_map;
map<string, int> coord_map;
map<pair<int, int>, int> coord_map2;
vector<int> OrigDests;
float avg_dest_intradist = 0;
float avg_dest_orig_dest = 0;
int iter_limit = 200;
set<int> all_destinations;
set<int> shortest_paths_determined;
int min_dist_to_dest = 10;
int max_dist_to_dest = 1;
bool random_origin_only = false;
int lambda_add = 0;
int max_x = 0;
int max_y = 0;
int file_node_counter = 0;
string input_filename = "";
map<int, int> dest_index;
int connectivity = 0;
bool map_display = false;
bool stocastic_go = false;
bool stocastic_go_observer = false;
bool stocastic_ficitious_play = false;
int strategy = 0;
double lambdaDest = 0;
double Budget = 0;
double Beta = 0;
bool acyclic=true;
bool grandparent_check=false;
vector<map<pair<int,float>,float > > lambda_obs;
vector<map<unsigned,float> > lambda_obs_dest_given_prefix;
//vector<vector<map<unsigned,float> > > prev_lambda_obs_dest_given_prefix;
vector<vector<float> > q_given_dest_and_prefix;
//vector<vector<vector<float> > > prev_q_given_dest_and_prefix;
vector<float> avg_q_prefix;
vector<float> phi_prefix;
set<unsigned> cutset_prefix;
//vector< set<unsigned> > prev_cutset_prefix;
vector< pair<unsigned, set<unsigned> > > prev_cutset_prefix;
//[iters][cut][repetitions,d0,d1,...dn]
vector<map<unsigned,vector<unsigned> > > prev_reward_per_prefix;    
map<unsigned,vector<unsigned> > merged_prev_cut_data;
vector<map<unsigned,vector<unsigned> > > prev_cut_data;

map<pair<int,float>,float> zeta_obs;
set<pair<int,float> > Active_AG;
set<pair<int,float> > X_Hat;
set<pair<int,float> > X_Crit;
map<pair<int,float>,int > X_Crit_dest;
multimap< pair<int, pair<int,float> >,vector<pair<int,float> > > subpaths;
bool Observer_Correct=false;
vector<float> prob_path;//Mu for stocastic observer
map<pair<unsigned, unsigned>, float> prefix_prob_path;//Mu for stocastic observer given a prefix 
vector<vector<float> > dest_prob_path;
map<int, int> DestinationsOrder1;
map<int, int> DestinationsOrder2;
vector<vector<float>> dest_prob_per_prefix;
//REMOVE IF NOT DOING DRAWINGS!!!
bitset<3> color_map[256][256]; //origin/destination/FinalPath
//bitset<3> color_map[10][10];//origin/destination/FinalPath

Graph main_graph;
set<int> main_graph_nodes;
//Graph main_graph_FO;
unique_ptr<DijkstraShortestPathAlg> main_dijkstra_alg;

vector<Graph> current_graph;
vector<Graph> current_graph_dest_to_dest;

vector<DijkstraShortestPathAlg> Dijkstra_algs;
vector<DijkstraShortestPathAlg> Dijkstra_algs_dest_to_dest;

map<int, pair<float, vector<std::shared_ptr<BasePath>>>> best_rel_cost_per_dest;

std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::duration<long int, std::ratio<1l, 1000000000l>>> start_time = chrono::high_resolution_clock::now();
std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::duration<long int, std::ratio<1l, 1000000000l>>> orig_start_time = chrono::high_resolution_clock::now();
std::chrono::_V2::system_clock::time_point main_start_time;
vector<vector<BasePath *>> orig_to_dest_shortest_paths;
vector<Link> lks;
int total_nodes = 0;
set<int> reachable_nodes; //For SaT Graphs
set<int> nodes_set;		  //For SaT Graphs
bool create_PO_graph = false;
int PO_level = 100;
unordered_set<int> FO_nodes;

vector<map<pair<pair<int, float>, pair<int, float> >, float> > Alpha;
map<pair<unsigned,unsigned>, float> AlphaPrefix;
map<pair<pair<int, float>, pair<int, float> >, float> C_LambdaAbs;
vector<float> barX;
float barQ=0.9;
int target_destination=0;
double epsilon=pow(0.1,20);
piecewise_constant_distribution<> paths_dist;
float initial_pred=0;
float savings=0;
float savings_adj_prefix=0;
size_t simulation_runs = 100000;
int saving_groups=10;
vector<float> avg_savings(saving_groups,0);
bool LoadPaths=false;
string input_paths_file = "input_paths.txt";
vector<vector<unsigned>> input_paths;
//std::map<std::pair<int,float>,std::pair<int,float> > obs_est_dest;
vector<pair<unsigned, float> > best_target_dest_choice;
vector<vector<pair<unsigned, float> > > prev_best_target_dest_choice;
std::shared_ptr<std::vector<std::set<unsigned int>>> paths_per_prefix;
std::shared_ptr<std::map<std::pair<unsigned int, unsigned int>, std::set<unsigned int>>> dest_paths_per_prefix;
vector<float> target_iterative_rewards;
size_t iterations=10000;
size_t current_iter=0;
float saturation=0;
int q_type = 0;
int total_cutsets = 0;
double normalization_param=0;



struct found_goal {}; // exception for termination

// visitor that terminates when we find the goal
template <class Vertex>
class astar_goal_visitor : public boost::default_astar_visitor
{
public:
  astar_goal_visitor(Vertex goal) : m_goal(goal) {}
  template <class Graph>
  void examine_vertex(Vertex u, Graph& g) {
    if(u == m_goal)
      throw found_goal();
  }
private:
  Vertex m_goal;
};

std::shared_ptr<std::vector<std::set<unsigned int>>> edges;
std::shared_ptr<std::map<int, std::vector<BaseVertex *>>> paths;
//std::map<int, map<int,float> > prob_cuts_per_path;
std::shared_ptr<std::vector<std::pair<unsigned int, unsigned int>>> prefixes;

template <typename S>
auto select_random(const S &s, size_t n)
{
    auto it = std::begin(s);
    // 'advance' the iterator n times
    std::advance(it, n);
    return it;
}
int get_manhattan_dist(int pos1, int pos2)
{
    int x_dist = abs(node_map[pos1].first - node_map[pos2].first);
    int y_dist = abs(node_map[pos1].second - node_map[pos2].second);

    return x_dist + y_dist;
}

void testDijkstraGraph()
{
    Graph *my_graph_pt = new Graph("data/test_1");
    DijkstraShortestPathAlg shortest_path_alg(my_graph_pt);
    BasePath *result =
        shortest_path_alg.get_shortest_path(
            my_graph_pt->get_vertex(0), my_graph_pt->get_vertex(5));
    result->PrintOut(std::cout);
}

void testYenAlg(const int &k,
                Link *lk,
                const int &size,
                const int &s,
                const int &d,
                const int &nodes,
                std::shared_ptr<PathMatrix> PthMat)
{
    auto start_time = std::chrono::high_resolution_clock::now();
    cout << "K:" << K << "," << "s:"<<s<<",d:"<<d<<",nodes:"<<nodes<<endl;
    int budget = INT_MAX;
    debug = true;

    Graph my_graph;
    auto all_paths = std::make_shared<set<vector<int>, vect_comp_int>>();

    my_graph.set_number_vertices(nodes);

    for (int i = 0; i < size; i++)
    {
        my_graph.add_link(lk[i].u, lk[i].v, lk[i].weight);
        //if(debug)
        //    std::cout<<"\tAdded Yen Graph node["<<lk[i].u<<","<<lk[i].v<<"],w="<<lk[ i ].weight<<endl;
    }
    if (debug)
        std::cout << "Yen Graph size:" << size << ",d:" << d << ",origin:" << s << endl;

    my_graph.setv();

    YenTopKShortestPathsAlg yenAlg(my_graph,
                                   my_graph.get_vertex(s),
                                   my_graph.get_vertex(d));

    //TESTING HACK
    //std::cout << "TESTING DIJKSTRA TO ALL PATHS" << endl;
    //yenAlg.get_shortest_distance_to_all_nodes(my_graph.get_vertex(s));
    //exit(1);
    //REMOVE ABOVE FROM HERE!!!!

    // Output the k-shortest paths
    int i = 0;
    while (yenAlg.has_next() && i < k)
    {
        ++i;
        BasePath current_path = *(yenAlg.next());
        if (optimal_distances[d] != 0)
        { //it is ok to initialize to 0 if unpopulated, on purpose!
            budget = optimal_distances[d] * optimal_limit;
            //std::cout<<"optimal_distance["<<d<<"]:"<<optimal_distances[d]<<",budget:"<<budget<<endl;
        }
        if (current_path.length() > budget)
        {
            //std::cout<<"paths_created for current destination:"<<d<<",finish, all future paths will have length bigger than:"<<current_path.length()<<",max_length allowed:"<<C<<endl;
            break;
        }
        if (budget == INT_MAX)
        { //first path is optimal for each destination
            optimal_distances[d] = current_path.length();
            std::cout << "optimal_distances[" << d << "]:" << optimal_distances[d] << endl;
        }

        //current_path.PrintOut(std::cout);
        current_path.add_paths_set(all_paths);
        //std::cout<<"\t added path"<<endl;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - orig_start_time;
    auto current_time = diff.count();
    cout << "Yen added:" << i << "paths" << "in,"<<current_time<<" secs"<<endl;
    PthMat->add_paths_set(all_paths, d);
}
int optimizing_observer(PathMatrix *PM, int nodes, int t_max)
{
    int t = 1;
    int disclosing_paths = 0;
    set<int> unsensed_nodes;
    PathMatrix PM_original = *PM;

    PathMatrix PM_temp;

    std::cout << "NOW OPTIMIZING OBSERVER NETWORK TO KEEP T_MAX LEVEL OF:" << t_max << endl;

    for (int w = 1; w < nodes; w++)
    {
        PM_temp = PM_original;
        unsensed_nodes.insert(w);
        PM_temp.set_W(&unsensed_nodes);

        for (t = 1; t < 10000; t++)
        {
            if (t > t_max)
            {
                std::cout << "\t\tFinished, this node is unsafe, t score raised by at least one level to:" << t << endl;
                break;
            }
            if (debug)
                std::cout << "\tWorking with t:" << t << ",candidate_node:" << w << endl;
            PM_temp.get_safe_nodes(t);
            if (PM_temp.erase_disclosing_paths(t, disclosing_paths))
            {
                if (debug)
                    std::cout << "\t\tFinished, erasing path would leave graph disconnected,final t score:" << t << " is lower than t_max:" << t_max << endl;
                break;
            }
        }

        if (t <= t_max)
        { //we can safely remove this node from observation grid
            std::cout << "\t t:" << t << ",t_max:" << t_max << ",safely removed " << w << " node from observation list" << endl;
        }
        else
        {
            std::cout << "\t t:" << t << ",t_max:" << t_max << ",unsafely removed " << w << " node , it needs to be in observation list" << endl;
            unsensed_nodes.erase(w);
        }
    }

    std::cout << "The following nodes can be safely ignored by the observer" << endl;
    for (auto w : unsensed_nodes)
        std::cout << w << ",";
    std::cout << endl;
    std::cout << "The following nodes need to be monitored by the observer to keep t_level untouched:" << endl;
    for (int i = 1; i < nodes; i++)
    {
        if (unsensed_nodes.find(i) != unsensed_nodes.end())
        {
            continue;
        }
        std::cout << i << ",";
    }
    std::cout << endl;
    return 0;
}
//optimizing_target_t_index optimizes in terms of length of path
//not in temrs of distance to goal
int optimizing_target_t_index(shared_ptr<PathMatrix> PM)
{
    int t = 1;
    int disclosing_paths = 0;
    for (t = 1; t < 1000; t++)
    {
        if (debug)
            std::cout << "Working with t:" << t << endl;
        PM->get_safe_nodes2(t, 2);
        //std::cout<<"after get_safe_nodes2"<<flush<<endl;
        if (PM->erase_disclosing_paths(t, disclosing_paths))
        {
            std::cout << "Finished, erasing path would leave graph disconnected,final T score:" << t << endl;
            break;
        }
        else if (debug)
        {
            std::cout << "AFTER ELIMINATING DISCLOSING PATHS(" << t << "), SURVIVING PATHS:" << endl;
            std::cout << "____________________________" << endl;
            PM->print_paths();
            std::cout << "____________________________" << endl;
        }
    }

    if (debug)
    {
        std::cout << "AFTER ELIMINATING DISCLOSING PATHS(" << t << "), SURVIVING PATHS:" << endl;
        std::cout << "____________________________" << endl;
        PM->print_paths();
        std::cout << "____________________________" << endl;
    }

    PM->print_final_results();
    return t;
}
//Greedy hill climb algorithm eliminating paths to
//reduce target's best distance to goal while remaining undisclosed
int optimizing_target_goal_dist(shared_ptr<PathMatrix> PM)
{
    start_time = std::chrono::system_clock::now();
    int t = 1;
    int disclosing_paths = 0;
    for (t = 1; t < 1000; t++)
    {
        if (debug)
            std::cout << "Working with t:" << t << endl;
        PM->get_safe_nodes2(t, 2);
        if (PM->erase_disclosing_paths(t, disclosing_paths))
        {
            std::cout << "Finished, erasing path would leave graph disconnected,final T score:" << t << endl;
            break;
        }
        else if (debug)
        {
            std::cout << "AFTER ELIMINATING DISCLOSING PATHS(" << t << "), SURVIVING PATHS:" << endl;
            std::cout << "____________________________" << endl;
            PM->print_paths();
            std::cout << "____________________________" << endl;
        }
    }

    if (debug)
    {
        std::cout << "AFTER ELIMINATING DISCLOSING PATHS(" << t << "), SURVIVING PATHS:" << endl;
        std::cout << "____________________________" << endl;
        PM->print_paths();
        std::cout << "____________________________" << endl;
    }

    PM->print_final_results();
    PM->set_InitialLambda(PM->get_max_lambda());
    PM->set_old_max_lambda(PM->get_max_lambda());
    PM->check_for_improvement_on_lambda_values();
    return t;
}
int optimizing_target_goal_dist_U_optimized(shared_ptr<PathMatrix> PM)
{
    start_time = std::chrono::system_clock::now();
    PM->populate_cover_matrix();
    PM->calculate_path_min_lambdas_from_cover_matrix();
    PM->record_best_solution();
    PM->set_old_max_lambda(PM->get_max_lambda());
    PM->set_InitialLambda(PM->get_max_lambda());
    for (int i = 0; i < iter_limit; i++)
    {
        //timings
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end_time - orig_start_time;
        auto current_time = diff.count();
        std::cout << ",Iteration:," << i << ",current_time:," << current_time;
        //elminate and reset
        //PM->elminate_max_lambda_paths_U();
        PM->remove_maximal_lambda_paths();
        //calling set_path in eliminatate_max_lambda_paths_U
        //PM->set_path_values_from_cover_paths();
        int max_lambda = PM->get_max_lambda();
        if (max_lambda == 0)
        {
            std::cout << "Finished, max_lambda cannot be improved from 0!" << endl;
            if (PM->is_lambda_connected() == Destinations.size())
            { //Need to ensure all destinations still reachable before updating best solution
                PM->set_old_max_lambda(0);
                PM->record_best_solution();
                PM->print_exit_statement();
            }
            exit(0);
        }

        int old_max_lambda = PM->get_old_max_lambda();
        if (max_lambda < old_max_lambda)
        {
            std::cout << "iteration:," << i << ",max_lambda reduced from " << max_lambda << " to old lambda" << old_max_lambda << ",continue" << endl;
            PM->record_best_solution();
            PM->set_old_max_lambda(max_lambda);
        }
        else if (max_lambda == old_max_lambda)
        {
            std::cout << "iteration:" << i << ",no improvement from old lambda,keep removing all paths with lambda==max_lambda:" << max_lambda << endl;
        }
        else
        {
            std::cout << "iteration:," << i << ",max_lambda increased from:," << old_max_lambda << ",to max_lambda:," << max_lambda << ",hence eliminating all paths whose lambda>old_lambda:," << old_max_lambda << endl;
        }
        PM->print_final_results();
    }
    return 0;
}

void create_only_grid_edges_from_file(vector<Link> *grid_edges)
{
    bool is_dest = false;
    int current_destination = 0;

    for (int x = 0; x < max_x; x++)
    {
        for (int y = 0; y < max_y; y++)
        {
            is_dest = false;
            pair<int, int> curr_pos = make_pair(x, y);
            if (coord_map2.find(curr_pos) == coord_map2.end())
            { //unreachable node
                //std::cout<<"\tnode is unreachable"<<endl;
                continue;
            }
            if (all_destinations.find(coord_map2[curr_pos]) != all_destinations.end())
            {	//Destinations have no outgoing edges
                //std::cout<<"\tnode is destination"<<endl;
                is_dest = true;
                current_destination = coord_map2[curr_pos];
            }
            pair<int, int> W_pos(x, y + 1);
            pair<int, int> E_pos(x + 1, y + 1);
            pair<int, int> D_pos(x + 1, y);
            pair<int, int> C_pos(x + 1, y - 1);
            pair<int, int> X_pos(x, y - 1);
            pair<int, int> Z_pos(x - 1, y - 1);
            pair<int, int> A_pos(x - 1, y);
            pair<int, int> Q_pos(x - 1, y + 1);
            //W-neighour
            if (y + 1 < max_y && coord_map2.find(W_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[W_pos], 1};
                //std::cout<<"\tW"<<coord_map2[curr_pos]<<"->"<<coord_map2[W_pos]<<flush<<endl;
                //std::cout<<"\t["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<W_pos.first<<","<<W_pos.second<<"]"<<endl;
                if (!is_dest)
                {
                    if (coord_map2[W_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges->push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[W_pos]);
                    }
                }
            }
            //E-neighour
            if (connectivity == 8 && y + 1 < max_y && x + 1 < max_x && coord_map2.find(E_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[E_pos], 1};
                //std::cout<<"\tE"<<coord_map2[curr_pos]<<"->"<<coord_map2[E_pos]<<flush<<endl;
                //std::cout<<"\t["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<E_pos.first<<","<<E_pos.second<<"]"<<endl;
                if (!is_dest)
                {
                    if (coord_map2[E_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges->push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[E_pos]);
                    }
                }
            }
            //D-neighour
            if (x + 1 < max_x && coord_map2.find(D_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[D_pos], 1};
                //std::cout<<"\tD["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<D_pos.first<<","<<D_pos.second<<"]"<<endl;
                if (!is_dest)
                {
                    if (coord_map2[D_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges->push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[D_pos]);
                    }
                }
            }
            //C-neighour
            if (connectivity == 8 && x + 1 < max_x && y - 1 > -1 && coord_map2.find(C_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[C_pos], 1};
                //std::cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[C_pos]<<endl;
                //std::cout<<"\tC["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<C_pos.first<<","<<C_pos.second<<"]"<<endl;
                if (!is_dest)
                {
                    if (coord_map2[C_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges->push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[C_pos]);
                    }
                }
            }
            //X-neighour
            if (y - 1 > -1 && coord_map2.find(X_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[X_pos], 1};
                //std::cout<<"X\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
                //std::cout<<"\tX["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<X_pos.first<<","<<X_pos.second<<"]"<<endl;
                if (!is_dest)
                {
                    if (coord_map2[X_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges->push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[X_pos]);
                    }
                }
            }
            //Z-neighour
            if (connectivity == 8 && y - 1 > -1 && coord_map2.find(Z_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[Z_pos], 1};
                //std::cout<<"\tZ["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<Z_pos.first<<","<<Z_pos.second<<"]"<<endl;
                //std::cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
                if (!is_dest)
                {
                    if (coord_map2[Z_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges->push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[Z_pos]);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.v);
                    }
                }
            }
            //A-neighour
            if (x - 1 > -1 && coord_map2.find(A_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[A_pos], 1};
                //std::cout<<"\tA["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<A_pos.first<<","<<A_pos.second<<"]"<<endl;
                //std::cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
                if (!is_dest)
                {
                    if (coord_map2[A_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges->push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[A_pos]);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.v);
                    }
                }
            }
            //Q-neighour
            if (connectivity == 8 && x - 1 > -1 && coord_map2.find(Q_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[Q_pos], 1};
                //std::cout<<"\tQ["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<Q_pos.first<<","<<Q_pos.second<<"]"<<endl;
                //std::cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[Q_pos]<<endl;
                if (!is_dest)
                {
                    if (coord_map2[Q_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges->push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[Q_pos]);
                    }
                }
            }
        }
    }
}
void create_grid_nodes_from_file()
{
    std::ifstream in(input_filename.c_str());
    int row = -4;
    int node_counter = 0;
    string str;
    //vector<string> vecOfStrs;
    char point{'.'};
    while (getline(in, str))
    {
        if (str.find("quarterly") != string::npos)
        {
            std::cout << "4-way Connectivity" << endl;
            connectivity = 4;
        }
        else if (str.find("octile") != string::npos)
        {
            std::cout << "8-way Connectivity" << endl;
            connectivity = 8;
        }

        max_x = max(size_t(max_x), str.length() - 1);

        // Line contains string of length > 0 then save it in vector
        //if(str.size() > 0)
        //  vecOfStrs.push_back(str);
        //create nodes
        //std::cout<<"row:"<<row<<",string:"<<str<<endl;
        if (row >= 0)
        { //skipping first 4 rows
            //std::cout<<str<<",size:,"<<str.size()<<endl;
            for (int i = 0; i < str.size(); i++)
            {
                if (str[i] != point)
                {
                    continue;
                }
                //string n1=to_string(row) + "_" + to_string(i);
                //coord_map[n1]=node_counter;
                node_map[node_counter] = make_pair(row, i);
                coord_map2[make_pair(row, i)] = node_counter++;
                //std::cout<<"node_map["<<node_counter-1<<"]:"<<node_map[node_counter-1].first<<","<<node_map[node_counter-1].second<<endl;
                //std::cout<<"coord_map2["<<row<<","<<i<<"]:,"<<coord_map2[make_pair(row,i)]<<endl;
            }
        }
        row++;
    }
    std::cout << "Connectivity:," << connectivity << endl;
    max_y = row;
    file_node_counter = node_counter - 1;
}
void random_orig_and_dest_placement_SaT()
{
    cout<<"hola random_orig_and_dest_placement_SaT"<<endl;
    Destinations.clear();
    all_destinations.clear();
    dest_index.clear();
    srand(random_seed);
    int r = rand() % nodes_set.size(); // not _really_ random
    auto n = *select_random(nodes_set, r);
    start = n;
    std::cout << "start:," << start << ",D:,";
    auto temp_nodes_set=nodes_set;
    temp_nodes_set.erase(start);
    while (Destinations.size() < random_positions)
    {
        if(temp_nodes_set.size()==0){
            cout <<endl<< "NO MORE NODES TO CHOOSE FOR DESTINATIONS!!!" << endl;
            exit(1);
        }
        int r = rand() % temp_nodes_set.size(); // not _really_ random
        auto n = *select_random(temp_nodes_set, r);
        if (n == start){ //skip origin as possible destination
            temp_nodes_set.erase(n);
            continue;
        }
        else if(all_destinations.count(n)>0){//Dest already chosen!!!
            //cout << "skip prev chosen dest:" << n << endl;
            temp_nodes_set.erase(n);
            continue;
        }
        else{
            temp_nodes_set.erase(n);
        }
        Destinations.push_back(n);
        all_destinations.insert(n);
        dest_index[n] = Destinations.size() - 1;
        std::cout << n << ",";
    }
    std::cout << endl;
}
void create_prefix_graph_from_SaT_file(){
    std::cout << "hello create_prefix,stocastic GO1" << flush << endl;
    std::ifstream in(input_filename.c_str());
    std::ifstream in2(input_filename.c_str());
    string str;
    //vector<vector<Link>> boundary_edges;
    //vector<vector<Link>> dest_to_dest_edges;
    vector<vector<Link>> lks_boundary;
    vector<vector<Link>> lks_dest;
    total_nodes = 0;
    int total_edges = 0;

    //read SaT_file to create main_graph
    //First pass to get list of all nodes
    while (getline(in2, str))
    {
        vector<string> results;
        boost::split(results, str, boost::is_any_of(","), boost::token_compress_on);
        auto ret = nodes_set.insert(stoi(results[0]));
        if (ret.second == true)
            total_nodes++;
        ret = nodes_set.insert(stoi(results[1]));
        if (ret.second == true)
            total_nodes++;
    }
    std::cout << "total_nodes=" << total_nodes << endl;
    
    //If random, choose from nodes_set randomly
    //if a destination proves unreachable, WE ARE FINISHED
    //SHOULD FIND A BETTER WAY
    if (random_destinations)
    { //Both origin and destination random selection
        random_orig_and_dest_placement_SaT();
    }
    
    //Now initialize all variables dependent on Destinations

    current_graph.resize(Destinations.size());
    current_graph_dest_to_dest.resize(Destinations.size());

    lks_boundary.resize(Destinations.size());
    lks_dest.resize(Destinations.size());

    while (getline(in, str))
    {
        vector<string> results;
        boost::split(results, str, boost::is_any_of(","), boost::token_compress_on);
        if (stoi(results[1]) == start)
        {
            continue; //skipping origin incoming edges
        }
        if (stoi(results[0]) == stoi(results[1]))
        { //no loops
            std::cout << "\t\tskipped loop,from," << stoi(results[0]) << ",to," << stoi(results[1]) << endl;
            continue;
        }
        if (all_destinations.find(stoi(results[0])) != all_destinations.end())
        {
            continue; //ignoring outgoing edges from destinations
        }
        else if (all_destinations.find(stoi(results[1])) != all_destinations.end())
        { //incoming destination edge
            //Adding reverse edge for backwards dest_to_dest and boundary graphs
            Link reverse_link{stoi(results[1]), stoi(results[0]), stof(results[2])};
            lks_boundary[dest_index[stoi(results[1])]].push_back(reverse_link);
            lks_dest[dest_index[stoi(results[1])]].push_back(reverse_link);
            std::cout << "Added2 reverse link for destination_index:," << dest_index[stoi(results[1])] << ",from:," << stoi(results[1]) << ",to:," << stoi(results[0]) << ",weight:," << stof(results[2]) << endl;
        }
        Link temp_link{stoi(results[0]), stoi(results[1]), 100*round2(stof(results[2]))};

        lks.push_back(temp_link);
        total_edges++;
    }
    std::cout << "main_graph, total_edges:," << total_edges << endl;
    file_node_counter = total_nodes;
    //Now create Dijkstra Graph
    main_graph.set_number_vertices(file_node_counter);
    set<int> DestinationsSet(Destinations.begin(), Destinations.end());
    DestAdjacentVec.assign(file_node_counter, false);
    for (long i = 0; i < total_edges; i++)
    {
        main_graph.add_link(lks[i].u, lks[i].v, lks[i].weight);
        if(DestinationsSet.count(lks[i].v)>0){
            //DestAdjacent.insert(lks[i].u);
            DestAdjacentVec[lks[i].u] = true;
            //cout << "\tnode:," << lks[i].u << ",is adjacent to destination," << lks[i].v << endl;
        }
    }

    main_graph.setv();
    std::cout << "total_nodes:" << total_nodes << ",main_graph:" << main_graph.get_number_vertices() << endl;
    main_dijkstra_alg = make_unique<DijkstraShortestPathAlg>(&main_graph);
    auto timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    std::cout << "before determine_all_shortest_paths,time:" << ctime(&timenow) << endl;
    main_dijkstra_alg->determine_all_shortest_paths(main_graph.get_vertex(start), -1);
    timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    std::cout << "after determine_all_shortest_paths,time:" << ctime(&timenow) << endl; 
    
    //STOCASTIC CODE
    //First get current cost for any node
    std::cout << "node costs_reg:" << endl;
    int reachable_nodes = 0;
    for (auto node : nodes_set)
    {
        if (!main_dijkstra_alg->is_node_reachable(node))
        {
            //std::cout<<"node "<<node<<" is not reachable"<<endl;
        }
        else
        {
            reachable_nodes++;
            //std::cout<<"start:"<<start<<",reachable node:"<<node<<",cost:"<<main_dijkstra_alg->get_cost(node)<<endl;
        }
    }
    std::cout << "Reachable nodes:" << reachable_nodes << endl;
    //CHECK destinations are reachable
    for (auto destination : Destinations)
    {
        if (!main_dijkstra_alg->is_node_reachable(destination))
        {
            cerr << "Destination:," << destination << ",not reachable from origin:," << start << ",stopping search" << endl;
            exit(1);
        }
        else
        {
            std::cout << "Destination:," << destination << ",is reachable from origin:," << start << ",continuing search" <<",cost:"<<main_dijkstra_alg->get_cost(destination)<< endl;
        }
    }
    std::cout << "hola12, Finished checking destinations are reachable" << endl;
    //States are not any more pair <node,g> but <node,prefix>
    //we call Λ^{d}_{ss'} the probability that the target transits from v to v' if its destination is d and has already incurred a cost x to reach v.
    //We need to determine if a state is feasible
    //That requires paths to destinations
    for (size_t dest = 0; dest < Destinations.size(); dest++)
    {
        //boundary graphs
        current_graph[dest].set_number_vertices(total_nodes); //skipping origin
        for (auto add : lks_boundary[dest])
        {
            //add reverse edges specific to boundary graphs
            if (!(main_dijkstra_alg->is_node_reachable(add.v) && main_dijkstra_alg->is_node_reachable(add.u)))
            {
                continue;
            }
            //std::cout<<"adding current_graph reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<",vertices:"<<file_node_counter - 1 <<endl;
            current_graph[dest].add_link(add.u, add.v, add.weight);
        }

        current_graph_dest_to_dest[dest].set_number_vertices(total_nodes); //skipping origin
        for (auto add : lks_dest[dest])
        {
            //add reverse edges specific to dest_to_dest graphs
            if (!(main_dijkstra_alg->is_node_reachable(add.v) && main_dijkstra_alg->is_node_reachable(add.u)))
            {
                continue;
            }
            //std::cout<<"adding current_graph_dest_to_dest reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<endl;
            current_graph_dest_to_dest[dest].add_link(add.u, add.v, add.weight);
        }
        for (long i = 0; i < total_edges; i++)
        {
            if (lks[i].u == start)
            {
                continue; //skipping edges leaving origin for boundary
            }
            else if (all_destinations.find(lks[i].v) != all_destinations.end())
            { //skipping edges into any destination for boundary graphs
                //but keeping them for dest to dest search
                /*current_graph_dest_to_dest[dest].add_link(lks[i].u, lks[ i ].v, lks[ i ].weight );
    std::cout<<"Added origin outgoing link for destination_index:,"<<dest<<",from:,"<<lks[i].u<<",to:,"<<lks[i].v<<",weight:,"<<lks[i].weight<<endl;*/
                continue;
            }
            //std::cout<<"adding current_graph_reverse_link edges["<<dest<<"]"<<",lks[i].v:"<<lks[i].v<<",lks[i].u:"<<lks[i].u<<"lks[i].weight:"<<lks[i].weight<<endl;
            current_graph[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
        }
        current_graph[dest].setv();
        current_graph[dest].set_number_vertices(total_nodes); //skipping origin
        std::cout << "total_nodes:" << total_nodes << ",current_graph[" << dest << "]:," << current_graph[dest].get_number_vertices() << endl;
        Dijkstra_algs.push_back(DijkstraShortestPathAlg(&(current_graph[dest])));

        std::cout << "adding bidirectional edges for dest_to_dest[" << dest << "]" << endl;
        for (long i = 0; i < total_edges; i++)
        {
            if (lks[i].u == start)
            { //skipping origin outgoing edges
                continue;
            }
            //FOR CALCULATING GEODETIC DISTANCE, WE CREATE A LINK IN BOTH DIRECTIONS AS LONG AS
            //ONE DIRECTION EXISTS
            //std::cout<<"adding bidirectional edges for dest_to_dest["<<dest<<"]"<<",lks[i].u:"<<lks[i].u<<",lks[i].v:"<<lks[i].v<<"lks[i].weight:"<<lks[i].weight<<endl;
            current_graph_dest_to_dest[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
            current_graph_dest_to_dest[dest].add_link(lks[i].u, lks[i].v, lks[i].weight);
        }
        current_graph_dest_to_dest[dest].setv();
        current_graph_dest_to_dest[dest].set_number_vertices(total_nodes); //skipping origin
        std::cout << "total_nodes:" << total_nodes << ",current_graph_dest_to_dest[" << dest << "]:," << current_graph_dest_to_dest[dest].get_number_vertices() << endl;

        /*current_graph_dest_to_dest[dest].clear();
    current_graph_dest_to_dest[dest].set_number_vertices(3);
    if(dest==0){
      Link temp_link{116622,116623,1.5};
      current_graph_dest_to_dest[0].add_link(temp_link.u,temp_link.v,temp_link.weight);
      Link temp_link2{116623,117023,1.5};
      current_graph_dest_to_dest[0].add_link(temp_link2.u,temp_link2.v,temp_link2.weight);
    }
    else{
      current_graph_dest_to_dest[1].add_link(117023,116623,1.5);
      current_graph_dest_to_dest[1].add_link(116623,116622,1.5);
    }
    current_graph_dest_to_dest[dest].dump_edges();*/
        Dijkstra_algs_dest_to_dest.push_back(DijkstraShortestPathAlg(&(current_graph_dest_to_dest[dest])));
        for (size_t dest2 = 0; dest2 < Destinations.size(); dest2++)
        {
            if (dest == dest2)
                continue;

            BasePath *result = Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
                current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(Destinations[dest2]));
            //BasePath* result = Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
            //	  current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(25157));
            //std::cout<<"\t"<<result->get_path_string()<<",Weight:"<<result->Weight()<<endl;
            std::cout << "\tdest_to_dest,dest1:," << Destinations[dest] << ",dest2:" << Destinations[dest2] << ",Weight:" << result->Weight() << ",length:" << result->length() << endl;
            //BasePath* main_result1 = main_dijkstra_alg->get_shortest_path(
            //	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]));
            //    std::cout<<"\torig_to_dest1,start:,"<<start<<",dest2:"<<Destinations[dest]<<",Weight:"<<main_result1->Weight()<<",length:"<<main_result1->length()<<endl;
            //  BasePath* main_result2 = main_dijkstra_alg->get_shortest_path(
            //	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest2]));
            //    std::cout<<"\torig_to_dest2,start:,"<<start<<",dest2:"<<Destinations[dest2]<<",Weight:"<<main_result2->Weight()<<",length:"<<main_result2->length()<<endl;
        }
        //Get all paths from dest
        //Get distance from node v to destination
    }

    for (size_t dest = 0; dest < Destinations.size(); dest++)
    {
            std::cout << "hola,destination:" << Destinations[dest] << endl;
            Dijkstra_algs_dest_to_dest[dest].determine_all_shortest_paths(current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), -1);
            for (auto node : nodes_set)
            {
                if (!main_dijkstra_alg->is_node_reachable(node))
                {
                    std::cout<<"node "<<node<<" is not reachable from start"<<endl;
                    continue;
                }
        }
    }
    //NOW CREATE AUGMENTED GRAPH
    //GIVEN BUDGET C
    //ORIGIN Note IS (Orig,0)
    main_dijkstra_alg->origin = main_graph.get_vertex(start);
    main_dijkstra_alg->set_Destinations_Order(&Destinations);
    for (auto dest : Destinations)
    {
        main_dijkstra_alg->Destinations.insert(dest);
    }
    if(LoadPaths){
        cout << "Loading paths from paths.txt for debugging" << endl;
        main_dijkstra_alg->load_paths(input_paths);
    }
    else{//remove paths.txt if it exists before appending to it
        const char* filename = "paths.txt";
        remove(filename);
        for (size_t dest = 0; dest < Destinations.size(); dest++)
        {
            //std::cout<<"start:"<<start<<",destination:"<<Destinations[dest]<<endl;
            //std::cout << "before finding paths:" << ctime(&timenow) << endl;
            std::cout<<"AGP,grandpatent_check="<<grandparent_check<<",acyclic:"<<acyclic<<",Destination["<<dest<<"]:"<<Destinations[dest]<<",budget:"<<Budget*main_dijkstra_alg->get_cost(Destinations[dest])<<endl;
            normalization_param=max(main_dijkstra_alg->createAGYen(main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]), Budget * main_dijkstra_alg->get_cost(Destinations[dest]),true),normalization_param);
            //main_dijkstra_alg->createAG(main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]), Budget*main_dijkstra_alg->get_cost(Destinations[dest]),acyclic,grandparent_check);
            //std::cout << "after finding paths:" << ctime(&timenow) << endl;
        }
    }
    //main_dijkstra_alg->printAGFile();
    main_dijkstra_alg->populateAGPnodes();
    
    //Now we can calculate strategies as in google doc summarizing version 7 of the paper
    if (strategy == 0){
        auto start_time = std::chrono::high_resolution_clock::now();
        //string method("min_path_time");
        calculate_min_path_strategy_prefixes(false);
        //elapsed_time(method, start_time);
    }
    else{
        cerr << "Strategy:" << strategy << " undefined!!!" << endl;
        exit(1);
    }

    //Now we calculate recursively phi and Cutset starting at origin
    auto start_time = std::chrono::high_resolution_clock::now();
    string method("calculate_q_recursive");
    calculate_q_recursive_prefix(0);
    //cout << "Iter:,"<<current_iter<<",cutset:," <<cutset_prefix<< endl;
    removeDominatedFromCutset();
    cout << "Iter:,"<<current_iter<<",clean_cutset:," <<cutset_prefix<< endl;
    cout << "Iter:," << current_iter<<",clean_cutset_endings:,"; print_cutset_endings(cutset_prefix); cout << endl;
    elapsed_time(method, start_time);
    if(stocastic_ficitious_play){
        //Need to calculate initial best target choice first
        best_avg_target_choice();
        fictitious_play2();
        exit(10);
    }


    auto initial_pred=dest_predictor_prefix(0);
    cout << "initial_pred=," << initial_pred.first << "," << initial_pred.second << endl;
    //Now run simulation
    //First create a weighted distribution of paths
    size_t startValue = 0;
    size_t endValue = prob_path.size()+1;
    size_t sizeOfVector = endValue - startValue;
    std::vector<size_t> path_index(sizeOfVector);
    generate (
            path_index.begin(),
            path_index.end(),
            [&](){
                return startValue++;
            });

    piecewise_constant_distribution<> temp_dist(std::begin(path_index),
                                                std::end(path_index),
                                                std::begin(prob_path));
    paths_dist = temp_dist;
    //paths_dist = make_shared<piecewise_constant_distribution>(paths_dist);
    //Make simulations trully random
    srand((int) time(0));
    start_time = std::chrono::high_resolution_clock::now();
    auto best_target_choice = best_avg_target_choice();
    elapsed_time("target_choice", start_time);
    //cout << "best_target_choice path[," << best_target_choice.first << "]:,";
    //main_dijkstra_alg->print_path(best_target_choice.first);
    //cout << ",reward:," << best_target_choice.second << endl;

    auto simulation_start_time = std::chrono::high_resolution_clock::now();
    for (unsigned i=0;i<simulation_runs;i++)
        simulation_observer_prefix();
    elapsed_time("simulation_observer_prefix", simulation_start_time);
    cout << "savings:,"<<savings/float(simulation_runs)<<",simulations:," << simulation_runs << ",savings:"<<savings<< endl;
    double sum = std::accumulate(avg_savings.begin(), avg_savings.end(), 0.0);
    double mean = sum / avg_savings.size();

    std::vector<double> diff(avg_savings.size());
    std::transform(avg_savings.begin(), avg_savings.end(), diff.begin(), [mean](double x) { return x - mean; });

    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / avg_savings.size());
    cout << "mean_savings:," << mean/(simulation_runs/saving_groups) << ",std_dev:," << stdev/(simulation_runs/saving_groups) << ",initial_savings:"<<initial_pred.second<<",dest_adj_savings:,"<<savings_adj_prefix/simulation_runs<<",avg_target_choice:,"<<best_target_choice<<",groups:," << saving_groups << endl;
    elapsed_time("main", main_start_time);
    exit(10);
}
void create_augmented_graph_from_SaT_file()
{ //
    std::cout << "hello create_augmented_graph,stocastic GO1" << flush << endl;
    std::ifstream in2(input_filename.c_str());
    string str;
    /*
    if(C>0){
        auto destinations_filename=input_filename.substr(0,input_filename.length()-4);
        destinations_filename+="_OD.csv";
        std::ifstream in(destinations_filename.c_str());
        std::cout<<"destinations_filename="<<destinations_filename<<endl;
        Destinations.clear();
        while (getline(in, str))
        {
            //str.erase(std::remove(str.begin(),str.end(),'\"'),str.end());
            vector<string> results;
            boost::split(results, str, boost::is_any_of(","), boost::token_compress_on);
            if(results[0]=="origin"){
                start=stoi(results[1]);
            }
            else{
                for(int i=1;i<results.size();i++){
                    Destinations.push_back(stoi(results[i]));
                }
            }
        }
    }*/
    lks.reserve(16000);
    total_nodes = 0;
    int total_edges = 0;

    //vector<vector<Link> > boundary_edges;
    //vector<vector<Link> >dest_to_dest_edges;
    vector<vector<Link>> lks_boundary;
    vector<vector<Link>> lks_dest;

    //std::cout<<"hello1"<<flush<<endl;
    //First pass to get list of all nodes
    std::ifstream in(input_filename.c_str());
    while (getline(in2, str))
    {
        vector<string> results;
        boost::split(results, str, boost::is_any_of(","), boost::token_compress_on);
        //std::cout << "results[1]" << results[1] << endl;
        if(!random_destinations){
            if (all_destinations.find(stoi(results[0])) != all_destinations.end()){
                continue;//skipping adding any nodes due to outgoing destination edges
            }
        }
        auto ret = nodes_set.insert(stoi(results[0]));
        if (ret.second == true)
            total_nodes++;
        ret = nodes_set.insert(stoi(results[1]));
        if (ret.second == true)
            total_nodes++;
    }
    std::cout << "total_nodes=" << total_nodes << endl;
    //If random, choose from nodes_set randomly
    //if a destination proves unreachable, WE ARE FINISHED
    //SHOULD FIND A BETTER WAY
    if (random_destinations)
    { //Both origin and destination random selection
        random_orig_and_dest_placement_SaT();
    }
    //Now initialize all variables dependent on Destinations
    //boundary_edges.resize(Destinations.size());
    //dest_to_dest_edges.resize(Destinations.size());

    current_graph.resize(Destinations.size());
    current_graph_dest_to_dest.resize(Destinations.size());

    lks_boundary.resize(Destinations.size());
    lks_dest.resize(Destinations.size());

    current_graph.resize(Destinations.size());
    current_graph_dest_to_dest.resize(Destinations.size());

    while (getline(in, str))
    {
        vector<string> results;
        boost::split(results, str, boost::is_any_of(","), boost::token_compress_on);
        //std::vector<std::string> results((std::istream_iterator<string>(iss)),
        //                             std::istream_iterator<string>());
        //std::cout<<results[0]<<flush<<endl;
        //std::cout<<results[1]<<flush<<endl;
        //std::cout<<results[2]<<flush<<endl;
        //std::cout<<"hello1a"<<flush<<endl;
        //std::cout<<","<<results[1]<<flush<<","<<results[2]<<flush<<endl;
        //std::cout<<"dest_index[results[1]]:"<<flush<<stoi(results[1])<<flush<<endl;
        if (stoi(results[1]) == start)
        {
            continue; //skipping origin incomming edges
        }
        if (stoi(results[0]) == stoi(results[1]))
        { //no loops
            std::cout << "\t\tskipped loop,from," << stoi(results[0]) << ",to," << stoi(results[1]) << endl;
            continue;
        }
        if (all_destinations.find(stoi(results[0])) != all_destinations.end())
        {
            std::cout<<"Ignoring outgoing edge from destination:,"<<stoi(results[0])<<",to,"<<stoi(results[1])<<endl;
            continue; //ignoring outgoing edges from destinations
        }
        else if (all_destinations.find(stoi(results[1])) != all_destinations.end())
        { //incoming destination edge
            //Adding reverse edge for backwards dest_to_dest and boundary graphs
            Link reverse_link{stoi(results[1]), stoi(results[0]), stof(results[2])};
            //std::cout << "Added reverse link for destination_index:," << dest_index[stoi(results[1])] << ",from:," << stoi(results[1]) << ",to:," << stoi(results[0]) << ",weight:," << stof(results[2]) << endl;
        }
        Link temp_link{stoi(results[0]), stoi(results[1]), 100*round2(stof(results[2]))};

        lks.push_back(temp_link);
        total_edges++;
    }
    std::cout << "main_graph, total_edges:," << total_edges << endl;

    //Remove any nodes who only exist due to outgoing edge from a destination-designed node

    file_node_counter = total_nodes;
    //Now create Dijkstra Graph
    main_graph.set_number_vertices(file_node_counter);
    for (long i = 0; i < total_edges; i++)
    {
        main_graph.add_link(lks[i].u, lks[i].v, lks[i].weight);
    }

    main_graph.setv();
    std::cout << "total_nodes:" << total_nodes << ",main_graph:" << main_graph.get_number_vertices() << endl;
    main_dijkstra_alg = make_unique<DijkstraShortestPathAlg>(&main_graph);

    //BasePath* result_path = main_dijkstra_alg->get_shortest_path(
    //main_graph.get_vertex(start), main_graph.get_vertex(Destinations[0]));
    //main_graph.dump_edges();
    //for (auto edge : lks){
    //  std::cout<<"from:,"<<edge.u<<",->,"<<edge.v<<",weight:,"<<edge.weight<<endl;
    //}
    auto timenow =
        chrono::system_clock::to_time_t(chrono::system_clock::now());
    std::cout << "before determine_all_shortest_paths,time:" << ctime(&timenow) << endl;
    main_dijkstra_alg->determine_all_shortest_paths(main_graph.get_vertex(start), -1);
    timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    std::cout << "after determine_all_shortest_paths,time:" << ctime(&timenow) << endl; //exit(1);

    //STOCASTIC CODE
    //First get current cost for any node
    std::cout << "node costs_reg:" << endl;
    int reachable_nodes = 0;
    for (auto node : nodes_set)
    {
        if (!main_dijkstra_alg->is_node_reachable(node))
        {
            //std::cout<<"node "<<node<<" is not reachable"<<endl;
        }
        else
        {
            reachable_nodes++;
            std::cout<<"start:"<<start<<",reachable node:"<<node<<",cost:"<<main_dijkstra_alg->get_cost(node)<<endl;
        }
    }
    std::cout << "Reachable nodes:" << reachable_nodes << endl;
    //CHECK destinations are reachable
    for (auto destination : Destinations)
    {
        if (!main_dijkstra_alg->is_node_reachable(destination))
        {
            cerr << "Destination:," << destination << ",not reachable from origin:," << start << ",stopping search" << endl;
            exit(1);
        }
        else
        {
            std::cout << "Destination:," << destination << ",is reachable from origin:," << start << ",continuing search" << endl;
        }
    }
    std::cout << "hola12" << endl;
    //States are pair <node,g>
    //we call Λ^{d}_{ss'} the probability that the target transits from v to v' if its destination is d and has already incurred a cost x to reach v.
    //vector<map<pair<int,int>,float> > Delta;
    //Delta.resize(Destinations.size());
    //We need to determine if a state is feasible
    //That requires paths to destinations
    for (size_t dest = 0; dest < Destinations.size(); dest++)
    {
        //boundary graphs
        current_graph[dest].set_number_vertices(total_nodes); //skipping origin
        for (auto add : lks_boundary[dest])
        {
            //add reverse edges specific to boundary graphs
            if (!(main_dijkstra_alg->is_node_reachable(add.v) && main_dijkstra_alg->is_node_reachable(add.u)))
            {
                continue;
            }
            //std::cout<<"adding current_graph reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<",vertices:"<<file_node_counter - 1 <<endl;
            current_graph[dest].add_link(add.u, add.v, add.weight);
        }

        current_graph_dest_to_dest[dest].set_number_vertices(total_nodes); //skipping origin
        for (auto add : lks_dest[dest])
        {
            //add reverse edges specific to dest_to_dest graphs
            if (!(main_dijkstra_alg->is_node_reachable(add.v) && main_dijkstra_alg->is_node_reachable(add.u)))
            {
                continue;
            }
            //std::cout<<"adding current_graph_dest_to_dest reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<endl;
            current_graph_dest_to_dest[dest].add_link(add.u, add.v, add.weight);
        }
        for (long i = 0; i < total_edges; i++)
        {
            if (lks[i].u == start)
            {
                continue; //skipping edges leaving origin for boundary
            }
            else if (all_destinations.find(lks[i].v) != all_destinations.end())
            { //skipping edges into any destination for boundary graphs
                //but keeping them for dest to dest search
                /*current_graph_dest_to_dest[dest].add_link(lks[i].u, lks[ i ].v, lks[ i ].weight );
    std::cout<<"Added origin outgoing link for destination_index:,"<<dest<<",from:,"<<lks[i].u<<",to:,"<<lks[i].v<<",weight:,"<<lks[i].weight<<endl;*/
                continue;
            }
            //std::cout<<"adding current_graph_reverse_link edges["<<dest<<"]"<<",lks[i].v:"<<lks[i].v<<",lks[i].u:"<<lks[i].u<<"lks[i].weight:"<<lks[i].weight<<endl;
            current_graph[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
        }
        current_graph[dest].setv();
        current_graph[dest].set_number_vertices(total_nodes); //skipping origin
        std::cout << "total_nodes:" << total_nodes << ",current_graph[" << dest << "]:," << current_graph[dest].get_number_vertices() << endl;
        Dijkstra_algs.push_back(DijkstraShortestPathAlg(&(current_graph[dest])));

        std::cout << "adding bidirectional edges for dest_to_dest[" << dest << "]" << endl;
        for (long i = 0; i < total_edges; i++)
        {
            if (lks[i].u == start)
            { //skipping origin outgoing edges
                continue;
            }
            //FOR CALCULATING GEODETIC DISTANCE, WE CREATE A LINK IN BOTH DIRECTIONS AS LONG AS
            //ONE DIRECTION EXISTS
            //std::cout<<"adding bidirectional edges for dest_to_dest["<<dest<<"]"<<",lks[i].u:"<<lks[i].u<<",lks[i].v:"<<lks[i].v<<"lks[i].weight:"<<lks[i].weight<<endl;
            current_graph_dest_to_dest[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
            current_graph_dest_to_dest[dest].add_link(lks[i].u, lks[i].v, lks[i].weight);
        }
        current_graph_dest_to_dest[dest].setv();
        current_graph_dest_to_dest[dest].set_number_vertices(total_nodes); //skipping origin
        std::cout << "total_nodes:" << total_nodes << ",current_graph_dest_to_dest[" << dest << "]:," << current_graph_dest_to_dest[dest].get_number_vertices() << endl;

        /*current_graph_dest_to_dest[dest].clear();
    current_graph_dest_to_dest[dest].set_number_vertices(3);
    if(dest==0){
      Link temp_link{116622,116623,1.5};
      current_graph_dest_to_dest[0].add_link(temp_link.u,temp_link.v,temp_link.weight);
      Link temp_link2{116623,117023,1.5};
      current_graph_dest_to_dest[0].add_link(temp_link2.u,temp_link2.v,temp_link2.weight);
    }
    else{
      current_graph_dest_to_dest[1].add_link(117023,116623,1.5);
      current_graph_dest_to_dest[1].add_link(116623,116622,1.5);
    }
    current_graph_dest_to_dest[dest].dump_edges();*/
        Dijkstra_algs_dest_to_dest.push_back(DijkstraShortestPathAlg(&(current_graph_dest_to_dest[dest])));
        for (size_t dest2 = 0; dest2 < Destinations.size(); dest2++)
        {
            if (dest == dest2)
                continue;

            BasePath *result = Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
                current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(Destinations[dest2]));
            //BasePath* result = Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
            //	  current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(25157));
            //std::cout<<"\t"<<result->get_path_string()<<",Weight:"<<result->Weight()<<endl;
            std::cout << "\tdest_to_dest,dest1:," << Destinations[dest] << ",dest2:" << Destinations[dest2] << ",Weight:" << result->Weight() << ",length:" << result->length() << endl;
            //BasePath* main_result1 = main_dijkstra_alg->get_shortest_path(
            //	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]));
            //    std::cout<<"\torig_to_dest1,start:,"<<start<<",dest2:"<<Destinations[dest]<<",Weight:"<<main_result1->Weight()<<",length:"<<main_result1->length()<<endl;
            //  BasePath* main_result2 = main_dijkstra_alg->get_shortest_path(
            //	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest2]));
            //    std::cout<<"\torig_to_dest2,start:,"<<start<<",dest2:"<<Destinations[dest2]<<",Weight:"<<main_result2->Weight()<<",length:"<<main_result2->length()<<endl;
        }
    }

    if (stocastic_go)
    {
        //STOCASTIC CODE Example for dest[0]
        //Get all paths from dest[0]
        //Get distance from node v to destination
        for (size_t dest = 0; dest < Destinations.size(); dest++)
        {
            std::cout << "hola,destination:" << Destinations[dest] << endl;
            Dijkstra_algs_dest_to_dest[dest].determine_all_shortest_paths(current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), -1);
            for (auto node : nodes_set)
            {
                if (!main_dijkstra_alg->is_node_reachable(node))
                {
                    std::cout<<"node "<<node<<" is not reachable from start"<<endl;
                    continue;
                }
                /*else
                {
                    std::cout << "\treachable node:" << node << ",cost to destination[" << Destinations[dest] << "]:" << Dijkstra_algs_dest_to_dest[dest].get_cost(node);
                    std::cout << ",distance to start:," << main_dijkstra_alg->get_cost(node) << endl;
                }*/
            }
        }

        //Dijkstra_algs[0].print_perim_paths();

        //NOW CREATE AUGMENTED GRAPH
        //GIVEN BUDGET C
        //ORIGIN Note IS (Orig,0)
        main_dijkstra_alg->origin = main_graph.get_vertex(start);
        main_dijkstra_alg->set_Destinations_Order(&Destinations);
        for (auto dest : Destinations)
        {
            main_dijkstra_alg->Destinations.insert(dest);
        }
        for (size_t dest = 0; dest < Destinations.size(); dest++)
        {
            //std::cout<<"start:"<<start<<",destination:"<<Destinations[dest]<<endl;
            std::cout<<"AG,grandpatent_check="<<grandparent_check<<endl;
            main_dijkstra_alg->createAG(main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]), Budget*main_dijkstra_alg->get_cost(Destinations[dest]),acyclic,grandparent_check);
        }
        main_dijkstra_alg->printAGFile();
        main_dijkstra_alg->populateAGnodes();

        if (strategy == 0)
        {
            calculate_min_path_strategy_AG();
        }
        else
        {
            //If we are doing a Budget-based strategy then the AG needs to be expanded
            //by adding all possible costs up to the budget
            //main_dijkstra_alg->expand_AG(Budget);
            calculate_Gibbs_strategy_AG(Beta, Budget);
        }
    }
    else if(stocastic_go_observer)
    {
        std::cout << "stocastic_go_observer" << endl;
        //Calculate furthest destination to graduate experiments
        float furthest_destination=0;
        for (size_t dest = 0; dest < Destinations.size(); dest++){
            furthest_destination=max(furthest_destination,float(main_dijkstra_alg->get_cost(Destinations[dest])));
        }
        //If C is provided, we use it to determine the budget.
        if(C>0){
            Budget=float(C)*furthest_destination+1;
            std::cout<<"C-based Budget:,"<<Budget<<endl;
        }
        //First calculate Lambda
        for (size_t dest = 0; dest < Destinations.size(); dest++)
        {
            Dijkstra_algs_dest_to_dest[dest].determine_all_shortest_paths(current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), -1);
            //Also populate barX for q calculations
            barX.push_back(0.5*main_dijkstra_alg->get_cost(Destinations[dest]));
            std::cout<<"barX["<<dest<<"]:"<<barX<<endl;
            for (auto node : nodes_set)
            {
                if (!main_dijkstra_alg->is_node_reachable(node))
                {
                    std::cout<<"node "<<node<<" is not reachable from destination"<<endl;
                    continue;
                }
                else
                {
                    //std::cout << "\treachable node:" << node << ",cost to destination[" << Destinations[dest] << "]:" << Dijkstra_algs_dest_to_dest[dest].get_cost(node);
                    //std::cout << ",distance to start:," << main_dijkstra_alg->get_cost(node) << endl;
                }
            }
        }

        //Dijkstra_algs[0].print_perim_paths();

        //NOW CREATE AUGMENTED GRAPH
        //GIVEN BUDGET C
        //ORIGIN Note IS (Orig,0)
        main_dijkstra_alg->origin = main_graph.get_vertex(start);
        main_dijkstra_alg->set_Destinations_Order(&Destinations);
        for (auto dest : Destinations)
        {
            main_dijkstra_alg->Destinations.insert(dest);
        }
        for (size_t dest = 0; dest < Destinations.size(); dest++)
        {
            std::cout<<"start:"<<start<<",destination:"<<Destinations[dest]<<endl;
            std::cout<<"grandparent_check="<<grandparent_check<<",acyclic_check:"<<acyclic<<endl;
            main_dijkstra_alg->createAG(main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]), Budget*main_dijkstra_alg->get_cost(Destinations[dest]),acyclic,grandparent_check);
        }
        main_dijkstra_alg->printAGFile();
        cout << "hola AGP" << endl;
        main_dijkstra_alg->populateAGnodes();

        std::cout << "before calculate_all_subpaths,time:" << ctime(&timenow) << endl;
        calculate_all_subpaths();
        std::cout << "after calculate_all_subpaths,time:" << ctime(&timenow) << flush<<endl;
        
        if (strategy == 0)
        {
            calculate_min_path_strategy_AG();
            //calculate_min_path_strategy_prefixes();
        }
        else
        {
            //If we are doing a Budget-based strategy then the AG needs to be expanded
            //by adding all possible costs up to the budget
            //main_dijkstra_alg->expand_AG(Budget);
            calculate_Gibbs_strategy_AG(Beta, Budget);
        }
        //std::cout<<"calling calculate_observer_lambda"<<flush<<endl;

        //For now the lambda(d) distribution is all destinations are equally probable
        vector<float> dest_dist;
        dest_dist.assign(Destinations.size(),1.0/float(Destinations.size()));
        //initialize lammbda_obs
        map<pair<int,float>,float > temp_map;
        lambda_obs.assign(Destinations.size(),temp_map);
        for(size_t dest=0;dest<Destinations.size();dest++){
            lambda_obs[dest][make_pair(start,0)]=dest_dist[dest];
        }
        auto timestart = std::chrono::high_resolution_clock::now();
        std::cout<<"calling calculate_observer_lambda"<<flush<<endl;
        calculate_observer_lambda(start,0);
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end_time - timestart;
        auto currentTime = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
        std::cout << "calculate_observer_lambda_runtime:,"<< currentTime <<endl;
        //For zeta, we recursively iterate backwards from each destination node
        std::cout<<"calling calculate_observer_zeta"<<flush<<endl;
        calculate_observer_zeta();
        //calculate target's real path
        std::cout<<"calling simulation_observer"<<flush<<endl;
        std::cout<<"target_destination:"<<target_destination<<endl;
        //timings
        end_time = std::chrono::high_resolution_clock::now();
        diff = end_time - orig_start_time;
        auto overall_time = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
        std::cout<<"overall_runtime:,"<<overall_time<<endl;

        auto temp_start_time = chrono::high_resolution_clock::now();
        exit(0);
        //Now calculate lambda for all states in the AG
    }
}
void calculate_observer_lambda(int current_node,float current_cost){
    if(all_destinations.count(current_node)>0){
        return;//End of recursive forward search
    }
    //int current_cost=0;
    auto AG_Edges = main_dijkstra_alg->getSGEdges();
    
    vector<float> lambda_initial;

    auto current_vertex = main_graph.get_vertex(current_node);

    //Numerator
    map<pair<int,int> ,float> numerator;
    map<int,float> denominator;
    for(size_t dest=0;dest<Destinations.size();dest++){
        for(auto dest_node : (*AG_Edges)[dest][current_node]){
            auto dest_vertex = main_graph.get_vertex(dest_node);
            float cost_child=current_cost + main_graph.get_edge_weight(current_vertex, dest_vertex);
            auto key = make_pair(make_pair(current_node, current_cost), make_pair(dest_node, cost_child));
            if(Alpha[dest].count(key)<1){
                continue;//skip 0 prob edges
            }
            else if(Alpha[dest][key]<epsilon){
                continue;//or quasi-0
            }
            //std::cout<<"dest_node:["<<Destinations[dest]<<"]["<<dest_node<<","<<cost_child<<"]"<<endl;
            denominator[dest_node]+=Alpha[dest][key]*lambda_obs[dest][make_pair(current_node,current_cost)];
            //numerator[make_pair(dest,dest_node)]=(Alpha[dest][key]*(*dest_dist)[dest]);
            numerator[make_pair(dest,dest_node)]=Alpha[dest][key]*lambda_obs[dest][make_pair(current_node,current_cost)];
        }
    }
    //Second passover once denominator is calculated
    Active_AG.insert(make_pair(start,0));
    for(size_t dest=0;dest<Destinations.size();dest++){
        for(auto dest_node : (*AG_Edges)[dest][current_node]){
            auto dest_vertex = main_graph.get_vertex(dest_node);
            float cost_child=current_cost + main_graph.get_edge_weight(current_vertex, dest_vertex);
            auto key = make_pair(make_pair(current_node, current_cost), make_pair(dest_node, cost_child));
            if(Alpha[dest].count(key)<1){
                continue;//skip 0 prob edges
            }
            else if(Alpha[dest][key]<epsilon){
                continue;//or quasi-0
            }
            //std::cout<<"dest_node:["<<Destinations[dest]<<"]["<<dest_node<<","<<cost_child<<"]"<<",";
            //std::cout<<"numerator:,"<<numerator[make_pair(dest,dest_node)]<<",denominator:,"<<denominator[dest_node]<<endl;
            lambda_obs[dest][make_pair(dest_node,cost_child)]=numerator[make_pair(dest,dest_node)]/denominator[dest_node];
            C_LambdaAbs[key]=denominator[dest_node];
            //std::cout<<"lambda,\u03BB["<<Destinations[dest]<<"]"<<"["<<dest_node<<","<<cost_child<<"]:,"<<lambda_obs[dest][make_pair(dest_node,cost_child)]<<endl;
            //std::cout<<"Clambda["<<current_node<<","<<current_cost<<"]["<<dest_node<<","<<cost_child<<"]="<<C_LambdaAbs[key]<<endl;
            Active_AG.insert(make_pair(dest_node,cost_child));
            //std::cout<<"Active_AG["<<dest_node<<","<<cost_child<<"]"<<endl;
            calculate_observer_lambda(dest_node,cost_child);
        }
    }
}
void calculate_observer_zeta(){
    auto AG_Edges = main_dijkstra_alg->getAGEdges();
    auto BW_AG_Edges = main_dijkstra_alg->getAGEdgesBw();
    auto destination_nodes=main_dijkstra_alg->getAGDestinationNodes();
    stack<pair<int,float> > pending;
    //std::cout<<"\tCalculate_observer_zeta"<<endl;
    //ForEach s in destinations add zeta value of 0 and add destinations to queue
    for(auto dest_node : *destination_nodes){
        auto s=make_pair(dest_node.first->getID(),dest_node.second);
        if(Active_AG.count(s)<1){
            continue;//Skipping destinations not reached by current strategy
        }
        zeta_obs[s]=0;
        pending.push(s);
        //std::cout<<"\tAdding destination ["<<s.first<<","<<s.second<<"] to stack"<<endl;
    }
    while(!pending.empty()){
        auto current_node=pending.top();
        //std::cout<<"\tcurrent_node:"<<current_node.first<<","<<current_node.second<<endl;
        //If z_s is known, then pop it from pending stack
        //And add AG parents as long as there is a non-zero lambda_obs value
        if(zeta_obs.find(current_node)!=zeta_obs.end()){
            //std::cout<<"\tCase 1, zeta is known for node:"<<current_node.first<<","<<current_node.second<<",adding children nodes:"<<endl;
            //std::cout<<"Removing1 node:["<<pending.top().first<<","<<pending.top().second<<"]"<<endl;
            pending.pop();
            //Add parent nodes to stack
            for (auto parent_node : (*BW_AG_Edges)[current_node]){
                if(Active_AG.count(current_node)<1){
                    continue;//Skipping nodes not reached by current strategy
                }
                if(zeta_obs.find(parent_node)==zeta_obs.end()){//parent node has no zeta so add to stack
                    pending.push(parent_node);
                    //std::cout<<"\t\tCase 2, zeta is missing for parent node:"<<parent_node.first<<","<<parent_node.second<<endl;
                    continue;
                }
            }
            //And go to next node in stack
            continue;
        }

        //So we need to calculate z(s)
        //Check if we know all zetas for all children(s)
        bool ready=true;
        for (auto dest_node : (*AG_Edges)[current_node.first]){
            auto current_vertex = main_graph.get_vertex(current_node.first);
            //std::cout<<"current_vertex:"<<current_vertex->getID()<<flush<<endl;
            auto dest_vertex = main_graph.get_vertex(dest_node);
            //std::cout<<"dest_vertex:"<<dest_vertex->getID()<<flush<<endl;
            float cost_child=current_node.second + main_graph.get_edge_weight(current_vertex, dest_vertex);
            if(Active_AG.count(make_pair(dest_node,cost_child))<1){
                //std::cout<<"node["<<dest_node<<","<<cost_child<<"] is missing from current AG, skipping"<<endl;
                continue;//Skipping nodes not reached by current strategy
            }
            //auto key = make_pair(make_pair(current_node, current_cost), make_pair(dest_node, cost_child));
            auto ag_dest_node=make_pair(dest_node,cost_child);
            if(zeta_obs.find(ag_dest_node)==zeta_obs.end()){//child node has no zeta so add to queue
                //std::cout<<"\t\tadding child node:["<<ag_dest_node.first<<","<<ag_dest_node.second<<"]"<<endl;
                pending.push(ag_dest_node);
                ready=false;//No break, we evaluate all children nodes so we can add any which have priority
            }
        }
        if(ready){
            //std::cout<<"Calculate Zeta,all children zetas available"<<endl;
            zeta_obs[current_node]=calculate_zeta(current_node);
            //std::cout<<"zeta["<<current_node.first<<","<<current_node.second<<"]:,"<<zeta_obs[current_node]<<endl;
            //std::cout<<"Removing2 node:["<<pending.top().first<<","<<pending.top().second<<"]"<<endl;
            pending.pop();
            //Add parent nodes with pending zeta calculation
            for (auto parent_node : (*BW_AG_Edges)[current_node]){
                //std::cout<<"\tevaluating for adding parent_node:["<<parent_node.first<<","<<parent_node.second<<"]"<<endl;
                if(Active_AG.count(parent_node)<1){
                    //std::cout<<"\tskipping parent_node:["<<parent_node.first<<","<<parent_node.second<<"]"<<endl;
                    continue;//Skipping nodes not reached by current strategy
                }
                if(zeta_obs.find(parent_node)==zeta_obs.end()){//parent node has no zeta so add to stack
                    pending.push(parent_node);
                    //std::cout<<"\t\tCase 2a, zeta is missing for parent node:"<<parent_node.first<<","<<parent_node.second<<endl;
                    continue;
                }

            }
            }
    }
    //Record All Relevant Data
    ofstream myfile;
    string filename = "OBS_V" + to_string(main_graph.getVertexNum()) + "_S" + to_string(random_seed) + ".txt";
    myfile.open(filename);
    
    ofstream myfile2;
    string filename2 = "../experiments/OBS_exp_V," + to_string(main_graph.getVertexNum()) + ",_S," + to_string(random_seed) + ",_Bu," + to_string(Budget) + ",_Be," + to_string(Beta) + ".txt";
    myfile2.open(filename2);
    auto key=make_pair(start,0.0);
    myfile2<<"zeta,\u03B6[,"<<key.first<<","<<key.second<<",]=,"<<zeta_obs[key]<<endl;
    for(auto node : zeta_obs){
        if(debug)
            myfile<<"zeta,\u03B6["<<node.first.first<<","<<node.first.second<<"]="<<calculate_term1(node.first)<<","<<calculate_term2(node.first)<<endl;
    }
    for(size_t dest=0;dest<Destinations.size();dest++){
        for(auto node : lambda_obs[dest]){
            if(debug)
                myfile<<"lambda,\u03BB|^{["<<node.first.first<<","<<node.first.second<<"]}"<<"("<<Destinations[dest]<<")="<<node.second<<endl;
        }
    }
    for(size_t dest=0;dest<Destinations.size();dest++){
        for(auto node : nodes_set){
            if(node==start){
                if(debug)
                    myfile<<"q(\u03B4("<<node<<","<<Destinations[dest]<<")="<<main_dijkstra_alg->get_cost(Destinations[dest])<<endl;
            }
            else{
                //NOTE:If we use the new q function, we need to print it out for each
                //AG node
                for(auto node : Active_AG){
                //myfile<<"q(\u03B4("<<node<<","<<Destinations[dest]<<")="<<Dijkstra_algs_dest_to_dest[dest].get_cost(node)<<endl;
                    if(debug)
                        myfile<<"q(\u03B4("<<node.first<<","<<node.second<<","<<Destinations[dest]<<")="<<calculate_q_avg(node,dest)<<endl;
                }
            }
        }
    }
    if(debug){
        for(auto node : C_LambdaAbs){
            myfile<<"Lambda,\u039B_{["<<node.first.first.first<<","<<node.first.first.second<<"]["<<node.first.second.first<<","<<node.first.second.second<<"]}="<<node.second<<endl;
        }
        for(auto node : X_Hat){
            myfile<<"X_Hat,X̂["<<node.first<<","<<node.second<<"]"<<endl;
        }
    }
    //Xcrit is calculated as 
    auto paths=main_dijkstra_alg->getAGpaths();
    for(auto path : *paths){
        float cost=0;
        if(path.size()==1){
            continue;
        }
        for(size_t i=1;i<path.size();i++){
            cost+=main_graph.get_edge_weight(path[i-1],path[i]);
            auto key=make_pair(path[i]->getID(),cost);
            if(X_Hat.count(key)>0){
                X_Crit.insert(key);
                break;
            }
        }
    }
    
    for(auto node : X_Crit){
        //Also calculate estimated AG node
        pair<int,float> chosen_AG_dest=make_pair(0,0.0);
        for (size_t dest = 0; dest < Destinations.size(); dest++){
            if(node.first==start){
                if(lambda_obs[dest][node]*calculate_q_avg(make_pair(start,0),dest)>=zeta_obs[node]){
                    X_Crit_dest[node]=Destinations[dest];
                    chosen_AG_dest=make_pair(Destinations[dest],node.second+main_dijkstra_alg->get_cost(Destinations[dest])>chosen_AG_dest.second);
                }
            }
            else{
                if(lambda_obs[dest][node]*calculate_q_avg(node,dest)>=zeta_obs[node]){
                    X_Crit_dest[node]=Destinations[dest];
                    chosen_AG_dest=make_pair(Destinations[dest],node.second+Dijkstra_algs_dest_to_dest[dest].get_cost(node.first));
                }
            }
        }
        if(debug)
            myfile<<"X_Crit,X̂^{crit}["<<node.first<<","<<node.second<<"],ChosenDest["<<chosen_AG_dest.first<<","<<chosen_AG_dest.second<<"]"<<endl;
        myfile2<<"X_Crit[,"<<node.first<<","<<node.second<<",],ChosenDest[,"<<chosen_AG_dest.first<<","<<chosen_AG_dest.second<<",]"<<endl;
        }

    myfile.close();
    myfile2.close();
}
float calculate_zeta(pair<int,float> s){
    auto AG_Edges = main_dijkstra_alg->getAGEdges();
    //First calculate inner maxD of eq. 14
    //std::cout<<"calculate_zeta for node["<<s.first<<","<<s.second<<endl;
    float maxD=0;
    for (size_t dest = 0; dest < Destinations.size(); dest++){
        if(lambda_obs[dest].find(s)==lambda_obs[dest].end()){
            continue;//Skipping nodes not in current AG-SubGraph
        }
        //std::cout<<"\tlambda_obs["<<dest<<"]["<<s.first<<","<<s.second<<"]:"<<lambda_obs[dest][s]<<endl;
        //std::cout<<"\tdelta["<<dest<<"]="<<Dijkstra_algs_dest_to_dest[dest].get_cost(s.first)<<endl;
        /*if(s.first==start){
            maxD=max(maxD,lambda_obs[dest][s]*float(main_dijkstra_alg->get_cost(Destinations[dest])));
        }
        else{
            maxD=max(maxD,lambda_obs[dest][s]*float(Dijkstra_algs_dest_to_dest[dest].get_cost(s.first)));
        }*/
        maxD=max(maxD,lambda_obs[dest][s]*calculate_q_avg(s,dest));
    }

    float zeta=maxD;
    //Second calculate the second term in eq. 14
    float sum=0;
    auto current_vertex = main_graph.get_vertex(s.first);
    for (auto dest_node : (*AG_Edges)[s.first]){
            auto dest_vertex = main_graph.get_vertex(dest_node);
            float cost_child=s.second + main_graph.get_edge_weight(current_vertex, dest_vertex);
            auto key=make_pair(s,make_pair(dest_node,cost_child));
            if(C_LambdaAbs.find(key)==C_LambdaAbs.end()){
                continue;//skip 0 probability edges
            }
            sum+=C_LambdaAbs[make_pair(s,make_pair(dest_node,cost_child))]*zeta_obs[make_pair(dest_node,cost_child)];
            //std::cout<<"s'["<<dest_node<<","<<cost_child<<"]"<<",C_Lambda_Abs="<<C_LambdaAbs[make_pair(s,make_pair(dest_node,cost_child))]<<",zeta_obs="<<zeta_obs[make_pair(dest_node,cost_child)]<<endl;
    }
    //X_Hat are those nodes where first term in eq 14 is dominant
    std::cout<<"node:["<<s.first<<","<<s.second<<"],maxD:"<<maxD<<",sum:"<<sum<<endl;
    if(maxD>=sum){
        X_Hat.insert(s);
    }
    //std::cout<<"term2:"<<sum<<endl;
    zeta=max(zeta,sum);
    return zeta;
}
float calculate_term1(pair<int,float> s){
    //Calculate term1 for eq. 14
    float maxD=0;
    for (size_t dest = 0; dest < Destinations.size(); dest++){
        if(lambda_obs[dest].find(s)==lambda_obs[dest].end()){
            continue;//Skipping nodes not in current AG-SubGraph
        }
        maxD=max(maxD,lambda_obs[dest][s]*calculate_q_avg(s,dest));
    }
    return maxD;
}
float calculate_term2(pair<int,float> s){
    auto AG_Edges = main_dijkstra_alg->getAGEdges();
    float sum=0;
    auto current_vertex = main_graph.get_vertex(s.first);
    for (auto dest_node : (*AG_Edges)[s.first]){
            auto dest_vertex = main_graph.get_vertex(dest_node);
            float cost_child=s.second + main_graph.get_edge_weight(current_vertex, dest_vertex);
            auto key=make_pair(s,make_pair(dest_node,cost_child));
            if(C_LambdaAbs.find(key)==C_LambdaAbs.end()){
                continue;//skip 0 probability edges
            }
            sum+=C_LambdaAbs[make_pair(s,make_pair(dest_node,cost_child))]*zeta_obs[make_pair(dest_node,cost_child)];
            //std::cout<<"s'["<<dest_node<<","<<cost_child<<"]"<<",C_Lambda_Abs="<<C_LambdaAbs[make_pair(s,make_pair(dest_node,cost_child))]<<",zeta_obs="<<zeta_obs[make_pair(dest_node,cost_child)]<<endl;
    }
    return sum;
}


//q(x) = x \times \frac{\bar{q}}{\bar{x}} if x \in [0, \bar{x}]
//q(x) = \bar{q} if x > \bar{x}
//where you can set \bar{q} = 0.9 and \bar{x} half of the distance between the origin and the destination.
float calculate_q(pair<int,float> node,int dest){
    //std::cout<<"calling calculate_q"<<endl;
    float x=0;
    if(node.first==start){
        x=main_dijkstra_alg->get_cost(start);
    }
    else{
        x=Dijkstra_algs_dest_to_dest[dest].get_cost(node.first);
    }

    if(x>barX[dest]){
        //std::cout<<"Far from dest("<<Destinations[dest]<<"),x="<<x<<",q["<<node.first<<","<<node.second<<"]=barQ="<<barQ<<endl;
        std::cout<<"final_q:"<<barQ<<endl;
        return barQ;
    }
    else{
        //std::cout<<"Close from dest("<<Destinations[dest]<<"),x="<<x<<",q["<<node.first<<","<<node.second<<"]=barQ="<<x*(barQ/barX[dest])<<endl;
        std::cout<<"final_q:"<<x*(barQ/barX[dest])<<endl;
        return x*(barQ/barX[dest]);
    }
}
float calculate_q_path(pair<int,float> node,int dest,float cost){
    //std::cout<<"calling calculate_q_path"<<endl;

    if(cost>barX[dest]){
        //std::cout<<"Far from dest("<<Destinations[dest]<<"),x="<<x<<",q["<<node.first<<","<<node.second<<"]=barQ="<<barQ<<endl;
        return barQ;
    }
    else{
        //std::cout<<"Close from dest("<<Destinations[dest]<<"),x="<<x<<",q["<<node.first<<","<<node.second<<"]=barQ="<<x*(barQ/barX[dest])<<endl;
        return cost*(barQ/barX[dest]);
    }
}
float calculate_q_avg(pair<int,float> node,int dest){//Eq 23 in paper,\bar{q}
    auto AG_Nodes = main_dijkstra_alg->getAGsubgraphs();
    auto key=make_pair(Destinations[dest],node);
    auto itr1 = subpaths.lower_bound(key);
    auto itr2 = subpaths.upper_bound(key);
    //std::cout<<"hola calculate_q_avg, subpaths from ["<<node.first<<","<<node.second<<"] to "<<"dest:"<<Destinations[dest]<<endl;
    float final_q=0;
    float total_prob=0;//debug check all AG nodes have combined prob=1.0
    //ForAll possible paths from s to d
    int counter=0;
    while (itr1 != itr2)    
    {
        if (itr1 -> first != key){
            break;
        }

        //std::cout<<"\tdest["<<counter++<<"]:"<<(itr1 -> first).first<<",s["<<(itr1 -> first).second.first<<",";
        //std::cout<<(itr1 -> first).second.second<<"]"<<endl;

        //loop through subpath to calculate partial q
        float path_prob=1.0;
        //recalculate barX, barQ is fixed at 0.9
        //vector<map<pair<pair<int, float>, pair<int, float>>, float>> Alpha;
        //barX.push_back(0.5*main_dijkstra_alg->get_cost(Destinations[dest]));
        //std::cout<<"\tpath_size:"<<itr1->second.size()<<endl;
        for(size_t i=0;i<itr1->second.size();i++){
            //std::cout<<"\t\tedge[:"<<itr1->second[i].first<<","<<itr1->second[i+1].first<<"]";
            auto prob1=Alpha[dest][make_pair(itr1->second[i],itr1->second[i+1])];
            path_prob=path_prob*prob1;
            //std::cout<<","<<path_prob<<endl;
        }
        //auto q1=calculate_q_path(itr1->second[0],dest,itr1->second.back().second-itr1->second.front().second);
        auto q1=calculate_q(itr1->second[0],dest);
        final_q+=q1*path_prob;
        //std::cout<<"\tpath_prob:,"<<path_prob<<",q1:,"<<q1<<endl;
        total_prob+=path_prob;
        itr1++;
        //std::cout << (itr1 -> first).first << "  "  << (itr1 -> first).second << "  "  << itr1 -> second << endl;        

    }
    //std::cout<<"total_prob:"<<total_prob<<",final_q:"<<final_q<<endl;
    if(final_q>1.0){
        cerr<<"Average_q:"<<final_q<<" cannot be bigger than 1, check code!"<<endl;
        exit(1);
    }
    else if(total_prob>1){
        cerr<<"total_prob:"<<total_prob<<" cannot be bigger than 1, check code!"<<endl;
        exit(1);

    }
    return final_q;

}
float calculate_q_recursive(pair<int,float> node,int dest){//Eq 25 in paper,\bar{q}
//q(s,d)=Sum_{s'\inX}Alpha^{d}_{ss'}{q(s',d)+W_{s,s'}}
//q((d,x),d)=0 for all xa
    //std::cout<<"hola1"<<flush<<endl;
    auto AG_SG_Edges = main_dijkstra_alg->getSGEdges();
    if(all_destinations.count(node.first)>0){//Destination node
        return 0;
    }
    float sum=0;
    auto orig_vertex = main_graph.get_vertex(node.first);
    for (auto child_node : (*AG_SG_Edges)[dest][node.first]){
        //std::cout<<"\tsum1="<<sum<<flush<<endl;
        auto dest_vertex = main_graph.get_vertex(child_node);
        float x = node.second + main_graph.get_edge_weight(orig_vertex, dest_vertex);
        auto key=make_pair(make_pair(node.first, node.second), make_pair(child_node, x));
        sum+=Alpha[dest][key]*(calculate_q_recursive(make_pair(child_node,x),dest)+ main_graph.get_edge_weight(orig_vertex, dest_vertex));
    }
    //std::cout<<"sum2="<<sum<<flush<<endl;
    return sum;
}

//We skip any code which does not need to be repeated when doing iterative fictitious play
void calculate_min_path_strategy_prefixes(bool skip=false){
    //When using only prefixes
    double lambdaDest = 1.0 / double(Destinations.size());
    //std::cout << "hola calculate_min_path_strategy_prefixes based on full paths" << endl;
    //Loop through all paths and count number of paths per detination
    //All destinations have equal likelyhood and all paths to same destination get an even dirstribution of that probability
    //First how many paths per destination
    paths = main_dijkstra_alg->getAGpaths2();
    vector<int> PathsPerDestination(Destinations.size(),0);
    map<int,int> PathsToDestination;
    vector<double> ProbPerDestinationPath;
    size_t counter = 0;
    if(!skip){
        for(auto d : Destinations){
            DestinationsOrder1[counter] = d;
            cout << "DestinationsOrder1["<<counter<<"]:" << d << endl;
            DestinationsOrder2[d] = counter++;
            cout << "DestinationsOrder2["<<d<<"]:" << counter-1 << endl;
        }
        cout << "DestinationsOrder1.size:" << DestinationsOrder1.size() << flush<<endl;
    }
    //Double loop to calculate prob per path
    counter = 0;
    if(!skip){
        PathsPerDestination.assign(Destinations.size(), 0);
    //Count paths per destination, eliminate any paths whose cost is suboptimal
        double PathCost = 0;
        set<int> PathsToEliminate;
        for(auto path : *paths){
            //cout << "\t working on path:" << counter++ <<",checking strategy compatibility"<< flush<<endl;
            //PathCost = 0;
            //for (size_t i=0;i<path.second.size()-1;i++){
            //    PathCost+= main_graph.get_edge_weight(path.second[i], path.second[i+1]);
            //}
            //double optimal_cost=main_dijkstra_alg->get_cost(path.second.back()->getID());
            //cout << "\t working on path:" << counter++ <<",SKIPPING strategy compatibility for DEBUGGING"<< endl;
            /*if(PathCost>optimal_cost){
                cout << "\t\teliminating path, path_cost:" << PathCost << ",optimalcost:" << optimal_cost << endl;
                PathsToEliminate.insert(counter - 1);
            continue;
        }*/
            PathsPerDestination[DestinationsOrder2[path.second.back()->getID()]] += 1;
        }
        //Now we fully eliminate from paths the sub-optimal paths
        //for (auto i : PathsToEliminate)
        //    paths->erase(i);
        cout << "PathsPerDestination.size:" << PathsPerDestination.size() << endl;
        //Now distribute probability equally per valid path per destination
        ProbPerDestinationPath.assign(Destinations.size(), 0.0);
        for(auto d : DestinationsOrder1){
            ProbPerDestinationPath[d.first] = lambdaDest / PathsPerDestination[d.first];
            cout << "ProbPerDestinationPath:" << ProbPerDestinationPath[d.first] << ",paths:" << PathsPerDestination[d.first] << ",lambdaDest:" << lambdaDest << flush << endl;
        }
        counter = 0;
        prob_path.assign(paths->size(), 0);
        dest_prob_path.assign(Destinations.size(), prob_path);
        for(auto path : *paths){
            auto destination=DestinationsOrder2[path.second.back()->getID()];
            if (PathsToEliminate.find(counter)!=PathsToEliminate.end()){
                prob_path[counter]=0.0;
                dest_prob_path[destination][counter] = 0;
            }
            else{
            prob_path[counter]=ProbPerDestinationPath[destination];
            dest_prob_path[destination][counter]=ProbPerDestinationPath[destination]/lambdaDest;
        }
            if(prob_path[counter]>0){
                //main_dijkstra_alg->print_path(counter);
                //cout << ",path_prob[" << counter << "]:" << prob_path[counter] <<", path_cost"<<main_dijkstra_alg->get_cost_to_dest(0,counter)<<",dest:,"<<destination<< endl;
                //cout << "dest_path_prob[" << destination<<"]["<< counter << "]:" <<dest_prob_path[destination][counter] << endl;
            }
            counter++;
        }

        float total_prob = 0;
        for (auto prob: prob_path){
            total_prob += prob;
        }
        cout <<"Aggregated path_prob:,"<< total_prob << endl;
    }
    //cout << "out of skip1" << endl;

    //Now calculate probabilty per prefix. Each prefix gets the probability of all possible paths it is part of.
    prefixes = main_dijkstra_alg->get_prefixes();
    vector<float> prob_per_prefix(prefixes->size(), 0.0);
    dest_prob_per_prefix.assign(Destinations.size(), prob_per_prefix);
    
    if(!skip){
        auto start_time = std::chrono::high_resolution_clock::now();
        string method("linkPrefixToPaths");
        paths_per_prefix=main_dijkstra_alg->linkPrefixToPaths();
        elapsed_time(method,start_time);
    }
    
    dest_paths_per_prefix = main_dijkstra_alg->get_dest_paths_per_prefix();

    start_time = std::chrono::high_resolution_clock::now();
    string method="dest_prob_per_prefix";
    counter=0;
    for (size_t i = 0; i < prefixes->size(); i++)
    {
        for(auto path_id : paths_per_prefix->at(i)){
            auto destination=DestinationsOrder2[paths->at(path_id).back()->getID()];
            prob_per_prefix[i] += prob_path[path_id];
            dest_prob_per_prefix[destination][i] += dest_prob_path[destination][path_id];
            //cout << "\t\tprefix_id:," << i << ",path_id:," << path_id << ",prob_per_prefix[]"<<prob_per_prefix[i]<<endl;
        }
        //cout << "\tprefix_id:,";
        //main_dijkstra_alg->print_prefix(i);cout <<"|"<< i << ",prob_per_prefix[]" << prob_per_prefix[i] << endl;
        //for (size_t d = 0; d < Destinations.size();d++)
        //{
        //    cout << "\t\tdestination:," << d << ",dest_prob_per_prefix:"<<dest_prob_per_prefix[d][i]<<endl;
        //}
    }
    elapsed_time(method, start_time);
    
    start_time = std::chrono::high_resolution_clock::now();
    method="prefix_prob_path";
    //Now calculate prob_per_path_given_prefix.
    prefix_prob_path.clear();//Make sure all previous probs are erased
    for (size_t i = 0; i < prefixes->size(); i++){
        for(auto path_id : paths_per_prefix->at(i)){
            if(prob_per_prefix[i]>0){
                auto id = make_pair(i, path_id);
                prefix_prob_path[id] = prob_path[path_id] / prob_per_prefix[i];
                //if(prefix_prob_path[id]>0){
                    //cout << "\t Iter:" << current_iter << ",prefix_prob_path[(," << i << "," << path_id << ",)]:," << prefix_prob_path[id] << ",pref:";
                    //main_dijkstra_alg->print_prefix(i);
                    //cout<< endl;
                //}
            }
        }
    }
    elapsed_time(method, start_time);

    //Now calculate prob_per_path_given_dest_and_prefix.
    
    vector<map<pair<unsigned,unsigned>,float> > prob_per_path_given_dest_and_prefix;
    map<pair<unsigned, unsigned>, float> temp;
    prob_per_path_given_dest_and_prefix.assign(Destinations.size(), temp);

    start_time = std::chrono::high_resolution_clock::now();
    method="prob_per_path_given_dest";
    for(size_t d=0;d<Destinations.size();d++){
        for(size_t pref=0;pref<prefixes->size();pref++){
            for(auto pth : paths_per_prefix->at(pref)){
                auto id = make_pair(pref, pth);
                if(dest_prob_per_prefix[d][pref]>0){
                    prob_per_path_given_dest_and_prefix[d][id]=dest_prob_path[d][pth]/dest_prob_per_prefix[d][pref];
                    //cout << "prob_per_path_given_dest_and_prefix[" << d << "][" << pref << "][" << pth << "]:," << prob_per_path_given_dest_and_prefix[d][id] << endl;
                }
            }
        }
    }
    elapsed_time(method, start_time);
    //Now we calculate lambda as a function of prefix and destination
    //For now the lambda(d) distribution is all destinations are equally probable
    //We also calculate q(d,prefix) since we need to use the same loop
    vector<float> dest_dist;
    dest_dist.assign(Destinations.size(),1.0/float(Destinations.size()));
    map<unsigned,float > temp_map;
    lambda_obs_dest_given_prefix.assign(Destinations.size(),temp_map);
    
    vector<float> temp_vec_float; temp_vec_float.assign(prefixes->size(), 0);
    q_given_dest_and_prefix.assign(Destinations.size(), temp_vec_float);
    
    start_time = std::chrono::high_resolution_clock::now();
    method="q_given_dest_and_prefix_time";
    for(size_t dest=0;dest<Destinations.size();dest++){
        for(size_t pref=0;pref<prefixes->size();pref++){
            //cout << "\tcalculating lambda and q per dest and prefix(" <<dest<<","<<pref<<")"<< flush << endl;
            auto id = make_pair(dest, pref);
            auto it = dest_paths_per_prefix->find(id);
            if(it==dest_paths_per_prefix->end())
                continue;
            for(auto pth : it->second){
                //cout << "pth:" << pth << flush << endl;
                auto id = make_pair(pref, pth);
                //skip quasi-0 probs to improve performance
                if(prob_per_path_given_dest_and_prefix[dest][id]<0.00001)
                    continue;
                lambda_obs_dest_given_prefix[dest][pref] += prefix_prob_path[id];
                q_given_dest_and_prefix[dest][pref] += prob_per_path_given_dest_and_prefix[dest][id]*(modified_q(main_dijkstra_alg->get_cost_to_dest(pref,pth)));
                //q_given_dest_and_prefix[dest][pref] = modified_q(q_given_dest_and_prefix[dest][pref]);
                //if(pref==0||pref==1||pref==175){
                //    cout << "\t\tIter:," << current_iter << ",pref:," << pref << ",d:," << dest << ",prob:,";
                //    cout<<prob_per_path_given_dest_and_prefix[dest][id] << ",cost:," << main_dijkstra_alg->get_cost_to_dest(pref, pth) << ",partial q_given_dest_and_prefix:," << prob_per_path_given_dest_and_prefix[dest][id] * main_dijkstra_alg->get_cost_to_dest(pref, pth) << endl;
                //}
            }
           //if(saturation>0){
           //     q_given_dest_and_prefix[dest][pref] = min((float) 1.0, q_given_dest_and_prefix[dest][pref]);
           //}
            //cout << "\t lambda_obs_dest_given_prefix[,"<<dest<<",][,"<<pref<<",]:," << lambda_obs_dest_given_prefix[dest][pref] << flush<<endl;
            
            //if(pref==0||pref==1||pref==175){
            //    cout << "\t q_given_dest_and_prefix[,"<<dest<<",][,"<<pref<<",]:," << q_given_dest_and_prefix[dest][pref] << flush<<endl;
            //}
        }
    }
    elapsed_time(method, start_time);
    //Now we can calculate avg_q per prefix as the max q for all dests * lambda(d,prefix).
    avg_q_prefix.assign(prefixes->size(), 0.0);
    phi_prefix.assign(prefixes->size(), 0.0);//Here is a good point to resize it as well.
    for(size_t pref=0;pref<prefixes->size();pref++){
        for(size_t dest=0;dest<Destinations.size();dest++){
            //avg_q_prefix[pref] = max((float) 2.0,max(avg_q_prefix[pref], lambda_obs_dest_given_prefix[dest][pref]*q_given_dest_and_prefix[dest][pref]));
            avg_q_prefix[pref] = max(avg_q_prefix[pref], lambda_obs_dest_given_prefix[dest][pref]*q_given_dest_and_prefix[dest][pref]);
            //if(saturation>0){
            //    avg_q_prefix[pref]=min((float)1.0,avg_q_prefix[pref]);
            //}
        }
        //if(pref==0||pref==1||pref==175){
        //    cout << "\t avg_q_prefix[,"<<pref<<",]:," << avg_q_prefix[pref] <<endl;
        //}
    }

    //Now we can calculate q_avg as the max estimated q reward for any destination for a given prefix
    

    //Now we calculate AlphaPrefix which is the transition probability between 2 adjacent prefixes
    start_time = std::chrono::high_resolution_clock::now();
    method = "AlphaPrefix";
    edges = main_dijkstra_alg->getAGPEdges();
    unsigned parent = 0;
    for(auto children : *edges){
        if(prob_per_prefix[parent]==0.0){
            parent++;
            continue;
        }
        for(auto child : children){
            auto edge = make_pair(parent, child);
            AlphaPrefix[edge] = prob_per_prefix[child] / prob_per_prefix[parent];
            //cout << "\tworking on parent prefix:";
            //main_dijkstra_alg->print_prefix(parent);
            //cout << "\tto child prefix:";
            //main_dijkstra_alg->print_prefix(child);
            //cout << "|||";
            //cout << "\tAlphaPrefix[" << edge.first<<","<<edge.second << "]:" << prob_per_prefix[child] << "/" << prob_per_prefix[parent] << "=" << AlphaPrefix[edge] << endl;
        }
        parent++;

    }
    elapsed_time(method, start_time);
}
void calculate_min_path_strategy_AG()
{
    std::cout << "hola calculate_min_path_strategy" << endl;
    //Initialize Alpha
    Alpha.clear();
    map<pair<pair<int, float>, pair<int, float>>, float> temp_map;
    Alpha.assign(Destinations.size(), temp_map);
    //auto orig_node=main_graph.get_vertex(start);
    ofstream myfile;
    string filename = "AG_V" + to_string(main_graph.getVertexNum()) + "_S" + to_string(random_seed) + ".txt";
    myfile.open(filename);

    auto AG_Nodes = main_dijkstra_alg->getAGsubgraphs();
    auto AG_SG_Edges = main_dijkstra_alg->getSGEdges();
    std::cout << "AG_Edges.size:" << AG_SG_Edges->size() << endl;

    for (size_t dest = 0; dest < Destinations.size(); dest++)
    {
        //std::cout << "\tworking on dest:," << Destinations[dest] << ",optimal distance:" << main_dijkstra_alg->get_cost(Destinations[dest]) << endl;
        /*if (Budget < main_dijkstra_alg->get_cost(Destinations[dest]))
        {
            //cerr << "Budget is smaller than optimal distance from origin to destination!!!" << endl;
            exit(1);
        }*/
        for (auto node : (*AG_Nodes)[dest])
        {
            double OptPaths = 0;
            auto orig_vertex = main_graph.get_vertex(node.first->getID());
            //if(all_destinations.count(node.first->getID())>0){
            //	continue;//skip destinations as origin nodes
            //}
            //std::cout << "\t\torig_vertex1:" << orig_vertex->getID();
            //std::cout << ",Edges size:" << (*AG_SG_Edges)[dest][orig_vertex->getID()].size() << endl;

            for (auto child_node : (*AG_SG_Edges)[dest][orig_vertex->getID()])
            {
                //std::cout << "\t\t\tchild_node:," << child_node << endl;
                auto dest_vertex = main_graph.get_vertex(child_node);
                //std::cout << "\t\tcost=" << node.second << "+" << main_graph.get_edge_weight(orig_vertex, dest_vertex);
                //std::cout << "+" << Dijkstra_algs_dest_to_dest[dest].get_cost(child_node) << endl;
                double x = node.second + main_graph.get_edge_weight(orig_vertex, dest_vertex);
                double cost = x + Dijkstra_algs_dest_to_dest[dest].get_cost(child_node);
                //Probability is 0 if path is not optimal
                //std::cout << "\t\t\tcost:" << cost << ",optimal_cost:" << main_dijkstra_alg->get_cost(Destinations[dest]) << endl;
                if (cost > main_dijkstra_alg->get_cost(Destinations[dest]))
                {
                    //std::cout << "\t\t\tpath is not optimal" << endl;
                    //Alpha[dest][make_pair(make_pair(node.first->getID(),node.second),
                    //make_pair(child_node,x))]=0;
                }
                else
                {
                    if((*AG_Nodes)[dest].count(make_pair(dest_vertex,x))==0){
                        if(debug){
                            std::cout<<"AGnode["<<child_node<<","<<x<<"] not in any valid path"<<endl;
                        }
                        continue;
                    }
                    else{
                        std::cout<<"AGnode["<<child_node<<","<<x<<"] belongs to valid path for destination:"<<Destinations[dest]<<endl;
                    }
                    //std::cout << "\t\t\tpath is optimal" << endl;
                    Alpha[dest][make_pair(make_pair(node.first->getID(), node.second),
                                          make_pair(child_node, x))] = 1.0;
                    OptPaths++;
                }
            }

            //Second pass to rationalize probabilities
            for (auto child_node : (*AG_SG_Edges)[dest][orig_vertex->getID()])
            {
                auto dest_vertex = main_graph.get_vertex(child_node);
                float x = node.second + main_graph.get_edge_weight(orig_vertex, dest_vertex);
                auto key=make_pair(make_pair(node.first->getID(), node.second),make_pair(child_node, x));
                if(Alpha[dest].count(key)<1){
                    continue;//skip 0 prob edges
                }

                Alpha[dest][key] = 1.0 / OptPaths;
                //std::cout << "Alpha[" << Destinations[dest] << ",(" << node.first->getID() << "," << node.second << ")->(" << child_node;
                //std::cout << "," << x << ")]:" << Alpha[dest][make_pair(make_pair(node.first->getID(), node.second), make_pair(child_node, x))] << endl;
                //std::cout << "," << x << ")]:" << Alpha[dest][key] << endl;
                //myfile << Destinations[dest] << "," << node.first->getID() << "," << node.second << "," << child_node;
                //myfile << "," << x << "," << Alpha[dest][make_pair(make_pair(node.first->getID(), node.second), make_pair(child_node, x))] << endl;
                myfile << "," << x << "," << Alpha[dest][key]<<endl;
            }
        }
    }
    myfile.close();
}
void calculate_Gibbs_strategy_AG(double beta = 0.5, double budget = 100)
{
    std::cout << "hola calculate_Gibbs_strategy,beta:" << beta << ",Budget:" << budget << endl;
    //Initialize Alpha
    Alpha.clear();
    map<pair<pair<int, float>, pair<int, float>>, float> temp_map;
    Alpha.assign(Destinations.size(), temp_map);
    //auto orig_node=main_graph.get_vertex(start);
    ofstream myfile;
    //string filename = "RNG_N" + to_string(main_graph.getVertexNum()) + "_S" + to_string(random_seed) + ".txt";
    string filename = "AG_V" + to_string(main_graph.getVertexNum()) + "_S" + to_string(random_seed) + ".txt";
    myfile.open(filename);

    auto AG_Nodes = main_dijkstra_alg->getAGsubgraphs();
    auto AG_SG_Edges = main_dijkstra_alg->getSGEdges();
    std::cout << "AG_Edges.size:" << AG_SG_Edges->size() << flush << endl;

    for (size_t dest = 0; dest < Destinations.size(); dest++)
    {
        std::cout << "\tGibbs,working on dest:," << Destinations[dest] << ",optimal distance:" << main_dijkstra_alg->get_cost(Destinations[dest]) << endl;
        double max_cost=budget*main_dijkstra_alg->get_cost(Destinations[dest]);
        /*if (Budget < main_dijkstra_alg->get_cost(Destinations[dest]))
        {
            cerr << "Budget is smaller than optimal distance from origin to destination!!!" << endl;
            exit(1);
        }*/
        for (auto node : (*AG_Nodes)[dest])
        {
            double NormConst = 0;
            auto orig_vertex = main_graph.get_vertex(node.first->getID());
            //std::cout << "\t\torig_vertex:" << orig_vertex->getID();
            //std::cout << ",Edges size:" << (*AG_SG_Edges)[dest][orig_vertex->getID()].size() << endl;
            if(debug)
                std::cout<<"working on Dest:"<<Destinations[dest]<<",parent Alpha node:["<<node.first->getID()<<","<<node.second<<"]"<<endl;

            for (auto child_node : (*AG_SG_Edges)[dest][orig_vertex->getID()])
            {
                //std::cout << "\t\t\tchild_node:," << child_node << endl;
                auto dest_vertex = main_graph.get_vertex(child_node);
                //std::cout << "\t\tcost=" << node.second << "+" << main_graph.get_edge_weight(orig_vertex, dest_vertex);
                //std::cout << "+" << Dijkstra_algs_dest_to_dest[dest].get_cost(child_node) << endl;
                double x = node.second + main_graph.get_edge_weight(orig_vertex, dest_vertex);
                double cost = x + Dijkstra_algs_dest_to_dest[dest].get_cost(child_node);
                //Probability is 0 if path is not feasible
                //std::cout << "\t\t\tcost:" << cost << ",optimal_cost:" << main_dijkstra_alg->get_cost(Destinations[dest]) << endl;
                if (cost >= max_cost)
                {
                    //std::cout << "\t\t\tpath is not feasible for budget:" << budget << endl;
                    //Alpha[dest][make_pair(make_pair(node.first->getID(),node.second), make_pair(child_node,x))]=0;
                    if(debug)
                        std::cout<<"unfeasible,Dest:"<<Destinations[dest]<<",skipping Alpha node:["<<node.first->getID()<<","<<node.second<<"]->["<<child_node<<","<<x<<"]"<<endl;
                }
                else
                {
                    if((*AG_Nodes)[dest].count(make_pair(dest_vertex,x))==0){
                        if(debug){
                            std::cout<<"AGnode["<<child_node<<","<<x<<"] not in any valid path"<<endl;
                        }
                        continue;
                    }
                    else if(debug){
                        std::cout<<"AGnode["<<child_node<<","<<x<<"] belongs to valid path for destination:"<<Destinations[dest]<<endl;
                    }
                    //std::cout << "\t\t\tpath is feasible" << endl;
                    //std::cout<<"Dest:"<<Destinations[dest]<<",feasible Alpha node:["<<node.first->getID()<<","<<node.second<<"]->["<<child_node<<","<<x<<"]"<<endl;
                    auto key = make_pair(make_pair(node.first->getID(), node.second),
                                         make_pair(child_node, x));
                    double exponent = main_graph.get_edge_weight(orig_vertex, dest_vertex) +
                                      Dijkstra_algs_dest_to_dest[dest].get_cost(child_node);
                    Alpha[dest][key] = exp(-beta * (exponent/1000.0));
                    //std::cout<<"Dest:"<<Destinations[dest]<<",feasible Alpha node:["<<node.first->getID()<<","<<node.second<<"]->["<<child_node<<","<<x<<"]"<<",exponent:"<<exponent<<endl;
                    if (Alpha[dest][key] > epsilon)
                        NormConst += Alpha[dest][key];
                }
            }

            //Second pass to rationalize probabilities
            for (auto child_node : (*AG_SG_Edges)[dest][orig_vertex->getID()])
            {
                auto dest_vertex = main_graph.get_vertex(child_node);
                double x = node.second + main_graph.get_edge_weight(orig_vertex, dest_vertex);
                auto key = make_pair(make_pair(node.first->getID(), node.second),
                                     make_pair(child_node, x));
                //std::cout<<"Dest:"<<Destinations[dest]<<",feasible Alpha node:["<<node.first->getID()<<","<<node.second<<"]->["<<child_node<<","<<x<<"]"<<",exp:"<<Alpha[dest][key]<<endl;

                if (Alpha[dest][key] / NormConst > epsilon)
                {
                    Alpha[dest][key] = Alpha[dest][key] / NormConst;
                    if(debug){
                        std::cout << "Alpha[" << Destinations[dest] << ",(" << node.first->getID() << "," << node.second << ")->(" << child_node;
                        std::cout << "," << x << ")]:" << Alpha[dest][make_pair(make_pair(node.first->getID(), node.second), make_pair(child_node, x))] << ",NormConst:"<<NormConst<<endl;
                    }
                    if(debug){
                        myfile <<Destinations[dest] << "," << node.first->getID() << "," << node.second << "," << child_node;
                        myfile << "," << x << "," << Alpha[dest][make_pair(make_pair(node.first->getID(), node.second), make_pair(child_node, x))] << endl;
                    }
                }
                else{
                    //std::cout<<"Erasing Alpha["<<dest<<"] key:["<<key.first.first<<","<<key.first.second<<"]["<<key.second.first<<","<<key.second.second<<"]"<<endl;
                    Alpha[dest].erase(key);
                }
            }
        }
        std::cout << "\tGibbs,Finished with dest:," << Destinations[dest] << ",optimal distance:" << main_dijkstra_alg->get_cost(Destinations[dest]) << endl;
    }
    myfile.close();
}

void simulation()
{
    std::cout<<"hello simulation"<<flush<<endl;
    auto AG_Edges = main_dijkstra_alg->getSGEdges();
    //1st choose Destination
    srand(random_seed);
    int dest = rand() % Destinations.size();
    target_destination=Destinations[dest];
    std::cout << endl << "Simulation, randomly selected destination:," << target_destination << ","<<"dest:,"<<dest<<endl;
    std::cout << "Simulated path:" << flush<<start;
    int current_node = start;
    float r;
    int counter = 0;
    float cost1 = 0;
    //Create output file for visualization
    ofstream myfile;
    string filename = "path_" + to_string(main_graph.getVertexNum()) + "_S" + to_string(random_seed) + ".txt";
    myfile.open(filename);
    myfile << start;
    while (true)
    {
        //std::cout<<"hola"<<fflush<<endl;
        r = ((float)rand() / (RAND_MAX));
        float select_prob = 0;
        for (auto dest_node : (*AG_Edges)[dest][current_node])
        {
            auto current_vertex = main_graph.get_vertex(current_node);
            auto dest_vertex = main_graph.get_vertex(dest_node);
            float cost2 = cost1 + main_graph.get_edge_weight(current_vertex, dest_vertex);
            auto key = make_pair(make_pair(current_node, cost1), make_pair(dest_node, cost2));
            select_prob += Alpha[dest][key];
            //std::cout<<"\tAlpha["<<Destinations[dest]<<",("<<current_node<<","<<cost1<<")->("<<dest_node;
            //std::cout<<","<<cost2<<")]:"<<Alpha[dest][key];
            //std::cout<<"\trandom_throw:,"<<r<<",select_prob:"<<select_prob<<",dest_node:"<<dest_node<<endl;
            if (r <= select_prob)
            { //
                current_node = dest_node;
                cost1 = cost2;
                //std::cout<<"\tchild is chosen,new current node:("<<current_node<<","<<cost2<<")"<<endl;
                std::cout << "," << current_node;
                myfile << "," << current_node;
                if (current_node == Destinations[dest])
                {
                    //std::cout<<"Finished, destination found"<<endl;
                    std::cout << endl;
                    myfile << endl;
                    return;
                }
                break;
            }
        }
    }
    myfile.close();
}
void simulation_observer()
{
    bool prediction_made=false;
    std::cout<<"hello simulation_observer"<<flush<<endl;
    auto AG_Edges = main_dijkstra_alg->getSGEdges();
    //1st choose Destination
    srand(random_seed);
    int dest = rand() % Destinations.size();
    target_destination=Destinations[dest];
    //std::cout << endl << "Simulation, randomly selected destination:," << target_destination << ","<<"dest:,"<<dest<<endl;
    //std::cout << "Simulated path:" << start;
    Dijkstra_algs[dest].print_paths();
    int current_node = start;
    float r;
    int counter = 0;
    float cost1 = 0;
    //Create output file for visualization
    ofstream myfile;
    string filename = "path_" + to_string(main_graph.getVertexNum()) + "_S" + to_string(random_seed) + ".txt";
    myfile.open(filename);
    myfile << start;
    int pred_dest=0;
    /*auto map1=Alpha[dest];
    map<pair<pair<int, float>, pair<int, float>>, float>::iterator it;
    for (it = map1.begin(); it != map1.end(); it++){
        std::cout<<"Alpha_key:"<<it->first.first.first;
        std::cout<<","<<it->first.first.second<<","<<it->first.second.first<<","<<it->first.second.second<<endl;
    }*/
while (true)
    {
        if(counter++>10000){
            cerr << "Counter:" << counter << endl;
            cerr<<"Simulation error, stuck, debug AG!"<<endl;
            break;
        }
        r = ((float)rand() / (RAND_MAX));
        //std::cout << "r:" << r << endl;
        float select_prob = 0;
        auto dest_nodes = (*AG_Edges)[dest][current_node];
        //std::cout << "available dest nodes:" << dest_nodes.size() << endl;
        for (auto dest_node : (*AG_Edges)[dest][current_node])
        {
            auto current_vertex = main_graph.get_vertex(current_node);
            //std::cout << "hola dest_node" << current_vertex->getID()<<endl;
            if(!prediction_made){
                //std::cout<<"checking if X_Crit node,"<<current_vertex->getID()<<",cost1:"<<cost1<<endl;
                if(X_Crit.count(make_pair(current_vertex->getID(),cost1))>0){
                    pred_dest=X_Crit_dest[make_pair(current_vertex->getID(),cost1)];
                    if(pred_dest==target_destination){
                        Observer_Correct=true;
                        //std::cout<<"Observer prediction is:"<<Observer_Correct<<", X_Crit_Dest==pred_dest:"<<pred_dest<<endl;
                    }
                    else{
                        //std::cout<<"Observer prediction is wrong, "<< pred_dest<<"!="<<target_destination<<endl;
                    }
                    prediction_made=true;
                    //Write down result of prediction, append to file
                    ofstream myfile2;
                    string filename2 = "../experiments/OBS_exp_V," + to_string(main_graph.getVertexNum()) + ",_S," + to_string(random_seed) + ",_Bu," + to_string(Budget) + ",_Be," + to_string(Beta) + ".txt";
                    myfile2.open(filename2,ios::app);
                    myfile2<<"prediction:,"<<Observer_Correct<<endl;
                    myfile2.close();
                    }
            }
            auto dest_vertex = main_graph.get_vertex(dest_node);
            float cost2 = cost1 + main_graph.get_edge_weight(current_vertex, dest_vertex);
            auto key = make_pair(make_pair(current_node, cost1), make_pair(dest_node, cost2));
            //Skip forbidden edges
            if(Alpha[dest].find(key)==Alpha[dest].end()){
                //std::cout<<"skipping Alpha["<<Destinations[dest]<<"] key:["<<key.first.first<<","<<key.first.second<<"]["<<key.second.first<<","<<key.second.second<<"]"<<endl;
                continue;
            }
            select_prob += Alpha[dest][key];
            //std::cout<<"\tAlpha["<<Destinations[dest]<<",("<<current_node<<","<<cost1<<")->("<<dest_node;
            //std::cout<<","<<cost2<<")]:"<<Alpha[dest][key];
            //std::cout<<"\trandom_throw:,"<<r<<",select_prob:"<<select_prob<<",dest_node:"<<dest_node<<endl;
            if(select_prob==0){
                cerr<<"Select_prob cannot be 0!"<<endl;
                exit(1);
            }
            if (r <= select_prob)
            { //
                current_node = dest_node;
                cost1 = cost2;
                //std::cout<<"\tchild is chosen,new current node:["<<current_node<<","<<cost2<<"]"<<endl;
                std::cout << "," << current_node;
                myfile << "," << current_node;
                if (current_node == Destinations[dest])
                {
                    //std::cout<<",Finished, destination found"<<endl;
                    std::cout << endl;
                    myfile << endl;
                    return;
                }
                break;
            }
        }
    }
    myfile.close();
}
void simulation_observer_prefix()
{
    bool prediction_made=false;
    //std::cout<<"hello simulation_observer"<<flush<<endl;
    auto AG_Edges = main_dijkstra_alg->getSGEdges();
    //Dijkstra_algs[dest].print_paths();
    size_t current_node = 0;
    float r;
    size_t counter = 0;
    float cost1 = 0;
    paths = main_dijkstra_alg->getAGpaths2();
    int temp_seed = rand() % (RAND_MAX);
    static int simulation_run=0;
    int saving_group = simulation_run / (simulation_runs / saving_groups);
    //cout << "saving_group:" << saving_group << ",simulation_run:," << simulation_run << endl;

    simulation_run++;
    std::mt19937 gen(temp_seed);  // seed as wanted
    // Create a vector of given size
    // Assigns values returned by successive
    // calls to lambda function to every element
    // of the vector
    //Choose a pathhs 
    unsigned path = static_cast<unsigned>(paths_dist(gen));
    //cout << "\t\tchosen_path_id:," << path << ", which had prob:" << prob_path[path] << endl;
    //cout << "chosen_path_vertices:";
    //main_dijkstra_alg->print_path(path);
    //cout << endl;
    //Now find first cut
    auto prefixes_per_path = main_dijkstra_alg->getPathsToPrefix();
    int cutset_pref = -1;

    //Check which is the first prefix adjacent to a destination
    auto dest_adj_pref = -1;

    for (auto pref : (*prefixes_per_path)[path]){
        //get last node in prefix
        if(DestAdjacentVec[main_dijkstra_alg->get_prefix_ending(pref)]){
            dest_adj_pref = pref;
            //cout << "pref[,"<<pref<<",]:";
            //main_dijkstra_alg->print_prefix(pref);
            //cout << " is adjacent to destination,actual_path:";
            //main_dijkstra_alg->print_path(path);
            //cout << endl;
            break;
        }
    }
    savings_adj_prefix+=dest_predictor_prefix(dest_adj_pref).second;

    //Now do normal cutset prediction
    for (auto pref : (*prefixes_per_path)[path]){
        //cout << "\t\t working on prefix[,"<<pref<<"]";
        //main_dijkstra_alg->print_prefix(pref);
        //cout << endl;
        if(cutset_prefix.count(pref)>0){
            //cout << "\tcutset node found at node:,";
            //main_dijkstra_alg->print_prefix(pref);
            cutset_pref = pref;
            //cout<< ",cutset_pref:," << pref << endl;
            break;
        }
    }
    //Check cutset wast found
    if(cutset_pref==-1){//We never foud the cutset!!!
        //THIS IS OK, REMAINING PATHS COULD BRING BETER REWARD IN AVERAGE
        //cout << ",no pred was made" << endl;
        //cout << "Debug me, cutset was not found!!!,path:";
        //main_dijkstra_alg->print_path(path);
        //cout << endl;
        //cout << "cutset:," << cutset_prefix << endl;
        //cerr << "Debug me, cutset was not found!!!" << endl;
        //exit(1);
    }
    else{
        //Now predict destination for cutset node
        auto pred=dest_predictor_prefix(cutset_pref);
        //cout << "Predicted dest:" << pred.first << ",q:" << pred.second;
        auto actual_dest = DestinationsOrder2[(*paths)[path].back()->getID()];
        if(pred.first==actual_dest){
            //cout << ",pred was right,full savings:" << pred.second<< endl;
            auto reward=main_dijkstra_alg->get_cost_to_dest(cutset_pref, path);
            savings += reward;
            avg_savings[saving_group] += reward;
        }
        else{
            //cout << ",pred was wrong, actual dest is:"<<actual_dest<<",savings lost:"<<pred.second << endl;
        }
    }

    //Create output file for visualization
/*    ofstream myfile;
    string filename = "path_" + to_string(main_graph.getVertexNum()) + "_S" + to_string(random_seed) + ".txt";
    myfile.open(filename);
    myfile << start;
    pair<unsigned, float> prediction;
    while (true)
    {
        if(counter++>10000){
            cerr << "Counter:" << counter << endl;
            cerr<<"Simulation error, stuck, debug AG!"<<endl;
            break;
        }
        r = ((float)rand() / (RAND_MAX));
        //std::cout << "r:" << r << endl;
        float select_prob = 0;
        auto current_prefix = main_graph.get_vertex(current_node);
        if(!prediction_made){//Check if node is in cutset
            if(cutset_prefix.count(current_node)>0){
                prediction=dest_predictor_prefix(current_node);
                prediction_made=true;
                //Write down result of prediction, append to file
                ofstream myfile2;
                string filename2 = "../experiments/OBS_exp_V," + to_string(main_graph.getVertexNum()) + ",_S," + to_string(random_seed) + ",_Bu," + to_string(Budget) + ",_Be," + to_string(Beta) + ".txt";
                myfile2.open(filename2,ios::app);
                myfile2<<"prediction:,"<<Observer_Correct<<endl;
                myfile2.close();
                }
            auto dest_vertex = main_graph.get_vertex(dest_node);
            float cost2 = cost1 + main_graph.get_edge_weight(current_vertex, dest_vertex);
            auto key = make_pair(make_pair(current_node, cost1), make_pair(dest_node, cost2));
            //Skip forbidden edges
            if(Alpha[dest].find(key)==Alpha[dest].end()){
                //std::cout<<"skipping Alpha["<<Destinations[dest]<<"] key:["<<key.first.first<<","<<key.first.second<<"]["<<key.second.first<<","<<key.second.second<<"]"<<endl;
                continue;
            }
            select_prob += Alpha[dest][key];
            //std::cout<<"\tAlpha["<<Destinations[dest]<<",("<<current_node<<","<<cost1<<")->("<<dest_node;
            //std::cout<<","<<cost2<<")]:"<<Alpha[dest][key];
            //std::cout<<"\trandom_throw:,"<<r<<",select_prob:"<<select_prob<<",dest_node:"<<dest_node<<endl;
            if(select_prob==0){
                cerr<<"Select_prob cannot be 0!"<<endl;
                exit(1);
            }
            if (r <= select_prob)
            { //
                current_node = dest_node;
                cost1 = cost2;
                //std::cout<<"\tchild is chosen,new current node:["<<current_node<<","<<cost2<<"]"<<endl;
                std::cout << "," << current_node;
                myfile << "," << current_node;
                if (current_node == Destinations[dest])
                {
                    //std::cout<<",Finished, destination found"<<endl;
                    std::cout << endl;
                    myfile << endl;
                    return;
                }
                break;
            }
        }
    }
    myfile.close();*/
}

void create_grid_nodes_and_edges_from_SaT_file()
{
    std::cout << "hello SaT" << flush << endl;
    std::ifstream in(input_filename.c_str());
    std::ifstream in2(input_filename.c_str());
    lks.reserve(16000);
    string str;
    total_nodes = 0;
    int total_edges = 0;

    vector<vector<Link>> boundary_edges;
    vector<vector<Link>> dest_to_dest_edges;
    vector<vector<Link>> lks_boundary;
    vector<vector<Link>> lks_dest;

    //std::cout<<"hello1"<<flush<<endl;
    //First pass to get list of all nodes
    while (getline(in2, str))
    {
        vector<string> results;
        boost::split(results, str, boost::is_any_of(","), boost::token_compress_on);
        auto ret = nodes_set.insert(stoi(results[0]));
        if (ret.second == true)
            total_nodes++;
        ret = nodes_set.insert(stoi(results[1]));
        if (ret.second == true)
            total_nodes++;
    }
    std::cout << "total_nodes=" << total_nodes << endl;
    //If random, choose from nodes_set randomly
    //if a destination proves unreachable, WE ARE FINISHED
    //SHOULD FIND A BETTER WAY
    if (random_destinations)
    { //Both origin and destination random selection
        random_orig_and_dest_placement_SaT();
    }
    //Now initialize all variables dependent on Destinations
    boundary_edges.resize(Destinations.size());
    dest_to_dest_edges.resize(Destinations.size());

    current_graph.resize(Destinations.size());
    current_graph_dest_to_dest.resize(Destinations.size());

    lks_boundary.resize(Destinations.size());
    lks_dest.resize(Destinations.size());

    while (getline(in, str))
    {
        vector<string> results;
        boost::split(results, str, boost::is_any_of(","), boost::token_compress_on);
        //std::vector<std::string> results((std::istream_iterator<string>(iss)),
        //                             std::istream_iterator<string>());
        //std::cout<<results[0]<<flush<<endl;
        //std::cout<<results[1]<<flush<<endl;
        //std::cout<<results[2]<<flush<<endl;
        //std::cout<<"hello1a"<<flush<<endl;
        //std::cout<<","<<results[1]<<flush<<","<<results[2]<<flush<<endl;
        //std::cout<<"dest_index[results[1]]:"<<flush<<stoi(results[1])<<flush<<endl;
        if (stoi(results[1]) == start)
        {
            continue; //skipping origin incoming edges
        }
        if (stoi(results[0]) == stoi(results[1]))
        { //no loops
            std::cout << "\t\tskipped loop,from," << stoi(results[0]) << ",to," << stoi(results[1]) << endl;
            continue;
        }
        if (all_destinations.find(stoi(results[0])) != all_destinations.end())
        {
            continue; //ignoring outgoing edges from destinations
        }
        else if (all_destinations.find(stoi(results[1])) != all_destinations.end())
        { //incoming destination edge
            //Adding reverse edge for backwards dest_to_dest and boundary graphs
            Link reverse_link{stoi(results[1]), stoi(results[0]), stof(results[2])};
            lks_boundary[dest_index[stoi(results[1])]].push_back(reverse_link);
            lks_dest[dest_index[stoi(results[1])]].push_back(reverse_link);
            std::cout << "Added2 reverse link for destination_index:," << dest_index[stoi(results[1])] << ",from:," << stoi(results[1]) << ",to:," << stoi(results[0]) << ",weight:," << stof(results[2]) << endl;
        }
        Link temp_link{stoi(results[0]), stoi(results[1]), 100*round2(stof(results[2]))};

        lks.push_back(temp_link);
        total_edges++;
    }
    std::cout << "main_graph, total_edges:," << total_edges << endl;
    file_node_counter = total_nodes;
    //Now create Dijkstra Graph
    main_graph.set_number_vertices(file_node_counter);
    for (long i = 0; i < total_edges; i++)
    {
        main_graph.add_link(lks[i].u, lks[i].v, lks[i].weight);
    }

    main_graph.setv();
    std::cout << "total_nodes:" << total_nodes << ",main_graph:" << main_graph.get_number_vertices() << endl;
    main_dijkstra_alg = make_unique<DijkstraShortestPathAlg>(&main_graph);

    //BasePath* result_path = main_dijkstra_alg->get_shortest_path(
    //main_graph.get_vertex(start), main_graph.get_vertex(Destinations[0]));
    //main_graph.dump_edges();
    //for (auto edge : lks){
    //  std::cout<<"from:,"<<edge.u<<",->,"<<edge.v<<",weight:,"<<edge.weight<<endl;
    //}
    auto timenow =
        chrono::system_clock::to_time_t(chrono::system_clock::now());
    std::cout << "before determine_all_shortest_paths,time:" << ctime(&timenow) << endl;
    main_dijkstra_alg->determine_all_shortest_paths(main_graph.get_vertex(start), -1);
    timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    std::cout << "after determine_all_shortest_paths,time:" << ctime(&timenow) << endl; //exit(1);
    std::cout << "node costs:" << endl;
    int reachable_nodes = 0;
    for (auto node : nodes_set)
    {
        if (!main_dijkstra_alg->is_node_reachable(node))
        {
            reachable_nodes++;
            std::cout<<"node "<<node<<" is not reachable"<<endl;
        }
        else
        {
            std::cout<<"reachable node:"<<node<<",cost:"<<main_dijkstra_alg->get_cost(node)<<endl;
        }
    }
    std::cout << "Reachable nodes:" << reachable_nodes << endl;

    //CHECK destinations are reachable
    for (auto destination : Destinations)
    {
        if (!main_dijkstra_alg->is_node_reachable(destination))
        {
            cerr << "Destination:," << destination << ",not reachable from origin:," << start << ",stopping search" << endl;
            exit(1);
        }
        else
        {
            std::cout << "Destination:," << destination << ",is reachable from origin:," << start << ",continuing search" << endl;
        }
    }
    //std::cout<<"\tdest_to_dest,origin:,"<<start<<",dest:"<<Destinations[0]<<",Weight:"<<result_path->Weight()<<",length:"<<result_path->length()<<endl;exit(1);

    //main_dijkstra_alg->dump_edges();
    //BasePath* main_result1 = main_dijkstra_alg->get_shortest_path(
    //main_graph.get_vertex(start), main_graph.get_vertex(Destinations[1]));
    //std::cout<<"\torig_to_dest1,start:,"<<start<<",dest2:"<<Destinations[0]<<",Weight:"<<main_result1->Weight()<<",length:"<<main_result1->length()<<endl;
    //BasePath* main_result2 = main_dijkstra_alg->get_shortest_path(
    //main_graph.get_vertex(start), main_graph.get_vertex(Destinations[1]));
    //std::cout<<"\torig_to_dest2,start:,"<<start<<",dest2:"<<Destinations[1]<<",Weight:"<<main_result2->Weight()<<",length:"<<main_result2->length()<<endl;

    //Now create dest_to_dest Dijkstra Graphs
    for (size_t dest = 0; dest < Destinations.size(); dest++)
    {
        //boundary graphs
        current_graph[dest].set_number_vertices(total_nodes); //skipping origin
        for (auto add : lks_boundary[dest])
        {
            //add reverse edges specific to boundary graphs
            if (!(main_dijkstra_alg->is_node_reachable(add.v) && main_dijkstra_alg->is_node_reachable(add.u)))
            {
                continue;
            }
            //std::cout<<"adding current_graph reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<",vertices:"<<file_node_counter - 1 <<endl;
            current_graph[dest].add_link(add.u, add.v, add.weight);
        }

        current_graph_dest_to_dest[dest].set_number_vertices(total_nodes); //skipping origin
        for (auto add : lks_dest[dest])
        {
            //add reverse edges specific to dest_to_dest graphs
            if (!(main_dijkstra_alg->is_node_reachable(add.v) && main_dijkstra_alg->is_node_reachable(add.u)))
            {
                continue;
            }
            //std::cout<<"adding current_graph_dest_to_dest reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<endl;
            current_graph_dest_to_dest[dest].add_link(add.u, add.v, add.weight);
        }
        for (long i = 0; i < total_edges; i++)
        {
            if (lks[i].u == start)
            {
                continue; //skipping edges leaving origin for boundary
            }
            else if (all_destinations.find(lks[i].v) != all_destinations.end())
            { //skipping edges into any destination for boundary graphs
                //but keeping them for dest to dest search
                /*current_graph_dest_to_dest[dest].add_link(lks[i].u, lks[ i ].v, lks[ i ].weight );
    std::cout<<"Added origin outgoing link for destination_index:,"<<dest<<",from:,"<<lks[i].u<<",to:,"<<lks[i].v<<",weight:,"<<lks[i].weight<<endl;*/
                continue;
            }
            //std::cout<<"adding current_graph_reverse_link edges["<<dest<<"]"<<",lks[i].v:"<<lks[i].v<<",lks[i].u:"<<lks[i].u<<"lks[i].weight:"<<lks[i].weight<<endl;
            current_graph[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
        }
        current_graph[dest].setv();
        current_graph[dest].set_number_vertices(total_nodes); //skipping origin
        std::cout << "total_nodes:" << total_nodes << ",current_graph[" << dest << "]:," << current_graph[dest].get_number_vertices() << endl;
        Dijkstra_algs.push_back(DijkstraShortestPathAlg(&(current_graph[dest])));

        std::cout << "adding bidirectional edges for dest_to_dest[" << dest << "]" << endl;
        for (long i = 0; i < total_edges; i++)
        {
            if (lks[i].u == start)
            { //skipping origin outgoing edges
                continue;
            }
            //FOR CALCULATING GEODETIC DISTANCE, WE CREATE A LINK IN BOTH DIRECTIONS AS LONG AS
            //ONE DIRECTION EXISTS
            //std::cout<<"adding bidirectional edges for dest_to_dest["<<dest<<"]"<<",lks[i].u:"<<lks[i].u<<",lks[i].v:"<<lks[i].v<<"lks[i].weight:"<<lks[i].weight<<endl;
            current_graph_dest_to_dest[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
            current_graph_dest_to_dest[dest].add_link(lks[i].u, lks[i].v, lks[i].weight);
        }
        current_graph_dest_to_dest[dest].setv();
        current_graph_dest_to_dest[dest].set_number_vertices(total_nodes); //skipping origin
        std::cout << "total_nodes:" << total_nodes << ",current_graph_dest_to_dest[" << dest << "]:," << current_graph_dest_to_dest[dest].get_number_vertices() << endl;

        /*current_graph_dest_to_dest[dest].clear();
    current_graph_dest_to_dest[dest].set_number_vertices(3);
    if(dest==0){
      Link temp_link{116622,116623,1.5};
      current_graph_dest_to_dest[0].add_link(temp_link.u,temp_link.v,temp_link.weight);
      Link temp_link2{116623,117023,1.5};
      current_graph_dest_to_dest[0].add_link(temp_link2.u,temp_link2.v,temp_link2.weight);
    }
    else{
      current_graph_dest_to_dest[1].add_link(117023,116623,1.5);
      current_graph_dest_to_dest[1].add_link(116623,116622,1.5);
    }
    current_graph_dest_to_dest[dest].dump_edges();*/
        Dijkstra_algs_dest_to_dest.push_back(DijkstraShortestPathAlg(&(current_graph_dest_to_dest[dest])));
        for (size_t dest2 = 0; dest2 < Destinations.size(); dest2++)
        {
            if (dest == dest2)
                continue;

            BasePath *result = Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
                current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(Destinations[dest2]));
            //BasePath* result = Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
            //current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(25157));
            //std::cout<<"\t"<<result->get_path_string()<<",Weight:"<<result->Weight()<<endl;
            std::cout << "\tdest_to_dest,dest1:," << Destinations[dest] << ",dest2:" << Destinations[dest2] << ",Weight:" << result->Weight() << ",length:" << result->length() << endl;
            //BasePath* main_result1 = main_dijkstra_alg->get_shortest_path(
            //	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]));
            //    std::cout<<"\torig_to_dest1,start:,"<<start<<",dest2:"<<Destinations[dest]<<",Weight:"<<main_result1->Weight()<<",length:"<<main_result1->length()<<endl;
            //  BasePath* main_result2 = main_dijkstra_alg->get_shortest_path(
            //	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest2]));
            //    std::cout<<"\torig_to_dest2,start:,"<<start<<",dest2:"<<Destinations[dest2]<<",Weight:"<<main_result2->Weight()<<",length:"<<main_result2->length()<<endl;
        }
        //exit(1);
    }
}
//The following is the algorithm to create a new graph, called NewGraph from now on, that is shaped to take into account the effects of partial observability on the connectivity of the nodes in the original graph.
//For every observable node O in the graph create a new subgraph where:
//Make O the origin node
// Every other observable node Xi losses its outgoing edges
//get optimal paths between O and every Xi( full Dijkstra search) in the subgraph.
//If there is a direct-edge in the original graph between O and Xi:
//NewGraph gets that edge if no alternative shorter path exists.
//Else (so no direct edge between O and Xi in the original graph):
//NewGraph gets a new edge from O to Xi, the weight is the cost of the optimal path.
void create_FO_graph()
{
    vector<Link> grid_edges;
    Graph main_graph_FO;
    Graph temp_graph;
    //IN CASE THIS IS THE FIRST RANDOM CHOICE
    srand(random_seed);

    std::cout << "Creating FO graph from Partially Ovservable graph,PO_level:," << PO_level << endl;
    auto start_FO_time = std::chrono::system_clock::now();
    //vector<Link> subgraph_edges;
    vector<Link> FO_edges;
    map<int, vector<Link>> visible_outgoing_edges;
    set<int> current_nodes;
    //unordered_map<pair<int,int>,float, hash_pair> best_known_distances;
    unordered_map<pair<int, int>, float, hash_pair> best_known_distances;
    //Graph temp_graph=main_graph;//we are going to remove all edges which are not active in any of the subgraphs
    //First create all edges as if everynode was observable
    create_only_grid_edges_from_file(&grid_edges);

    //First select who is observable from original graph
    FO_nodes.insert(start); //origin must be visible!!!
    for (auto dest : Destinations)
    {
        FO_nodes.insert(dest); //origin must be visible!!!
    }
    for (auto node : main_graph_nodes)
    {
        if (rand() % 100 < PO_level)
        {
            FO_nodes.insert(node);
        }
    }
    //std::cout<<"FO_nodes:,"; for (auto node : FO_nodes) std::cout<<node<<",";
    //std::cout<<endl;

    //first create regular graph as if everything was observable
    temp_graph.set_number_vertices(main_graph_nodes.size());
    for (long i = 0; i < grid_edges.size(); i++)
    {
        temp_graph.add_link(grid_edges[i].u, grid_edges[i].v, grid_edges[i].weight);
        /*if(grid_edges[i].u==1795){
      std::cout<<"1795->"<<grid_edges[i].v<<",w:,"<<grid_edges[i].weight<<endl;
      if(FO_nodes.find(grid_edges[i].v)!=FO_nodes.end()){
    std::cout<<",visible"<<endl;
      }
      else
    std::cout<<",invisible"<<endl;
    }*/
    }
    temp_graph.setv();

    //Second iterate over each observable node
    //to create an initial graph where all visible nodes loose their outgoing edges
    unordered_set<int> NC_nodes; //Nodes whose connections are changed because they are connected to a partially visible node
    //set<int> FO_nodes_final;
    for (auto edge : grid_edges)
    {
        if (FO_nodes.find(edge.u) != FO_nodes.end())
        { //it is an outgoing invisible node
            //if(edge.u==1795){
            //std::cout<<"grid_edge:,"<<edge.u<<",->,"<<edge.v<<",origin visible"<<endl;
            //}
            //edge is between 2 observable nodes so skip it but copy it to FO graph
            if (FO_nodes.find(edge.v) != FO_nodes.end())
            {
                //std::cout<<"\tFO_edge:,"<<edge.u<<",->,"<<edge.v<<",destination visible"<<endl;
                FO_edges.push_back(edge);
                best_known_distances[make_pair<int, int>(int(edge.u), int(edge.v))] = edge.weight;
                temp_graph.remove_edge(make_pair<int, int>(int(edge.u), int(edge.v)));
                visible_outgoing_edges[edge.u].push_back(edge);
                continue;
            }
            //std::cout<<"FO_edge:,"<<edge.u<<",->,"<<edge.v<<",destination invisible"<<endl;
            temp_graph.remove_edge(make_pair<int, int>(int(edge.u), int(edge.v)));
            NC_nodes.insert(edge.u);
            visible_outgoing_edges[edge.u].push_back(edge);
            continue;
        }
        else
        {
            //std::cout<<"grid_edge:,"<<edge.u<<",->,"<<edge.v<<",origin invisible"<<endl;
            current_nodes.insert(edge.u);
            current_nodes.insert(edge.v);
        }
    }

    //Third, create a graph for each visible node as an origin
    //to found the costs associated to reaching other nodes
    std::cout << "Nodes which require new connections:,NC_nodes:," << NC_nodes.size() << ",FO_nodes:," << FO_nodes.size() << endl;
    int counter = 0;
    auto start_time = chrono::high_resolution_clock::now();
    //std::chrono::duration<double> aggregated_dijkstra_time(0);
    double aggregated_dijkstra_time = 0;

    for (auto origin : NC_nodes)
    {
        counter++;
        if (counter % 10 == 0)
        {
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end_time - start_time;
            auto current_time = diff.count();
            aggregated_dijkstra_time += current_time;
            std::cout << "\t\tNC_counter:" << counter << ",total:," << NC_nodes.size() << ",time:," << current_time << ",aggregated_dijkstra_time:," << aggregated_dijkstra_time << endl;
            start_time = std::chrono::high_resolution_clock::now();
            //aggregated_dijkstra_time=std::chrono::milliseconds(0);
        }

        //auto subgraph_edges=temp_edges;
        for (auto edge : visible_outgoing_edges[origin])
        {
            temp_graph.recover_removed_edge(make_pair<int, int>(int(edge.u), int(edge.v)));
            //std::cout<<"recovered edge from "<<edge.u<<" to "<<edge.v<<endl;
            current_nodes.insert(edge.u);
            current_nodes.insert(edge.v);
        }
        //Fourth,do a Dijkstra and store as a single edge all optimal distances
        //from the origin node to all reachable nodes
        //auto start_dijkstra_time=std::chrono::high_resolution_clock::now();
        auto temp_dijkstra_alg = make_unique<DijkstraShortestPathAlg>(&temp_graph);

        temp_dijkstra_alg->determine_all_shortest_paths(temp_graph.get_vertex(origin), INT_MAX);
        //if(origin==52){
        //  auto path=temp_dijkstra_alg->recover_shortest_perim_path(82);
        // std::cout<<"Path[52,82]:,W:,"<<path->Weight()<<",";path->PrintOut(std::cout);
        //}
        //auto end_dijkstra_time = std::chrono::high_resolution_clock::now();aggregated_dijkstra_time += end_dijkstra_time-start_dijkstra_time;
        //Now
        //std::cout<<"\t\tNC_counter:,"<<counter<<",dijkstra_time:,"<<(end_dijkstra_time-start_dijkstra_time).count()<<",size:"<<temp_dijkstra_alg->get_perim_size()<<endl;
        //start_dijkstra_time=std::chrono::high_resolution_clock::now();
        temp_dijkstra_alg->improve_distances(&FO_nodes, &FO_edges, &best_known_distances);
        //end_dijkstra_time = std::chrono::high_resolution_clock::now();aggregated_dijkstra_time += end_dijkstra_time-start_dijkstra_time;
        //std::cout<<"\t\tNC_counter2:,"<<counter<<",improve_distance_time:,"<<(end_dijkstra_time-start_dijkstra_time).count()<<",size2:"<<temp_dijkstra_alg->get_perim_size()<<endl;
        //Now
        //exit(1);
        //Remove current origin nodes and edges for next subgraph
        for (auto edge : visible_outgoing_edges[origin])
        {
            current_nodes.erase(edge.u);
            current_nodes.erase(edge.v);
            //subgraph_edges.pop_back();
            temp_graph.remove_edge(make_pair<int, int>(int(edge.u), int(edge.v)));
        }
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_FO_time;
    auto current_time = diff.count();
    std::cout << "finished building FD_graph,time:," << current_time << ",visible_nodes:," << FO_nodes.size() << ",FO_edges:" << FO_edges.size() << endl;

    //Fifth, create all secondary graphs for dest_to_dest and boundaries
    //Now initialize all variables dependent on Destinations
    vector<vector<Link>> boundary_edges;
    vector<vector<Link>> dest_to_dest_edges;
    vector<vector<Link>> lks_boundary;
    vector<vector<Link>> lks_dest;
    boundary_edges.resize(Destinations.size());
    dest_to_dest_edges.resize(Destinations.size());

    current_graph.resize(Destinations.size());
    current_graph_dest_to_dest.resize(Destinations.size());

    lks_boundary.resize(Destinations.size());
    lks_dest.resize(Destinations.size());

    int total_edges = 0;
    for (auto edge : FO_edges)
    {
        //if(edge.u==1795)
        //std::cout<<"cleaning,FO_edge:,"<<edge.u<<",->,"<<edge.v<<endl;
        if (edge.v == start)
        {
            continue;
        }
        if (edge.u == edge.v)
        { //no loops
            std::cout << "\t\tskipped loop,from," << edge.u << ",to," << edge.v << endl;
            continue;
        }
        if (all_destinations.find(edge.u) != all_destinations.end())
        {
            continue; //ignoring outgoing edges from destinations
        }
        else if (all_destinations.find(edge.v) != all_destinations.end())
        { //incoming destination edge
            //Adding reverse edge for backwards dest_to_dest and boundary graphs
            Link reverse_link{edge.v, edge.u, edge.weight};
            lks_boundary[dest_index[edge.v]].push_back(reverse_link);
            lks_dest[dest_index[edge.v]].push_back(reverse_link);
            //std::cout<<"Added reverse link for destination_index:,"<<dest_index[edge.v]<<",from:,"<<edge.v<<",to:,"<<edge.u<<",weight:,"<<edge.weight<<endl;
        }
        Link temp_link{edge.u, edge.v, edge.weight};

        lks.push_back(temp_link);
        total_edges++;
    }

    std::cout << "main_graph, total_edges:," << total_edges << ",grid_edges:," << grid_edges.size() << endl;
    //Check if any nodes FO nodes are unreachable
    unordered_set<int> final_FO_nodes;
    for (auto edge : FO_edges)
    {
        //if(edge.u==1795)
        //std::cout<<"final_add,FO_edge:,"<<edge.u<<",->,"<<edge.v<<endl;
        final_FO_nodes.insert(edge.u);
        final_FO_nodes.insert(edge.v);
    }

    for (auto node : FO_nodes)
    {
        if (final_FO_nodes.find(node) == final_FO_nodes.end())
        {
            std::cout << "\t\t FO_node:," << node << ",not reachable in final_FO_nodes, is this right?" << endl;
        }
    }

    //Now create Dijkstra Graph
    main_graph.clear();
    //main_graph.set_number_vertices( file_node_counter );
    main_graph.set_number_vertices(main_graph_nodes.size());
    for (auto edge : FO_edges)
    {
        main_graph.add_link(edge.u, edge.v, edge.weight);
        //std::cout<<"main.u:,"<<edge.u<<",main.v:,"<<edge.v<<",main.w:,"<<edge.weight<<",code:,"<<main_graph.get_edge_code(edge.u,edge.v)<<endl;
        //std::cout<<"main.u:,"<<edge.u<<",main.v:,"<<edge.v<<",main.w:,"<<edge.weight<<endl;
    }
    main_graph.setv();
    main_graph.set_number_vertices(main_graph_nodes.size());
    total_edges = FO_edges.size();
    std::cout << "total_FO_nodes:" << FO_nodes.size() << ",main_graph:" << main_graph.get_number_vertices() << ",final_FO_nodes:," << final_FO_nodes.size() << endl;
    main_dijkstra_alg = make_unique<DijkstraShortestPathAlg>(&main_graph);

    auto timenow =
        chrono::system_clock::to_time_t(chrono::system_clock::now());
    std::cout << "before determine_all_shortest_paths,time:" << ctime(&timenow) << endl;
    main_dijkstra_alg->determine_all_shortest_paths(main_graph.get_vertex(start), -1);

    /*auto path=main_dijkstra_alg->recover_shortest_perim_path(84);
  std::cout<<"Path["<<start<<",84]:,W:,"<<path->Weight()<<",";path->PrintOut(std::cout);
  auto path2=main_dijkstra_alg->recover_shortest_perim_path(20);
  std::cout<<"Path["<<start<<",20]:,W:,"<<path2->Weight()<<",";path2->PrintOut(std::cout);
  auto path3=main_dijkstra_alg->recover_shortest_perim_path(52);
  std::cout<<"Path["<<start<<",52]:,W:,"<<path3->Weight()<<",";path3->PrintOut(std::cout);*/

    std::cout << "perim_size:," << main_dijkstra_alg->get_perim_size() << endl;
    timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
    std::cout << "after determine_all_shortest_paths,time:" << ctime(&timenow) << endl; //exit(1);

    //CHECK destinations are reachable
    for (auto destination : Destinations)
    {
        if (!main_dijkstra_alg->is_node_reachable(destination))
        {
            cerr << "Destination:," << destination << ",not reachable from origin:," << start << ",stopping search" << endl;
            exit(1);
        }
        else
        {
            std::cout << "Destination:," << destination << ",is reachable from origin:," << start << ",continuing search" << endl;
        }
    }

    total_nodes = main_graph.get_number_vertices();
    //Now create dest_to_dest Dijkstra Graphs
    for (size_t dest = 0; dest < Destinations.size(); dest++)
    {
        //boundary graphs
        //current_graph[dest].set_number_vertices(total_nodes);//skipping origin
        current_graph[dest].set_number_vertices(main_graph_nodes.size());
        for (auto add : lks_boundary[dest])
        {
            //add reverse edges specific to boundary graphs
            if (!(main_dijkstra_alg->is_node_reachable(add.v) && main_dijkstra_alg->is_node_reachable(add.u)))
            {
                continue;
            }
            //std::cout<<"adding current_graph reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<",vertices:"<<file_node_counter - 1 <<endl;
            current_graph[dest].add_link(add.u, add.v, add.weight);
        }

        current_graph_dest_to_dest[dest].set_number_vertices(main_graph_nodes.size()); //skipping origin
        for (auto add : lks_dest[dest])
        {
            //add reverse edges specific to dest_to_dest graphs
            if (!(main_dijkstra_alg->is_node_reachable(add.v) && main_dijkstra_alg->is_node_reachable(add.u)))
            {
                continue;
            }
            //std::cout<<"adding current_graph_dest_to_dest reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<endl;
            current_graph_dest_to_dest[dest].add_link(add.u, add.v, add.weight);
        }
        for (long i = 0; i < total_edges; i++)
        {
            if (lks[i].u == start)
            {
                continue; //skipping edges leaving origin for boundary
            }
            else if (all_destinations.find(lks[i].v) != all_destinations.end())
            { //skipping edges into any destination for boundary graphs
                //but keeping them for dest to dest search
                /*current_graph_dest_to_dest[dest].add_link(lks[i].u, lks[ i ].v, lks[ i ].weight );
    std::cout<<"Added origin outgoing link for destination_index:,"<<dest<<",from:,"<<lks[i].u<<",to:,"<<lks[i].v<<",weight:,"<<lks[i].weight<<endl;*/
                continue;
            }
            //std::cout<<"adding current_graph_reverse_link edges["<<dest<<"]"<<",lks[i].v:"<<lks[i].v<<",lks[i].u:"<<lks[i].u<<"lks[i].weight:"<<lks[i].weight<<endl;
            current_graph[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
        }
        current_graph[dest].setv();
        current_graph[dest].set_number_vertices(main_graph_nodes.size());
        std::cout << "total_nodes:" << total_nodes << ",current_graph[" << dest << "]:," << current_graph[dest].get_number_vertices() << endl;
        Dijkstra_algs.push_back(DijkstraShortestPathAlg(&(current_graph[dest])));

        std::cout << "adding bidirectional edges for dest_to_dest[" << dest << "]" << endl;
        for (long i = 0; i < total_edges; i++)
        {
            if (lks[i].u == start)
            { //skipping origin outgoing edges
                continue;
            }
            //FOR CALCULATING GEODETIC DISTANCE, WE CREATE A LINK IN BOTH DIRECTIONS AS LONG AS
            //ONE DIRECTION EXISTS
            //std::cout<<"adding bidirectional edges for dest_to_dest["<<dest<<"]"<<",lks[i].u:"<<lks[i].u<<",lks[i].v:"<<lks[i].v<<"lks[i].weight:"<<lks[i].weight<<endl;
            current_graph_dest_to_dest[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
            current_graph_dest_to_dest[dest].add_link(lks[i].u, lks[i].v, lks[i].weight);
        }
        current_graph_dest_to_dest[dest].setv();
        current_graph_dest_to_dest[dest].set_number_vertices(main_graph_nodes.size());
        std::cout << "total_nodes:" << total_nodes << ",current_graph_dest_to_dest[" << dest << "]:," << current_graph_dest_to_dest[dest].get_number_vertices() << endl;

        /*current_graph_dest_to_dest[dest].clear();
    current_graph_dest_to_dest[dest].set_number_vertices(3);
    if(dest==0){
      Link temp_link{116622,116623,1.5};
      current_graph_dest_to_dest[0].add_link(temp_link.u,temp_link.v,temp_link.weight);
      Link temp_link2{116623,117023,1.5};
      current_graph_dest_to_dest[0].add_link(temp_link2.u,temp_link2.v,temp_link2.weight);
    }
    else{
      current_graph_dest_to_dest[1].add_link(117023,116623,1.5);
      current_graph_dest_to_dest[1].add_link(116623,116622,1.5);
    }
    current_graph_dest_to_dest[dest].dump_edges();*/
        Dijkstra_algs_dest_to_dest.push_back(DijkstraShortestPathAlg(&(current_graph_dest_to_dest[dest])));
        for (size_t dest2 = 0; dest2 < Destinations.size(); dest2++)
        {
            if (dest == dest2)
                continue;

            //BasePath* result = Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
            //  current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(Destinations[dest2]));
            //BasePath* result = Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
            //	  current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(25157));
            //std::cout<<"\t"<<result->get_path_string()<<",Weight:"<<result->Weight()<<endl;
            //std::cout<<"\tdest_to_dest,dest1:,"<<Destinations[dest]<<",dest2:"<<Destinations[dest2]<<",Weight:"<<result->Weight()<<",length:"<<result->length()<<endl;
            //BasePath* main_result1 = main_dijkstra_alg->get_shortest_path(
            //	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]));
            //    std::cout<<"\torig_to_dest1,start:,"<<start<<",dest2:"<<Destinations[dest]<<",Weight:"<<main_result1->Weight()<<",length:"<<main_result1->length()<<endl;
            //  BasePath* main_result2 = main_dijkstra_alg->get_shortest_path(
            //	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest2]));
            //    std::cout<<"\torig_to_dest2,start:,"<<start<<",dest2:"<<Destinations[dest2]<<",Weight:"<<main_result2->Weight()<<",length:"<<main_result2->length()<<endl;
        }
    }
}

void create_grid_edges_from_file()
{
    std::cout << "max_x:" << max_x << ",max_y:" << max_y << endl;
    vector<Link> grid_edges;
    vector<vector<Link>> boundary_edges;
    vector<vector<Link>> dest_to_dest_edges;

    //Check connectivity
    if (connectivity != 4 && connectivity != 8)
    {
        cerr << "Connectivity must be either 4 or 8, file not valid!" << endl;
        exit(1);
    }

    boundary_edges.resize(Destinations.size());
    dest_to_dest_edges.resize(Destinations.size());
    vector<set<int>> dest_to_dest_nodes;
    vector<set<int>> boundary_edges_nodes;
    dest_to_dest_nodes.resize(Destinations.size());
    boundary_edges_nodes.resize(Destinations.size());

    std::cout << "Destinations:";
    for (auto dest : Destinations)
    {
        std::cout << "," << dest;
        std::cout << ",pos:," << node_map[dest].first << "," << node_map[dest].second << endl;
    }
    //for(auto it : all_destinations) std::cout<<","<<it;
    std::cout << endl;
    pair<int, int> start_pos = node_map[start];

    for (size_t i = 0; i < Destinations.size(); i++)
    {
        dest_index[Destinations[i]] = i;
    }

    int counter = 0;
    vector<Graph> dest_to_dest_grid;
    bool is_dest = false;
    int current_destination = 0;
    for (int x = 0; x < max_x; x++)
    {
        for (int y = 0; y < max_y; y++)
        {
            is_dest = false;
            pair<int, int> curr_pos = make_pair(x, y);
            if (coord_map2.find(curr_pos) == coord_map2.end())
            { //unreachable node
                //std::cout<<"\tnode is unreachable"<<endl;
                continue;
            }
            if (all_destinations.find(coord_map2[curr_pos]) != all_destinations.end())
            {	//Destinations have no outgoing edges
                //std::cout<<"\tnode is destination"<<endl;
                is_dest = true;
                current_destination = coord_map2[curr_pos];
            }
            pair<int, int> W_pos(x, y + 1);
            pair<int, int> E_pos(x + 1, y + 1);
            pair<int, int> D_pos(x + 1, y);
            pair<int, int> C_pos(x + 1, y - 1);
            pair<int, int> X_pos(x, y - 1);
            pair<int, int> Z_pos(x - 1, y - 1);
            pair<int, int> A_pos(x - 1, y);
            pair<int, int> Q_pos(x - 1, y + 1);
            //W-neighour
            if (y + 1 < max_y && coord_map2.find(W_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[W_pos], 1};
                //std::cout<<"\tW"<<coord_map2[curr_pos]<<"->"<<coord_map2[W_pos]<<flush<<endl;
                //std::cout<<"\t["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<W_pos.first<<","<<W_pos.second<<"]"<<endl;
                if (!is_dest)
                {
                    if (coord_map2[W_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges.push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[W_pos]);

                        for (int i = 0; i < Destinations.size(); i++)
                        { //Destinations unreachable from boundary graph
                            if (all_destinations.find(coord_map2[W_pos]) == all_destinations.end())
                                boundary_edges[i].push_back(temp_link);
                            //boundary_edges_nodes[i].insert(temp_link.u);
                            //boundary_edges_nodes[i].insert(temp_link.v);
                        }
                    }
                    for (int i = 0; i < Destinations.size(); i++)
                    {
                        dest_to_dest_edges[i].push_back(temp_link);
                        dest_to_dest_nodes[i].insert(coord_map2[curr_pos]);
                        dest_to_dest_nodes[i].insert(coord_map2[W_pos]);
                    }
                }
                else
                { //add only to corresponding dest_to_dest grid
                    dest_to_dest_edges[dest_index[current_destination]].push_back(temp_link);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[curr_pos]);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[W_pos]);
                    if (all_destinations.find(coord_map2[W_pos]) == all_destinations.end())
                        boundary_edges[dest_index[current_destination]].push_back(temp_link);
                    //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.u);
                    //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.v);
                }
            }
            //E-neighour
            if (connectivity == 8 && y + 1 < max_y && x + 1 < max_x && coord_map2.find(E_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[E_pos], 1};
                //std::cout<<"\tE"<<coord_map2[curr_pos]<<"->"<<coord_map2[E_pos]<<flush<<endl;
                //std::cout<<"\t["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<E_pos.first<<","<<E_pos.second<<"]"<<endl;
                if (!is_dest)
                {
                    if (coord_map2[E_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges.push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[E_pos]);
                        for (int i = 0; i < Destinations.size(); i++)
                        {
                            if (all_destinations.find(coord_map2[E_pos]) == all_destinations.end())
                                boundary_edges[i].push_back(temp_link);
                            //boundary_edges_nodes[i].insert(temp_link.u);
                            //boundary_edges_nodes[i].insert(temp_link.v);
                        }
                    }
                    for (int i = 0; i < Destinations.size(); i++)
                    {
                        dest_to_dest_edges[i].push_back(temp_link);
                        dest_to_dest_nodes[i].insert(coord_map2[curr_pos]);
                        dest_to_dest_nodes[i].insert(coord_map2[E_pos]);
                    }
                }
                else
                { //add only to corresponding dest_to_dest grid
                    dest_to_dest_edges[dest_index[current_destination]].push_back(temp_link);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[curr_pos]);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[E_pos]);
                    if (all_destinations.find(coord_map2[E_pos]) == all_destinations.end())
                        boundary_edges[dest_index[current_destination]].push_back(temp_link);
                    //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.u);
                    //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.v);
                }
            }
            //D-neighour
            if (x + 1 < max_x && coord_map2.find(D_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[D_pos], 1};
                //std::cout<<"\tD["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<D_pos.first<<","<<D_pos.second<<"]"<<endl;
                if (!is_dest)
                {
                    if (coord_map2[D_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges.push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[D_pos]);
                        for (int i = 0; i < Destinations.size(); i++)
                            if (all_destinations.find(coord_map2[D_pos]) == all_destinations.end())
                            {
                                boundary_edges[i].push_back(temp_link);
                                //boundary_edges_nodes[i].insert(temp_link.u);
                                //boundary_edges_nodes[i].insert(temp_link.v);
                            }
                    }
                    for (int i = 0; i < Destinations.size(); i++)
                    {
                        dest_to_dest_edges[i].push_back(temp_link);
                        dest_to_dest_nodes[i].insert(temp_link.u);
                        dest_to_dest_nodes[i].insert(temp_link.v);
                    }
                }
                else
                { //add only to corresponding dest_to_dest grid
                    dest_to_dest_edges[dest_index[current_destination]].push_back(temp_link);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(temp_link.u);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(temp_link.v);
                    if (all_destinations.find(coord_map2[D_pos]) == all_destinations.end())
                    {
                        boundary_edges[dest_index[current_destination]].push_back(temp_link);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.u);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.v);
                    }
                }
            }
            //C-neighour
            if (connectivity == 8 && x + 1 < max_x && y - 1 > -1 && coord_map2.find(C_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[C_pos], 1};
                //std::cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[C_pos]<<endl;
                //std::cout<<"\tC["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<C_pos.first<<","<<C_pos.second<<"]"<<endl;
                if (!is_dest)
                {
                    if (coord_map2[C_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges.push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[C_pos]);
                        for (int i = 0; i < Destinations.size(); i++)
                            if (all_destinations.find(coord_map2[C_pos]) == all_destinations.end())
                            {
                                boundary_edges[i].push_back(temp_link);
                                //boundary_edges_nodes[i].insert(temp_link.u);
                                //boundary_edges_nodes[i].insert(temp_link.v);
                            }
                    }
                    for (int i = 0; i < Destinations.size(); i++)
                    {
                        dest_to_dest_edges[i].push_back(temp_link);
                        dest_to_dest_nodes[i].insert(temp_link.u);
                        dest_to_dest_nodes[i].insert(temp_link.v);
                    }
                }
                else
                { //add only to corresponding dest_to_dest grid
                    dest_to_dest_edges[dest_index[current_destination]].push_back(temp_link);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[curr_pos]);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[C_pos]);
                    if (all_destinations.find(coord_map2[C_pos]) == all_destinations.end())
                    {
                        boundary_edges[dest_index[current_destination]].push_back(temp_link);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.u);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.v);
                    }
                }
            }
            //X-neighour
            if (y - 1 > -1 && coord_map2.find(X_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[X_pos], 1};
                //std::cout<<"X\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
                //std::cout<<"\tX["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<X_pos.first<<","<<X_pos.second<<"]"<<endl;
                if (!is_dest)
                {
                    if (coord_map2[X_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges.push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[X_pos]);
                        for (int i = 0; i < Destinations.size(); i++)
                            if (all_destinations.find(coord_map2[X_pos]) == all_destinations.end())
                            {
                                boundary_edges[i].push_back(temp_link);
                                //boundary_edges_nodes[i].insert(temp_link.u);
                                //boundary_edges_nodes[i].insert(temp_link.v);
                            }
                    }
                    for (int i = 0; i < Destinations.size(); i++)
                    {
                        dest_to_dest_edges[i].push_back(temp_link);
                        dest_to_dest_nodes[i].insert(coord_map2[curr_pos]);
                        dest_to_dest_nodes[i].insert(coord_map2[X_pos]);
                    }
                }
                else
                { //add only to corresponding dest_to_dest grid
                    dest_to_dest_edges[dest_index[current_destination]].push_back(temp_link);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[curr_pos]);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[X_pos]);
                    if (all_destinations.find(coord_map2[X_pos]) == all_destinations.end())
                    {
                        boundary_edges[dest_index[current_destination]].push_back(temp_link);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.u);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.v);
                    }
                }
            }
            //Z-neighour
            if (connectivity == 8 && y - 1 > -1 && coord_map2.find(Z_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[Z_pos], 1};
                //std::cout<<"\tZ["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<Z_pos.first<<","<<Z_pos.second<<"]"<<endl;
                //std::cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
                if (!is_dest)
                {
                    if (coord_map2[Z_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges.push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[Z_pos]);
                        for (int i = 0; i < Destinations.size(); i++)
                            if (all_destinations.find(coord_map2[Z_pos]) == all_destinations.end())
                            {
                                boundary_edges[i].push_back(temp_link);
                                //boundary_edges_nodes[i].insert(temp_link.u);
                                //boundary_edges_nodes[i].insert(temp_link.v);
                            }
                    }
                    for (int i = 0; i < Destinations.size(); i++)
                    {
                        dest_to_dest_edges[i].push_back(temp_link);
                        dest_to_dest_nodes[i].insert(coord_map2[curr_pos]);
                        dest_to_dest_nodes[i].insert(coord_map2[Z_pos]);
                    }
                }
                else
                { //add only to corresponding dest_to_dest grid
                    dest_to_dest_edges[dest_index[current_destination]].push_back(temp_link);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[curr_pos]);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[Z_pos]);
                    if (all_destinations.find(coord_map2[Z_pos]) == all_destinations.end())
                    {
                        boundary_edges[dest_index[current_destination]].push_back(temp_link);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.u);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.v);
                    }
                }
            }
            //A-neighour
            if (x - 1 > -1 && coord_map2.find(A_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[A_pos], 1};
                //std::cout<<"\tA["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<A_pos.first<<","<<A_pos.second<<"]"<<endl;
                //std::cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
                if (!is_dest)
                {
                    if (coord_map2[A_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges.push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[A_pos]);
                        for (int i = 0; i < Destinations.size(); i++)
                            if (all_destinations.find(coord_map2[A_pos]) == all_destinations.end())
                            {
                                boundary_edges[i].push_back(temp_link);
                                //boundary_edges_nodes[i].insert(temp_link.u);
                                //boundary_edges_nodes[i].insert(temp_link.v);
                            }
                    }
                    for (int i = 0; i < Destinations.size(); i++)
                    {
                        dest_to_dest_edges[i].push_back(temp_link);
                        dest_to_dest_nodes[i].insert(coord_map2[curr_pos]);
                        dest_to_dest_nodes[i].insert(coord_map2[A_pos]);
                    }
                }
                else
                { //add only to corresponding dest_to_dest grid
                    dest_to_dest_edges[dest_index[current_destination]].push_back(temp_link);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[curr_pos]);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[A_pos]);
                    if (all_destinations.find(coord_map2[A_pos]) == all_destinations.end())
                    {
                        boundary_edges[dest_index[current_destination]].push_back(temp_link);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.u);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.v);
                    }
                }
            }
            //Q-neighour
            if (connectivity == 8 && x - 1 > -1 && coord_map2.find(Q_pos) != coord_map2.end())
            {
                Link temp_link{coord_map2[curr_pos], coord_map2[Q_pos], 1};
                //std::cout<<"\tQ["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<Q_pos.first<<","<<Q_pos.second<<"]"<<endl;
                //std::cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[Q_pos]<<endl;
                if (!is_dest)
                {
                    if (coord_map2[Q_pos] != start)
                    { //no imcoming edges into origin of main graph
                        grid_edges.push_back(temp_link);
                        main_graph_nodes.insert(coord_map2[curr_pos]);
                        main_graph_nodes.insert(coord_map2[Q_pos]);
                        for (int i = 0; i < Destinations.size(); i++)
                        {
                            boundary_edges[i].push_back(temp_link);
                            //boundary_edges_nodes[i].insert(temp_link.u);
                            //boundary_edges_nodes[i].insert(temp_link.v);
                        }
                    }
                    for (int i = 0; i < Destinations.size(); i++)
                    {
                        dest_to_dest_edges[i].push_back(temp_link);
                        dest_to_dest_nodes[i].insert(temp_link.u);
                        dest_to_dest_nodes[i].insert(temp_link.v);
                    }
                }
                else
                { //add only to corresponding dest_to_dest grid
                    dest_to_dest_edges[dest_index[current_destination]].push_back(temp_link);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[curr_pos]);
                    dest_to_dest_nodes[dest_index[current_destination]].insert(coord_map2[Q_pos]);
                    if (all_destinations.find(coord_map2[Q_pos]) == all_destinations.end())
                    {
                        boundary_edges[dest_index[current_destination]].push_back(temp_link);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.u);
                        //boundary_edges_nodes[dest_index[current_destination]].insert(temp_link.v);
                    }
                }
            }
        }
    }
    std::cout << "start_node:," << start << ",start_pos:" << start_pos.first << "," << start_pos.second << ",grid nodes:," << file_node_counter << ",grid_edges:," << grid_edges.size() << endl;

    int edges_counter = 0;
    for (auto link : grid_edges)
    {
        //std::cout<<"\t"<<link.u<<"->"<<link.v<<",w:"<<link.weight<<endl;
        edges_counter++;
    }
    std::cout << "edges_counter:" << edges_counter << endl;

    //Now create Dijkstra Graph
    main_graph.set_number_vertices(file_node_counter + 1);
    for (long i = 0; i < grid_edges.size(); i++)
    {
        main_graph.add_link(grid_edges[i].u, grid_edges[i].v, grid_edges[i].weight);
    }
    main_graph.setv();
    main_graph.set_number_vertices(file_node_counter + 1);

    std::cout << "maing_graph_nodes_counter:" << main_graph_nodes.size() << ",main_graph:," << main_graph.getVertexNum() << endl;

    main_dijkstra_alg = make_unique<DijkstraShortestPathAlg>(&main_graph);
    //Now create dest_to_dest Dijkstra Graphs
    current_graph.resize(Destinations.size());
    current_graph_dest_to_dest.resize(Destinations.size());

    for (size_t dest = 0; dest < Destinations.size(); dest++)
    {
        //boundary graphs
        //Note, count was off for boundary_edges in N150x150 map, NEED TO INVESTIGATE!!!
        //For now we just re-count to make sure they are right
        boundary_edges_nodes[dest].clear();
        for (long i = 0; i < boundary_edges[dest].size(); i++)
        {
            boundary_edges_nodes[dest].insert(boundary_edges[dest][i].u);
            boundary_edges_nodes[dest].insert(boundary_edges[dest][i].v);
        }

        current_graph[dest].set_number_vertices(file_node_counter + 1);
        //current_graph[dest].set_number_vertices(main_graph.getVertexNum());
        for (long i = 0; i < boundary_edges[dest].size(); i++)
        {
            current_graph[dest].add_link(boundary_edges[dest][i].u, boundary_edges[dest][i].v, boundary_edges[dest][i].weight);
        }
        current_graph[dest].setv();
        current_graph[dest].set_number_vertices(file_node_counter + 1);

        std::cout << "boundary_edges_nodes[:" << dest << "[:," << boundary_edges_nodes[dest].size() << ",graph:," << current_graph[dest].getVertexNum() << endl;
        //std::cout<<"boundary_edges_nodes[:"<<dest<<"[:,"<<main_graph.getVertexNum()<<",graph:,"<<current_graph[dest].getVertexNum()<<endl;

        Dijkstra_algs.push_back(DijkstraShortestPathAlg(&(current_graph[dest])));
        //dest_to_dest graphs
        current_graph_dest_to_dest[dest].set_number_vertices(file_node_counter + 1);
        for (long i = 0; i < dest_to_dest_edges[dest].size(); i++)
        {
            current_graph_dest_to_dest[dest].add_link(dest_to_dest_edges[dest][i].u, dest_to_dest_edges[dest][i].v, dest_to_dest_edges[dest][i].weight);
        }
        current_graph_dest_to_dest[dest].setv();
        current_graph_dest_to_dest[dest].set_number_vertices(file_node_counter + 1);

        Dijkstra_algs_dest_to_dest.push_back(DijkstraShortestPathAlg(&(current_graph_dest_to_dest[dest])));
        for (size_t dest2 = 0; dest2 < Destinations.size(); dest2++)
        {
            if (dest == dest2)
                continue;

            BasePath *result = Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
                current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(Destinations[dest2]));
            //BasePath* result = Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
            //	  current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(25157));
            //std::cout<<"\t"<<result->get_path_string()<<",Weight:"<<result->Weight()<<endl;
            std::cout << "\tdest_to_dest,dest1:," << Destinations[dest] << ",dest2:" << Destinations[dest2] << ",Weight:" << result->Weight() << ",length:" << result->length() << endl;
            if (result->length() == 0)
            {
                cerr << "Cannot find path from dest1:," << Destinations[dest] << ",to dest:," << Destinations[dest2] << endl;
                exit(1);
            }
        }
    }
    //exit(1);
    //std::cout<<"getting shortest destination path to all nodes from origin in main graph"<<endl;
    //main_dijkstra_alg->determine_all_shortest_paths(main_graph.get_vertex(start),-1);
}

void create_grid_nodes(int n)
{
    //ASSUMING ORIG IS NODE 0 at x=0,y=0, FUTURE FIX THIS FOR ANY ORIG!!!!!
    int counter = 0;
    //First map coordinate to integer named cells
    for (int y = 0; y < n; y++)
    {
        for (int x = 0; x < n; x++)
        {
            string n1 = to_string(x) + "_" + to_string(y);
            coord_map[n1] = counter;
            node_map[counter] = make_pair(x, y);
            //std::cout<<n1<<"->"<<counter<<endl;
            //  if(all_destinations.find(counter)!=all_destinations.end()){
            //	std::cout<<FORERED<<"D "<<RESETTEXT;
            //    }
            //  else if(counter==start){
            //	std::cout<<FOREGRN<<"O "<<RESETTEXT;
            //    }
            //  else{
            //if(debug)
            //std::cout<<counter<<"\t";
            //std::cout<<"_ ";
            //}
            //std::cout<<RESETTEXT;
            counter++;
        }
        //if(debug)
        //std::cout<<endl;
    }
    //exit(1);
}

//dest_to_dest graph generates extra reverse paths from origin destination
//it also genertes extra reverse paths for the actual origin node in case the optimal path uses the origin node
void create_grid_edges_dest_to_dest(int n, vector<Link> &grid_edges, int from_dest)
{
    string dest_to_dest_filename = "d";
    dest_to_dest_filename += to_string(from_dest);
    dest_to_dest_filename += "_geodetic.dot";
    std::ofstream outfile(dest_to_dest_filename);

    std::cout << "Destination:";
    for (auto it : all_destinations)
        std::cout << "," << it;
    std::cout << endl;

    int counter = 0;
    int edges = 2 * (2 * n * n - 2 * n); //edges in a square graph
    //Now create edges between cells
    for (int x = 0; x < n; x++)
    {
        for (int y = 0; y < n; y++)
        {
            string n1 = std::to_string(x) + "_" + to_string(y);
            if (all_destinations.find(coord_map[n1]) != all_destinations.end() &&
                coord_map[n1] != from_dest)
            { //No edges leaving non-origin destination
                continue;
            }
            //RIGHT NEIGHOUR
            if (x + 1 < n)
            {
                string n2 = std::to_string(x + 1) + "_" + to_string(y);
                if (coord_map[n1] != from_dest || all_destinations.find(coord_map[n2]) == all_destinations.end())
                { //destinations cannot be connected directly!
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    counter++;
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
            //LEFT NEIGHOUR
            if (x - 1 > -1)
            {
                string n2 = std::to_string(x - 1) + "_" + to_string(y);
                if (coord_map[n1] != from_dest || all_destinations.find(coord_map[n2]) == all_destinations.end())
                { //destinations cannot be connected directly!
                    counter++;
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
            //TOP NEIGHOUR
            if (y + 1 < n)
            {
                string n2 = std::to_string(x) + "_" + to_string(y + 1);
                if (coord_map[n1] != from_dest || all_destinations.find(coord_map[n2]) == all_destinations.end())
                { //destinations cannot be connected directly!
                    counter++;
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
            //BOTTOM NEIGHOUR
            if (y - 1 > -1)
            {
                string n2 = std::to_string(x) + "_" + to_string(y - 1);
                if (coord_map[n1] != from_dest || all_destinations.find(coord_map[n2]) == all_destinations.end())
                { //destinations cannot be connected directly!
                    counter++;
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
        }
    }
    outfile.close();
    //std::cout<<"nodes:,"<<n*n<<"total_edges:,"<<counter<<",max_edges:,"<<edges<<endl;
}

void create_grid_edges_boundary(int n, vector<Link> &grid_edges, int from_dest)
{
    string dest_to_dest_filename = "d";
    dest_to_dest_filename += to_string(from_dest);
    dest_to_dest_filename += "_boundary.dot";
    std::ofstream outfile(dest_to_dest_filename);

    std::cout << "Destination:";
    for (auto it : all_destinations)
        std::cout << "," << it;
    std::cout << endl;

    int counter = 0;
    int edges = 2 * (2 * n * n - 2 * n); //edges in a square graph
    //Now create edges between cells
    for (int x = 0; x < n; x++)
    {
        for (int y = 0; y < n; y++)
        {
            string n1 = std::to_string(x) + "_" + to_string(y);
            if (all_destinations.find(coord_map[n1]) != all_destinations.end() &&
                coord_map[n1] != from_dest)
            { //No edges leaving non-origin destination
                continue;
            }
            if (coord_map[n1] == start)
            { //no outgoing edges from origin
                continue;
            }
            //RIGHT NEIGHOUR
            if (x + 1 < n)
            {
                string n2 = std::to_string(x + 1) + "_" + to_string(y);
                if (all_destinations.find(coord_map[n2]) == all_destinations.end())
                { //destinations cannot be reached in boundary graph
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    counter++;
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
            //LEFT NEIGHOUR
            if (x - 1 > -1)
            {
                string n2 = std::to_string(x - 1) + "_" + to_string(y);
                if (all_destinations.find(coord_map[n2]) == all_destinations.end())
                { //destinations cannot be reached in boundary graph
                    counter++;
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
            //TOP NEIGHOUR
            if (y + 1 < n)
            {
                string n2 = std::to_string(x) + "_" + to_string(y + 1);
                if (all_destinations.find(coord_map[n2]) == all_destinations.end())
                { //destinations cannot be reached in boundary graph
                    counter++;
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
            //BOTTOM NEIGHOUR
            if (y - 1 > -1)
            {
                string n2 = std::to_string(x) + "_" + to_string(y - 1);
                if (all_destinations.find(coord_map[n2]) == all_destinations.end())
                { //destinations cannot be reached in boundary graph
                    counter++;
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
        }
    }
    outfile.close();
    //std::cout<<"nodes:,"<<n*n<<"total_edges:,"<<counter<<",max_edges:,"<<edges<<endl;
}

void create_grid_edges(int n, vector<Link> &grid_edges)
{
    std::ofstream outfile("LastGraph.dot");

    set<int> dest_set;
    std::cout << "Destinations:";
    for (auto it : all_destinations)
        std::cout << "," << it;
    std::cout << endl;

    //for(auto it : coord_map) std::cout<<it.first<<"->"<<it.second<<endl;

    //string orig_str=std::to_string(node_map[start].first) + "_" + to_string(node_map[start].second);
    int counter = 0;
    int edges = 2 * (2 * n * n - 2 * n); //edges in a square graph
    //Now create edges between cells
    for (int y = 0; y < n; y++)
    {
        for (int x = 0; x < n; x++)
        {
            string n1 = std::to_string(x) + "_" + to_string(y);
            if (all_destinations.find(coord_map[n1]) != all_destinations.end())
            { //No edges leaving destinations
                continue;
            }
            //RIGHT NEIGHOUR
            if (x + 1 < n)
            {
                string n2 = std::to_string(x + 1) + "_" + to_string(y);
                if (coord_map[n2] != start)
                { //no imcoming edges into origin
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    counter++;
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
            //LEFT NEIGHOUR
            if (x - 1 > -1)
            {
                string n2 = std::to_string(x - 1) + "_" + to_string(y);
                if (coord_map[n2] != start)
                {
                    counter++;
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
            //TOP NEIGHOUR
            if (y + 1 < n)
            {
                string n2 = std::to_string(x) + "_" + to_string(y + 1);
                if (coord_map[n2] != start)
                {
                    counter++;
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
            //BOTTOM NEIGHOUR
            if (y - 1 > -1)
            {
                string n2 = std::to_string(x) + "_" + to_string(y - 1);
                if (coord_map[n2] != start)
                {
                    counter++;
                    Link temp_link{coord_map[n1], coord_map[n2], 1};
                    grid_edges.push_back(temp_link);
                    //std::cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
                    if (debug)
                        outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
                }
            }
        }
    }
    outfile.close();
    //std::cout<<"nodes:,"<<n*n<<"total_edges:,"<<counter<<",max_edges:,"<<edges<<endl;
}

//Choose origin and destinations according to min_dist_orig,max_dist_orig,min_dest_dist,max_dest_dist
void random_dest_and_origin_placement()
{
    //std::cout<<"hola random_dest_and_origin_placement"<<flush<<endl;
    srand(random_seed);

    int min_dist_orig = 1;
    int max_dist_orig = 10000;
    int min_dest_dist = 1;
    int max_dest_dist = 10000;

    //if(start==0)//not manually chosen
    start = rand() % node_map.size();
    if (random_origin_only)
    { //positions always the same, using srand(2)
        random_seed = random_positions + N;
        srand(random_seed);
    }
    //std::cout<<"hola00"<<flush<<endl;
    set<int> available_pos;

    for (int i = 0; i < node_map.size(); i++)
    {
        if (i == start && !random_origin_only)
            continue; //origin is not available
        int dist = get_manhattan_dist(start, i);
        if (!random_origin_only && (dist < min_dist_orig || dist > max_dist_orig))
            continue;
        available_pos.insert(i);
    }
    if (debug)
    {
        std::cout << "available_pos for first_dest:";
        for (auto x : available_pos)
            std::cout << x << ",";
        std::cout << endl;
    }
    //First choose a random destination
    Destinations.clear();
    auto first_pos_index = rand() % available_pos.size();
    auto first_dest = *select_random(available_pos, first_pos_index);
    OrigDests.push_back(get_manhattan_dist(start, first_dest));
    avg_dest_orig_dest += OrigDests.back();
    //std::cout<<"first_dest:"<<first_dest<<endl;
    Destinations.push_back(first_dest);
    available_pos.erase(first_dest);
    //Now calculate available positions given chosen destinations
    int counter = 1;
    while (Destinations.size() < random_positions)
    {
        vector<int> candidate_available_pos;
        for (auto cand_node : available_pos)
        {
            for (auto node_dest : Destinations)
            {
                int dist = get_manhattan_dist(node_dest, cand_node);
                if ((dist >= min_dest_dist && dist <= max_dest_dist))
                {
                    candidate_available_pos.push_back(cand_node);
                    //std::cout<<"\t\tcand_node:,"<<cand_node<<",node_dest:,"<<node_dest<<",dist:,"<<dist<<",position added,counter:,"<<counter++<<endl;
                    break;
                }
            }
        }
        //Exit if we have no available choice
        if (candidate_available_pos.size() == 0)
        { //there are no more possible random nodes under existing constraints
            std::cout << "random_positions unachivable, exiting" << endl;
            exit(1);
        }
        //Add to destinations random legal choice
        int next_dest = 0;
        int counter_select = 0;
        while (true)
        {
            auto next_pos_index = rand() % available_pos.size();
            next_dest = *select_random(available_pos, next_pos_index);
            if (Destinations.size() == 0)
                break; //first destination can be at any distance
            int min_intra_dist = INT_MAX;
            for (auto node : Destinations)
            {
                min_intra_dist = min(min_intra_dist, get_manhattan_dist(next_dest, node));
            }
            if (min_intra_dist <= min_dist_to_dest && min_intra_dist >= max_dist_to_dest)
            {
                std::cout << "\t\t dest:" << next_dest << "is " << min_intra_dist << " from closest destination" << endl;
                break; //
            }
            //else if(debug){
            //	std::cout<<"\t\t dest:"<<next_dest<<"is "<<min_intra_dist<<" from closest destination,D="<<Destinations<<endl;
            //     }
            if (counter_select++ > 100000)
            {
                cerr << "cannot find random destination at max destination distance=" << min_dist_to_dest << ",min:" << max_dist_to_dest << ",N=" << N << ",R=" << random_positions << ",already_chosen:" << Destinations << endl;
                exit(25);
            }
        }
        Destinations.push_back(next_dest);
        available_pos.erase(next_dest);

        //std::cout<<"candidate_available_pos:"<<candidate_available_pos<<",chosen dests:"<<Destinations<<endl;
        OrigDests.push_back(get_manhattan_dist(start, next_dest));
        avg_dest_orig_dest += OrigDests.back();
    }
    avg_dest_orig_dest = avg_dest_orig_dest / float(OrigDests.size());
    for (auto node1 : Destinations)
    {
        int counter = 0;
        float current_intradist = 0;
        for (auto node2 : Destinations)
        {
            if (node1 == node2)
                continue;
            current_intradist += (get_manhattan_dist(node1, node2));
            counter++;
        }
        avg_dest_intradist += current_intradist / float(counter);
    }
    avg_dest_intradist = avg_dest_intradist / float(Destinations.size());

    set<int> temp_dest_set;
    for (auto it : Destinations)
        temp_dest_set.insert(it);
    if (temp_dest_set.find(start) != temp_dest_set.end())
    {
        cerr << "chosen destination as origin, exiting" << endl;
        exit(26);
    }

    std::cout << "\n\n\n"
         << endl;
    //Printout current configuration
    for (auto it : Destinations)
        all_destinations.insert(it);
    counter = 0;

    for (int y = 0; y < N; y++)
    {
        for (int x = 0; x < N; x++)
        {
            if (all_destinations.find(counter) != all_destinations.end())
            {
                std::cout << FORERED << "D " << RESETTEXT;
            }
            else if (counter == start)
            {
                std::cout << FOREGRN << "O " << RESETTEXT;
            }
            else
            {
                if ((y == 0) || (y == N - 1))
                {
                    std::cout << ". ";
                }
                else if ((x == 0) || (x == N - 1))
                {
                    std::cout << ". ";
                }
                else
                {
                    std::cout << "  ";
                }
            }
            counter++;
        }
        std::cout << endl;
    }
    std::cout << "\n\nstart_pos:," << start << ",Random_Destinations:" << Destinations << ",random_seed:" << random_seed << ",add_lambda:" << lambda_add << endl;
    //exit(1);
}
void from_file_random_dest_and_origin_selection()
{
    std::cout<<"hola random_dest_and_origin1"<<endl;
    srand(random_seed);
    start = rand() % node_map.size();
    if (random_origin_only)
    { //positions always the same, using srand(2)
        random_seed = random_positions + N;
        srand(random_seed);
    }
    Destinations.clear();
    auto item = node_map.begin();
    while (Destinations.size() < random_positions)
    {
        auto candidate_pos = rand() % node_map.size();
        if (candidate_pos == start) //origin cannot be destination
            continue;
        if (all_destinations.find(candidate_pos) == all_destinations.end())
        {
            all_destinations.insert(candidate_pos);
            Destinations.push_back(candidate_pos);
        }
    }

    //bitset<3> color_map[256][256];//origin/destination/FinalPath

    std::cout << "\n\nstart_pos:," << start << ",Random_Destinations:" << Destinations << ",random_seed:" << random_seed << ",add_lambda:" << lambda_add << endl;
}

int main(int argc, char *argv[])
{

    main_start_time = std::chrono::high_resolution_clock::now();
    int array_size = 0;
    int nodes = 0;
    vector<Link> lks;

    int bifurcated_destinations = 0;
    int total_pairs = 0;
    int total_triples = 0;
    int perimeter_nodes = 0;
    int neighourhood_nodes = 0;
    //time variables
    auto end_time = chrono::high_resolution_clock::now();
    auto diff = end_time - start_time;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
    auto ms_bifur_backward = ms;
    auto ms_pairs = ms;
    auto ms_triples = ms;
    //cpu_timer timer;
    //std::cout<<"timer:"<<timer.format()<<endl;
    //int edges=2*(2*n*n-2*n);//edges in a square graph
    bool optimize_target = false;
    bool optimize_target_u = false;
    bool optimize_observer = false;
    bool optimize_budget = false;
    bool SaT_map = false;
    bool AugGraph = false;
    LoadPaths = false;

    vector<vector<Link>> grid_vect_orig_dest;
    vector<vector<Link>> grid_vect_dest_to_dest;

    for (int i = 0; i < argc; i++)
    {
        std::string action(argv[i]);
        if (action == "OnlineAG")
        {
            //using main2 as the main for online AG until I merge Dijkstra and A* variants
            main2(argc,argv); 
            exit(0);
        }
        if (action == "AG")
        {
            AugGraph = true;
        }
        if (action == "load_paths")
        {
            LoadPaths = true;
            read_paths_from_file();
        }
        if (action == "debug")
        {
            debug = true;
        }
        if (action == "M")
        {
            std::string action(argv[i]);
            map_display = true;
        }
        if (action == "random_origin_only")
        {
            std::cout << "random_origin_only=true" << endl;
            random_origin_only = true;
        }
        if (action == "optimize_budget")
        {
            optimize_budget = true;
        }
        if (action == "stocastic_go")
        {
            optimize_budget = true;
            stocastic_go = true;
        }
        if (action == "stocastic_go_observer")
        {
            optimize_budget = true;
            stocastic_go_observer = true;
        }
        if (action == "stocastic_fictitious_play")
        {
            optimize_budget = true;
            stocastic_go_observer = true;//for initial s
            stocastic_ficitious_play = true;
        }
        if (action == "min_cal")
        {
            min_cal = true;
        }
        if (action == "optimize_target")
        {
            optimize_target = true;
        }
        if (action == "optimize_lambda_U")
        {
            optimize_target_u = true;
        }
        if (action == "optimize_observer")
        {
            optimize_observer = true;
        }
        if (action.find("filename=") != string::npos)
        {
            string temp = action.substr(9, action.length());
            input_filename = temp;
            std::cout << "input_filename:" << input_filename << endl;
        }
        if (action.find("filename_SaT=") != string::npos)
        {
            string temp = action.substr(13, action.length());
            input_filename = temp;
            std::cout << "input_SaT_filename:" << input_filename << endl;
            SaT_map = true;
        }
        if (action.find("LA=") != string::npos)
        {
            string temp = action.substr(3, action.length());
            lambda_add = stoi(temp);
        }
        //Partial observable graph choose obsevable nodes from 1 to 99% chance
        if (action.find("Z=") != string::npos)
        {
            create_PO_graph = true;
            string temp = action.substr(2, action.length());
            PO_level = stoi(temp);
            if (PO_level < 1 || PO_level > 100)
            {
                cerr << "PO_level must be between 1 and 99 to randomly select observable nodes." << endl;
                exit(1);
            }
            std::cout << "PO_level:" << PO_level << endl;
        }
        if (action.find("IL=") != string::npos)
        {
            string temp = action.substr(3, action.length());
            input_lambda = stof(temp);
            if (debug)
                std::cout << "IL:" << input_lambda << endl;
        }
        if (action.find("LA=") != string::npos)
        {
            string temp = action.substr(3, action.length());
            lambda_add = stoi(temp);
            if (debug)
                std::cout << "LA:" << lambda_add << endl;
        }
        if (action.find("MDD1=") != string::npos)
        {
            string temp = action.substr(5, action.length());
            min_dist_to_dest = stof(temp);
            if (debug)
                std::cout << "min_dist_to_dest:" << min_dist_to_dest << endl;
        }
        if (action.find("MDD2=") != string::npos)
        {
            string temp = action.substr(5, action.length());
            max_dist_to_dest = stof(temp);
            if (debug)
                std::cout << "max_dist_to_dest:" << max_dist_to_dest << endl;
        }
        if (action.find("C=") != string::npos)
        {
            string temp = action.substr(2, action.length());
            C = stoi(temp);
            if (debug)
                std::cout << "C:" << C << endl;
        }
        if (action.find("O=") != string::npos)
        {
            string temp = action.substr(2, action.length());
            optimal_limit = stof(temp);
            std::cout << "O:" << optimal_limit << endl;
        }
        if (action.find("sat=") != string::npos)
        {
            string temp = action.substr(4, action.length());
            //cout << "temp:" << temp << flush << endl;
            saturation = stof(temp);
            std::cout << "sat:," << saturation << endl;
        }
        if (action.find("type=") != string::npos)
        {
            string temp = action.substr(5, action.length());
            //cout << "temp:" << temp << flush << endl;
            q_type = stoi(temp);
            std::cout << "q_type:," << q_type << endl;
        }
        if (action.find("K=") != string::npos)
        {
            string temp = action.substr(2, action.length());
            K = stoi(temp);
            if (debug)
                std::cout << "K:" << C << endl;
        }
        if (action.find("D=") != string::npos)
        {
            Destinations.clear();
            string temp = action.substr(2, action.length());
            std::stringstream ss(temp);
            int i;
            while (ss >> i)
            {
                Destinations.push_back(i);
                all_destinations.insert(i);

                if (ss.peek() == ',')
                    ss.ignore();
            }
            for (size_t i = 0; i < Destinations.size(); i++)
            {
                dest_index[Destinations[i]] = i;
            }
        }
        if (action.find("OG=") != string::npos)
        {
            string temp = action.substr(3, action.length());
            start = stoi(temp);
            std::cout << "start=" << start << endl;
        }
        if (action.find("R=") != string::npos)
        {
            //if (debug)
                std::cout << "Using random destinations" << endl;
            random_destinations = true;
            string temp = action.substr(2, action.length());
            random_positions = stoi(temp);
        }
        if (action.find("S=") != string::npos)
        {
            string temp = action.substr(2, action.length());
            random_seed = stoi(temp);
            std::cout << "Using random seed=" << random_seed << endl;
        }
        if (action.find("I=") != string::npos)
        {
            string temp = action.substr(2, action.length());
            iterations = stoi(temp);
            std::cout << "Random Simulation Iterations=" << iterations << endl;
        }
        if (action.find("N=") != string::npos)
        {
            Destinations.clear();
            string temp = action.substr(2, action.length());
            N = stoi(temp);
            std::cout << "N:" << N << endl;
        }
        if (action.find("Stg=") != string::npos)
        {
            string temp = action.substr(4, action.length());
            strategy = stoi(temp);
            std::cout << "strategy:" << strategy << endl;
        }
        if (action.find("Beta=") != string::npos)
        {
            string temp = action.substr(5, action.length());
            Beta = stof(temp);
            std::cout << "Beta:" << Beta << endl;
        }
        if (action.find("Budget=") != string::npos)
        {
            string temp = action.substr(7, action.length());
            Budget = stof(temp);
            std::cout << "Budget:" << Budget << endl;
        }
        if (action == "cyclic")
        {
            acyclic = false;
        }
        if (action == "grandparent_check")
        {
            acyclic = false;
            grandparent_check = true;
            std::cout<<"grandparent_check="<<grandparent_check<<endl;
        }
        if (action.find("i=") != string::npos){
            string temp = action.substr(2, action.length());
            iterations = stoi(temp);
            std::cout << "iterations:" << iterations << endl;
        }
    }
    //Check input file exists
    ifstream f(input_filename);
    if (input_filename.length() > 0){
        if(!f.good()){
            cerr<<"graph input file:"<<input_filename<<" does not exist!"<<endl;
            exit(31);
        }
    }
    int positions = N * N - 1;
    if (optimize_target == false &&
        optimize_observer == false &&
        optimize_target_u == false &&
        optimize_budget == false)
    {
        std::cout << "Choose either optimize_target, optimize_observer or optimize_budget" << endl;
        exit(1);
    }
    else if (optimize_target && (optimize_observer || optimize_budget) ||
             optimize_budget && (optimize_target || optimize_observer))
    {
        std::cout << "Choose only optimize_target or optimize_observer" << endl;
        exit(1);
    }

    //Create grid nodes, necessary for random dest and origin placement
    if (input_filename.length() == 0)
    {
        create_grid_nodes(N);
        if (random_destinations)
        {
            random_dest_and_origin_placement();
        }
    }
    else if (SaT_map == true)
    {
        if (AugGraph)
        {
            //create_augmented_graph_from_SaT_file();
            create_prefix_graph_from_SaT_file();
        }
/*		else if(OnlineAugGraph){
            create_online_augmented_graph_from_SaT_file();
        }*/
        else
        {
            create_grid_nodes_and_edges_from_SaT_file();
        }
        //exit(1);
    }
    else
    {
        create_grid_nodes_from_file();
        if (random_destinations)
            from_file_random_dest_and_origin_selection();
        if (PO_level == 100)
            create_grid_edges_from_file();
        else
        {
            create_FO_graph();
            SaT_map = true; //Hack, treat the same way because dest_to_dest graphs are done
        }
    }
    //First, no dests
    if (input_filename.length() == 0)
    {
        vector<Link> grid_vect;
        std::cout << "Creating main grid with original origin" << endl;
        create_grid_edges(N, grid_vect);

        lks.resize(grid_vect.size());
        cout << "lks_size:" << lks.size() << endl;
        for (size_t i = 0; i < grid_vect.size(); i++)
        {
            lks[i] = grid_vect[i];
        }
        array_size = sizeof(lks) / sizeof(lks[0]);
        nodes = N * N;
        //Create main graph
        main_graph.set_number_vertices(nodes);
        std::cout << "main_graph numberof nodes:" << nodes << endl;

        for (int i = 0; i < array_size; i++)
        {
            //std::cout<<"Adding link from:"<< lks[ i ].u<<",to:"<<lks[i].v<<",weight:"<<lks[i].weight<<endl;
            main_graph.add_link(lks[i].u, lks[i].v, lks[i].weight);
        }
        main_graph.setv();
        main_dijkstra_alg = make_unique<DijkstraShortestPathAlg>(&main_graph);
        //Second, we include outer links for destinations

        grid_vect_orig_dest.resize(Destinations.size());
        grid_vect_dest_to_dest.resize(Destinations.size());
        for (size_t dest = 0; dest < Destinations.size(); dest++)
        {
            vector<int> mod_dests;
            for (size_t temp_dest = 0; temp_dest < Destinations.size(); temp_dest++)
            {
                if (dest != temp_dest)
                {
                    mod_dests.push_back(Destinations[temp_dest]);
                }
            }
            std::cout << "Creating grid with origin at destination:" << Destinations[dest] << endl;
            create_grid_edges_boundary(N, grid_vect_orig_dest[dest], Destinations[dest]);		 //blocking leaving edges from destination
            create_grid_edges_dest_to_dest(N, grid_vect_dest_to_dest[dest], Destinations[dest]); //blocking leaving edges from destination
        }
    }

    /*link lks[] = {  {0, 1},
        {0, 12},
        {0, 13},
        {1, 2},
        {2, 1},
        {2, 3},
        {2, 5},
        {4, 3},
        {4, 5},
        {4, 6},
        {5, 2},
        {5, 4},
        {5, 7},
        {6, 4},
        {6, 8},
        {7, 5},
        {7, 10},
        {8, 6},
        {8, 10},
        {9, 8},
        {9, 11},
        {10, 7},
        {10, 8},
        {10, 11},
        {11, 9},
        {11, 10},
        {11, 12},
        {11, 13},
        {13, 11},
        {13, 14}};*/

    if (K == -1) //Default value
        K = 2000;
    //std::cout<<"YenNodes:"<<nodes<<endl;
    auto PM = std::make_shared<PathMatrix>();
    //set<int> unsensed_nodes;
    //unsensed_nodes->insert(12);
    //PM->set_W(unsensed_nodes);
    PM->set_C(C);
    PM->set_min_cal(min_cal);
    //std::cout<<"C:,"<<C<<",min_cal:,"<<min_cal<<endl;
    std::cout << "hola Destinations:" << Destinations << endl;

    Link *lks_array = &lks[0];
    optimize_budget = false;
    if (!optimize_budget)
    {
        for (auto end : Destinations)
        {
            testYenAlg(K,
                       lks_array,
                       lks.size(),
                       start,
                       end,
                       nodes,
                       PM);

            auto end_time = chrono::high_resolution_clock::now();
            auto diff = end_time - start_time;
            start_time = end_time; //restart
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
            std::cout << "Time(ms) for producing yen paths in secs:" << ms << endl;
            YenTime = ms;
        }
        exit(1);
    }
    else
    { //optimize_budget!
        if (debug)
            //std::cout<<FORERED<<RESETTEXT<<"First get distances for each destination up to input_lambda="<<input_lambda<<RESETTEXT<<endl;
            std::cout << "First get distances for each destination up to input_lambda=" << input_lambda << endl;
        //Zero, Create set of initial DIJKSTRA Graphs where each destination is an origin with outgoing links (the rest remain as sinks)

        double min_rho = INT_MAX;
        double max_rho = 0;
        //vector<unique_ptr<Graph> > graph_dests;
        vector<Graph *> graph_dests;

        vector<set<int>> perimeter_nodes_vect(Destinations.size());
        vector<set<int>> neighourhood_nodes_vect(Destinations.size());
        if (!SaT_map)
        {
            current_graph.resize(Destinations.size());
            current_graph_dest_to_dest.resize(Destinations.size());
        }
        for (size_t dest = 0; dest < Destinations.size(); dest++)
        {
            if (input_filename.length() == 0)
            {
                //perimeter_graph_nodes keep track of node count for graph getting the correct weights
                set<int> all_perimeter_graph_nodes;
                set<int> all_dest_to_dest_graph_nodes;
                Link lks2[grid_vect_orig_dest[dest].size()];
                Link lks3[grid_vect_dest_to_dest[dest].size()];
                for (size_t i = 0; i < grid_vect_orig_dest[dest].size(); i++)
                {
                    lks2[i] = grid_vect_orig_dest[dest][i];
                    all_perimeter_graph_nodes.insert(grid_vect_orig_dest[dest][i].u);
                    all_perimeter_graph_nodes.insert(grid_vect_orig_dest[dest][i].v);
                }
                for (size_t i = 0; i < grid_vect_dest_to_dest[dest].size(); i++)
                {
                    lks3[i] = grid_vect_dest_to_dest[dest][i];
                    all_dest_to_dest_graph_nodes.insert(grid_vect_dest_to_dest[dest][i].u);
                    all_dest_to_dest_graph_nodes.insert(grid_vect_dest_to_dest[dest][i].v);
                }
                int array_size2 = sizeof(lks2) / sizeof(lks2[0]);
                int array_size3 = sizeof(lks3) / sizeof(lks3[0]);
                current_graph[dest].set_number_vertices(all_perimeter_graph_nodes.size());
                current_graph_dest_to_dest[dest].set_number_vertices(all_dest_to_dest_graph_nodes.size());
                for (int i = 0; i < array_size2; i++)
                {
                    current_graph[dest].add_link(lks2[i].u, lks2[i].v, lks2[i].weight);
                    if (debug)
                        std::cout << "\tAdded DJKSTR Graph node[" << lks2[i].u << "," << lks2[i].v << "],w=" << lks2[i].weight << endl;
                }
                for (int i = 0; i < array_size3; i++)
                {
                    current_graph_dest_to_dest[dest].add_link(lks3[i].u, lks3[i].v, lks3[i].weight);
                }
                if (debug)
                    std::cout << "DJKSTR Graph size:" << array_size2 << ",origin:" << Destinations[dest] << endl;
                current_graph[dest].setv();
                current_graph_dest_to_dest[dest].setv();
                //graph_dests.push_back(make_unique<Graph>(current_graph));
                //graph_dests.push_back(&current_graph);
                DijkstraShortestPathAlg dijkstra_alg(&(current_graph[dest]));
                Dijkstra_algs.push_back(dijkstra_alg);
                DijkstraShortestPathAlg dijkstra_dest_to_dest(&(current_graph_dest_to_dest[dest]));
                Dijkstra_algs_dest_to_dest.push_back(dijkstra_alg);
            }
            //Calculate all minimum distances, give us best possible undisclosing factor
            //if(debug){
            std::cout << "Paths from dest:" << Destinations[dest] << " to alternative destinations:" << endl;
            std::cout << "____________________________" << endl;
            //}

            for (size_t dest2 = 0; dest2 < Destinations.size(); dest2++)
            {
                if (dest2 == dest)
                    continue;
                BasePath *result =
                    Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
                        current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(Destinations[dest2]));
                if (debug)
                    std::cout << "\t" << result->get_path_string() << ",Weight:" << result->Weight() << endl;
                else
                    std::cout << "," << result->Weight();
                min_rho = min(min_rho, result->Weight());
                max_rho = max(max_rho, result->Weight());
            }
            std::cout << endl
                 << "____________________________" << endl;
        }

        Dijkstra_algs_dest_to_dest.clear();
        current_graph_dest_to_dest.clear();
        std::cout << "Destinations calculated" << endl;

        float upper_disc_dist = (min_rho / 2.0);
        if (debug)
        {
            //std::cout<<FORERED<<RESETTEXT<<"First finished,"<<FOREGRN<<RESETTEXT<<"rho_min=,"<<min_rho-1<<",UPPER DISCLOSING DISTANCE:"<<(min_rho/2)-1<<RESETTEXT<<endl;
            std::cout << "FIRST finished,"
                 << "rho_min=," << min_rho << ",UPPER DISCLOSING DISTANCE:" << (min_rho / 2) << endl;
        }
        else
        {
            //std::cout<<FORERED<<RESETTEXT<<"Minimum geodetic distance between destinations(rho_min)=,"<<min_rho-1<<",UPPER DISCLOSING DISTANCE:"<<upper_disc_dist<<RESETTEXT<<endl;
            std::cout << "FIRST Minimum geodetic distance between destinations(rho_min)=," << min_rho << ",UPPER DISCLOSING DISTANCE:," << upper_disc_dist << ",max_rho=" << max_rho << endl;
        }
        //std::cout<<FORERED<<RESETTEXT<<"Second, get P(d) for all Destinations:"<<Destinations<<", for boundary:"<<upper_disc_dist+1<<RESETTEXT<<endl;
        std::cout << "____________________________" << endl;
        if (input_lambda == 0)
        {
            input_lambda = upper_disc_dist + 1.0 + lambda_add;
        }
        else if (input_lambda < upper_disc_dist)
        {
            cerr << "lambda cannot be smaller than best possible lambda:" << upper_disc_dist << endl;
            exit(1);
        }
        else
        {
            std::cout << "New input_lambda=," << input_lambda << endl;
        }
        std::cout << "SECOND, get P(d) for all Destinations:" << Destinations << ", for boundary:" << input_lambda << ",time:," << ms << endl;

        //map<int,shared_ptr<BasePath> > current_perim_paths;
        //vector<map<int,shared_ptr<BasePath> > > perim_paths;
        map<int, DijkstraShortestPathAlg> Dijkstra_bifurcations_map;
        for (size_t dest = 0; dest < Destinations.size(); dest++)
        {
            Dijkstra_algs[dest].determine_all_shortest_paths(current_graph[dest].get_vertex(Destinations[dest]), input_lambda);
            Dijkstra_algs[dest].get_perim_nodes(&perimeter_nodes_vect[dest], input_lambda);
            Dijkstra_algs[dest].get_perim_nodes(&neighourhood_nodes_vect[dest], 0);
            //STOCASTIC CODE
            int reachable_nodes = 0;
            for (auto node : nodes_set)
            {
                if (!main_dijkstra_alg->is_node_reachable(node))
                {
                    //std::cout<<"node "<<node<<" is not reachable from start"<<endl;
                    continue;
                }
                else
                {
                    reachable_nodes++;
                    auto cost_to_dest = Dijkstra_algs[0].get_cost(node);
                    if (cost_to_dest == 0)
                    {
                        continue;
                    }
                    //std::cout<<"reachable node:"<<node<<",cost to destination:"<<Dijkstra_algs[0].get_cost(node)<<endl;
                    //Feasible analysis, for now initial cost of all nodes is 0
                    //auto budget=40;

                    //std::cout<<"\t feasible for budget=40"<<endl;
                }
            }
            std::cout << "Reachable nodes:" << reachable_nodes << endl;
            //std::cout<<FOREGRN<<RESETTEXT<<"Dest["<<Destinations[dest]<<"],Nd:";
            std::cout << "Dest[" << Destinations[dest] << "],Nd_size:" << perimeter_nodes_vect[dest].size() << endl;
            /*for(auto it : perimeter_nodes_vect[dest]){
          std::cout<<","<<it;
        }*/
            if (perimeter_nodes_vect[dest].size() == 0)
            {
                cerr << "Aborting, destination:" << Destinations[dest] << " has no perimeter nodes!!!" << endl;
                exit(101);
            }
            perimeter_nodes += perimeter_nodes_vect[dest].size();
            neighourhood_nodes += neighourhood_nodes_vect[dest].size();
            //std::cout<<FOREGRN<<RESETTEXT<<"\nDest["<<Destinations[dest]<<"],Nd|:";
            std::cout << "\nDest[" << Destinations[dest] << "],Nd|:"
                 << ",Nd|_size:" << neighourhood_nodes_vect[dest].size() << endl;
            //for(auto it : neighourhood_nodes_vect[dest]){
            //  std::cout<<","<<it;
            //}
            //std::cout<<endl<<RESETTEXT;

            //std::cout<<"hola6"<<flush<<endl;
            //Dijkstra_algs[dest].get_perim_paths(current_perim_paths);
            //perim_paths.push_back(current_perim_paths);
            std::cout << "____________________________" << endl;
            if (debug)
                Dijkstra_algs[dest].print_paths();
        }
        //std::cout<<"hola"<<endl;
        /*for(auto it : perim_paths){
        for (auto it2 : it){
          std::cout<<"\t";(*it2.second).PrintOut(std::cout);
        }
      }*/

        if (!SaT_map)
        { //Already done in create_grid_nodes_and_edges_from_SaT_file
            std::cout << "getting shortest destination path to all nodes from origin in main graph" << endl;
            main_dijkstra_alg->determine_all_shortest_paths(main_graph.get_vertex(start), INT_MAX);
            std::cout << "finished getting shortest paths to all nodes from origin in the main graph." << endl;
        }
        //Calculate max_lambda
        int max_lambda = 0;
        vector<int> Optimal_path_dests;

        /*for(auto i=1;i<100;i++){
        if(FO_nodes.find(i)==FO_nodes.end()){
          continue;
        }
        auto path=main_dijkstra_alg->recover_shortest_perim_path(i);
        std::cout<<"Cost["<<i<<"]="<<path->Weight()<<endl;
        std::cout<<"Path["<<i<<"]="; path->PrintOut(std::cout);
      }
        
      auto path20=main_dijkstra_alg->recover_shortest_perim_path(20);
      std::cout<<"Cost[20]="<<path20->Weight()<<endl;
      auto path60=main_dijkstra_alg->recover_shortest_perim_path(60);
      std::cout<<"Cost[60]="<<path60->Weight()<<endl;
      auto path91=main_dijkstra_alg->recover_shortest_perim_path(91);
      std::cout<<"Cost[91]="<<path91->Weight()<<endl;
      auto path96=main_dijkstra_alg->recover_shortest_perim_path(96);
      std::cout<<"Cost[96]="<<path96->Weight()<<endl;
      exit(1);*/

        for (auto it : all_destinations)
        {
            auto path = main_dijkstra_alg->recover_shortest_perim_path(it);
            std::cout << "Destination:," << it << ",optimal_path_length:," << path->length() << ",Cost:" << path->Weight() << endl;
            path->PrintOut(std::cout);
            Optimal_path_dests.push_back(path->Weight());
            if (max_lambda < path->Weight())
            {
                max_lambda = path->Weight();
            }
        }
        if (input_lambda > max_lambda)
        {
            cerr << "lambda cannot be higher than max_cost_path:" << max_lambda << endl;
            exit(1);
        }

        shortest_paths_determined.insert(start);
        set<shared_ptr<BasePath>> Final_paths;

        std::cout << "____________Pd________________" << endl;
        for (size_t dest = 0; dest < Destinations.size(); dest++)
        {
            //Dijkstra_algs[dest].print_paths();

            //std::cout<<FOREGRN<<RESETTEXT<<"P(d)["<<Destinations[dest]<<"]:"<<endl;
            std::cout << "P(d)[" << Destinations[dest] << "]:" << endl;

            //std::map<int,BasePath*> *current_perim_pahts;
            //std::cout<<FOREGRN<<RESETTEXT;Dijkstra_algs[dest].make_full_paths(main_dijkstra_alg->get_pt_paths());
            //Dijkstra_algs[dest].make_full_paths(main_dijkstra_alg->get_pt_paths());
            Dijkstra_algs[dest].print_perim_paths();
            //std::cout<<RESETTEXT;
            //Dijkstra_algs[dest].populate_final_paths(&Final_paths);
            //final_path1.append_to_path(tail1);
        }
        std::cout << "____________________________" << endl;
        //Time up to Pd, including generating the full list of paths from origin
        end_time = chrono::high_resolution_clock::now();
        diff = end_time - orig_start_time;
        ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();

        auto temp_start_time = chrono::high_resolution_clock::now();
        //Thirdly, Find crossings & also check if no destination has crossings
        set<int> paired_destinations;
        set<int> unpaired_destinations;
        //std::cout<<FORERED<<RESETTEXT<<"THIRD, get Nd,d' for min_rho+1"<<RESETTEXT<<endl;
        std::cout << "THIRD, get Nd,d' for min_rho,time:," << ms << endl;
        map<pair<int, int>, set<int>> Ndd;
        map<pair<int, int>, vector<pair<shared_ptr<BasePath>, shared_ptr<BasePath>>>> Pdd;
        map<shared_ptr<BasePath>, set<shared_ptr<BasePath>>> path_connections;
        for (size_t dest1 = 0; dest1 < Destinations.size(); dest1++)
        {
            int best_bifurc_node = start;
            for (size_t dest2 = 0; dest2 < Destinations.size(); dest2++)
            {
                if (dest2 <= dest1)
                    continue;
                set<int> intersect;
                set_intersection(neighourhood_nodes_vect[dest1].begin(), neighourhood_nodes_vect[dest1].end(),
                                 neighourhood_nodes_vect[dest2].begin(), neighourhood_nodes_vect[dest2].end(),
                                 std::inserter(intersect, intersect.begin()));
                if (intersect.size() > 0)
                {
                    bifurcated_destinations++;
                    //std::cout<<FOREGRN<<RESETTEXT<<"N["<<Destinations[dest1]<<","<<Destinations[dest2]<< "] has the following intersections:";
                    std::cout << "N[" << Destinations[dest1] << "," << Destinations[dest2] << "] has the following intersections:";
                    for (auto it : intersect)
                        std::cout << "," << it; //std::cout<<RESETTEXT<<endl;
                    std::cout << endl;
                    paired_destinations.insert(Destinations[dest1]);
                    paired_destinations.insert(Destinations[dest2]);
                    Ndd[make_pair(Destinations[dest1], Destinations[dest2])] = intersect;
                    total_pairs += intersect.size();
                }
                else
                {
                    //std::cout<<FOREGRN<<RESETTEXT<<"N(d,d')["<<Destinations[dest1]<<","<<Destinations[dest2]<< "] has no intersections."<<endl;
                    std::cout << "N(d,d')[" << Destinations[dest1] << "," << Destinations[dest2] << "] has no intersections." << endl;
                }
                auto end_time = chrono::high_resolution_clock::now();
                auto diff = end_time - temp_start_time;
                ms_pairs += std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();

                int shortest_bifurc_dist = INT_MAX;
                BasePath *best_bifur_path;

                //DijkstraShortestPathAlg temp_dijkstra_alg(&main_graph);
                temp_start_time = chrono::high_resolution_clock::now();
                for (auto it_intersect : intersect)
                {
                    if (shortest_paths_determined.find(it_intersect) == shortest_paths_determined.end())
                    { //in case bifurcation is not new
                        Dijkstra_bifurcations_map.insert(make_pair(it_intersect, DijkstraShortestPathAlg(&main_graph)));
                        if (debug)
                            std::cout << "getting shortest destination path to all nodes from bifurcation:" << it_intersect << endl;
                        Dijkstra_bifurcations_map[it_intersect].determine_all_shortest_paths(main_graph.get_vertex(it_intersect), INT_MAX);
                        shortest_paths_determined.insert(it_intersect);
                    }
                }
                end_time = chrono::high_resolution_clock::now();
                diff = end_time - temp_start_time;
                ms_bifur_backward += std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
                //Final_paths.push_back(best_bifur_path);
                //int bifurcation_counter=0;
                for (auto it_intersect : intersect)
                {
                    //std::cout<<"\tshortest path from origin node:,"<<start<<",to bifurcation node:"<<it_intersect<<",for destination pair:"<<Destinations[dest1]<<","<<Destinations[dest2]<<endl;
                    //main_dijkstra_alg->clear();
                    shared_ptr<BasePath> head_path =
                        main_dijkstra_alg->recover_shortest_perim_path(it_intersect);
                    //std::cout<<"head_path:"<<head_path->get_path_string()<<endl;
                    shared_ptr<BasePath> tail1 =
                        Dijkstra_algs[dest1].recover_shortest_perim_path(it_intersect);
                    //std::cout<<"tail1:"<<tail1->get_path_string()<<endl;
                    shared_ptr<BasePath> tail2 =
                        Dijkstra_algs[dest2].recover_shortest_perim_path(it_intersect);
                    //std::cout<<"tail2:"<<tail1->get_path_string()<<endl;

                    BasePath path1 = *head_path;
                    path1.append_to_path(tail1);

                    BasePath path2 = *head_path;
                    path2.append_to_path(tail2);
                    //std::cout<<"path1:"<<path1.get_path_string()<<endl;
                    //std::cout<<"path2:"<<path2.get_path_string()<<endl;

                    //std::cout<<"\t";head_path->PrintOut(std::cout);
                    //std::cout<<"\t";path1.PrintOut(std::cout);
                    //std::cout<<"\t";path2.PrintOut(std::cout);
                    auto path1_ptr = make_shared<BasePath>(path1);
                    auto path2_ptr = make_shared<BasePath>(path2);
                    Pdd[make_pair(Destinations[dest1], Destinations[dest2])].push_back(make_pair(path1_ptr, path2_ptr));
                    /*path_connections[path1_ptr].insert(path2_ptr);
        path_connections[path2_ptr].insert(path1_ptr);*/
                    /*if(shortest_bifurc_dist>result->length()){
          shortest_bifurc_dist=result->length();
          best_bifur_path=result;
          best_bifurc_node=it_intersect;
        }*/
                }

                /*if(intersect.size()>0){
        std::cout<<endl<<"Shortest path from origin node:,"<<start<<",to bifurcation node:"<<best_bifurc_node<<",for destination pair:"<<Destinations[dest1]<<","<<Destinations[dest2]<<"is:";
        best_bifur_path->PrintOut(std::cout);std::cout<<endl;
        std::cout<<FORERED<<RESETTEXT"FOURTH, concatenation of origin-Pd,d'-PairedDestinations, final pair of paths:"<<RESETTEXT<<endl;
        //First paired path
        //shared_ptr<BasePath> tail1 = make_
        //  main_dijkstra_alg->get_shortest_path(
        //      main_graph.get_vertex(best_bifurc_node),main_graph.get_vertex(Destinations[dest1]));
        //BasePath final_path1=*best_bifur_path;
        //final_path1.append_to_path(tail1);
        //std::cout<<FOREGRN<<RESETTEXT"\t\t First paired path:";final_path1.PrintOut(std::cout);std::cout<<RESETTEXT;
        //Final_paths.push_back(final_path1);
        //Second paired path
        //BasePath* tail2 = 
        //  main_dijkstra_alg->get_shortest_path(
        //      main_graph.get_vertex(best_bifurc_node),main_graph.get_vertex(Destinations[dest2]));
        //BasePath final_path2=*best_bifur_path;
        //final_path2.append_to_path(tail2);
        //std::cout<<FOREGRN<<RESETTEXT<<"\t\t Second paired path:";final_path2.PrintOut(std::cout);std::cout<<RESETTEXT;
        //Final_paths.push_back(final_path2);
          }*/
            }
        }

        std::cout << "____________________________" << endl;
        for (auto it : Pdd)
        {
            for (auto it2 : it.second)
            {
                if (debug)
                {
                    std::cout << "P(d,d')[" << it.first.first << "," << it.first.second << "]" << endl;
                    //std::cout<<FOREGRN<<RESETTEXT<<"\tP(d,d')["<<it.first.first<<","<<it.first.second<<"],";
                    //need to get biggest path to find bifurcation node
                    /*if(it2.first->length()>=it2.second->length()){
          std::cout<<"bifurcation_node:"<<it2.first->get_path_string()<<flush;std::cout<<","<<it2.first->get_node(min(it2.first->length(),input_lambda));
        }
        else{
          std::cout<<"bifurcation_node:"<<it2.second->get_path_string()<<flush;std::cout<<","<<it2.first->get_node(min(it2.second->length(),input_lambda));
        }*/
                    //std::cout<<",path1:"<<it2.first->get_path_string()<<",path2:"<<it2.second->get_path_string()<<endl<<RESETTEXT;
                    if (debug)
                        std::cout << endl
                             << "\tpath1:" << flush << it2.first->get_path_string() << ",path2:" << it2.second->get_path_string() << endl;
                }
                //Final_paths.insert(it2.first);
                //Final_paths.insert(it2.second);

                //Get best pair per destination in terms of max_rel_cost for any path in the pair

                shared_ptr<BasePath> opt_path1 = main_dijkstra_alg->recover_shortest_perim_path(it.first.first);
                float opt_cost1 = opt_path1->Weight();
                float current_cost1 = it2.first->Weight();
                float rel_cost1 = current_cost1 / opt_cost1;

                shared_ptr<BasePath> opt_path2 = main_dijkstra_alg->recover_shortest_perim_path(it.first.second);
                float opt_cost2 = opt_path2->Weight();
                float current_cost2 = it2.second->Weight();
                float rel_cost2 = current_cost2 / opt_cost2;

                if (best_rel_cost_per_dest.find(it.first.first) == best_rel_cost_per_dest.end())
                {
                    best_rel_cost_per_dest[it.first.first].first = INT_MAX;
                }
                if (max(rel_cost1, rel_cost2) < best_rel_cost_per_dest[it.first.first].first)
                {
                    best_rel_cost_per_dest[it.first.first].first = max(rel_cost1, rel_cost2);
                    best_rel_cost_per_dest[it.first.first].second.clear();
                    best_rel_cost_per_dest[it.first.first].second.push_back(it2.first);
                    best_rel_cost_per_dest[it.first.first].second.push_back(it2.second);
                    if (debug)
                    {
                        std::cout << "added paired path:" << it2.first->get_path_string() << " for destination:," << it.first.first << ",new_best_val:," << max(rel_cost1, rel_cost2) << endl;
                        std::cout << "added paired path:" << it2.second->get_path_string() << " for destination:," << it.first.first << ",new_best_val:," << max(rel_cost1, rel_cost2) << endl;
                        std::cout << "size of paths for destination:," << it.first.first << ",is:," << best_rel_cost_per_dest[it.first.first].second.size();
                    }
                }

                if (best_rel_cost_per_dest.find(it.first.second) == best_rel_cost_per_dest.end())
                {
                    best_rel_cost_per_dest[it.first.second].first = INT_MAX;
                }

                if (max(rel_cost1, rel_cost2) < best_rel_cost_per_dest[it.first.second].first)
                {
                    best_rel_cost_per_dest[it.first.second].first = max(rel_cost1, rel_cost2);
                    best_rel_cost_per_dest[it.first.second].second.clear();
                    best_rel_cost_per_dest[it.first.second].second.push_back(it2.first);
                    best_rel_cost_per_dest[it.first.second].second.push_back(it2.second);
                    if (debug)
                    {
                        std::cout << "added paired path:" << it2.first->get_path_string() << " for destination:," << it.first.second << ",new_best_val:," << max(rel_cost1, rel_cost2) << endl;
                        std::cout << "added paired path:" << it2.second->get_path_string() << " for destination:," << it.first.second << ",new_best_val:," << max(rel_cost1, rel_cost2) << endl;
                        std::cout << "size of paths for destination:" << it.first.second << ",is:," << best_rel_cost_per_dest[it.first.second].second.size();
                    }
                }
            }
        }
        std::cout << "____________________________" << endl;
        //std::cout<<FORERED<<RESETTEXT<<"FOURTH, get P(d,d',d'')"<<RESETTEXT<<endl;
        end_time = chrono::high_resolution_clock::now();
        diff = end_time - orig_start_time;
        ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
        std::cout << "FOURTH, get P(d,d',d''),time:," << ms << endl;
        if (debug)
        {
            std::cout << "After P(d,d'), current Final paths:" << endl;
            int counter = 0;
            for (auto it : Final_paths)
            {
                std::cout << "P_lambda[" << counter++ << "]->" << it->get_path_string() << endl;
            }
        }
        //timings
        auto temp_end_time = chrono::high_resolution_clock::now();
        auto diff = temp_end_time - temp_start_time;
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
        std::cout << "Time(ms) for PDD:," << ms << endl;

        temp_start_time = chrono::high_resolution_clock::now();
        //Calculate shortest path to common nodes
        vector<BasePath *> Budget_optimized_final_paths;
        //Find shortest paths between singled out destinations and paired destinations
        set_difference(all_destinations.begin(), all_destinations.end(), paired_destinations.begin(),
                       paired_destinations.end(), std::inserter(unpaired_destinations, unpaired_destinations.begin()));
        map<pair<int, pair<int, int>>, vector<tuple<shared_ptr<BasePath>, shared_ptr<BasePath>, shared_ptr<BasePath>>>> Pddd;
        for (auto it_all_dests : all_destinations)
        {
            //if(debug)
            std::cout << "Working on dest:" << it_all_dests << endl; //std::cout<<",neighourhood_nodes:"<<neighourhood_nodes_vect<<endl;
            shared_ptr<BasePath> opt_dest_path =
                main_dijkstra_alg->recover_shortest_perim_path(it_all_dests);
            if (opt_dest_path->Weight() <= input_lambda)
            {
                std::cout << "skipping destination:" << it_all_dests << ",optimal_path_weight:" << opt_dest_path->Weight() << ",lambda:" << input_lambda << "from d in tripple d,d',d'' because it does not need cover" << endl;
                //So best path for destination is any optimal path
                best_rel_cost_per_dest[it_all_dests].first = 1.0;
                best_rel_cost_per_dest[it_all_dests].second.clear();
                best_rel_cost_per_dest[it_all_dests].second.push_back(opt_dest_path);
                continue;
            }

            for (auto it : perimeter_nodes_vect[dest_index[it_all_dests]])
            {
                if (!main_dijkstra_alg->is_node_reachable(it))
                {
                    std::cout << "skipping perimeter node " << it << " because is is not reachable from origin." << endl;
                    continue;
                }

                int h1 = it;
                if (h1 == start)
                    continue;
                //NEED TO ENSURE NEITHER H1 OR H2 ARE DESTINATIONS
                //CAN HAPPEN WHEN WE ALLOW CYCLES
                if (all_destinations.find(h1) != all_destinations.end())
                {
                    std::cout << "\tskipping destination_node:" << h1 << endl;
                    continue;
                }

                if (debug)
                    std::cout << "\th1:" << h1 << endl;
                shared_ptr<BasePath> head_path =
                    main_dijkstra_alg->recover_shortest_perim_path(it);
                if (debug)
                    std::cout << "\thead_path:" << head_path->get_path_string() << endl;
                BasePath FinalPath0 = *head_path;
                FinalPath0.append_to_path(Dijkstra_algs[dest_index[it_all_dests]].recover_shortest_perim_path(h1));
                auto path_ptr0 = make_shared<BasePath>(FinalPath0);
                if (debug)
                    std::cout << "\t\t\t\tFinalPath1" << FinalPath0.get_path_string() << endl;

                for (auto it2 : Ndd)
                {
                    if (it2.first.first == it_all_dests || it2.first.second == it_all_dests)
                        continue; //d neq d' nor d''
                    if (debug)
                        std::cout << "\t\t connecting to Ndd(" << it2.first.first << "," << it2.first.second << ")" << endl;
                    for (auto it3 : it2.second)
                    {
                        int h2 = it3;
                        //Origin cannot be h2, no incoming edges!
                        if (h2 == start)
                            continue;
                        //NEED TO ENSURE NEITHER H1 OR H2 ARE DESTINATIONS
                        //CAN HAPPEN WHEN WE ALLOW CYCLES
                        if (all_destinations.find(h2) != all_destinations.end())
                            continue;
                        if (debug)
                            std::cout << "\t\t\t"
                                 << ",h1:," << h1 << ",h2:," << h2 << ",h1->h2 path:";
                        if (!Dijkstra_bifurcations_map[h2].is_node_reachable(h1))
                        {
                            //std::cout<<"\tskipping h2node "<<h2<<" because is is not reachable from h1"<<h1<<endl;
                            continue;
                        }
                        shared_ptr<BasePath> h1_h2_path;
                        if (h1 != h2)
                        {
                            //if(h1!=start){
                            h1_h2_path = Dijkstra_bifurcations_map[h2].recover_shortest_perim_path(h1);
                            if (debug)
                                std::cout << h1_h2_path->get_path_string() << endl;
                        }
                        //}
                        BasePath FinalPath1 = *head_path;
                        //if(h1!=h2&&h1!=start){
                        if (h1 != h2)
                        {
                            FinalPath1.append_to_path(h1_h2_path);
                        }
                        FinalPath1.append_to_path(Dijkstra_algs[dest_index[it2.first.first]].recover_shortest_perim_path(h2));
                        /*}
            else{
              FinalPath1.append_to_path(Dijkstra_bifurcations_map[it2.first.first].recover_shortest_perim_path(destination));
              FinalPath1.append_to_path(main_dijkstra_alg->recover_shortest_perim_path(h2));
            }*/

                        if (debug)
                            std::cout << "\t\t\t\tFinalPath2" << FinalPath1.get_path_string() << endl;
                        BasePath FinalPath2 = *head_path;
                        //if(h1!=h2&&h1!=start){
                        if (h1 != h2)
                        {
                            FinalPath2.append_to_path(h1_h2_path);
                        }
                        FinalPath2.append_to_path(Dijkstra_algs[dest_index[it2.first.second]].recover_shortest_perim_path(h2));
                        /*else{
              FinalPath2.append_to_path(Dijkstra_bifurcations_map[it2.first.second].recover_shortest_perim_path(destination));
            }*/
                        if (debug)
                            std::cout << "\t\t\t\tFinalPath3" << FinalPath2.get_path_string() << endl;
                        pair<int, int> dest_pair(it2.first.first, it2.first.second);

                        auto path_ptr1 = make_shared<BasePath>(FinalPath1);
                        auto path_ptr2 = make_shared<BasePath>(FinalPath2);

                        tuple<shared_ptr<BasePath>, shared_ptr<BasePath>, shared_ptr<BasePath>> Final_path_triple = make_tuple(path_ptr0, path_ptr1, path_ptr2);
                        //std::cout<<"\t FinalPath0 inserted:"<<path_ptr0->get_path_string()<<endl;
                        //std::cout<<"\t FinalPath1 inserted:"<<path_ptr1->get_path_string()<<endl;
                        //std::cout<<"\t FinalPath2 inserted:"<<path_ptr2->get_path_string()<<endl;

                        Pddd[make_pair(it_all_dests, dest_pair)].push_back(Final_path_triple);
                        total_triples++;
                        /*path_connections[path_ptr0].insert(path_ptr1);
            path_connections[path_ptr0].insert(path_ptr2);
            path_connections[path_ptr1].insert(path_ptr0);
            path_connections[path_ptr1].insert(path_ptr2);
            path_connections[path_ptr2].insert(path_ptr0);
            path_connections[path_ptr2].insert(path_ptr1);*/

                        //Get best pair per destination in terms of max_rel_cost for any path in the pair
                        shared_ptr<BasePath> opt_path1 = main_dijkstra_alg->recover_shortest_perim_path(it_all_dests);
                        float opt_cost1 = opt_path1->Weight();
                        float current_cost1 = path_ptr0->Weight();
                        float rel_cost1 = current_cost1 / opt_cost1;

                        shared_ptr<BasePath> opt_path2 = main_dijkstra_alg->recover_shortest_perim_path(dest_pair.first);
                        float opt_cost2 = opt_path2->Weight();
                        float current_cost2 = path_ptr1->Weight();
                        float rel_cost2 = current_cost2 / opt_cost2;

                        float max_cost = max(rel_cost1, rel_cost2);

                        shared_ptr<BasePath> opt_path3 = main_dijkstra_alg->recover_shortest_perim_path(dest_pair.second);
                        float opt_cost3 = opt_path3->Weight();
                        float current_cost3 = path_ptr2->Weight();
                        float rel_cost3 = current_cost3 / opt_cost3;

                        max_cost = max(rel_cost3, max_cost);

                        if (best_rel_cost_per_dest.find(it_all_dests) == best_rel_cost_per_dest.end())
                        {
                            best_rel_cost_per_dest[it_all_dests].first = INT_MAX;
                        }

                        if (max_cost < best_rel_cost_per_dest[it_all_dests].first ||
                            best_rel_cost_per_dest[it_all_dests].first == 0)
                        {
                            best_rel_cost_per_dest[it_all_dests].first = max_cost; //update best value
                            best_rel_cost_per_dest[it_all_dests].second.clear();
                            best_rel_cost_per_dest[it_all_dests].second.push_back(path_ptr0);
                            best_rel_cost_per_dest[it_all_dests].second.push_back(path_ptr1);
                            best_rel_cost_per_dest[it_all_dests].second.push_back(path_ptr2);
                            if (debug)
                            {
                                std::cout << "added tripple paths:" << path_ptr0->get_path_string() << " for single destination:," << it_all_dests << ",new_best_val:," << max_cost << endl;
                                std::cout << "added tripple paths:" << path_ptr1->get_path_string() << " for single destination:," << it_all_dests << ",new_best_val:," << max_cost << endl;
                                std::cout << "added tripple paths:" << path_ptr2->get_path_string() << " for single destination:," << it_all_dests << ",new_best_val:," << max_cost << endl;
                                std::cout << "size of paths for destination:" << it_all_dests << ",is:," << best_rel_cost_per_dest[it_all_dests].second.size();
                            }
                        }

                        /*if(max_cost<best_rel_cost_per_dest[dest_pair.first].first){
              best_rel_cost_per_dest[dest_pair.first].first=max_cost;//update best value
              best_rel_cost_per_dest[dest_pair.first].second.clear();
              best_rel_cost_per_dest[dest_pair.first].second.push_back(path_ptr0);
              best_rel_cost_per_dest[dest_pair.first].second.push_back(path_ptr1);
              best_rel_cost_per_dest[dest_pair.first].second.push_back(path_ptr2);
              //if(debug){
            //std::cout<<"added tripple paths:"<<path_ptr0->get_path_string()<<" for paired destination:,"<<dest_pair.first<<",new_best_val:,"<<max_cost<<endl;
            //std::cout<<"added tripple paths:"<<path_ptr1->get_path_string()<<" for paired destination:,"<<dest_pair.first<<",new_best_val:,"<<max_cost<<endl;
            //std::cout<<"added tripple paths:"<<path_ptr1->get_path_string()<<" for paired destination:,"<<dest_pair.first<<",new_best_val:,"<<max_cost<<endl;
              //}
            }
            
            if(max_cost<best_rel_cost_per_dest[dest_pair.second].first){
              best_rel_cost_per_dest[dest_pair.second].first=max_cost;//update best value
              best_rel_cost_per_dest[dest_pair.second].second.clear();
              best_rel_cost_per_dest[dest_pair.second].second.push_back(path_ptr0);
              best_rel_cost_per_dest[dest_pair.second].second.push_back(path_ptr1);
              best_rel_cost_per_dest[dest_pair.second].second.push_back(path_ptr2);
              //if(debug){
            //std::cout<<"added tripple paths:"<<path_ptr0->get_path_string()<<" for paired destination:,"<<dest_pair.second<<",new_best_val:,"<<max_cost<<endl;
            //std::cout<<"added tripple paths:"<<path_ptr1->get_path_string()<<" for paired destination:,"<<dest_pair.second<<",new_best_val:,"<<max_cost<<endl;
            //std::cout<<"added tripple paths:"<<path_ptr1->get_path_string()<<" for paired destination:,"<<dest_pair.second<<",new_best_val:,"<<max_cost<<endl;
              //}
            }*/
                    }
                }
            }
        }
        //timings
        temp_end_time = chrono::high_resolution_clock::now();
        diff = temp_end_time - temp_start_time;
        ms_triples += std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
        std::cout << "Time(ms) for PDDD:," << ms << ",triples:" << total_triples << endl;

        //Now simply combine all paths in best_rel_cost_per_destination
        float max_rel_cost = 1.0;
        for (auto it : best_rel_cost_per_dest)
        {
            max_rel_cost = max(max_rel_cost, it.second.first);
            for (auto it2 : it.second.second)
            {
                Final_paths.insert(it2);
            }
        }

        int counter = 0;
        for (auto it : Final_paths)
        {
            shared_ptr<BasePath> opt_path = main_dijkstra_alg->recover_shortest_perim_path(it->get_node(0));
            float opt_cost = opt_path->Weight();
            float current_cost = it->Weight();
            float rel_cost = current_cost / opt_cost;
            std::cout << ",rel_cost:," << rel_cost << ",P_lambda[0," << counter++ << "]->" << it->get_path_string() << endl;
            for (size_t i = 0; i < it->length(); i++)
            {
                int id = it->get_node(i);
                auto pos = node_map[id];
                if (map_display)
                {
                    color_map[pos.first][pos.second][2] = 1;
                }
            }
        }

        counter = 0;
        pair<int, int> old_pos; /*
      for(auto it : Final_paths){
        std::cout<<"final_path_pos["<<counter++<<"]";
        for(size_t i=0;i<it->length();i++){
          int id=it->get_node(i);
          auto pos=node_map[id];
          std::cout<<"["<<pos.first<<","<<pos.second<<"],";
          if(i>1){
        if(abs(pos.first-old_pos.first)+abs(pos.second-old_pos.second)>1){
          std::cout<<"pos:"<<pos.first<<","<<pos.second<<",old_pos:"<<old_pos.first<<","<<old_pos.second<<endl;
          exit(1);
        }
          }
          old_pos=pos;
        }
        std::cout<<endl;
      }*/
        //timings
        end_time = chrono::high_resolution_clock::now();
        diff = end_time - start_time;
        start_time = end_time; //restart
        ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();

        std::cout << "N:," << N << ",R:," << random_positions << ",min_dist_to_dest:," << min_dist_to_dest << ",max_dist_to_dest:," << max_dist_to_dest << ",max_rel_cost:," << max_rel_cost << ",";
        std::cout << "final_paths:," << Final_paths.size() << ",overall_time:," << ms << ",ms_bifur_backwards:," << ms_bifur_backward << ",ms_pairs:," << ms_pairs << ",ms_triples:," << ms_triples;
        std::cout << ",boundary_nodes:," << perimeter_nodes << ",bifurcations:," << bifurcated_destinations << ",total_perim_nodes:," << perimeter_nodes << ",total_neighourhood_nodes:," << neighourhood_nodes << ",total_pairs:," << total_pairs << ",total_triples:" << total_triples;
        std::cout << ",labmda_star:," << upper_disc_dist << ",input_lambda:," << input_lambda << ",max_lambda:," << max_lambda;
        std::cout << ",Destination OP lengths:," << Optimal_path_dests << ",connectivity:," << connectivity << ",start:," << start;
        std::cout << ",Total_nodes:," << nodes_set.size() << ",reachable_nodes:," << main_dijkstra_alg->get_perim_size() << endl;

        if (!map_display)
            exit(0);

        //Destination and origin map positions in case no random selection
        for (auto dest : Destinations)
        {
            color_map[node_map[dest].first][node_map[dest].second][1] = 1;
        }
        color_map[node_map[start].first][node_map[start].second][0] = 1;

        int size = max_x * max_y;
        for (int i = 0; i < max_x; i++)
        {
            for (int j = 0; j < max_y; j++)
            {
                pair<int, int> current_pos = make_pair(i, j);
                if (coord_map2.find(current_pos) == coord_map2.end())
                { //not passable
                    //std::cout<<"@@";
                    std::cout << "@";
                }
                else
                {
                    //auto current_node=coord_map2[current_pos];
                    if (color_map[i][j][0] == 1)
                    { //origin node;
                        //std::cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
                        std::cout << FOREGRN << "\u25A0" << RESETTEXT;
                        //std::cout<<FOREGRN<<"."<<RESETTEXT;
                        continue;
                    }
                    else
                    { //origin and neighourhood
                        //ORIGIN
                        if (size > 100000)
                        {
                            if (i < max_x - 1)
                            {
                                if (j < max_y - 1 && color_map[i + 1][j + 1][0] == 1)
                                { //origin node;
                                    //std::cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FOREGRN << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                                else if (j > 0 && color_map[i + 1][j - 1][0] == 1)
                                { //origin node;
                                    //std::cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FOREGRN << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                                else if (color_map[i + 1][j][0] == 1)
                                { //origin node;
                                    //std::cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FOREGRN << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                            }
                            if (i > 0)
                            {
                                if (color_map[i - 1][j][0] == 1)
                                { //origin node;
                                    //std::cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FOREGRN << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                                else if (j > 0 && color_map[i - 1][j - 1][0] == 1)
                                { //origin node;
                                    //std::cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FOREGRN << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                                else if (j < max_y - 1 && color_map[i - 1][j + 1][0] == 1)
                                { //origin node;
                                    //std::cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FOREGRN << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                            }
                        }

                        //DESTINATIONS
                        if (color_map[i][j][1] == 1)
                        { //destination node
                            //std::cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
                            std::cout << FORERED << "\u25A0" << RESETTEXT;
                            //std::cout<<FORERED<<"."<<RESETTEXT;
                            continue;
                        }
                        if (size > 100000)
                        {
                            if (j > 0)
                            {
                                if (color_map[i][j - 1][0] == 1)
                                { //destinaiton node
                                    //std::cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FOREGRN << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                            }
                            if (j < max_y - 1)
                            {
                                if (color_map[i][j + 1][0] == 1)
                                { //destination node;
                                    //std::cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FOREGRN << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                            }
                            if (i < max_x - 1)
                            {
                                if (j < max_y - 1 && color_map[i + 1][j + 1][1] == 1)
                                { //destination node;
                                    //std::cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FORERED << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                                else if (j > 0 && color_map[i + 1][j - 1][1] == 1)
                                { //destination node;
                                    //std::cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FORERED << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                                else if (color_map[i + 1][j][1] == 1)
                                { //destination node;
                                    //std::cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FORERED << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                            }
                            if (i > 0)
                            {
                                if (color_map[i - 1][j][1] == 1)
                                { //destination node;
                                    //std::cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FORERED << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                                else if (j > 0 && color_map[i - 1][j - 1][1] == 1)
                                { //destination node;
                                    //std::cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FORERED << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                                else if (j < max_y - 1 && color_map[i - 1][j + 1][1] == 1)
                                { //destination node;
                                    //std::cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FORERED << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                            }
                            if (j > 0)
                            {
                                if (color_map[i][j - 1][1] == 1)
                                { //destinaiton node
                                    //std::cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FORERED << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                            }
                            if (j < max_y - 1)
                            {
                                if (color_map[i][j + 1][1] == 1)
                                { //destination node;
                                    //std::cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
                                    std::cout << FORERED << "\u25A0" << RESETTEXT;
                                    continue;
                                }
                            }
                        }
                        //FINAL PATHS
                        if (color_map[i][j][2] == 1)
                        { //path node
                            //std::cout<<FOREBLU<<"\u25A0\u25A0"<<RESETTEXT;
                            //std::cout<<FOREBLU<<"\u25A0"<<RESETTEXT;
                            std::cout << FOREBLU << "*" << RESETTEXT;
                            continue;
                        }
                        //std::cout<<"..";
                        std::cout << ".";
                    }
                }
            }
            std::cout << endl;
        }
        exit(0);
    }
    //Need to make copy of paths before erasing any of them
    PathMatrix PM_orig = *PM;
    int t_max = optimizing_target_goal_dist(PM);
    optimizing_observer(&PM_orig, nodes, t_max);
}

/////////////////////////
// auxiliary types for A*
struct location
{
  float y, x; // lat, long
  location(float x1,float y1) 
        {
            x=x1;
            y=y1;
        }
};
typedef float cost;

template <class Name, class LocMap>
class city_writer {
public:
  city_writer(Name n, LocMap l, float _minx, float _maxx,
              float _miny, float _maxy,
              unsigned int _ptx, unsigned int _pty)
    : name(n), loc(l), minx(_minx), maxx(_maxx), miny(_miny),
      maxy(_maxy), ptx(_ptx), pty(_pty) {}
  template <class Vertex>
  void operator()(ostream& out, const Vertex& v) const {
    float px = 1 - (loc[v].x - minx) / (maxx - minx);
    float py = (loc[v].y - miny) / (maxy - miny);
    out << "[label=\"" << name[v] << "\", pos=\""
        << static_cast<unsigned int>(ptx * px) << ","
        << static_cast<unsigned int>(pty * py)
        << "\", fontsize=\"11\"]";
  }
private:
  Name name;
  LocMap loc;
  float minx, maxx, miny, maxy;
  unsigned int ptx, pty;
};

template <class WeightMap>
class time_writer {
public:
  time_writer(WeightMap w) : wm(w) {}
  template <class Edge>
  void operator()(ostream &out, const Edge& e) const {
    out << "[label=\"" << wm[e] << "\", fontsize=\"11\"]";
  }
private:
  WeightMap wm;
};


// euclidean distance heuristic
template <class Graph, class CostType, class LocMap>
class distance_heuristic : public boost::astar_heuristic<Graph, CostType>
{
public:
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  distance_heuristic(LocMap l, Vertex goal)
    : m_location(l), m_goal(goal) {}
  CostType operator()(Vertex u)
  {
    CostType dx = m_location[m_goal].x - m_location[u].x;
    CostType dy = m_location[m_goal].y - m_location[u].y;
    return ::sqrt(dx * dx + dy * dy);
  }
private:
  LocMap m_location;
  Vertex m_goal;
};
/////ASTAR declarations finished
void read_Astar_Graph(vector<pair<int,int> > &edges,vector<cost> &weights,vector<location> &locations)
{ 
    typedef std::pair<int, int> edge;
    std::cout << "hello read_Astar_Graph" << flush << endl;
    std::ifstream in(input_filename.c_str());//in is for graph edges and weights
    auto input_filename2=input_filename;
    input_filename2.insert(input_filename.length()-4,"_Labels");//extension should always be csv
    //std::cout<<input_filename2<<endl;
    std::ifstream in2(input_filename2.c_str());//in2 is for coordinates file
    
    if(!in.good()){
        cerr<<"file:"<<input_filename<<" does not exist!"<<endl;
        exit(1);
    }
    else if(!in2.good()){
        cerr<<"file:"<<input_filename2<<" does not exist!"<<endl;
        exit(1);
    }

    string str;
    total_nodes = 0;
    int total_edges = 0;

    while (getline(in2, str))
    {
        vector<string> results;
        boost::split(results, str, boost::is_any_of(","), boost::token_compress_on);
        locations.push_back(location(stof(results[1]),stof(results[2])));

        total_nodes++;
    }
    std::cout << "total_nodes=" << total_nodes << endl;

    while (getline(in, str))
    {
        vector<string> results;
        boost::split(results, str, boost::is_any_of(","), boost::token_compress_on);
        //std::vector<std::string> results((std::istream_iterator<string>(iss)),
        //                             std::istream_iterator<string>());
        //std::cout<<results[0]<<flush<<endl;
        //std::cout<<results[1]<<flush<<endl;
        //std::cout<<results[2]<<flush<<endl;
        //std::cout<<"hello1a"<<flush<<endl;
        //std::cout<<","<<results[1]<<flush<<","<<results[2]<<flush<<endl;
        //std::cout<<"dest_index[results[1]]:"<<flush<<stoi(results[1])<<flush<<endl;
        if (stoi(results[1]) == start)
        {
            std::cout<<"skipping incoming edges to start node:"<<stoi(results[0])<<endl;
            continue; //skipping origin incomming edges
        }
        if (stoi(results[0]) == stoi(results[1]))
        { //no loops
            std::cout << "\t\tskipped loop,from," << stoi(results[0]) << ",to," << stoi(results[1]) << endl;
            continue;
        }
        if (all_destinations.find(stoi(results[0])) != all_destinations.end())
        {
            std::cout<<"skipping outgoing edges from destination:"<<stoi(results[0])<<endl;
            continue; //ignoring outgoing edges from destinations
        }
        edges.push_back(make_pair(stoi(results[0]), stoi(results[1])));
        weights.push_back(stof(results[2]));
        total_edges++;
    }
    std::cout << "main_graph, total_edges:," << total_edges << endl;
    //std::cout << "total_nodes:" << total_nodes << ",main_graph:" << main_graph.get_number_vertices() << endl;
}

int main2(int argc, char **argv)
{
  using namespace boost;
  // specify some types
  typedef adjacency_list<listS, vecS, directedS, no_property,
    property<edge_weight_t, cost> > mygraph_t;
  typedef property_map<mygraph_t, edge_weight_t>::type WeightMap;
  typedef mygraph_t::vertex_descriptor vertex;
  typedef mygraph_t::edge_descriptor edge_descriptor;
  typedef mygraph_t::vertex_iterator vertex_iterator;
  typedef std::pair<int, int> edge;
  
  
  // specify data
  /*enum nodes {
    Troy, LakePlacid, Plattsburgh, Massena, Watertown, Utica,
    Syracuse, Rochester, Buffalo, Ithaca, Binghamton, Woodstock,
    NewYork, N
  };
  const char *name[] = {
    "Troy", "Lake Placid", "Plattsburgh", "Massena",
    "Watertown", "Utica", "Syracuse", "Rochester", "Buffalo",
    "Ithaca", "Binghamton", "Woodstock", "New York"
  };
  location locations[] = { // lat/long
    {42.73, 73.68}, {44.28, 73.99}, {44.70, 73.46},
    {44.93, 74.89}, {43.97, 75.91}, {43.10, 75.23},
    {43.04, 76.14}, {43.17, 77.61}, {42.89, 78.86},
    {42.44, 76.50}, {42.10, 75.91}, {42.04, 74.11},
    {40.67, 73.94}
  };*/
  /*
  edge edge_array[] = {
    edge(Troy,Utica), edge(Troy,LakePlacid),
    edge(Troy,Plattsburgh), edge(LakePlacid,Plattsburgh),
    edge(Plattsburgh,Massena), edge(LakePlacid,Massena),
    edge(Massena,Watertown), edge(Watertown,Utica),
    edge(Watertown,Syracuse), edge(Utica,Syracuse),
    edge(Syracuse,Rochester), edge(Rochester,Buffalo),
    edge(Syracuse,Ithaca), edge(Ithaca,Binghamton),
    edge(Ithaca,Rochester), edge(Binghamton,Troy),
    edge(Binghamton,Woodstock), edge(Binghamton,NewYork),
    edge(Syracuse,Binghamton), edge(Woodstock,Troy),
    edge(Woodstock,NewYork)
  };*/
  vector<edge> edge_array;
  //edge_array.push_back(make_pair(0,1));
  //edge_array.push_back(make_pair(1,2));
  //edge_array.push_back(make_pair(2,3));
  /*cost weights[] = { // estimated travel time (mins)
    96, 134, 143, 65, 115, 133, 117, 116, 74, 56,
    84, 73, 69, 70, 116, 147, 173, 183, 74, 71, 124
  };*/
  vector<cost> weights;
  //weights.assign(3,1.0);
  vector<location> locations;
  //locations.push_back({0,1});
  //locations.push_back({1,1});
  //locations.push_back({2,1});

  read_Astar_Graph(edge_array,weights,locations);
  //unsigned int num_edges = sizeof(edge_array) / sizeof(edge);
  unsigned int num_edges = edge_array.size();
  std::cout<<"edge_array:";
  for(auto element : edge_array){
      std::cout<<element.first<<","<<element.second<<endl;
  }
  /*location locations[] = { // lat/long
    {42.73, 73.68}, {44.28, 73.99}, {44.70, 73.46}};*/
  
  
  // create graph
  mygraph_t g(N);
  WeightMap weightmap = get(edge_weight, g);
  for(std::size_t j = 0; j < num_edges; ++j) {
    edge_descriptor e; bool inserted;
    boost::tie(e, inserted) = add_edge(edge_array[j].first,
                                       edge_array[j].second, g);
    weightmap[e] = weights[j];
  }
  
  
  // pick random start/goal
  boost::mt19937 gen(time(0));
  //vertex start = random_vertex(g, gen);
  //vertex goal = random_vertex(g, gen);
  vertex goal = Destinations[0];
  
  
  std::cout << "Start vertex: " << start << endl;
  std::cout << "Goal vertex: " << goal << endl;
  std::cout<<"Num_vertices:"<<num_vertices(g)<<",num_edges:"<<num_edges<<endl;
  
  /*ofstream dotfile;
  dotfile.open("test-astar-cities.dot");
  write_graphviz(dotfile, g,
                 city_writer<const char **, location*>
                  (name, locations, 73.46, 78.86, 40.67, 44.93,
                   480, 400),
                 time_writer<WeightMap>(weightmap));*/
  
  
  vector<mygraph_t::vertex_descriptor> p(num_vertices(g));
  vector<cost> d(num_vertices(g));
  try {
    // call astar named parameter interface
    astar_search_tree
      (g, start,
       distance_heuristic<mygraph_t, cost, vector<location> >
        (locations, goal),
       predecessor_map(&p[0]).distance_map(&d[0]).
       visitor(astar_goal_visitor<vertex>(goal)));
  
  
  } catch(found_goal fg) { // found a path to the goal
    list<vertex> shortest_path;
    for(vertex v = goal;; v = p[v]) {
      shortest_path.push_front(v);
      if(p[v] == v)
        break;
    }
    std::cout << "Shortest path from " << start << " to "
         << goal << ": ";
    list<vertex>::iterator spi = shortest_path.begin();
    std::cout << start;
    for(++spi; spi != shortest_path.end(); ++spi)
      std::cout << " -> " << *spi;
    std::cout << endl << "Total travel cost: " << d[goal] << endl;
    return 0;
  }
  
  std::cout << "Didn't find a path from " << start << "to"
       << goal << "!" << endl;
  return 0;
  
}

float round2(float n){
    return roundf(n * 100) / 100;
}

void calculate_all_subpaths(){
    //Paths is made of all vectorized AG paths
    //Key is destination node (end of path) and AG origin node s
    //We scan all possible paths to generate subpaths list
    auto paths=main_dijkstra_alg->getAGpaths();
    vector<pair<int,float> >costed_path;
    int counter=0;
    //int path_counter=0;
    for(auto path : *paths){
        //path_counter++;
        //std::cout<<"path_counter:"<<path_counter<<",size:"<<path.size()<<flush<<endl;
        float cost=0;
        if(path.size()==1){
            continue;
        }

        //Forward scan to get costed path
        costed_path.clear();
        for(size_t i=0;i<path.size()-1;i++){
            auto key=make_pair(path[i]->getID(),cost);
            //This node has a subpath to d
            costed_path.push_back(key);
            cost+=main_graph.get_edge_weight(path[i],path[i+1]);
            //std::cout<<"["<<path[i]->getID()<<","<<cost<<"],"<<flush<<endl;
        }

        auto last = costed_path.end();
        for(size_t i=0;i<path.size()-1;i++){
            auto first = costed_path.begin() + i;
            vector<pair<int,float> > subpath(first, last);
            pair<int,pair<int,float> > key2=make_pair(path.back()->getID(),costed_path[i]);
            subpaths.insert(make_pair(key2,subpath));
            //std::cout<<"dest:"<<path.back()->getID()<<",path:";
            //for(auto node : path){
            //	std::cout<< node->getID()<<",";
            //}
            //std::cout<<endl;
            //std::cout<<"subpath,parent_node:["<<costed_path[i].first<<","<<costed_path[i].second<<"]"<<endl;
            //for(auto node : subpath){
            //	std::cout<< node.first<<",";
            //}
            //std::cout<<endl;
            counter++;
        }
    }
    //std::cout<<"paths:"<<(*paths).size()<<",subpathcounter:,"<<counter<<endl;
}
float calculate_q_recursive_prefix(unsigned node){
    //auto edges = main_dijkstra_alg->getAGPEdges();
    //auto paths = main_dijkstra_alg->getAGpaths2();
    //auto prefixes = main_dijkstra_alg->get_prefixes();
    //cout << "hola1" << endl;
    if((*prefixes)[node].second==(*paths)[(*prefixes)[node].first].size()-1){
        //cout << "\tnode is destination node" << endl;
        return 0;
    }
    float sum=0;
    //cout<<"hola2"<<flush<<endl;
    auto orig_vertex = (*paths)[(*prefixes)[node].first][(*prefixes)[node].second];
    //cout << "\torig_vertex:" << orig_vertex->getID() << endl;
    for(auto child : (*edges)[node]){
        //cout<<"\tsum1="<<sum<<flush<<endl;
        //auto dest_vertex = (*paths)[(*prefixes)[child].first][(*prefixes)[child].second];
        //cout << "\t\tchild_vertex:," << dest_vertex->getID()<< flush << endl;
        auto key = make_pair(node, child);
        //sum+=Alpha[dest][key]*(calculate_q_recursive(make_pair(child_node,x),dest)+ main_graph.get_edge_weight(orig_vertex, dest_vertex));
        //sum += AlphaPrefix[key] * (calculate_q_recursive_prefix(child)+ main_graph.get_edge_weight(orig_vertex, dest_vertex));
        sum += AlphaPrefix[key] * (calculate_q_recursive_prefix(child) );
    }
    //phi_prefix[node] = sum;
    //cout<<"phi_prefix[,"<<node<<",]:"<<phi_prefix[node]<<",avg_q_prefix:"<<avg_q_prefix[node]<<endl;
    if(avg_q_prefix[node]>=0.98*sum){//To account for rounding errors
        //cout << "\t prefix[,"<<node<<",]:";main_dijkstra_alg->print_prefix(node);
        //cout << " belongs to cutset,sum_prefix:" << sum << ",avg_q_prefix:" << avg_q_prefix[node];
        //auto pred=dest_predictor_prefix(node);
        //cout << ",pred:," << pred.first<<",|,"<<pred.second<<endl;
        cutset_prefix.insert(node);
    }
    //else{
        //cout << "\t prefix[,"<<node<<",]:";main_dijkstra_alg->print_prefix(node);
        //cout << " does2 not belong to cutset,sum_prefix:" << sum << ",avg_q_prefix:" << avg_q_prefix[node] << endl;
        //cout<<"phi_prefix[,"<<node<<",]:"<<phi_prefix[node]<<",avg_q_prefix:"<<avg_q_prefix[node]<<endl;
        //auto pred=dest_predictor_prefix(node);
        //cout << ",pred:," << pred.first<<",|,"<<pred.second<<endl;
        //phi_prefix[node] = max(phi_prefix[node], avg_q_prefix[node]);
        //cout << "\tphi second term is lower than avg_q_prefix,so maximizing:" << node << endl;
    //}
    //phi_prefix[node] = max((float) 2.0,max(sum, avg_q_prefix[node]));
    phi_prefix[node] = max(sum, avg_q_prefix[node]);
    return phi_prefix[node];
}
//given a node (in cutset!), predict which destination to defend
/*pair<unsigned, float> prev_dest_predictor_prefix(unsigned node, unsigned iter){
    size_t chosen_dest = 0;
    float max_q = 0;
    float avg_q = 0;
    for (size_t d = 0; d < Destinations.size();d++){
        if(prev_lambda_obs_dest_given_prefix[iter][d].find(node)==prev_lambda_obs_dest_given_prefix[iter][d].end()){
            continue;//destination prob is 0 for this path
        }
        //cout << "hola3" << flush << endl;
        avg_q = prev_lambda_obs_dest_given_prefix[iter][d][node] * prev_q_given_dest_and_prefix[iter][d][node];
        if(node==5){
            cout<<"\t\td:"<<d<<",lambda:"<<prev_lambda_obs_dest_given_prefix[iter][d][node]<<",q:"<<prev_q_given_dest_and_prefix[iter][d][node]<<",curr_q:"<<avg_q<<endl;
            cout<<"\t\td:"<<d<<",lambda:"<<lambda_obs_dest_given_prefix[d][node]<<",q:"<<q_given_dest_and_prefix[d][node]<<",curr_q:"<<avg_q<<endl;
        }
        if(avg_q>max_q){
            chosen_dest=d;
            max_q = avg_q;
        }
    }
    cout << "\t iter:" << iter << ",pref:" << node << ",chosen_dest:" << chosen_dest << ",reward:" << max_q << endl;
    return make_pair(chosen_dest,max_q);
}*/
//given a node (in cutset!), predict which destination to defend
pair<unsigned, float> dest_predictor_prefix(unsigned node){
    size_t chosen_dest = 0;
    float max_q = 0;
    float avg_q = 0;
    for (size_t d = 0; d < Destinations.size();d++){
        if(lambda_obs_dest_given_prefix[d].find(node)==lambda_obs_dest_given_prefix[d].end()){
            continue;//destination prob is 0 for this path
        }
        avg_q = lambda_obs_dest_given_prefix[d][node] * q_given_dest_and_prefix[d][node];
        if(node==0||node==1||node==175){
            //cout<<"\td:"<<d<<",node:,"<<node<<",lambda:"<<lambda_obs_dest_given_prefix[d][node]<<",q:"<<q_given_dest_and_prefix[d][node]<<",curr_q:"<<avg_q<<endl;
        }
        if(avg_q>max_q){
            chosen_dest=d;
            max_q = avg_q;
        }
    }
    //if(node==0||node==1||node==175){
        //cout<<"\tchosen_dest:,"<<chosen_dest<<",node:,"<<node<<"max_q:,"<<max_q<<endl;
    //}
    return make_pair(chosen_dest,max_q);
}
float best_avg_target_choice(){
    auto paths = main_dijkstra_alg->getAGpaths2();
    auto prefixes_per_path = main_dijkstra_alg->getPathsToPrefix();
    bool cutset_found = false;
    float reward = 0;
    cout<<"hola best_avg,paths:"<<paths->size()<<endl;

    pair<unsigned, float> temp_pair(0, INT_MAX);
    best_target_dest_choice.assign(Destinations.size(), temp_pair);

    //pair<unsigned,float> best_target_choice(0,INT_MAX);

    for(size_t path=0;path<paths->size();path++){
        //cout<<"\t\thola path:"<<path<<endl;
        auto actual_dest = DestinationsOrder2[(*paths)[path].back()->getID()];
        cutset_found = false;
        reward = 0;
        int pred_dest=0;
        unsigned cutset_pref = 0;
        for (auto pref : (*prefixes_per_path)[path]){
            if(cutset_prefix.count(pref)>0){
                /*if(cutset_found){
                    cout << "cutset was already found" << endl;
                }
                cout << "path:,";
                main_dijkstra_alg->print_path(path);
                cout << endl
                     << "prefix:,";
                main_dijkstra_alg->print_prefix(pref);
                cout << endl;*/
                pred_dest = dest_predictor_prefix(pref).first;
                if (pred_dest == actual_dest){
                    //reward = main_dijkstra_alg->get_cost_to_dest(pref, path);
                    reward = modified_q(main_dijkstra_alg->get_cost_to_dest(pref, path));
                    cutset_pref = pref;
                }
                //cout<<"\tFound " <<reward<<" reward for target, pred_dest:"<<pred_dest;
                //cout<<",actual_dest:"<<actual_dest<<",pref,:"<<pref<<",path:,"<<path<<endl;
                //WE ONLY DO ONE SHOT CHOICE AND REWARD
                /*else{//observer will have to wait to next cutset_node
                    continue;
                }*/
                cutset_found = true;
            }
            if(cutset_found){
                break;
            }
        }
        if(!cutset_found){
            cout << "PATH:," << path << ", HAS NO CUTSET, DEBUG ME!!!";
        }
        //check if we found a cheaper path for target to pick for actual dest and overall
        //cout << "\t\t\t\treward:" << reward << ",dest:" << actual_dest << endl;
        if(best_target_dest_choice[actual_dest].second>=(reward*0.98)){//Add equal to allow for choosing shorter paths
            //preference of shorter path if reward is basically the same
            if(best_target_dest_choice[actual_dest].second>=(reward*0.98)&& best_target_dest_choice[actual_dest].second<(reward*1.02)){
                if(main_dijkstra_alg->path_size(path)<main_dijkstra_alg->path_size(best_target_dest_choice[actual_dest].first)){//New path is shorter and same reward because of rounding errors
                    //cout << "\t\tIter:," << prev_cutset_prefix.size();
                    best_target_dest_choice[actual_dest].second=reward;
                    best_target_dest_choice[actual_dest].first=path;
                    //cout << ",updating partial better choice, path["<<path<<"]:,";
                    main_dijkstra_alg->print_path(path);
                    //cout << ",reward:," << reward << endl;

                }
            }
            else{
                //cout << "\t\t\tpref:";
                //main_dijkstra_alg->print_prefix(pref);
                //cout << "\t\tIter:," << prev_cutset_prefix.size();
                best_target_dest_choice[actual_dest].second=reward;
                best_target_dest_choice[actual_dest].first=path;
                //cout << ",skipping partial better choice, path["<<path<<"]:,";
                //main_dijkstra_alg->print_path(path);
                //cout << ",reward:," << reward <<",best_reward"<<best_target_dest_choice[actual_dest].second<<endl;
                }
        }
        //else{
            //cout << "\t\t\tpref:";
            //main_dijkstra_alg->print_prefix(pref);
            //cout << "\t\tIter:," << prev_cutset_prefix.size();
            //cout << ",skipping partial worse choice, path["<<path<<"]:,";
            //main_dijkstra_alg->print_path(path);
            //cout << ",reward:," << reward <<",best_reward"<<best_target_dest_choice[actual_dest].second<<endl;
        //}
    }
    //Now calculate best_average_choice
    float avg_reward = 0;
    size_t counter = 0;
    for(auto choice : best_target_dest_choice){
        if(debug){
            cout << "Iter:," << current_iter;
            cout << ",1Found better choice, path["<<choice.first<<"]:,";
            main_dijkstra_alg->print_path(choice.first);
            cout << ",reward:," << choice.second << endl;
        }
        //cout << "\tDest:," << counter++ << ",best_reward_for_Dest:," << choice.second << endl;
        avg_reward += choice.second;
    }
    avg_reward = avg_reward / float(Destinations.size());
    if(debug){
        cout << "Iter:," << current_iter<<",avg_reward for target choice:," << avg_reward << endl;
    }
    return avg_reward;
}

void elapsed_time(string method, std::chrono::_V2::system_clock::time_point start_time){
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;
    auto current_time = diff.count();
    cout << method << ",time:,"<<current_time<< endl;
};
void store_prev_strategies(){
    cout << "Storing Prev Strategy" << endl;
    //First store the components to estimate a destination for the observer
    //prev_lambda_obs_dest_given_prefix.push_back(lambda_obs_dest_given_prefix);
    //prev_q_given_dest_and_prefix.push_back(q_given_dest_and_prefix);
    //auto iter = prev_lambda_obs_dest_given_prefix.size() - 1;
    //cout << "\tstored d:0,lambda:" << lambda_obs_dest_given_prefix[0][5] << ",q:" << q_given_dest_and_prefix[0][5] << endl;
    //Then store cutset as well, now we have (Cx,Dx)
    store_cutset(cutset_prefix);
    cutset_prefix.clear();
    //Now store previous Q set of paths from the target
    prev_best_target_dest_choice.push_back(best_target_dest_choice);
    cout<<"size of best_target_dest_choice:"<<best_target_dest_choice.size()<<endl;

}
void calculate_next_observer_strategy(){
    //Calculate new cutset based on the average of previous mus. We create an average mu policy based on all previous mus.
    auto paths = main_dijkstra_alg->getAGpaths2();
    auto prefixes_per_path = main_dijkstra_alg->getPathsToPrefix();
    
    auto previousTime = std::chrono::high_resolution_clock::now();
    //Reset the necessary calculations
    //for (size_t d = 0; d < Destinations.size();d++){
    //    std::fill(dest_prob_per_prefix[d].begin(), dest_prob_per_prefix[d].end(), 0.0);//fastest reset
   // }
    prob_path.assign(paths->size(), 0.0);
    dest_prob_path.assign(Destinations.size(), prob_path);
    //Average previous mus, calculate the frequency of all chosen paths in previous mus
    map<int, int> chosen_times;
    vector<int> paths_per_dest(Destinations.size(),0);
    for(auto mu : prev_best_target_dest_choice){
        for (size_t d = 0; d < Destinations.size();d++){
            //cout<<"mu size:"<<mu.size()<<flush<<endl;
            chosen_times[mu[d].first]++;
            paths_per_dest[d]++;
        }
    }
    for(auto p : chosen_times){
        //cout<<"Iter:,"<<current_iter<<",path:,"<<p.first<<",times:,"<<p.second;
        //auto dest=paths
        auto d = DestinationsOrder2[(*paths)[p.first].back()->getID()];
        //cout<<",d:,"<<d;
        float weight = float(p.second) / float(iterations);
        //prof per path given dest is number of times path is chosen vs number of iterations
        dest_prob_path[d][p.first] = weight;
        //cout<<",dest_prob_path[,"<<d<<",][,"<<p.first<<",]:,"<<weight<<",lambdaDest:,"<<lambdaDest;
        prob_path[p.first]=weight*lambdaDest;
        //cout<<",prob_path["<<p.first<<"],"<< prob_path[p.first]<<endl;
    }
    elapsed_time("\t\tSetup_observer_aggregated_target_paths", previousTime);
    calculate_min_path_strategy_prefixes(true);
    //Now we calculate recursively phi and Cutset starting at origin
    auto start_time = std::chrono::high_resolution_clock::now();
    string method("calculate_q_recursive");
    float cutset_reward=calculate_q_recursive_prefix(0);
    removeDominatedFromCutset();
    cout << "Iter:,"<<current_iter<<",cutset_reward:,"<<cutset_reward<<",clean_cutset:," <<cutset_prefix<< endl;
    cout << "Iter:," << current_iter<<",clean_cutset_endings:,"; print_cutset_endings(cutset_prefix); cout << endl;
    //if(debug){
        for (auto pref : cutset_prefix){
                cout << "\t cutset_pref[" << pref<<"]:,";
                main_dijkstra_alg->print_prefix(pref);
                //auto pred_dest = prev_dest_predictor_prefix(pref,prev_cutset_prefix.size()-1).first;
                auto pred_dest = dest_predictor_prefix(pref).first;
                cout <<",pred_dest:,"<<pred_dest<<",:,"<<DestinationsOrder1[pred_dest]<<endl;
        }
    //}
    elapsed_time(method, start_time);
}
void calculate_next_target_strategy(){
    cout << "Calling calculate_next_target_strategy,Iter:," << current_iter<<endl;
    auto paths = main_dijkstra_alg->getAGpaths2();
    auto prefixes_per_path = main_dijkstra_alg->getPathsToPrefix();
    float reward = 0;
    pair<unsigned, float> temp_pair(0, INT_MAX);
    best_target_dest_choice.assign(Destinations.size(), temp_pair);
    vector<bool> cutsets_found(prev_cutset_prefix.size(),false);
    unsigned cutsets_counter = 0;
    unsigned cutsets_found_counter = 0;
    //HACK TO TEST OBSERVER REACTION TO PATHSET RATIOS IN EXAMPLE3
    /*int path1=0;int path2=0;
    auto n = rand() % 100;
    if(n<60){
        path1 = 0;
        path2 = 1;
    }
    else{
        path1=2;
        path2=3;
    }*/
    //Calculate new target pathset based on the average of previous cutsets
    //Iterate over each path and calculate target best response ()
    auto previousTime = std::chrono::high_resolution_clock::now();

    for(size_t path=0;path<paths->size();path++){
        //auto currentTime = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - previousTime);
        //cout << "path:," << path <<",time:"<<duration.count()<<endl;
        //previousTime = std::chrono::high_resolution_clock::now();

        auto actual_dest = DestinationsOrder2[(*paths)[path].back()->getID()];
        reward = 0;
        cutsets_found.assign(prev_cutset_prefix.size(),false);
        cutsets_found_counter = 0;
        //Check if current prefix node is in any of the previous cutsets, if not rho is 0.
        for (auto& pref : (*prefixes_per_path)[path]){
            //cout << "pref:," << pref <<flush<< endl;
            for (size_t iter = 0;iter<prev_cutset_prefix.size();iter++)
            {
                //cout << "iter:," << iter <<flush<< endl;
                if(cutsets_found[iter]){
                    //cutsets_counter++;
                    continue;
                }
                pair<unsigned, set<unsigned> > & temp_cutset = prev_cutset_prefix[iter];

                if(temp_cutset.second.count(pref)>0){
                    //NEED TO MAKE SURE THE dest_predictor_prefix FUNCTION IS CORRECT IN THIS CONTEXT
                    //ALSO NEED TO ENSURE THAT THERE IS NO DOUBLE COUNT OF REWARD
                    cutsets_found[iter]=true;
                    cutsets_found_counter++;
                    //auto pred_dest = prev_dest_predictor_prefix(pref,iter).first;
                    //cout<<"\tcutset:,"<<pref<<"[";main_dijkstra_alg->print_prefix(pref);
                    //cout<<"],for path:,"<<path<<",found,pred_dest:,"<<pred_dest<<",actual_dest:,"<<actual_dest<<",repetitions:,"<<temp_cutset.first<<endl;
                    //auto weight = main_dijkstra_alg->get_cost_to_dest(pref, path);
                    auto weight = modified_q(main_dijkstra_alg->get_cost_to_dest(pref, path));
                    //SQUARE WEIGHT HACK
                    //weight = weight * weight;
                    //if(saturation>0){
                    //    weight = weight / saturation;
                    //    weight = min(weight, (float)1.0);
                    //}
                    //weight = modified_q(weight);
                    if(prev_cut_data[iter][pref][actual_dest]>0){
                        reward += (weight * float(prev_cut_data[iter][pref][actual_dest]));
                    }
                    //cout<<"\tIter:,"<<current_iter<<",cutset:,"<<pref<<"[";main_dijkstra_alg->print_prefix(pref);
                    //cout<<"],for path:,"<<path<<",found,actual_dest:,"<<actual_dest<<",repetitions:,"<<temp_cutset.first<<",prev_cut_data:"<<prev_cut_data[iter][pref];
                    //cout << ",reward:," << (weight * float(prev_cut_data[iter][pref][actual_dest])) << ",max_reward:," << weight << endl;
                }
                //cutsets_counter++;
            }
            //Reward per node is a function of how many cutsets if temp_cutset is chosen then probability of making the choice is 100%
            //If all cutsets were found,we are finished with the path.
            if(cutsets_found_counter==prev_cutset_prefix.size()){
                //cout << "\t\t found all cutsets for pref:"<<pref<<"[";
                //main_dijkstra_alg->print_prefix(pref);
                //cout << "]" << endl;
                break;
            }
        }
        reward = (reward / float(current_iter+1));
        //if(saturation>0)
        //    reward = min(reward, (float)1.0);
        /*if(current_iter==140){
            cout << "\t\titer:," << current_iter << ",avg_reward:," << reward << ",for path:[" << path<<"]:,";
            main_dijkstra_alg->print_path(path);
            cout << endl;
        }*/

        //Now update the best_path for the target for this destination if necessary
        //if(path==path1||path==path2)
        if(best_target_dest_choice[actual_dest].second>=(reward*0.98)){//Add equal to allow for choosing shorter paths
            //preference of shorter path if reward is basically the same
            if(best_target_dest_choice[actual_dest].second>=(reward*0.98)&& best_target_dest_choice[actual_dest].second<(reward*1.02)){
                if(main_dijkstra_alg->path_size(path)<main_dijkstra_alg->path_size(best_target_dest_choice[actual_dest].first)){//New path is shorter and same reward because of rounding errors
                    best_target_dest_choice[actual_dest].second=reward;
                    best_target_dest_choice[actual_dest].first=path;
                }
            }
            else{
                //cout << "\t\t\tpref:";
                //main_dijkstra_alg->print_prefix(pref);
                //cout << "\tIter:," << prev_cutset_prefix.size();
                best_target_dest_choice[actual_dest].second=reward;
                best_target_dest_choice[actual_dest].first=path;
                //cout << ",found better choice, path["<<path<<"]:,";
                //main_dijkstra_alg->print_path(path);
                //cout << ",reward:," << reward << endl;
                }
        }
    }
    float avg_reward = 0;
    int counter = 0;
    set<int> chosen_paths;
    for (auto choice : best_target_dest_choice){
        chosen_paths.insert(choice.first);
        cout << "Iter:," << current_iter;
        cout << ",2Found better choice, path["<<choice.first<<"]:,";
        main_dijkstra_alg->print_path(choice.first);
        cout << ",reward:," << choice.second << endl;
        avg_reward += choice.second;
        counter++;
    }
    avg_reward = avg_reward / float(Destinations.size());
    target_iterative_rewards.push_back(avg_reward);
    float avg_iterative_target_reward = 0;
    for(auto reward: target_iterative_rewards)
        avg_iterative_target_reward += reward;
    avg_iterative_target_reward = avg_iterative_target_reward / float(target_iterative_rewards.size());

    cout << "Iter:," << current_iter << ",last_avg_reward for target choice:," << avg_reward <<",avg_iterative_reward:,"<<avg_iterative_target_reward<< ",chosen_paths:," << chosen_paths<<endl;
}
void calculate_next_target_strategy2(){
    cout << "Calling calculate_next_target_strategy2,Iter:," << current_iter<<endl;
    auto paths = main_dijkstra_alg->getAGpaths2();
    auto prefixes_per_path = main_dijkstra_alg->getPathsToPrefix();
    double reward = 0;
    pair<unsigned, float> temp_pair(0, INT_MAX);
    best_target_dest_choice.assign(Destinations.size(), temp_pair);
    vector<bool> cutsets_found(prev_cutset_prefix.size(),false);
    unsigned cutsets_counter = 0;
    unsigned cutsets_found_counter = 0;
    //HACK TO TEST OBSERVER REACTION TO PATHSET RATIOS IN EXAMPLE3
    /*int path1=0;int path2=0;
    auto n = rand() % 100;
    if(n<60){
        path1 = 0;
        path2 = 1;
    }
    else{
        path1=2;
        path2=3;
    }*/
    //Calculate new target pathset based on the average of previous cutsets
    //Iterate over each path and calculate target best response ()
    //auto previousTime = std::chrono::high_resolution_clock::now();

    for(size_t path=0;path<paths->size();path++){
        //auto currentTime = std::chrono::high_resolution_clock::now();
        //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - previousTime);
        //cout << "path:," << path <<",time:"<<duration.count()<<endl;
        //previousTime = std::chrono::high_resolution_clock::now();

        auto actual_dest = DestinationsOrder2[(*paths)[path].back()->getID()];
        reward = 0;
        cutsets_found_counter = 0;
        //Check if current prefix node is in any of the previous cutsets, if not rho is 0.
        for (auto& pref : (*prefixes_per_path)[path]){
            if(merged_prev_cut_data.find(pref)!=merged_prev_cut_data.end()){
                /*if(pref==45){
                    cout << "\t\tpref:," << pref << ",freq:," << merged_prev_cut_data[pref][actual_dest] << ",actual_dest:," << actual_dest << ",path:,";
                    main_dijkstra_alg->print_path(path);
                    cout << endl;
                }*/
                if(merged_prev_cut_data[pref][actual_dest]>0){
                    cutsets_found_counter += merged_prev_cut_data[pref][actual_dest];
                    auto weight = modified_q(main_dijkstra_alg->get_cost_to_dest(pref, path));
                    reward += weight * merged_prev_cut_data[pref][actual_dest];
                }
            }
            else{
                continue;
            }
        }
        reward = (reward / float(current_iter+1));
        //if(saturation>0)
        //    reward = min(reward, (float)1.0);
        if(current_iter==140){
            cout << "\t\titer:," << current_iter << ",avg_reward:," << reward << ",for path:[" << path<<"]:,";
            main_dijkstra_alg->print_path(path);
            cout << endl;
        }

        //Now update the best_path for the target for this destination if necessary
        //if(path==path1||path==path2)
        if(best_target_dest_choice[actual_dest].second>=(reward*0.98)){//Add equal to allow for choosing shorter paths
            //preference of shorter path if reward is basically the same
            if(best_target_dest_choice[actual_dest].second>=(reward*0.98)&& best_target_dest_choice[actual_dest].second<(reward*1.02)){
                if(main_dijkstra_alg->path_size(path)<main_dijkstra_alg->path_size(best_target_dest_choice[actual_dest].first)){//New path is shorter and same reward because of rounding errors
                    best_target_dest_choice[actual_dest].second=reward;
                    best_target_dest_choice[actual_dest].first=path;
                }
            }
            else{
                //cout << "\t\t\tpref:";
                //main_dijkstra_alg->print_prefix(pref);
                //cout << "\tIter:," << prev_cutset_prefix.size();
                best_target_dest_choice[actual_dest].second=reward;
                best_target_dest_choice[actual_dest].first=path;
                //cout << ",found better choice, path["<<path<<"]:,";
                //main_dijkstra_alg->print_path(path);
                //cout << ",reward:," << reward << endl;
                }
        }
    }
    float avg_reward = 0;
    int counter = 0;
    set<int> chosen_paths;
    for (auto choice : best_target_dest_choice){
        chosen_paths.insert(choice.first);
        cout << "Iter:," << current_iter;
        cout << ",2Found better choice, path["<<choice.first<<"]:,";
        main_dijkstra_alg->print_path(choice.first);
        cout << ",reward:," << choice.second << endl;
        avg_reward += choice.second;
        counter++;
    }
    avg_reward = avg_reward / float(Destinations.size());
    target_iterative_rewards.push_back(avg_reward);
    float avg_iterative_target_reward = 0;
    for(auto reward: target_iterative_rewards)
        avg_iterative_target_reward += reward;
    avg_iterative_target_reward = avg_iterative_target_reward / float(target_iterative_rewards.size());

    cout << "Iter:," << current_iter << ",last_avg_reward for target choice:," << avg_reward <<",avg_iterative_reward:,"<<avg_iterative_target_reward<< ",chosen_paths:," << chosen_paths<<endl;
}
void fictitious_play2(){
//    iterations = 1000;
    lambdaDest = 1.0 / float(Destinations.size());
    for (current_iter = 0;current_iter<iterations;current_iter++){
        cout<<"iter:,"<<current_iter<<",calling store_prev_strategies"<<endl;
        /*for (auto pref : cutset_prefix){
            cout<<"Santi,Iter:,"<<current_iter<<",pref:,"<<pref<<",dest:,"<<dest_predictor_prefix(pref).first<<endl;
        }*/
        store_prev_strategies();
        auto start_time = std::chrono::high_resolution_clock::now(); string method("calculate_next_target_strategy2_");
        method += to_string(current_iter);
        calculate_next_target_strategy2();
        elapsed_time(method, start_time);
        start_time = std::chrono::high_resolution_clock::now(); method="calculate_next_observer_strategy";
        method += to_string(current_iter);
        calculate_next_observer_strategy();
        elapsed_time(method, start_time);

    }

    for (auto pref : cutset_prefix){
        cout << "\t cutset_pref[" << pref<<"]:,";
        main_dijkstra_alg->print_prefix(pref);
        //auto pred_dest = prev_dest_predictor_prefix(pref,prev_cutset_prefix.size()-1).first;
        auto pred_dest = dest_predictor_prefix(pref).first;
        cout <<",pred_dest:,"<<pred_dest<<",:,"<<DestinationsOrder1[pred_dest]<<endl;
    }
    calculate_mu_rho_reward();

    // Calculate percentage change for last 20 iterations
    // calculate_percentage_difference(target_iterative_rewards,160);
    // vector<float> target_iterative_rewards_capped(target_iterative_rewards.begin()+150 , target_iterative_rewards.end());
    vector<float> target_iterative_rewards_capped(target_iterative_rewards.begin() , target_iterative_rewards.end());
    calculate_statistics(target_iterative_rewards_capped);
}
void fictitious_play(){
    //For every iteration recalculate best_target_response based on the opposite player strategy
    size_t iterations = 10;
    lambdaDest = 1.0 / float(Destinations.size());
    float total_prob = 0;
    edges = main_dijkstra_alg->getAGPEdges();
    prefixes = main_dijkstra_alg->get_prefixes();
    paths = main_dijkstra_alg->getAGpaths2();
    vector<vector<float> > prev_prob;
    //vector<float> avg_dest_prob_path;
    float temp_avg = 0;
    prev_prob.push_back(prob_path);

    // Now we need to recalculate mu and rho
    // First recalculate cutset based on the new path_probs

    start_time = std::chrono::high_resolution_clock::now();
    auto best_target_choice = best_avg_target_choice();
    elapsed_time("target_choice", start_time);

    for(size_t iter=0;iter<iterations;iter++){
        cout << "fictitious_play, iteration:" << iter << endl;
        for (size_t d = 0; d < Destinations.size();d++){
            std::fill(dest_prob_per_prefix[d].begin(), dest_prob_per_prefix[d].end(), 0.0);//fastest reset
        }
        for (size_t i= 0; i < prob_path.size();i++){
            temp_avg = 0;
            for (size_t j = 0; j < iter;j++){
                cout << "j:" << j << "i:" << i << flush << endl;
                temp_avg += prev_prob[j][i];
            }
            auto dest = DestinationsOrder2[(*paths)[i].back()->getID()];
            cout << "hola1" << flush << endl;
            if(i==best_target_dest_choice[dest].first){//This path prob gets average with 100% new prob for last mu
                temp_avg += lambdaDest;
            }
            else{//
                temp_avg += 0;
            }
            cout << "hola2" << flush << endl;
            prob_path[i] /= temp_avg/float(iter + 2);
            for (size_t d = 0; d < Destinations.size();d++){
                dest_prob_path[d][i]=temp_avg/lambdaDest;
            }
        }
            //cout << "prev best_path for d:" << d << ",is:" << best_target_dest_choice[d].first << endl;
            //prob_path[best_target_dest_choice[d].first] = 1.0 * lambdaDest;
            //dest_prob_path[d][best_target_dest_choice[d].first]=prob_path[best_target_dest_choice[d].first]/lambdaDest;
            //cout << "new dest_prob_path[" << d << "][" << best_target_dest_choice[d].first << "]:" << dest_prob_path[d][best_target_dest_choice[d].first]<<endl;
        //size_t counter = 0;
        //for (auto prob: prob_path){
        //    cout<<"path["<<counter<<"]:"<<prob_path[counter++]<<endl;
        //    total_prob += prob;
        //}

        cout << "starting with calculate_min_path_strategy_prefix,iter:," << iter << endl;
        calculate_min_path_strategy_prefixes(true);
        cout << "finished with calculate_min_path_strategy_prefix,iter:," << iter << endl;

        auto start_time = std::chrono::high_resolution_clock::now();
        string method("calculate_q_recursive");
        method += to_string(iter);
        calculate_q_recursive_prefix(0);
        elapsed_time(method, start_time);
        exit(1);

        //Now Calculate the new target_strategy
        start_time = std::chrono::high_resolution_clock::now();
        method += to_string(iter);
        auto best_target_choice = best_avg_target_choice();
        elapsed_time(method, start_time);
        //NEED TO CALCULATE NEW PROB_PATH TO ADD TO PREV_PROB
        cout<<"NEED TO CALCULATE NEW PROB_PATH TO ADD TO PREV_PROB"<<endl;
        exit(1);
    }
    
}

pair<float,float> calculate_percentage_difference(std::vector<float> values, int starting_iteration){

    if (values.size() < 2) {
        std::cout << "Not enough values to calculate percentage differences." << std::endl;
        std::cerr << "Not enough values to calculate percentage differences." << std::endl;
        exit(1);
    }
    vector<float> differences;

    for (size_t i = starting_iteration; i < values.size(); ++i) {
        double previousValue = values[i - 1];
        double currentValue = values[i];
        
        // Handle the case where the previous value is zero to avoid division by zero
        if (previousValue == 0) {
            std::cout << "Undefined (previous value is zero)." << std::endl;
        } else {
            float percentageDifference = ((currentValue - previousValue) / previousValue) * 100;
            differences.push_back(percentageDifference);
            std::cout << "Percentage difference between value " << i << " and value " << i - 1 << ": "
                      << std::abs(percentageDifference) << "%" << std::endl;
        }
    }
    //vector<float> differences_capped;
    //differences_capped.assign(differences.begin() + starting_iteration, differences.end());
    cout << "differences:," <<differences<< endl;

    return calculate_statistics(differences);
}
pair<float,float> calculate_statistics(const std::vector<float>& values) {
    size_t n = values.size();
    if (n == 0) {
        cout << "No values available for calculate_statistics." << std::endl;
        cerr << "No values available for calculate_statistics." << std::endl;
        exit(1);
        return make_pair(0.0,0.0);
    }

    // Calculate the average
    float mean = std::accumulate(values.begin(), values.end(), 0.0f) / n;
    
    // Calculate the standard deviation
    float variance = 0.0f;
    for (float value : values) {
        variance += (value - mean) * (value - mean);
    }
    variance /= n;
    float stddev = std::sqrt(variance);

    // Output the results
    std::cout << "Mean:," << mean << ",stdev:," << stddev << ",CV:, " << stddev/mean << ",entries:"<< n << endl;
    return make_pair(mean, stddev);
}
void store_cutset(set<unsigned>& cutset_prefix){
    total_cutsets++;
    cout << "Calling store_cutset" << endl;
    //First create 2D cutset_map where [cut][repetitions,d0,d1,...dn]
    map<unsigned, unsigned> cutset_data;
    for(auto& elem : cutset_prefix){
        auto data=dest_predictor_prefix(elem);
        cutset_data[elem] = data.first;
    }



    //Check if set previously seen
    bool found = false;
    //for(auto& elem : prev_cutset_prefix)
    for(size_t i=0;i<prev_cutset_prefix.size();i++){
        //Found the cutset
        if(cutset_prefix == prev_cutset_prefix[i].second){
            //vector<unsigned> D(Destinations.size() + 1, 0);
            //D[0] = 1;
            //update D vector where 1st element is number of repetitions of cutset and the rest are how many times each d is chosen per cut. 
            prev_cutset_prefix[i].first++;
            found = true;
            cout << "Found cutset, prev_cutset_prefix remains to:," << prev_cutset_prefix.size() <<", repetitions:,"<<prev_cutset_prefix[i].first<<flush<<endl;
            //Now update destination count for the cutset
            for(auto& elem : prev_cutset_prefix[i].second)//Iterate over prefixes
            {
                prev_cut_data[i][elem][cutset_data[elem]]++;
                merged_prev_cut_data[elem][cutset_data[elem]]++;
                //cout << "\tprev_cut_data[" << i << "]" << "["<<elem<<"]["<<prev_cut_data[i][elem]<<endl;
                //cout << "\tmerged_prev_cut_data["<<elem<<"]:"<<merged_prev_cut_data[elem]<<endl;
            }
            break;
        }
    }
    //If not found, add as a new cutset
    if(!found){
        prev_cutset_prefix.push_back(make_pair(1, cutset_prefix));
        cout << "New cutset:,"<<cutset_prefix<<",prev_cutset_prefix size grows to:," << prev_cutset_prefix.size() << endl;

        //Now update destination count for the cutset
        map<unsigned, vector<unsigned>> temp_map;
        prev_cut_data.push_back(temp_map);
        vector<unsigned> temp_count(Destinations.size(), 0);
        for(auto& elem : prev_cutset_prefix.back().second){//Iterate over prefixes
            prev_cut_data.back()[elem] = temp_count;
            prev_cut_data.back()[elem][cutset_data[elem]]++;
            //cout << "\tprev_cut_data[" << prev_cutset_prefix.size() - 1 << "]" << "["<<elem<<"]["<<prev_cut_data.back()[elem]<<endl;
            merged_prev_cut_data[elem] = temp_count;
            merged_prev_cut_data[elem][cutset_data[elem]]++;
            //cout << "\tmerged_prev_cut_data["<<elem<<"]:"<<merged_prev_cut_data[elem]<<endl;
        }
    }
    if(current_iter%100==99){
        for(size_t i=0;i<prev_cutset_prefix.size();i++){
            cout << "Freq,iter:," << current_iter << ",";
            print_cutset_endings(prev_cutset_prefix[i].second);
            cout << "," << float(prev_cutset_prefix[i].first) / float(current_iter + 1) << endl;
        }
    }
}
void removeDominatedFromCutset(){
    set<int> nodesToRemove;
    //auto edges = main_dijkstra_alg->getAGPEdges();
    
    for (const auto& node : cutset_prefix) {
        exploreDescendants(node, edges, nodesToRemove);
    }
    
    for (const auto& node : nodesToRemove) {
        cutset_prefix.erase(node);
    }
}

void exploreDescendants(int node, std::shared_ptr<std::vector<std::set<unsigned int>>> edges,
                        std::set<int>& nodesToRemove) {
    for (const auto& child : (*edges)[node]) {
        nodesToRemove.insert(child);
        //cout << "Removing prefix[" << child << "]";
        //main_dijkstra_alg->print_prefix(child);
        //cout << ",child of:,";
        //main_dijkstra_alg->print_prefix(node);
        //cout << endl;
        exploreDescendants(child, edges, nodesToRemove);
    }
}
void print_cutset_endings(set<unsigned>& cutset){
    set<unsigned> ordered_endings;
    cout << "[";
    for(auto pref : cutset)
    {
        ordered_endings.insert(main_dijkstra_alg->get_prefix_ending(pref));
    }
    for(auto ending : ordered_endings){
        cout<<","<<ending;
    }
    cout << ",]";
}
void read_paths_from_file(){
    ifstream file(input_paths_file);

    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + input_paths_file);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<unsigned> vec;
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            vec.push_back(std::stoi(value));
        }

        input_paths.push_back(vec);
    }

    cout << "input paths:;";
    for(auto path : input_paths){
        cout << path << ";";
    }
    file.close();
}
/*void calculate_reward_last_cutset(){
    //iterate over all prefix paths to destinations
    //add the pred_destination reward(s) and average out with the failures
    float cutset_average_reward = 0;
    for (auto pref : cutset_prefix){
        for(size_t dest=0;dest<Destinations.size();dest++){
            auto pred_dest = dest_predictor_prefix(pref).first;

            for (auto path : dest_paths_per_prefix){
                float positive_reward = main_dijkstra_alg->get_cost_to_dest(pref, path);
            }
        }
    }
}*/
void calculate_mu_rho_reward(){
    //for all possible paths
    //for prefix k=0 to k=length of path
    //calculate average chance prefix is in target path
    //We got it in prob_path from last observer calculation.
    //calculate average chance observer predicts destination correctly.
    //Note: prediction chance is 0 if prefix is not yet inside cut.
    //mulitply chances by q(weight) of the corresponding suffix
    auto prefixes_per_path = main_dijkstra_alg->getPathsToPrefix();
    float mu = 0;//mu prob
    float q = 0;// q(w(gamma(k+)))
    float cutsets_found_counter = 0;
    //WE ADAPT EQUATION 4 (phi(mu,rho)) so that we iterate over all existing cutsets per path
    //Phi as a function of all paths for the aggregated observer and target strategies
    float phi_mu_rho = 0;
    float rho = 0;
    for(size_t path=0;path<paths->size();path++){
        //mu prob is prob_path
        mu = prob_path[path];
        if(mu<0.000001)
            continue;
        //now we calculate observer probability it gets destination right
        auto actual_dest = DestinationsOrder2[(*paths)[path].back()->getID()];
        //iterate over k = 0;
        for (auto& pref : (*prefixes_per_path)[path]){
            //cout << "\t\tprefix:,"; main_dijkstra_alg->print_prefix(pref);
            cutsets_found_counter=0;
            //create reference to prev_cutset_prefix
            //pair<unsigned, set<unsigned> > & temp_cutset = prev_cutset_prefix[iter];
            if(merged_prev_cut_data.find(pref)!=merged_prev_cut_data.end()){
                if(merged_prev_cut_data[pref][actual_dest]>0){
                    cutsets_found_counter += merged_prev_cut_data[pref][actual_dest];
                }
            }
            if(cutsets_found_counter==0)
                continue;
            rho = cutsets_found_counter;
            rho=rho/((float)total_cutsets);
            q = main_dijkstra_alg->get_cost_to_dest(pref, path);
            cout << "\t\t q before adjust:," << q << endl;
            //SQUARE WEIGHT HACK
            //q = q * q;
            q = modified_q(q);
            phi_mu_rho += rho * mu * q;
            cout<<"\t\t q after adjust:,"<<q<<",mu:"<<mu<<",rho:"<<rho<<",phi_mu_rho:"<<phi_mu_rho<<endl;
            if(debug){
            cout << "\t\t pref:,"; main_dijkstra_alg->print_prefix(pref);
            cout<< ",dest:," << actual_dest << ",cutsets_found:," << cutsets_found_counter << ",total_cutsets:," << total_cutsets << endl;
            
            cout << "\tpref:,";
            main_dijkstra_alg->print_prefix(pref);
            cout<< ",path_reward:," << mu * rho * q << ",mu:," << mu << ",rho:," << rho << ",q:" << q << endl;}
        }
        cout << "rho:,"<< rho<<",for path:,";main_dijkstra_alg->print_path(path);cout<<",partial phi_mu_rho:,"<<phi_mu_rho<<endl;
        //cout << "," << endl;
        //cout << "reward:,"<< rho*mu<<"for path:,";main_dijkstra_alg->print_path(path);
        //cout << "," << endl;

    }
    cout << "phi_mu_rho:," << phi_mu_rho << endl;
}
/*float calculate_all_cut_prob()
{
    int cutset_repetitions = 0;
    //First, Is there a cut for this prefix
    for (size_t i = 0;i<prev_cutset_prefix.size();i++){
            if(prev_cutset_prefix[i].second.count(pref)>0){
                cutset_repetitions += prev_cutset_prefix[i].first;
            }
    }
    return float(cutset_repetitions)/
}*/
float modified_q(float q){
    if(q_type==0){//Unmodifier,return normalized result
        return q/normalization_param;
    }
    else if(q_type==1){//Division
        //cout << "\t\tq:," << q << ",sat:," << saturation << ",modified:,"<<q / saturation << endl;
        return q / saturation;
        //return min(q / saturation, float(1.0));
    }
    else if(q_type==2){//Min of Division vs 1
        //cout<<"q_0:,"<<q<<",q/sat:,"<<q/saturation<<"max:,"<<max(q / saturation, float(1.0))<<endl;
        return min(q / saturation, float(1.0));
    }
    else if(q_type==3){//Max of Division vs 1
        //cout<<"q_0:,"<<q<<",q/sat:,"<<q/saturation<<"max:,"<<max(q / saturation, float(1.0))<<endl;
        return max(q / saturation, float(1.0));
    }
    else if(q_type==4){//Squared and divided
        return (q*q) / saturation;
    }
    else{
        cerr<<"q_type:"<<q_type<<" is undefined"<<endl;
        exit(1);
    }
}