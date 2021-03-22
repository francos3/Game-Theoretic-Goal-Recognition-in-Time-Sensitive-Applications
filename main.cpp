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
#include "YenTopKShortestPathsAlg.h"
#include <boost/numeric/ublas/matrix.hpp>
#include "colors.h"
#include <tuple>
#include <boost/algorithm/string.hpp>
#include <ctime>

//#include<graphics.h>
//#include <curses.h>
//#include <boost/timer/timer.hpp>
//using namespace boost::timer;
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
int strategy = 0;
double Budget = 0;
double Beta = 0;
bool acyclic=true;
bool grandparent_check=false;
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
vector<vector<BasePath *>> orig_to_dest_shortest_paths;
vector<Link> lks;
int total_nodes = 0;
set<int> reachable_nodes; //For SaT Graphs
set<int> nodes_set;		  //For SaT Graphs
bool create_PO_graph = false;
int PO_level = 100;
unordered_set<int> FO_nodes;

vector<map<pair<pair<int, float>, pair<int, float>>, float>> Alpha;

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
	result->PrintOut(cout);
}

void testYenAlg(const int &k,
				Link *lk,
				const int &size,
				const int &s,
				const int &d,
				const int &nodes,
				std::shared_ptr<PathMatrix> PthMat)
{
	int budget = INT_MAX;

	Graph my_graph;
	auto all_paths = std::make_shared<set<vector<int>, vect_comp_int>>();

	my_graph.set_number_vertices(nodes);

	for (int i = 0; i < size; i++)
	{
		my_graph.add_link(lk[i].u, lk[i].v, lk[i].weight);
		//if(debug)
		//cout<<"\tAdded Yen Graph node["<<lk[i].u<<","<<lk[i].v<<"],w="<<lk[ i ].weight<<endl;
	}
	if (debug)
		cout << "Yen Graph size:" << size << ",d:" << d << ",origin:" << s << endl;

	my_graph.setv();

	YenTopKShortestPathsAlg yenAlg(my_graph,
								   my_graph.get_vertex(s),
								   my_graph.get_vertex(d));

	//TESTING HACK
	cout << "TESTING DIJKSTRA TO ALL PATHS" << endl;
	yenAlg.get_shortest_distance_to_all_nodes(my_graph.get_vertex(s));
	exit(1);
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
			//cout<<"optimal_distance["<<d<<"]:"<<optimal_distances[d]<<",budget:"<<budget<<endl;
		}
		if (current_path.length() > budget)
		{
			//cout<<"paths_created for current destination:"<<d<<",finish, all future paths will have length bigger than:"<<current_path.length()<<",max_length allowed:"<<C<<endl;
			break;
		}
		if (budget == INT_MAX)
		{ //first path is optimal for each destination
			optimal_distances[d] = current_path.length();
			cout << "optimal_distances[" << d << "]:" << optimal_distances[d] << endl;
		}

		//current_path.PrintOut(cout);
		current_path.add_paths_set(all_paths);
		//cout<<"\t added path"<<endl;
	}
	PthMat->add_paths_set(all_paths, d);
}
int optimizing_observer(PathMatrix *PM, int nodes, int t_max)
{
	int t = 1;
	int disclosing_paths = 0;
	set<int> unsensed_nodes;
	PathMatrix PM_original = *PM;

	PathMatrix PM_temp;

	cout << "NOW OPTIMIZING OBSERVER NETWORK TO KEEP T_MAX LEVEL OF:" << t_max << endl;

	for (int w = 1; w < nodes; w++)
	{
		PM_temp = PM_original;
		unsensed_nodes.insert(w);
		PM_temp.set_W(&unsensed_nodes);

		for (t = 1; t < 10000; t++)
		{
			if (t > t_max)
			{
				cout << "\t\tFinished, this node is unsafe, t score raised by at least one level to:" << t << endl;
				break;
			}
			if (debug)
				cout << "\tWorking with t:" << t << ",candidate_node:" << w << endl;
			PM_temp.get_safe_nodes(t);
			if (PM_temp.erase_disclosing_paths(t, disclosing_paths))
			{
				if (debug)
					cout << "\t\tFinished, erasing path would leave graph disconnected,final t score:" << t << " is lower than t_max:" << t_max << endl;
				break;
			}
		}

		if (t <= t_max)
		{ //we can safely remove this node from observation grid
			cout << "\t t:" << t << ",t_max:" << t_max << ",safely removed " << w << " node from observation list" << endl;
		}
		else
		{
			cout << "\t t:" << t << ",t_max:" << t_max << ",unsafely removed " << w << " node , it needs to be in observation list" << endl;
			unsensed_nodes.erase(w);
		}
	}

	cout << "The following nodes can be safely ignored by the observer" << endl;
	for (auto w : unsensed_nodes)
		cout << w << ",";
	cout << endl;
	cout << "The following nodes need to be monitored by the observer to keep t_level untouched:" << endl;
	for (int i = 1; i < nodes; i++)
	{
		if (unsensed_nodes.find(i) != unsensed_nodes.end())
		{
			continue;
		}
		cout << i << ",";
	}
	cout << endl;
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
			cout << "Working with t:" << t << endl;
		PM->get_safe_nodes2(t, 2);
		//cout<<"after get_safe_nodes2"<<flush<<endl;
		if (PM->erase_disclosing_paths(t, disclosing_paths))
		{
			cout << "Finished, erasing path would leave graph disconnected,final T score:" << t << endl;
			break;
		}
		else if (debug)
		{
			cout << "AFTER ELIMINATING DISCLOSING PATHS(" << t << "), SURVIVING PATHS:" << endl;
			cout << "____________________________" << endl;
			PM->print_paths();
			cout << "____________________________" << endl;
		}
	}

	if (debug)
	{
		cout << "AFTER ELIMINATING DISCLOSING PATHS(" << t << "), SURVIVING PATHS:" << endl;
		cout << "____________________________" << endl;
		PM->print_paths();
		cout << "____________________________" << endl;
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
			cout << "Working with t:" << t << endl;
		PM->get_safe_nodes2(t, 2);
		if (PM->erase_disclosing_paths(t, disclosing_paths))
		{
			cout << "Finished, erasing path would leave graph disconnected,final T score:" << t << endl;
			break;
		}
		else if (debug)
		{
			cout << "AFTER ELIMINATING DISCLOSING PATHS(" << t << "), SURVIVING PATHS:" << endl;
			cout << "____________________________" << endl;
			PM->print_paths();
			cout << "____________________________" << endl;
		}
	}

	if (debug)
	{
		cout << "AFTER ELIMINATING DISCLOSING PATHS(" << t << "), SURVIVING PATHS:" << endl;
		cout << "____________________________" << endl;
		PM->print_paths();
		cout << "____________________________" << endl;
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
		cout << ",Iteration:," << i << ",current_time:," << current_time;
		//elminate and reset
		//PM->elminate_max_lambda_paths_U();
		PM->remove_maximal_lambda_paths();
		//calling set_path in eliminatate_max_lambda_paths_U
		//PM->set_path_values_from_cover_paths();
		int max_lambda = PM->get_max_lambda();
		if (max_lambda == 0)
		{
			cout << "Finished, max_lambda cannot be improved from 0!" << endl;
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
			cout << "iteration:," << i << ",max_lambda reduced from " << max_lambda << " to old lambda" << old_max_lambda << ",continue" << endl;
			PM->record_best_solution();
			PM->set_old_max_lambda(max_lambda);
		}
		else if (max_lambda == old_max_lambda)
		{
			cout << "iteration:" << i << ",no improvement from old lambda,keep removing all paths with lambda==max_lambda:" << max_lambda << endl;
		}
		else
		{
			cout << "iteration:," << i << ",max_lambda increased from:," << old_max_lambda << ",to max_lambda:," << max_lambda << ",hence eliminating all paths whose lambda>old_lambda:," << old_max_lambda << endl;
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
				//cout<<"\tnode is unreachable"<<endl;
				continue;
			}
			if (all_destinations.find(coord_map2[curr_pos]) != all_destinations.end())
			{	//Destinations have no outgoing edges
				//cout<<"\tnode is destination"<<endl;
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
				//cout<<"\tW"<<coord_map2[curr_pos]<<"->"<<coord_map2[W_pos]<<flush<<endl;
				//cout<<"\t["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<W_pos.first<<","<<W_pos.second<<"]"<<endl;
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
				//cout<<"\tE"<<coord_map2[curr_pos]<<"->"<<coord_map2[E_pos]<<flush<<endl;
				//cout<<"\t["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<E_pos.first<<","<<E_pos.second<<"]"<<endl;
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
				//cout<<"\tD["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<D_pos.first<<","<<D_pos.second<<"]"<<endl;
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
				//cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[C_pos]<<endl;
				//cout<<"\tC["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<C_pos.first<<","<<C_pos.second<<"]"<<endl;
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
				//cout<<"X\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
				//cout<<"\tX["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<X_pos.first<<","<<X_pos.second<<"]"<<endl;
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
				//cout<<"\tZ["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<Z_pos.first<<","<<Z_pos.second<<"]"<<endl;
				//cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
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
				//cout<<"\tA["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<A_pos.first<<","<<A_pos.second<<"]"<<endl;
				//cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
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
				//cout<<"\tQ["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<Q_pos.first<<","<<Q_pos.second<<"]"<<endl;
				//cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[Q_pos]<<endl;
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
			cout << "4-way Connectivity" << endl;
			connectivity = 4;
		}
		else if (str.find("octile") != string::npos)
		{
			cout << "8-way Connectivity" << endl;
			connectivity = 8;
		}

		max_x = max(size_t(max_x), str.length() - 1);

		// Line contains string of length > 0 then save it in vector
		//if(str.size() > 0)
		//  vecOfStrs.push_back(str);
		//create nodes
		//cout<<"row:"<<row<<",string:"<<str<<endl;
		if (row >= 0)
		{ //skipping first 4 rows
			//cout<<str<<",size:,"<<str.size()<<endl;
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
				//cout<<"node_map["<<node_counter-1<<"]:"<<node_map[node_counter-1].first<<","<<node_map[node_counter-1].second<<endl;
				//cout<<"coord_map2["<<row<<","<<i<<"]:,"<<coord_map2[make_pair(row,i)]<<endl;
			}
		}
		row++;
	}
	cout << "Connectivity:," << connectivity << endl;
	max_y = row;
	file_node_counter = node_counter - 1;
}
void random_orig_and_dest_placement_SaT()
{
	Destinations.clear();
	all_destinations.clear();
	dest_index.clear();
	srand(random_seed);
	int r = rand() % nodes_set.size(); // not _really_ random
	auto n = *select_random(nodes_set, r);
	start = n;
	cout << "start:," << start << ",D:,";
	while (Destinations.size() < random_positions)
	{
		int r = rand() % nodes_set.size(); // not _really_ random
		auto n = *select_random(nodes_set, r);
		if (n == start) //skip origin as possible destination
			continue;
		Destinations.push_back(n);
		all_destinations.insert(n);
		dest_index[n] = Destinations.size() - 1;
		cout << n << ",";
	}
	cout << endl;
}
void create_augmented_graph_from_SaT_file()
{ //
	cout << "hello create_augmented_graph,stocastic GO" << flush << endl;
	std::ifstream in(input_filename.c_str());
	std::ifstream in2(input_filename.c_str());
	lks.reserve(16000);
	string str;
	total_nodes = 0;
	int total_edges = 0;

	//vector<vector<Link> > boundary_edges;
	//vector<vector<Link> >dest_to_dest_edges;
	vector<vector<Link>> lks_boundary;
	vector<vector<Link>> lks_dest;

	//cout<<"hello1"<<flush<<endl;
	//First pass to get list of all nodes
	while (getline(in2, str))
	{
		vector<string> results;
		boost::split(results, str, boost::is_any_of(","), boost::token_compress_on);
		cout << "results[1]" << results[1] << endl;
		auto ret = nodes_set.insert(stoi(results[0]));
		if (ret.second == true)
			total_nodes++;
		ret = nodes_set.insert(stoi(results[1]));
		if (ret.second == true)
			total_nodes++;
		cout << "total_nodes=" << total_nodes << endl;
	}
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
		//cout<<results[0]<<flush<<endl;
		//cout<<results[1]<<flush<<endl;
		//cout<<results[2]<<flush<<endl;
		//cout<<"hello1a"<<flush<<endl;
		//cout<<","<<results[1]<<flush<<","<<results[2]<<flush<<endl;
		//cout<<"dest_index[results[1]]:"<<flush<<stoi(results[1])<<flush<<endl;
		if (stoi(results[1]) == start)
		{
			continue; //skipping origin incomming edges
		}
		if (stoi(results[0]) == stoi(results[1]))
		{ //no loops
			cout << "\t\tskipped loop,from," << stoi(results[0]) << ",to," << stoi(results[1]) << endl;
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
			cout << "Added reverse link for destination_index:," << dest_index[stoi(results[1])] << ",from:," << stoi(results[1]) << ",to:," << stoi(results[0]) << ",weight:," << stof(results[2]) << endl;
		}
		Link temp_link{stoi(results[0]), stoi(results[1]), stof(results[2])};

		lks.push_back(temp_link);
		total_edges++;
	}
	cout << "main_graph, total_edges:," << total_edges << endl;
	file_node_counter = total_nodes;
	//Now create Dijkstra Graph
	main_graph.set_number_vertices(file_node_counter);
	for (long i = 0; i < total_edges; i++)
	{
		main_graph.add_link(lks[i].u, lks[i].v, lks[i].weight);
	}

	main_graph.setv();
	cout << "total_nodes:" << total_nodes << ",main_graph:" << main_graph.get_number_vertices() << endl;
	main_dijkstra_alg = make_unique<DijkstraShortestPathAlg>(&main_graph);

	//BasePath* result_path = main_dijkstra_alg->get_shortest_path(
	//main_graph.get_vertex(start), main_graph.get_vertex(Destinations[0]));
	//main_graph.dump_edges();
	//for (auto edge : lks){
	//  cout<<"from:,"<<edge.u<<",->,"<<edge.v<<",weight:,"<<edge.weight<<endl;
	//}
	auto timenow =
		chrono::system_clock::to_time_t(chrono::system_clock::now());
	cout << "before determine_all_shortest_paths,time:" << ctime(&timenow) << endl;
	main_dijkstra_alg->determine_all_shortest_paths(main_graph.get_vertex(start), -1);
	timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
	cout << "after determine_all_shortest_paths,time:" << ctime(&timenow) << endl; //exit(1);

	//STOCASTIC CODE
	//First get current cost for any node
	cout << "node costs_reg:" << endl;
	int reachable_nodes = 0;
	for (auto node : nodes_set)
	{
		if (!main_dijkstra_alg->is_node_reachable(node))
		{
			reachable_nodes++;
			//cout<<"node "<<node<<" is not reachable"<<endl;
		}
		else
		{
			//cout<<"reachable node:"<<node<<",cost:"<<main_dijkstra_alg->get_cost(node)<<endl;
		}
	}
	cout << "Reachable nodes:" << reachable_nodes << endl;
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
			cout << "Destination:," << destination << ",is reachable from origin:," << start << ",continuing search" << endl;
		}
	}
	cout << "hola12" << endl;
	//States are pair <node,g>
	//we call Λ^{d}_{ss'} the probability that the target transits from v to v' if its destination is d  and has already incurred a cost x to reach v.
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
			//cout<<"adding current_graph reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<",vertices:"<<file_node_counter - 1 <<endl;
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
			//cout<<"adding current_graph_dest_to_dest reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<endl;
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
	cout<<"Added origin outgoing link for destination_index:,"<<dest<<",from:,"<<lks[i].u<<",to:,"<<lks[i].v<<",weight:,"<<lks[i].weight<<endl;*/
				continue;
			}
			//cout<<"adding current_graph_reverse_link edges["<<dest<<"]"<<",lks[i].v:"<<lks[i].v<<",lks[i].u:"<<lks[i].u<<"lks[i].weight:"<<lks[i].weight<<endl;
			current_graph[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
		}
		current_graph[dest].setv();
		current_graph[dest].set_number_vertices(total_nodes); //skipping origin
		cout << "total_nodes:" << total_nodes << ",current_graph[" << dest << "]:," << current_graph[dest].get_number_vertices() << endl;
		Dijkstra_algs.push_back(DijkstraShortestPathAlg(&(current_graph[dest])));

		cout << "adding bidirectional edges for dest_to_dest[" << dest << "]" << endl;
		for (long i = 0; i < total_edges; i++)
		{
			if (lks[i].u == start)
			{ //skipping origin outgoing edges
				continue;
			}
			//FOR CALCULATING GEODETIC DISTANCE, WE CREATE A LINK IN BOTH DIRECTIONS AS LONG AS
			//ONE DIRECTION EXISTS
			//cout<<"adding bidirectional edges for dest_to_dest["<<dest<<"]"<<",lks[i].u:"<<lks[i].u<<",lks[i].v:"<<lks[i].v<<"lks[i].weight:"<<lks[i].weight<<endl;
			current_graph_dest_to_dest[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
			current_graph_dest_to_dest[dest].add_link(lks[i].u, lks[i].v, lks[i].weight);
		}
		current_graph_dest_to_dest[dest].setv();
		current_graph_dest_to_dest[dest].set_number_vertices(total_nodes); //skipping origin
		cout << "total_nodes:" << total_nodes << ",current_graph_dest_to_dest[" << dest << "]:," << current_graph_dest_to_dest[dest].get_number_vertices() << endl;

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
			//cout<<"\t"<<result->get_path_string()<<",Weight:"<<result->Weight()<<endl;
			cout << "\tdest_to_dest,dest1:," << Destinations[dest] << ",dest2:" << Destinations[dest2] << ",Weight:" << result->Weight() << ",length:" << result->length() << endl;
			//BasePath* main_result1 = main_dijkstra_alg->get_shortest_path(
			//	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]));
			//    cout<<"\torig_to_dest1,start:,"<<start<<",dest2:"<<Destinations[dest]<<",Weight:"<<main_result1->Weight()<<",length:"<<main_result1->length()<<endl;
			//  BasePath* main_result2 = main_dijkstra_alg->get_shortest_path(
			//	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest2]));
			//    cout<<"\torig_to_dest2,start:,"<<start<<",dest2:"<<Destinations[dest2]<<",Weight:"<<main_result2->Weight()<<",length:"<<main_result2->length()<<endl;
		}
	}

	if (stocastic_go)
	{
		//STOCASTIC CODE Example for dest[0]
		//Get all paths from dest[0]
		//Get distance from node v to destination
		for (size_t dest = 0; dest < Destinations.size(); dest++)
		{
			cout << "hola,destination:" << Destinations[dest] << endl;
			Dijkstra_algs_dest_to_dest[dest].determine_all_shortest_paths(current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), -1);
			for (auto node : nodes_set)
			{
				if (!main_dijkstra_alg->is_node_reachable(node))
				{
					//cout<<"node "<<node<<" is not reachable from start"<<endl;
					continue;
				}
				else
				{
					cout << "\treachable node:" << node << ",cost to destination[" << Destinations[dest] << "]:" << Dijkstra_algs_dest_to_dest[dest].get_cost(node);
					cout << ",distance to start:," << main_dijkstra_alg->get_cost(node) << endl;
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
			//cout<<"start:"<<start<<",destination:"<<Destinations[dest]<<endl;
			cout<<"grandpatent_check="<<grandparent_check<<endl;
			main_dijkstra_alg->createAG(main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]), Budget,acyclic,grandparent_check);
		}
		main_dijkstra_alg->printAGFile();
		main_dijkstra_alg->populateAGnodes(Budget);

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
		simulation();
		exit(0);
	}
}
void calculate_min_path_strategy_AG()
{
	cout << "hola calculate_min_path_strategy" << endl;
	//Initialize Alpha
	Alpha.clear();
	map<pair<pair<int, float>, pair<int, float>>, float> temp_map;
	Alpha.assign(Destinations.size(), temp_map);
	//auto orig_node=main_graph.get_vertex(start);
	ofstream myfile;
	string filename = "AG" + to_string(main_graph.getVertexNum()) + "_S" + to_string(random_seed) + ".txt";
	myfile.open(filename);

	auto AG_Nodes = main_dijkstra_alg->getAGsubgraphs();
	auto AG_SG_Edges = main_dijkstra_alg->getSGEdges();
	cout << "AG_Edges.size:" << AG_SG_Edges->size() << endl;

	for (size_t dest = 0; dest < Destinations.size(); dest++)
	{
		cout << "\tworking on dest:," << Destinations[dest] << ",optimal distance:" << main_dijkstra_alg->get_cost(Destinations[dest]) << endl;
		if (Budget < main_dijkstra_alg->get_cost(Destinations[dest]))
		{
			cerr << "Budget is smaller than optimal distance from origin to destination!!!" << endl;
			exit(1);
		}
		for (auto node : (*AG_Nodes)[dest])
		{
			double OptPaths = 0;
			auto orig_vertex = main_graph.get_vertex(node.first->getID());
			cout << "\t\torig_vertex:" << orig_vertex->getID();
			cout << ",Edges size:" << (*AG_SG_Edges)[dest][orig_vertex->getID()].size() << endl;

			for (auto child_node : (*AG_SG_Edges)[dest][orig_vertex->getID()])
			{
				cout << "\t\t\tchild_node:," << child_node << endl;
				auto dest_vertex = main_graph.get_vertex(child_node);
				cout << "\t\tcost=" << node.second << "+" << main_graph.get_edge_weight(orig_vertex, dest_vertex);
				cout << "+" << Dijkstra_algs_dest_to_dest[dest].get_cost(child_node) << endl;
				double x = node.second + main_graph.get_edge_weight(orig_vertex, dest_vertex);
				double cost = x + Dijkstra_algs_dest_to_dest[dest].get_cost(child_node);
				//Probability is 0 if path is not optimal
				cout << "\t\t\tcost:" << cost << ",optimal_cost:" << main_dijkstra_alg->get_cost(Destinations[dest]) << endl;
				if (cost > main_dijkstra_alg->get_cost(Destinations[dest]))
				{
					cout << "\t\t\tpath is not optimal" << endl;
					//Alpha[dest][make_pair(make_pair(node.first->getID(),node.second),
					//make_pair(child_node,x))]=0;
				}
				else
				{
					cout << "\t\t\tpath is optimal" << endl;
					Alpha[dest][make_pair(make_pair(node.first->getID(), node.second),
										  make_pair(child_node, x))] = 1.0;
					OptPaths++;
				}
			}

			//Second pass to rationalize probabilities
			for (auto child_node : (*AG_SG_Edges)[dest][orig_vertex->getID()])
			{
				auto dest_vertex = main_graph.get_vertex(child_node);
				double x = node.second + main_graph.get_edge_weight(orig_vertex, dest_vertex);

				if (Alpha[dest][make_pair(make_pair(node.first->getID(), node.second),
										  make_pair(child_node, x))] > 0)
				{
					Alpha[dest][make_pair(make_pair(node.first->getID(), node.second),
										  make_pair(child_node, x))] = 1.0 / OptPaths;
				}
				cout << "Alpha[" << Destinations[dest] << ",(" << node.first->getID() << "," << node.second << ")->(" << child_node;
				cout << "," << x << ")]:" << Alpha[dest][make_pair(make_pair(node.first->getID(), node.second), make_pair(child_node, x))] << endl;
				myfile << Destinations[dest] << "," << node.first->getID() << "," << node.second << "," << child_node;
				myfile << "," << x << "," << Alpha[dest][make_pair(make_pair(node.first->getID(), node.second), make_pair(child_node, x))] << endl;
			}
		}
	}
	myfile.close();
}
void calculate_Gibbs_strategy_AG(double beta = 0.5, double budget = 100)
{
	cout << "hola calculate_Gibbs_strategy,beta:" << beta << ",Budget:" << budget << endl;
	//Initialize Alpha
	Alpha.clear();
	map<pair<pair<int, float>, pair<int, float>>, float> temp_map;
	Alpha.assign(Destinations.size(), temp_map);
	//auto orig_node=main_graph.get_vertex(start);
	ofstream myfile;
	string filename = "RNG_N" + to_string(main_graph.getVertexNum()) + "_S" + to_string(random_seed) + ".txt";
	myfile.open(filename);

	auto AG_Nodes = main_dijkstra_alg->getAGsubgraphs();
	auto AG_SG_Edges = main_dijkstra_alg->getSGEdges();
	cout << "AG_Edges.size:" << AG_SG_Edges->size() << flush << endl;

	for (size_t dest = 0; dest < Destinations.size(); dest++)
	{
		cout << "\tworking on dest:," << Destinations[dest] << ",optimal distance:" << main_dijkstra_alg->get_cost(Destinations[dest]) << endl;
		if (Budget < main_dijkstra_alg->get_cost(Destinations[dest]))
		{
			cerr << "Budget is smaller than optimal distance from origin to destination!!!" << endl;
			exit(1);
		}
		for (auto node : (*AG_Nodes)[dest])
		{
			double NormConst = 0;
			auto orig_vertex = main_graph.get_vertex(node.first->getID());
			cout << "\t\torig_vertex:" << orig_vertex->getID();
			cout << ",Edges size:" << (*AG_SG_Edges)[dest][orig_vertex->getID()].size() << endl;

			for (auto child_node : (*AG_SG_Edges)[dest][orig_vertex->getID()])
			{
				cout << "\t\t\tchild_node:," << child_node << endl;
				auto dest_vertex = main_graph.get_vertex(child_node);
				cout << "\t\tcost=" << node.second << "+" << main_graph.get_edge_weight(orig_vertex, dest_vertex);
				cout << "+" << Dijkstra_algs_dest_to_dest[dest].get_cost(child_node) << endl;
				double x = node.second + main_graph.get_edge_weight(orig_vertex, dest_vertex);
				double cost = x + Dijkstra_algs_dest_to_dest[dest].get_cost(child_node);
				//Probability is 0 if path is not feasible
				cout << "\t\t\tcost:" << cost << ",optimal_cost:" << main_dijkstra_alg->get_cost(Destinations[dest]) << endl;
				if (cost >= budget)
				{
					cout << "\t\t\tpath is not feasible for budget:" << budget << endl;
					//Alpha[dest][make_pair(make_pair(node.first->getID(),node.second),
					//make_pair(child_node,x))]=0;
				}
				else
				{
					cout << "\t\t\tpath is feasible" << endl;
					auto key = make_pair(make_pair(node.first->getID(), node.second),
										 make_pair(child_node, x));
					double exponent = main_graph.get_edge_weight(orig_vertex, dest_vertex) +
									  Dijkstra_algs_dest_to_dest[dest].get_cost(child_node);
					Alpha[dest][key] = exp(-beta * (exponent));
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

				if (Alpha[dest][key] > 0)
				{
					Alpha[dest][key] = Alpha[dest][key] / NormConst;
					cout << "Alpha[" << Destinations[dest] << ",(" << node.first->getID() << "," << node.second << ")->(" << child_node;
					cout << "," << x << ")]:" << Alpha[dest][make_pair(make_pair(node.first->getID(), node.second), make_pair(child_node, x))] << endl;
					myfile << Destinations[dest] << "," << node.first->getID() << "," << node.second << "," << child_node;
					myfile << "," << x << "," << Alpha[dest][make_pair(make_pair(node.first->getID(), node.second), make_pair(child_node, x))] << endl;
				}
			}
		}
	}
	myfile.close();
}

void simulation()
{
	auto AG_Edges = main_dijkstra_alg->getSGEdges();
	//1st choose Destination
	srand(random_seed);
	int dest = rand() % Destinations.size();
	int d = Destinations[dest];
	cout << endl
		 << "Simulation, randomly selected destination:" << d << endl;
	cout << "Simulated path:" << start;
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
		r = ((float)rand() / (RAND_MAX));
		float select_prob = 0;
		for (auto dest_node : (*AG_Edges)[dest][current_node])
		{
			auto current_vertex = main_graph.get_vertex(current_node);
			auto dest_vertex = main_graph.get_vertex(dest_node);
			float cost2 = cost1 + main_graph.get_edge_weight(current_vertex, dest_vertex);
			auto key = make_pair(make_pair(current_node, cost1), make_pair(dest_node, cost2));
			select_prob += Alpha[dest][key];
			//cout<<"\tAlpha["<<Destinations[dest]<<",("<<current_node<<","<<cost1<<")->("<<dest_node;
			//cout<<","<<cost2<<")]:"<<Alpha[dest][key];
			//cout<<"\trandom_throw:,"<<r<<",select_prob:"<<select_prob<<",dest_node:"<<dest_node<<endl;
			if (r <= select_prob)
			{ //
				current_node = dest_node;
				cost1 = cost2;
				//cout<<"\tchild is chosen,new current node:("<<current_node<<","<<cost2<<")"<<endl;
				cout << "," << current_node;
				myfile << "," << current_node;
				if (current_node == Destinations[dest])
				{
					//cout<<"Finished, destination found"<<endl;
					cout << endl;
					myfile << endl;
					return;
				}
				break;
			}
		}
	}
	myfile.close();
}

void create_grid_nodes_and_edges_from_SaT_file()
{
	cout << "hello SaT" << flush << endl;
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

	//cout<<"hello1"<<flush<<endl;
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
	cout << "total_nodes=" << total_nodes << endl;
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
		//cout<<results[0]<<flush<<endl;
		//cout<<results[1]<<flush<<endl;
		//cout<<results[2]<<flush<<endl;
		//cout<<"hello1a"<<flush<<endl;
		//cout<<","<<results[1]<<flush<<","<<results[2]<<flush<<endl;
		//cout<<"dest_index[results[1]]:"<<flush<<stoi(results[1])<<flush<<endl;
		if (stoi(results[1]) == start)
		{
			continue; //skipping origin incomming edges
		}
		if (stoi(results[0]) == stoi(results[1]))
		{ //no loops
			cout << "\t\tskipped loop,from," << stoi(results[0]) << ",to," << stoi(results[1]) << endl;
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
			cout << "Added reverse link for destination_index:," << dest_index[stoi(results[1])] << ",from:," << stoi(results[1]) << ",to:," << stoi(results[0]) << ",weight:," << stof(results[2]) << endl;
		}
		Link temp_link{stoi(results[0]), stoi(results[1]), stof(results[2])};

		lks.push_back(temp_link);
		total_edges++;
	}
	cout << "main_graph, total_edges:," << total_edges << endl;
	file_node_counter = total_nodes;
	//Now create Dijkstra Graph
	main_graph.set_number_vertices(file_node_counter);
	for (long i = 0; i < total_edges; i++)
	{
		main_graph.add_link(lks[i].u, lks[i].v, lks[i].weight);
	}

	main_graph.setv();
	cout << "total_nodes:" << total_nodes << ",main_graph:" << main_graph.get_number_vertices() << endl;
	main_dijkstra_alg = make_unique<DijkstraShortestPathAlg>(&main_graph);

	//BasePath* result_path = main_dijkstra_alg->get_shortest_path(
	//main_graph.get_vertex(start), main_graph.get_vertex(Destinations[0]));
	//main_graph.dump_edges();
	//for (auto edge : lks){
	//  cout<<"from:,"<<edge.u<<",->,"<<edge.v<<",weight:,"<<edge.weight<<endl;
	//}
	auto timenow =
		chrono::system_clock::to_time_t(chrono::system_clock::now());
	cout << "before determine_all_shortest_paths,time:" << ctime(&timenow) << endl;
	main_dijkstra_alg->determine_all_shortest_paths(main_graph.get_vertex(start), -1);
	timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
	cout << "after determine_all_shortest_paths,time:" << ctime(&timenow) << endl; //exit(1);
	cout << "node costs:" << endl;
	int reachable_nodes = 0;
	for (auto node : nodes_set)
	{
		if (!main_dijkstra_alg->is_node_reachable(node))
		{
			reachable_nodes++;
			//cout<<"node "<<node<<" is not reachable"<<endl;
		}
		else
		{
			//cout<<"reachable node:"<<node<<",cost:"<<main_dijkstra_alg->get_cost(node)<<endl;
		}
	}
	cout << "Reachable nodes:" << reachable_nodes << endl;

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
			cout << "Destination:," << destination << ",is reachable from origin:," << start << ",continuing search" << endl;
		}
	}
	//cout<<"\tdest_to_dest,origin:,"<<start<<",dest:"<<Destinations[0]<<",Weight:"<<result_path->Weight()<<",length:"<<result_path->length()<<endl;exit(1);

	//main_dijkstra_alg->dump_edges();
	//BasePath* main_result1 = main_dijkstra_alg->get_shortest_path(
	//main_graph.get_vertex(start), main_graph.get_vertex(Destinations[1]));
	//cout<<"\torig_to_dest1,start:,"<<start<<",dest2:"<<Destinations[0]<<",Weight:"<<main_result1->Weight()<<",length:"<<main_result1->length()<<endl;
	//BasePath* main_result2 = main_dijkstra_alg->get_shortest_path(
	//main_graph.get_vertex(start), main_graph.get_vertex(Destinations[1]));
	//cout<<"\torig_to_dest2,start:,"<<start<<",dest2:"<<Destinations[1]<<",Weight:"<<main_result2->Weight()<<",length:"<<main_result2->length()<<endl;

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
			//cout<<"adding current_graph reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<",vertices:"<<file_node_counter - 1 <<endl;
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
			//cout<<"adding current_graph_dest_to_dest reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<endl;
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
	cout<<"Added origin outgoing link for destination_index:,"<<dest<<",from:,"<<lks[i].u<<",to:,"<<lks[i].v<<",weight:,"<<lks[i].weight<<endl;*/
				continue;
			}
			//cout<<"adding current_graph_reverse_link edges["<<dest<<"]"<<",lks[i].v:"<<lks[i].v<<",lks[i].u:"<<lks[i].u<<"lks[i].weight:"<<lks[i].weight<<endl;
			current_graph[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
		}
		current_graph[dest].setv();
		current_graph[dest].set_number_vertices(total_nodes); //skipping origin
		cout << "total_nodes:" << total_nodes << ",current_graph[" << dest << "]:," << current_graph[dest].get_number_vertices() << endl;
		Dijkstra_algs.push_back(DijkstraShortestPathAlg(&(current_graph[dest])));

		cout << "adding bidirectional edges for dest_to_dest[" << dest << "]" << endl;
		for (long i = 0; i < total_edges; i++)
		{
			if (lks[i].u == start)
			{ //skipping origin outgoing edges
				continue;
			}
			//FOR CALCULATING GEODETIC DISTANCE, WE CREATE A LINK IN BOTH DIRECTIONS AS LONG AS
			//ONE DIRECTION EXISTS
			//cout<<"adding bidirectional edges for dest_to_dest["<<dest<<"]"<<",lks[i].u:"<<lks[i].u<<",lks[i].v:"<<lks[i].v<<"lks[i].weight:"<<lks[i].weight<<endl;
			current_graph_dest_to_dest[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
			current_graph_dest_to_dest[dest].add_link(lks[i].u, lks[i].v, lks[i].weight);
		}
		current_graph_dest_to_dest[dest].setv();
		current_graph_dest_to_dest[dest].set_number_vertices(total_nodes); //skipping origin
		cout << "total_nodes:" << total_nodes << ",current_graph_dest_to_dest[" << dest << "]:," << current_graph_dest_to_dest[dest].get_number_vertices() << endl;

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
			//cout<<"\t"<<result->get_path_string()<<",Weight:"<<result->Weight()<<endl;
			cout << "\tdest_to_dest,dest1:," << Destinations[dest] << ",dest2:" << Destinations[dest2] << ",Weight:" << result->Weight() << ",length:" << result->length() << endl;
			//BasePath* main_result1 = main_dijkstra_alg->get_shortest_path(
			//	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]));
			//    cout<<"\torig_to_dest1,start:,"<<start<<",dest2:"<<Destinations[dest]<<",Weight:"<<main_result1->Weight()<<",length:"<<main_result1->length()<<endl;
			//  BasePath* main_result2 = main_dijkstra_alg->get_shortest_path(
			//	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest2]));
			//    cout<<"\torig_to_dest2,start:,"<<start<<",dest2:"<<Destinations[dest2]<<",Weight:"<<main_result2->Weight()<<",length:"<<main_result2->length()<<endl;
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

	cout << "Creating FO graph from Partially Ovservable graph,PO_level:," << PO_level << endl;
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
	//cout<<"FO_nodes:,"; for (auto node : FO_nodes) cout<<node<<",";
	//cout<<endl;

	//first create regular graph as if everything was observable
	temp_graph.set_number_vertices(main_graph_nodes.size());
	for (long i = 0; i < grid_edges.size(); i++)
	{
		temp_graph.add_link(grid_edges[i].u, grid_edges[i].v, grid_edges[i].weight);
		/*if(grid_edges[i].u==1795){
      cout<<"1795->"<<grid_edges[i].v<<",w:,"<<grid_edges[i].weight<<endl;
      if(FO_nodes.find(grid_edges[i].v)!=FO_nodes.end()){
	cout<<",visible"<<endl;
      }
      else
	cout<<",invisible"<<endl;
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
			//cout<<"grid_edge:,"<<edge.u<<",->,"<<edge.v<<",origin visible"<<endl;
			//}
			//edge is between 2 observable nodes so skip it but copy it to FO graph
			if (FO_nodes.find(edge.v) != FO_nodes.end())
			{
				//cout<<"\tFO_edge:,"<<edge.u<<",->,"<<edge.v<<",destination visible"<<endl;
				FO_edges.push_back(edge);
				best_known_distances[make_pair<int, int>(int(edge.u), int(edge.v))] = edge.weight;
				temp_graph.remove_edge(make_pair<int, int>(int(edge.u), int(edge.v)));
				visible_outgoing_edges[edge.u].push_back(edge);
				continue;
			}
			//cout<<"FO_edge:,"<<edge.u<<",->,"<<edge.v<<",destination invisible"<<endl;
			temp_graph.remove_edge(make_pair<int, int>(int(edge.u), int(edge.v)));
			NC_nodes.insert(edge.u);
			visible_outgoing_edges[edge.u].push_back(edge);
			continue;
		}
		else
		{
			//cout<<"grid_edge:,"<<edge.u<<",->,"<<edge.v<<",origin invisible"<<endl;
			current_nodes.insert(edge.u);
			current_nodes.insert(edge.v);
		}
	}

	//Third, create a graph for each visible node as an origin
	//to found the costs associated to reaching other nodes
	cout << "Nodes which require new connections:,NC_nodes:," << NC_nodes.size() << ",FO_nodes:," << FO_nodes.size() << endl;
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
			cout << "\t\tNC_counter:" << counter << ",total:," << NC_nodes.size() << ",time:," << current_time << ",aggregated_dijkstra_time:," << aggregated_dijkstra_time << endl;
			start_time = std::chrono::high_resolution_clock::now();
			//aggregated_dijkstra_time=std::chrono::milliseconds(0);
		}

		//auto subgraph_edges=temp_edges;
		for (auto edge : visible_outgoing_edges[origin])
		{
			temp_graph.recover_removed_edge(make_pair<int, int>(int(edge.u), int(edge.v)));
			//cout<<"recovered edge from "<<edge.u<<" to "<<edge.v<<endl;
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
		// cout<<"Path[52,82]:,W:,"<<path->Weight()<<",";path->PrintOut(cout);
		//}
		//auto end_dijkstra_time = std::chrono::high_resolution_clock::now();aggregated_dijkstra_time += end_dijkstra_time-start_dijkstra_time;
		//Now
		//cout<<"\t\tNC_counter:,"<<counter<<",dijkstra_time:,"<<(end_dijkstra_time-start_dijkstra_time).count()<<",size:"<<temp_dijkstra_alg->get_perim_size()<<endl;
		//start_dijkstra_time=std::chrono::high_resolution_clock::now();
		temp_dijkstra_alg->improve_distances(&FO_nodes, &FO_edges, &best_known_distances);
		//end_dijkstra_time = std::chrono::high_resolution_clock::now();aggregated_dijkstra_time += end_dijkstra_time-start_dijkstra_time;
		//cout<<"\t\tNC_counter2:,"<<counter<<",improve_distance_time:,"<<(end_dijkstra_time-start_dijkstra_time).count()<<",size2:"<<temp_dijkstra_alg->get_perim_size()<<endl;
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
	cout << "finished building FD_graph,time:," << current_time << ",visible_nodes:," << FO_nodes.size() << ",FO_edges:" << FO_edges.size() << endl;

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
		//cout<<"cleaning,FO_edge:,"<<edge.u<<",->,"<<edge.v<<endl;
		if (edge.v == start)
		{
			continue;
		}
		if (edge.u == edge.v)
		{ //no loops
			cout << "\t\tskipped loop,from," << edge.u << ",to," << edge.v << endl;
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
			//cout<<"Added reverse link for destination_index:,"<<dest_index[edge.v]<<",from:,"<<edge.v<<",to:,"<<edge.u<<",weight:,"<<edge.weight<<endl;
		}
		Link temp_link{edge.u, edge.v, edge.weight};

		lks.push_back(temp_link);
		total_edges++;
	}

	cout << "main_graph, total_edges:," << total_edges << ",grid_edges:," << grid_edges.size() << endl;
	//Check if any nodes FO nodes are unreachable
	unordered_set<int> final_FO_nodes;
	for (auto edge : FO_edges)
	{
		//if(edge.u==1795)
		//cout<<"final_add,FO_edge:,"<<edge.u<<",->,"<<edge.v<<endl;
		final_FO_nodes.insert(edge.u);
		final_FO_nodes.insert(edge.v);
	}

	for (auto node : FO_nodes)
	{
		if (final_FO_nodes.find(node) == final_FO_nodes.end())
		{
			cout << "\t\t FO_node:," << node << ",not reachable in final_FO_nodes, is this right?" << endl;
		}
	}

	//Now create Dijkstra Graph
	main_graph.clear();
	//main_graph.set_number_vertices( file_node_counter );
	main_graph.set_number_vertices(main_graph_nodes.size());
	for (auto edge : FO_edges)
	{
		main_graph.add_link(edge.u, edge.v, edge.weight);
		//cout<<"main.u:,"<<edge.u<<",main.v:,"<<edge.v<<",main.w:,"<<edge.weight<<",code:,"<<main_graph.get_edge_code(edge.u,edge.v)<<endl;
		//cout<<"main.u:,"<<edge.u<<",main.v:,"<<edge.v<<",main.w:,"<<edge.weight<<endl;
	}
	main_graph.setv();
	main_graph.set_number_vertices(main_graph_nodes.size());
	total_edges = FO_edges.size();
	cout << "total_FO_nodes:" << FO_nodes.size() << ",main_graph:" << main_graph.get_number_vertices() << ",final_FO_nodes:," << final_FO_nodes.size() << endl;
	main_dijkstra_alg = make_unique<DijkstraShortestPathAlg>(&main_graph);

	auto timenow =
		chrono::system_clock::to_time_t(chrono::system_clock::now());
	cout << "before determine_all_shortest_paths,time:" << ctime(&timenow) << endl;
	main_dijkstra_alg->determine_all_shortest_paths(main_graph.get_vertex(start), -1);

	/*auto path=main_dijkstra_alg->recover_shortest_perim_path(84);
  cout<<"Path["<<start<<",84]:,W:,"<<path->Weight()<<",";path->PrintOut(cout);
  auto path2=main_dijkstra_alg->recover_shortest_perim_path(20);
  cout<<"Path["<<start<<",20]:,W:,"<<path2->Weight()<<",";path2->PrintOut(cout);
  auto path3=main_dijkstra_alg->recover_shortest_perim_path(52);
  cout<<"Path["<<start<<",52]:,W:,"<<path3->Weight()<<",";path3->PrintOut(cout);*/

	cout << "perim_size:," << main_dijkstra_alg->get_perim_size() << endl;
	timenow = chrono::system_clock::to_time_t(chrono::system_clock::now());
	cout << "after determine_all_shortest_paths,time:" << ctime(&timenow) << endl; //exit(1);

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
			cout << "Destination:," << destination << ",is reachable from origin:," << start << ",continuing search" << endl;
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
			//cout<<"adding current_graph reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<",vertices:"<<file_node_counter - 1 <<endl;
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
			//cout<<"adding current_graph_dest_to_dest reverse edges["<<dest<<"]"<<",add.u:"<<add.u<<",add.v:"<<add.v<<"add.weight:"<<add.weight<<endl;
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
	cout<<"Added origin outgoing link for destination_index:,"<<dest<<",from:,"<<lks[i].u<<",to:,"<<lks[i].v<<",weight:,"<<lks[i].weight<<endl;*/
				continue;
			}
			//cout<<"adding current_graph_reverse_link edges["<<dest<<"]"<<",lks[i].v:"<<lks[i].v<<",lks[i].u:"<<lks[i].u<<"lks[i].weight:"<<lks[i].weight<<endl;
			current_graph[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
		}
		current_graph[dest].setv();
		current_graph[dest].set_number_vertices(main_graph_nodes.size());
		cout << "total_nodes:" << total_nodes << ",current_graph[" << dest << "]:," << current_graph[dest].get_number_vertices() << endl;
		Dijkstra_algs.push_back(DijkstraShortestPathAlg(&(current_graph[dest])));

		cout << "adding bidirectional edges for dest_to_dest[" << dest << "]" << endl;
		for (long i = 0; i < total_edges; i++)
		{
			if (lks[i].u == start)
			{ //skipping origin outgoing edges
				continue;
			}
			//FOR CALCULATING GEODETIC DISTANCE, WE CREATE A LINK IN BOTH DIRECTIONS AS LONG AS
			//ONE DIRECTION EXISTS
			//cout<<"adding bidirectional edges for dest_to_dest["<<dest<<"]"<<",lks[i].u:"<<lks[i].u<<",lks[i].v:"<<lks[i].v<<"lks[i].weight:"<<lks[i].weight<<endl;
			current_graph_dest_to_dest[dest].add_link(lks[i].v, lks[i].u, lks[i].weight);
			current_graph_dest_to_dest[dest].add_link(lks[i].u, lks[i].v, lks[i].weight);
		}
		current_graph_dest_to_dest[dest].setv();
		current_graph_dest_to_dest[dest].set_number_vertices(main_graph_nodes.size());
		cout << "total_nodes:" << total_nodes << ",current_graph_dest_to_dest[" << dest << "]:," << current_graph_dest_to_dest[dest].get_number_vertices() << endl;

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
			//cout<<"\t"<<result->get_path_string()<<",Weight:"<<result->Weight()<<endl;
			//cout<<"\tdest_to_dest,dest1:,"<<Destinations[dest]<<",dest2:"<<Destinations[dest2]<<",Weight:"<<result->Weight()<<",length:"<<result->length()<<endl;
			//BasePath* main_result1 = main_dijkstra_alg->get_shortest_path(
			//	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest]));
			//    cout<<"\torig_to_dest1,start:,"<<start<<",dest2:"<<Destinations[dest]<<",Weight:"<<main_result1->Weight()<<",length:"<<main_result1->length()<<endl;
			//  BasePath* main_result2 = main_dijkstra_alg->get_shortest_path(
			//	  main_graph.get_vertex(start), main_graph.get_vertex(Destinations[dest2]));
			//    cout<<"\torig_to_dest2,start:,"<<start<<",dest2:"<<Destinations[dest2]<<",Weight:"<<main_result2->Weight()<<",length:"<<main_result2->length()<<endl;
		}
	}
}

void create_grid_edges_from_file()
{
	cout << "max_x:" << max_x << ",max_y:" << max_y << endl;
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

	cout << "Destinations:";
	for (auto dest : Destinations)
	{
		cout << "," << dest;
		cout << ",pos:," << node_map[dest].first << "," << node_map[dest].second << endl;
	}
	//for(auto it : all_destinations) cout<<","<<it;
	cout << endl;
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
				//cout<<"\tnode is unreachable"<<endl;
				continue;
			}
			if (all_destinations.find(coord_map2[curr_pos]) != all_destinations.end())
			{	//Destinations have no outgoing edges
				//cout<<"\tnode is destination"<<endl;
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
				//cout<<"\tW"<<coord_map2[curr_pos]<<"->"<<coord_map2[W_pos]<<flush<<endl;
				//cout<<"\t["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<W_pos.first<<","<<W_pos.second<<"]"<<endl;
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
				//cout<<"\tE"<<coord_map2[curr_pos]<<"->"<<coord_map2[E_pos]<<flush<<endl;
				//cout<<"\t["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<E_pos.first<<","<<E_pos.second<<"]"<<endl;
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
				//cout<<"\tD["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<D_pos.first<<","<<D_pos.second<<"]"<<endl;
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
				//cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[C_pos]<<endl;
				//cout<<"\tC["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<C_pos.first<<","<<C_pos.second<<"]"<<endl;
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
				//cout<<"X\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
				//cout<<"\tX["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<X_pos.first<<","<<X_pos.second<<"]"<<endl;
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
				//cout<<"\tZ["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<Z_pos.first<<","<<Z_pos.second<<"]"<<endl;
				//cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
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
				//cout<<"\tA["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<A_pos.first<<","<<A_pos.second<<"]"<<endl;
				//cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[X_pos]<<endl;
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
				//cout<<"\tQ["<<curr_pos.first<<","<<curr_pos.second<<"]->["<<Q_pos.first<<","<<Q_pos.second<<"]"<<endl;
				//cout<<"\t"<<coord_map2[curr_pos]<<"->"<<coord_map2[Q_pos]<<endl;
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
	cout << "start_node:," << start << ",start_pos:" << start_pos.first << "," << start_pos.second << ",grid nodes:," << file_node_counter << ",grid_edges:," << grid_edges.size() << endl;

	int edges_counter = 0;
	for (auto link : grid_edges)
	{
		//cout<<"\t"<<link.u<<"->"<<link.v<<",w:"<<link.weight<<endl;
		edges_counter++;
	}
	cout << "edges_counter:" << edges_counter << endl;

	//Now create Dijkstra Graph
	main_graph.set_number_vertices(file_node_counter + 1);
	for (long i = 0; i < grid_edges.size(); i++)
	{
		main_graph.add_link(grid_edges[i].u, grid_edges[i].v, grid_edges[i].weight);
	}
	main_graph.setv();
	main_graph.set_number_vertices(file_node_counter + 1);

	cout << "maing_graph_nodes_counter:" << main_graph_nodes.size() << ",main_graph:," << main_graph.getVertexNum() << endl;

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

		cout << "boundary_edges_nodes[:" << dest << "[:," << boundary_edges_nodes[dest].size() << ",graph:," << current_graph[dest].getVertexNum() << endl;
		//cout<<"boundary_edges_nodes[:"<<dest<<"[:,"<<main_graph.getVertexNum()<<",graph:,"<<current_graph[dest].getVertexNum()<<endl;

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
			//cout<<"\t"<<result->get_path_string()<<",Weight:"<<result->Weight()<<endl;
			cout << "\tdest_to_dest,dest1:," << Destinations[dest] << ",dest2:" << Destinations[dest2] << ",Weight:" << result->Weight() << ",length:" << result->length() << endl;
			if (result->length() == 0)
			{
				cerr << "Cannot find path from dest1:," << Destinations[dest] << ",to dest:," << Destinations[dest2] << endl;
				exit(1);
			}
		}
	}
	//exit(1);
	//cout<<"getting shortest destination path to all nodes from origin in main graph"<<endl;
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
			//cout<<n1<<"->"<<counter<<endl;
			//  if(all_destinations.find(counter)!=all_destinations.end()){
			//	cout<<FORERED<<"D "<<RESETTEXT;
			//    }
			//  else if(counter==start){
			//	cout<<FOREGRN<<"O "<<RESETTEXT;
			//    }
			//  else{
			//if(debug)
			//cout<<counter<<"\t";
			//cout<<"_ ";
			//}
			//cout<<RESETTEXT;
			counter++;
		}
		//if(debug)
		//cout<<endl;
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

	cout << "Destination:";
	for (auto it : all_destinations)
		cout << "," << it;
	cout << endl;

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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
					if (debug)
						outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
				}
			}
		}
	}
	outfile.close();
	//cout<<"nodes:,"<<n*n<<"total_edges:,"<<counter<<",max_edges:,"<<edges<<endl;
}

void create_grid_edges_boundary(int n, vector<Link> &grid_edges, int from_dest)
{
	string dest_to_dest_filename = "d";
	dest_to_dest_filename += to_string(from_dest);
	dest_to_dest_filename += "_boundary.dot";
	std::ofstream outfile(dest_to_dest_filename);

	cout << "Destination:";
	for (auto it : all_destinations)
		cout << "," << it;
	cout << endl;

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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
					if (debug)
						outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
				}
			}
		}
	}
	outfile.close();
	//cout<<"nodes:,"<<n*n<<"total_edges:,"<<counter<<",max_edges:,"<<edges<<endl;
}

void create_grid_edges(int n, vector<Link> &grid_edges)
{
	std::ofstream outfile("LastGraph.dot");

	set<int> dest_set;
	cout << "Destinations:";
	for (auto it : all_destinations)
		cout << "," << it;
	cout << endl;

	//for(auto it : coord_map) cout<<it.first<<"->"<<it.second<<endl;

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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
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
					//cout<<"\tn1:"<<n1<<",n2:"<<n2<<",temp_link:"<<coord_map[n1]<<"->"<<coord_map[n2]<<endl;
					if (debug)
						outfile << "\"" << coord_map[n1] << "\"->\"" << coord_map[n2] << "\"" << endl;
				}
			}
		}
	}
	outfile.close();
	//cout<<"nodes:,"<<n*n<<"total_edges:,"<<counter<<",max_edges:,"<<edges<<endl;
}

//Choose origin and destinations according to min_dist_orig,max_dist_orig,min_dest_dist,max_dest_dist
void random_dest_and_origin_placement()
{
	//cout<<"hola random_dest_and_origin_placement"<<flush<<endl;
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
	//cout<<"hola00"<<flush<<endl;
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
		cout << "available_pos for first_dest:";
		for (auto x : available_pos)
			cout << x << ",";
		cout << endl;
	}
	//First choose a random destination
	Destinations.clear();
	auto first_pos_index = rand() % available_pos.size();
	auto first_dest = *select_random(available_pos, first_pos_index);
	OrigDests.push_back(get_manhattan_dist(start, first_dest));
	avg_dest_orig_dest += OrigDests.back();
	//cout<<"first_dest:"<<first_dest<<endl;
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
					//cout<<"\t\tcand_node:,"<<cand_node<<",node_dest:,"<<node_dest<<",dist:,"<<dist<<",position added,counter:,"<<counter++<<endl;
					break;
				}
			}
		}
		//Exit if we have no available choice
		if (candidate_available_pos.size() == 0)
		{ //there are no more possible random nodes under existing constraints
			cout << "random_positions unachivable, exiting" << endl;
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
				cout << "\t\t dest:" << next_dest << "is " << min_intra_dist << " from closest destination" << endl;
				break; //
			}
			//else if(debug){
			//	cout<<"\t\t dest:"<<next_dest<<"is "<<min_intra_dist<<" from closest destination,D="<<Destinations<<endl;
			//     }
			if (counter_select++ > 100000)
			{
				cerr << "cannot find random destination at max destination distance=" << min_dist_to_dest << ",min:" << max_dist_to_dest << ",N=" << N << ",R=" << random_positions << ",already_chosen:" << Destinations << endl;
				exit(25);
			}
		}
		Destinations.push_back(next_dest);
		available_pos.erase(next_dest);

		//cout<<"candidate_available_pos:"<<candidate_available_pos<<",chosen dests:"<<Destinations<<endl;
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

	cout << "\n\n\n"
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
				cout << FORERED << "D " << RESETTEXT;
			}
			else if (counter == start)
			{
				cout << FOREGRN << "O " << RESETTEXT;
			}
			else
			{
				if ((y == 0) || (y == N - 1))
				{
					cout << ". ";
				}
				else if ((x == 0) || (x == N - 1))
				{
					cout << ". ";
				}
				else
				{
					cout << "  ";
				}
			}
			counter++;
		}
		cout << endl;
	}
	cout << "\n\nstart_pos:," << start << ",Random_Destinations:" << Destinations << ",random_seed:" << random_seed << ",add_lambda:" << lambda_add << endl;
	exit(1);
}
void from_file_random_dest_and_origin_selection()
{
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

	cout << "\n\nstart_pos:," << start << ",Random_Destinations:" << Destinations << ",random_seed:" << random_seed << ",add_lambda:" << lambda_add << endl;
}

int main(int argc, char *argv[])
{

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
	//cout<<"timer:"<<timer.format()<<endl;
	//int edges=2*(2*n*n-2*n);//edges in a square graph
	bool optimize_target = false;
	bool optimize_target_u = false;
	bool optimize_observer = false;
	bool optimize_budget = false;
	bool SaT_map = false;
	bool AugGraph = false;

	vector<vector<Link>> grid_vect_orig_dest;
	vector<vector<Link>> grid_vect_dest_to_dest;

	for (int i = 0; i < argc; i++)
	{
		std::string action(argv[i]);
		if (action == "AG")
		{
			AugGraph = true;
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
			cout << "random_origin_only=true" << endl;
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
			cout << "input_filename:" << input_filename << endl;
		}
		if (action.find("filename_SaT=") != string::npos)
		{
			string temp = action.substr(13, action.length());
			input_filename = temp;
			cout << "input_SaT_filename:" << input_filename << endl;
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
			cout << "PO_level:" << PO_level << endl;
		}
		if (action.find("IL=") != string::npos)
		{
			string temp = action.substr(3, action.length());
			input_lambda = stof(temp);
			if (debug)
				cout << "IL:" << input_lambda << endl;
		}
		if (action.find("LA=") != string::npos)
		{
			string temp = action.substr(3, action.length());
			lambda_add = stoi(temp);
			if (debug)
				cout << "LA:" << lambda_add << endl;
		}
		if (action.find("MDD1=") != string::npos)
		{
			string temp = action.substr(5, action.length());
			min_dist_to_dest = stof(temp);
			if (debug)
				cout << "min_dist_to_dest:" << min_dist_to_dest << endl;
		}
		if (action.find("MDD2=") != string::npos)
		{
			string temp = action.substr(5, action.length());
			max_dist_to_dest = stof(temp);
			if (debug)
				cout << "max_dist_to_dest:" << max_dist_to_dest << endl;
		}
		if (action.find("C=") != string::npos)
		{
			string temp = action.substr(2, action.length());
			C = stoi(temp);
			if (debug)
				cout << "C:" << C << endl;
		}
		if (action.find("O=") != string::npos)
		{
			string temp = action.substr(2, action.length());
			optimal_limit = stof(temp);
			cout << "O:" << optimal_limit << endl;
		}
		if (action.find("K=") != string::npos)
		{
			string temp = action.substr(2, action.length());
			K = stoi(temp);
			if (debug)
				cout << "K:" << C << endl;
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
			cout << "start=" << start << endl;
		}
		if (action.find("R=") != string::npos)
		{
			if (debug)
				cout << "Using random destinations" << endl;
			random_destinations = true;
			string temp = action.substr(2, action.length());
			random_positions = stoi(temp);
		}
		if (action.find("S=") != string::npos)
		{
			string temp = action.substr(2, action.length());
			random_seed = stoi(temp);
			cout << "Using random seed=" << random_seed << endl;
		}
		if (action.find("N=") != string::npos)
		{
			Destinations.clear();
			string temp = action.substr(2, action.length());
			N = stoi(temp);
			cout << "N:" << N << endl;
		}
		if (action.find("Stg=") != string::npos)
		{
			string temp = action.substr(4, action.length());
			strategy = stoi(temp);
			cout << "strategy:" << strategy << endl;
		}
		if (action.find("Beta=") != string::npos)
		{
			string temp = action.substr(5, action.length());
			Beta = stof(temp);
			cout << "Beta:" << Beta << endl;
		}
		if (action.find("Budget=") != string::npos)
		{
			string temp = action.substr(7, action.length());
			Budget = stof(temp);
			cout << "Budget:" << Budget << endl;
		}
		if (action == "cyclic")
		{
			acyclic = false;
		}
		if (action == "grandparent_check")
		{
			acyclic = false;
			grandparent_check = true;
			cout<<"grandparent_check="<<grandparent_check<<endl;
		}
	}
	//Check input file exists
	ifstream f(input_filename);
	if(!f.good()){
		cerr<<"graph input file:"<<input_filename<<" does not exist!"<<endl;
		exit(31);
	}
	int positions = N * N - 1;
	if (optimize_target == false &&
		optimize_observer == false &&
		optimize_target_u == false &&
		optimize_budget == false)
	{
		cout << "Choose either optimize_target, optimize_observer or optimize_budget" << endl;
		exit(1);
	}
	else if (optimize_target && (optimize_observer || optimize_budget) ||
			 optimize_budget && (optimize_target || optimize_observer))
	{
		cout << "Choose only optimize_target or optimize_observer" << endl;
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
			create_augmented_graph_from_SaT_file();
		}
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
		cout << "Creating main grid with original origin" << endl;
		create_grid_edges(N, grid_vect);

		lks.resize(grid_vect.size());
		for (size_t i = 0; i < grid_vect.size(); i++)
		{
			lks[i] = grid_vect[i];
		}
		array_size = sizeof(lks) / sizeof(lks[0]);
		nodes = N * N;
		//Create main graph
		main_graph.set_number_vertices(nodes);
		cout << "main_graph numberof nodes:" << nodes << endl;

		for (int i = 0; i < array_size; i++)
		{
			//cout<<"Adding link from:"<< lks[ i ].u<<",to:"<<lks[i].v<<",weight:"<<lks[i].weight<<endl;
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
			cout << "Creating grid with origin at destination:" << Destinations[dest] << endl;
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
	//cout<<"YenNodes:"<<nodes<<endl;
	auto PM = std::make_shared<PathMatrix>();
	//set<int> unsensed_nodes;
	//unsensed_nodes->insert(12);
	//PM->set_W(unsensed_nodes);
	PM->set_C(C);
	PM->set_min_cal(min_cal);
	//cout<<"C:,"<<C<<",min_cal:,"<<min_cal<<endl;
	cout << "Destinations:" << Destinations << endl;

	Link *lks_array = &lks[0];
	if (!optimize_budget)
	{
		for (auto end : Destinations)
		{
			testYenAlg(K,
					   lks_array,
					   array_size,
					   start,
					   end,
					   nodes,
					   PM);

			auto end_time = chrono::high_resolution_clock::now();
			auto diff = end_time - start_time;
			start_time = end_time; //restart
			auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
			cout << "Time(ms) for producing yen paths in secs:" << ms << endl;
			YenTime = ms;
		}
	}
	else
	{ //optimize_budget!
		if (debug)
			//cout<<FORERED<<RESETTEXT<<"First get distances for each destination up to input_lambda="<<input_lambda<<RESETTEXT<<endl;
			cout << "First get distances for each destination up to input_lambda=" << input_lambda << endl;
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
						cout << "\tAdded DJKSTR Graph node[" << lks2[i].u << "," << lks2[i].v << "],w=" << lks2[i].weight << endl;
				}
				for (int i = 0; i < array_size3; i++)
				{
					current_graph_dest_to_dest[dest].add_link(lks3[i].u, lks3[i].v, lks3[i].weight);
				}
				if (debug)
					cout << "DJKSTR Graph size:" << array_size2 << ",origin:" << Destinations[dest] << endl;
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
			cout << "Paths from dest:" << Destinations[dest] << " to alternative destinations:" << endl;
			cout << "____________________________" << endl;
			//}

			for (size_t dest2 = 0; dest2 < Destinations.size(); dest2++)
			{
				if (dest2 == dest)
					continue;
				BasePath *result =
					Dijkstra_algs_dest_to_dest[dest].get_shortest_path(
						current_graph_dest_to_dest[dest].get_vertex(Destinations[dest]), current_graph_dest_to_dest[dest].get_vertex(Destinations[dest2]));
				if (debug)
					cout << "\t" << result->get_path_string() << ",Weight:" << result->Weight() << endl;
				else
					cout << "," << result->Weight();
				min_rho = min(min_rho, result->Weight());
				max_rho = max(max_rho, result->Weight());
			}
			cout << endl
				 << "____________________________" << endl;
		}

		Dijkstra_algs_dest_to_dest.clear();
		current_graph_dest_to_dest.clear();
		cout << "Destinations calculated" << endl;

		float upper_disc_dist = (min_rho / 2.0);
		if (debug)
		{
			//cout<<FORERED<<RESETTEXT<<"First finished,"<<FOREGRN<<RESETTEXT<<"rho_min=,"<<min_rho-1<<",UPPER DISCLOSING DISTANCE:"<<(min_rho/2)-1<<RESETTEXT<<endl;
			cout << "FIRST finished,"
				 << "rho_min=," << min_rho << ",UPPER DISCLOSING DISTANCE:" << (min_rho / 2) << endl;
		}
		else
		{
			//cout<<FORERED<<RESETTEXT<<"Minimum geodetic distance between destinations(rho_min)=,"<<min_rho-1<<",UPPER DISCLOSING DISTANCE:"<<upper_disc_dist<<RESETTEXT<<endl;
			cout << "FIRST Minimum geodetic distance between destinations(rho_min)=," << min_rho << ",UPPER DISCLOSING DISTANCE:," << upper_disc_dist << ",max_rho=" << max_rho << endl;
		}
		//cout<<FORERED<<RESETTEXT<<"Second, get P(d) for all Destinations:"<<Destinations<<", for boundary:"<<upper_disc_dist+1<<RESETTEXT<<endl;
		cout << "____________________________" << endl;
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
			cout << "New input_lambda=," << input_lambda << endl;
		}
		cout << "SECOND, get P(d) for all Destinations:" << Destinations << ", for boundary:" << input_lambda << ",time:," << ms << endl;

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
					//cout<<"node "<<node<<" is not reachable from start"<<endl;
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
					//cout<<"reachable node:"<<node<<",cost to destination:"<<Dijkstra_algs[0].get_cost(node)<<endl;
					//Feasible analysis, for now initial cost of all nodes is 0
					//auto budget=40;

					//cout<<"\t feasible for budget=40"<<endl;
				}
			}
			cout << "Reachable nodes:" << reachable_nodes << endl;
			exit(1);
			//cout<<FOREGRN<<RESETTEXT<<"Dest["<<Destinations[dest]<<"],Nd:";
			cout << "Dest[" << Destinations[dest] << "],Nd_size:" << perimeter_nodes_vect[dest].size() << endl;
			/*for(auto it : perimeter_nodes_vect[dest]){
	      cout<<","<<it;
	    }*/
			if (perimeter_nodes_vect[dest].size() == 0)
			{
				cerr << "Aborting, destination:" << Destinations[dest] << " has no perimeter nodes!!!" << endl;
				exit(101);
			}
			perimeter_nodes += perimeter_nodes_vect[dest].size();
			neighourhood_nodes += neighourhood_nodes_vect[dest].size();
			//cout<<FOREGRN<<RESETTEXT<<"\nDest["<<Destinations[dest]<<"],Nd|:";
			cout << "\nDest[" << Destinations[dest] << "],Nd|:"
				 << ",Nd|_size:" << neighourhood_nodes_vect[dest].size() << endl;
			//for(auto it : neighourhood_nodes_vect[dest]){
			//  cout<<","<<it;
			//}
			//cout<<endl<<RESETTEXT;

			//cout<<"hola6"<<flush<<endl;
			//Dijkstra_algs[dest].get_perim_paths(current_perim_paths);
			//perim_paths.push_back(current_perim_paths);
			cout << "____________________________" << endl;
			if (debug)
				Dijkstra_algs[dest].print_paths();
		}
		//cout<<"hola"<<endl;
		/*for(auto it : perim_paths){
	    for (auto it2 : it){
	      cout<<"\t";(*it2.second).PrintOut(cout);
	    }
	  }*/

		if (!SaT_map)
		{ //Already done in create_grid_nodes_and_edges_from_SaT_file
			cout << "getting shortest destination path to all nodes from origin in main graph" << endl;
			main_dijkstra_alg->determine_all_shortest_paths(main_graph.get_vertex(start), INT_MAX);
			cout << "finished getting shortest paths to all nodes from origin in the main graph." << endl;
		}
		//Calculate max_lambda
		int max_lambda = 0;
		vector<int> Optimal_path_dests;

		/*for(auto i=1;i<100;i++){
	    if(FO_nodes.find(i)==FO_nodes.end()){
	      continue;
	    }
	    auto path=main_dijkstra_alg->recover_shortest_perim_path(i);
	    cout<<"Cost["<<i<<"]="<<path->Weight()<<endl;
	    cout<<"Path["<<i<<"]="; path->PrintOut(cout);
	  }
	    
	  auto path20=main_dijkstra_alg->recover_shortest_perim_path(20);
	  cout<<"Cost[20]="<<path20->Weight()<<endl;
	  auto path60=main_dijkstra_alg->recover_shortest_perim_path(60);
	  cout<<"Cost[60]="<<path60->Weight()<<endl;
	  auto path91=main_dijkstra_alg->recover_shortest_perim_path(91);
	  cout<<"Cost[91]="<<path91->Weight()<<endl;
	  auto path96=main_dijkstra_alg->recover_shortest_perim_path(96);
	  cout<<"Cost[96]="<<path96->Weight()<<endl;
	  exit(1);*/

		for (auto it : all_destinations)
		{
			auto path = main_dijkstra_alg->recover_shortest_perim_path(it);
			cout << "Destination:," << it << ",optimal_path_length:," << path->length() << ",Cost:" << path->Weight() << endl;
			path->PrintOut(cout);
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

		cout << "____________Pd________________" << endl;
		for (size_t dest = 0; dest < Destinations.size(); dest++)
		{
			//Dijkstra_algs[dest].print_paths();

			//cout<<FOREGRN<<RESETTEXT<<"P(d)["<<Destinations[dest]<<"]:"<<endl;
			cout << "P(d)[" << Destinations[dest] << "]:" << endl;

			//std::map<int,BasePath*> *current_perim_pahts;
			//cout<<FOREGRN<<RESETTEXT;Dijkstra_algs[dest].make_full_paths(main_dijkstra_alg->get_pt_paths());
			//Dijkstra_algs[dest].make_full_paths(main_dijkstra_alg->get_pt_paths());
			Dijkstra_algs[dest].print_perim_paths();
			//cout<<RESETTEXT;
			//Dijkstra_algs[dest].populate_final_paths(&Final_paths);
			//final_path1.append_to_path(tail1);
		}
		cout << "____________________________" << endl;
		//Time up to Pd, including generating the full list of paths from origin
		end_time = chrono::high_resolution_clock::now();
		diff = end_time - orig_start_time;
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();

		auto temp_start_time = chrono::high_resolution_clock::now();
		//Thirdly, Find crossings & also check if no destination has crossings
		set<int> paired_destinations;
		set<int> unpaired_destinations;
		//cout<<FORERED<<RESETTEXT<<"THIRD, get Nd,d' for min_rho+1"<<RESETTEXT<<endl;
		cout << "THIRD, get Nd,d' for min_rho,time:," << ms << endl;
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
					//cout<<FOREGRN<<RESETTEXT<<"N["<<Destinations[dest1]<<","<<Destinations[dest2]<< "] has the following intersections:";
					cout << "N[" << Destinations[dest1] << "," << Destinations[dest2] << "] has the following intersections:";
					for (auto it : intersect)
						cout << "," << it; //cout<<RESETTEXT<<endl;
					cout << endl;
					paired_destinations.insert(Destinations[dest1]);
					paired_destinations.insert(Destinations[dest2]);
					Ndd[make_pair(Destinations[dest1], Destinations[dest2])] = intersect;
					total_pairs += intersect.size();
				}
				else
				{
					//cout<<FOREGRN<<RESETTEXT<<"N(d,d')["<<Destinations[dest1]<<","<<Destinations[dest2]<< "] has no intersections."<<endl;
					cout << "N(d,d')[" << Destinations[dest1] << "," << Destinations[dest2] << "] has no intersections." << endl;
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
							cout << "getting shortest destination path to all nodes from bifurcation:" << it_intersect << endl;
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
					//cout<<"\tshortest path from origin node:,"<<start<<",to bifurcation node:"<<it_intersect<<",for destination pair:"<<Destinations[dest1]<<","<<Destinations[dest2]<<endl;
					//main_dijkstra_alg->clear();
					shared_ptr<BasePath> head_path =
						main_dijkstra_alg->recover_shortest_perim_path(it_intersect);
					//cout<<"head_path:"<<head_path->get_path_string()<<endl;
					shared_ptr<BasePath> tail1 =
						Dijkstra_algs[dest1].recover_shortest_perim_path(it_intersect);
					//cout<<"tail1:"<<tail1->get_path_string()<<endl;
					shared_ptr<BasePath> tail2 =
						Dijkstra_algs[dest2].recover_shortest_perim_path(it_intersect);
					//cout<<"tail2:"<<tail1->get_path_string()<<endl;

					BasePath path1 = *head_path;
					path1.append_to_path(tail1);

					BasePath path2 = *head_path;
					path2.append_to_path(tail2);
					//cout<<"path1:"<<path1.get_path_string()<<endl;
					//cout<<"path2:"<<path2.get_path_string()<<endl;

					//cout<<"\t";head_path->PrintOut(cout);
					//cout<<"\t";path1.PrintOut(cout);
					//cout<<"\t";path2.PrintOut(cout);
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
		cout<<endl<<"Shortest path from origin node:,"<<start<<",to bifurcation node:"<<best_bifurc_node<<",for destination pair:"<<Destinations[dest1]<<","<<Destinations[dest2]<<"is:";
		best_bifur_path->PrintOut(cout);cout<<endl;
		cout<<FORERED<<RESETTEXT"FOURTH, concatenation of origin-Pd,d'-PairedDestinations, final pair of paths:"<<RESETTEXT<<endl;
		//First paired path
		//shared_ptr<BasePath> tail1 = make_
		//  main_dijkstra_alg->get_shortest_path(
		//      main_graph.get_vertex(best_bifurc_node),main_graph.get_vertex(Destinations[dest1]));
		//BasePath final_path1=*best_bifur_path;
		//final_path1.append_to_path(tail1);
		//cout<<FOREGRN<<RESETTEXT"\t\t First paired path:";final_path1.PrintOut(cout);cout<<RESETTEXT;
		//Final_paths.push_back(final_path1);
		//Second paired path
		//BasePath* tail2 = 
		//  main_dijkstra_alg->get_shortest_path(
		//      main_graph.get_vertex(best_bifurc_node),main_graph.get_vertex(Destinations[dest2]));
		//BasePath final_path2=*best_bifur_path;
		//final_path2.append_to_path(tail2);
		//cout<<FOREGRN<<RESETTEXT<<"\t\t Second paired path:";final_path2.PrintOut(cout);cout<<RESETTEXT;
		//Final_paths.push_back(final_path2);
	      }*/
			}
		}

		cout << "____________________________" << endl;
		for (auto it : Pdd)
		{
			for (auto it2 : it.second)
			{
				if (debug)
				{
					cout << "P(d,d')[" << it.first.first << "," << it.first.second << "]" << endl;
					//cout<<FOREGRN<<RESETTEXT<<"\tP(d,d')["<<it.first.first<<","<<it.first.second<<"],";
					//need to get biggest path to find bifurcation node
					/*if(it2.first->length()>=it2.second->length()){
		  cout<<"bifurcation_node:"<<it2.first->get_path_string()<<flush;cout<<","<<it2.first->get_node(min(it2.first->length(),input_lambda));
		}
		else{
		  cout<<"bifurcation_node:"<<it2.second->get_path_string()<<flush;cout<<","<<it2.first->get_node(min(it2.second->length(),input_lambda));
		}*/
					//cout<<",path1:"<<it2.first->get_path_string()<<",path2:"<<it2.second->get_path_string()<<endl<<RESETTEXT;
					if (debug)
						cout << endl
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
						cout << "added paired path:" << it2.first->get_path_string() << " for destination:," << it.first.first << ",new_best_val:," << max(rel_cost1, rel_cost2) << endl;
						cout << "added paired path:" << it2.second->get_path_string() << " for destination:," << it.first.first << ",new_best_val:," << max(rel_cost1, rel_cost2) << endl;
						cout << "size of paths for destination:," << it.first.first << ",is:," << best_rel_cost_per_dest[it.first.first].second.size();
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
						cout << "added paired path:" << it2.first->get_path_string() << " for destination:," << it.first.second << ",new_best_val:," << max(rel_cost1, rel_cost2) << endl;
						cout << "added paired path:" << it2.second->get_path_string() << " for destination:," << it.first.second << ",new_best_val:," << max(rel_cost1, rel_cost2) << endl;
						cout << "size of paths for destination:" << it.first.second << ",is:," << best_rel_cost_per_dest[it.first.second].second.size();
					}
				}
			}
		}
		cout << "____________________________" << endl;
		//cout<<FORERED<<RESETTEXT<<"FOURTH, get P(d,d',d'')"<<RESETTEXT<<endl;
		end_time = chrono::high_resolution_clock::now();
		diff = end_time - orig_start_time;
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
		cout << "FOURTH, get P(d,d',d''),time:," << ms << endl;
		if (debug)
		{
			cout << "After P(d,d'), current Final paths:" << endl;
			int counter = 0;
			for (auto it : Final_paths)
			{
				cout << "P_lambda[" << counter++ << "]->" << it->get_path_string() << endl;
			}
		}
		//timings
		auto temp_end_time = chrono::high_resolution_clock::now();
		auto diff = temp_end_time - temp_start_time;
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
		cout << "Time(ms) for PDD:," << ms << endl;

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
			cout << "Working on dest:" << it_all_dests << endl; //cout<<",neighourhood_nodes:"<<neighourhood_nodes_vect<<endl;
			shared_ptr<BasePath> opt_dest_path =
				main_dijkstra_alg->recover_shortest_perim_path(it_all_dests);
			if (opt_dest_path->Weight() <= input_lambda)
			{
				cout << "skipping destination:" << it_all_dests << ",optimal_path_weight:" << opt_dest_path->Weight() << ",lambda:" << input_lambda << "from d in tripple d,d',d'' because it does not need cover" << endl;
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
					cout << "skipping perimeter node " << it << " because is is not reachable from origin." << endl;
					continue;
				}

				int h1 = it;
				if (h1 == start)
					continue;
				//NEED TO ENSURE NEITHER H1 OR H2 ARE DESTINATIONS
				//CAN HAPPEN WHEN WE ALLOW CYCLES
				if (all_destinations.find(h1) != all_destinations.end())
				{
					cout << "\tskipping destination_node:" << h1 << endl;
					continue;
				}

				if (debug)
					cout << "\th1:" << h1 << endl;
				shared_ptr<BasePath> head_path =
					main_dijkstra_alg->recover_shortest_perim_path(it);
				if (debug)
					cout << "\thead_path:" << head_path->get_path_string() << endl;
				BasePath FinalPath0 = *head_path;
				FinalPath0.append_to_path(Dijkstra_algs[dest_index[it_all_dests]].recover_shortest_perim_path(h1));
				auto path_ptr0 = make_shared<BasePath>(FinalPath0);
				if (debug)
					cout << "\t\t\t\tFinalPath1" << FinalPath0.get_path_string() << endl;

				for (auto it2 : Ndd)
				{
					if (it2.first.first == it_all_dests || it2.first.second == it_all_dests)
						continue; //d neq d' nor d''
					if (debug)
						cout << "\t\t connecting to Ndd(" << it2.first.first << "," << it2.first.second << ")" << endl;
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
							cout << "\t\t\t"
								 << ",h1:," << h1 << ",h2:," << h2 << ",h1->h2 path:";
						if (!Dijkstra_bifurcations_map[h2].is_node_reachable(h1))
						{
							//cout<<"\tskipping h2node "<<h2<<" because is is not reachable from h1"<<h1<<endl;
							continue;
						}
						shared_ptr<BasePath> h1_h2_path;
						if (h1 != h2)
						{
							//if(h1!=start){
							h1_h2_path = Dijkstra_bifurcations_map[h2].recover_shortest_perim_path(h1);
							if (debug)
								cout << h1_h2_path->get_path_string() << endl;
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
							cout << "\t\t\t\tFinalPath2" << FinalPath1.get_path_string() << endl;
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
							cout << "\t\t\t\tFinalPath3" << FinalPath2.get_path_string() << endl;
						pair<int, int> dest_pair(it2.first.first, it2.first.second);

						auto path_ptr1 = make_shared<BasePath>(FinalPath1);
						auto path_ptr2 = make_shared<BasePath>(FinalPath2);

						tuple<shared_ptr<BasePath>, shared_ptr<BasePath>, shared_ptr<BasePath>> Final_path_triple = make_tuple(path_ptr0, path_ptr1, path_ptr2);
						//cout<<"\t FinalPath0 inserted:"<<path_ptr0->get_path_string()<<endl;
						//cout<<"\t FinalPath1 inserted:"<<path_ptr1->get_path_string()<<endl;
						//cout<<"\t FinalPath2 inserted:"<<path_ptr2->get_path_string()<<endl;

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
								cout << "added tripple paths:" << path_ptr0->get_path_string() << " for single destination:," << it_all_dests << ",new_best_val:," << max_cost << endl;
								cout << "added tripple paths:" << path_ptr1->get_path_string() << " for single destination:," << it_all_dests << ",new_best_val:," << max_cost << endl;
								cout << "added tripple paths:" << path_ptr2->get_path_string() << " for single destination:," << it_all_dests << ",new_best_val:," << max_cost << endl;
								cout << "size of paths for destination:" << it_all_dests << ",is:," << best_rel_cost_per_dest[it_all_dests].second.size();
							}
						}

						/*if(max_cost<best_rel_cost_per_dest[dest_pair.first].first){
		      best_rel_cost_per_dest[dest_pair.first].first=max_cost;//update best value
		      best_rel_cost_per_dest[dest_pair.first].second.clear();
		      best_rel_cost_per_dest[dest_pair.first].second.push_back(path_ptr0);
		      best_rel_cost_per_dest[dest_pair.first].second.push_back(path_ptr1);
		      best_rel_cost_per_dest[dest_pair.first].second.push_back(path_ptr2);
		      //if(debug){
			//cout<<"added tripple paths:"<<path_ptr0->get_path_string()<<" for paired destination:,"<<dest_pair.first<<",new_best_val:,"<<max_cost<<endl;
			//cout<<"added tripple paths:"<<path_ptr1->get_path_string()<<" for paired destination:,"<<dest_pair.first<<",new_best_val:,"<<max_cost<<endl;
			//cout<<"added tripple paths:"<<path_ptr1->get_path_string()<<" for paired destination:,"<<dest_pair.first<<",new_best_val:,"<<max_cost<<endl;
		      //}
		    }
		    
		    if(max_cost<best_rel_cost_per_dest[dest_pair.second].first){
		      best_rel_cost_per_dest[dest_pair.second].first=max_cost;//update best value
		      best_rel_cost_per_dest[dest_pair.second].second.clear();
		      best_rel_cost_per_dest[dest_pair.second].second.push_back(path_ptr0);
		      best_rel_cost_per_dest[dest_pair.second].second.push_back(path_ptr1);
		      best_rel_cost_per_dest[dest_pair.second].second.push_back(path_ptr2);
		      //if(debug){
			//cout<<"added tripple paths:"<<path_ptr0->get_path_string()<<" for paired destination:,"<<dest_pair.second<<",new_best_val:,"<<max_cost<<endl;
			//cout<<"added tripple paths:"<<path_ptr1->get_path_string()<<" for paired destination:,"<<dest_pair.second<<",new_best_val:,"<<max_cost<<endl;
			//cout<<"added tripple paths:"<<path_ptr1->get_path_string()<<" for paired destination:,"<<dest_pair.second<<",new_best_val:,"<<max_cost<<endl;
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
		cout << "Time(ms) for PDDD:," << ms << ",triples:" << total_triples << endl;

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
			cout << ",rel_cost:," << rel_cost << ",P_lambda[0," << counter++ << "]->" << it->get_path_string() << endl;
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
	    cout<<"final_path_pos["<<counter++<<"]";
	    for(size_t i=0;i<it->length();i++){
	      int id=it->get_node(i);
	      auto pos=node_map[id];
	      cout<<"["<<pos.first<<","<<pos.second<<"],";
	      if(i>1){
		if(abs(pos.first-old_pos.first)+abs(pos.second-old_pos.second)>1){
		  cout<<"pos:"<<pos.first<<","<<pos.second<<",old_pos:"<<old_pos.first<<","<<old_pos.second<<endl;
		  exit(1);
		}
	      }
	      old_pos=pos;
	    }
	    cout<<endl;
	  }*/
		//timings
		end_time = chrono::high_resolution_clock::now();
		diff = end_time - start_time;
		start_time = end_time; //restart
		ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();

		cout << "N:," << N << ",R:," << random_positions << ",min_dist_to_dest:," << min_dist_to_dest << ",max_dist_to_dest:," << max_dist_to_dest << ",max_rel_cost:," << max_rel_cost << ",";
		cout << "final_paths:," << Final_paths.size() << ",overall_time:," << ms << ",ms_bifur_backwards:," << ms_bifur_backward << ",ms_pairs:," << ms_pairs << ",ms_triples:," << ms_triples;
		cout << ",boundary_nodes:," << perimeter_nodes << ",bifurcations:," << bifurcated_destinations << ",total_perim_nodes:," << perimeter_nodes << ",total_neighourhood_nodes:," << neighourhood_nodes << ",total_pairs:," << total_pairs << ",total_triples:" << total_triples;
		cout << ",labmda_star:," << upper_disc_dist << ",input_lambda:," << input_lambda << ",max_lambda:," << max_lambda;
		cout << ",Destination OP lengths:," << Optimal_path_dests << ",connectivity:," << connectivity << ",start:," << start;
		cout << ",Total_nodes:," << nodes_set.size() << ",reachable_nodes:," << main_dijkstra_alg->get_perim_size() << endl;

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
					//cout<<"@@";
					cout << "@";
				}
				else
				{
					//auto current_node=coord_map2[current_pos];
					if (color_map[i][j][0] == 1)
					{ //origin node;
						//cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
						cout << FOREGRN << "\u25A0" << RESETTEXT;
						//cout<<FOREGRN<<"."<<RESETTEXT;
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
									//cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FOREGRN << "\u25A0" << RESETTEXT;
									continue;
								}
								else if (j > 0 && color_map[i + 1][j - 1][0] == 1)
								{ //origin node;
									//cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FOREGRN << "\u25A0" << RESETTEXT;
									continue;
								}
								else if (color_map[i + 1][j][0] == 1)
								{ //origin node;
									//cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FOREGRN << "\u25A0" << RESETTEXT;
									continue;
								}
							}
							if (i > 0)
							{
								if (color_map[i - 1][j][0] == 1)
								{ //origin node;
									//cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FOREGRN << "\u25A0" << RESETTEXT;
									continue;
								}
								else if (j > 0 && color_map[i - 1][j - 1][0] == 1)
								{ //origin node;
									//cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FOREGRN << "\u25A0" << RESETTEXT;
									continue;
								}
								else if (j < max_y - 1 && color_map[i - 1][j + 1][0] == 1)
								{ //origin node;
									//cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FOREGRN << "\u25A0" << RESETTEXT;
									continue;
								}
							}
						}

						//DESTINATIONS
						if (color_map[i][j][1] == 1)
						{ //destination node
							//cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
							cout << FORERED << "\u25A0" << RESETTEXT;
							//cout<<FORERED<<"."<<RESETTEXT;
							continue;
						}
						if (size > 100000)
						{
							if (j > 0)
							{
								if (color_map[i][j - 1][0] == 1)
								{ //destinaiton node
									//cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FOREGRN << "\u25A0" << RESETTEXT;
									continue;
								}
							}
							if (j < max_y - 1)
							{
								if (color_map[i][j + 1][0] == 1)
								{ //destination node;
									//cout<<FOREGRN<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FOREGRN << "\u25A0" << RESETTEXT;
									continue;
								}
							}
							if (i < max_x - 1)
							{
								if (j < max_y - 1 && color_map[i + 1][j + 1][1] == 1)
								{ //destination node;
									//cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FORERED << "\u25A0" << RESETTEXT;
									continue;
								}
								else if (j > 0 && color_map[i + 1][j - 1][1] == 1)
								{ //destination node;
									//cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FORERED << "\u25A0" << RESETTEXT;
									continue;
								}
								else if (color_map[i + 1][j][1] == 1)
								{ //destination node;
									//cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FORERED << "\u25A0" << RESETTEXT;
									continue;
								}
							}
							if (i > 0)
							{
								if (color_map[i - 1][j][1] == 1)
								{ //destination node;
									//cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FORERED << "\u25A0" << RESETTEXT;
									continue;
								}
								else if (j > 0 && color_map[i - 1][j - 1][1] == 1)
								{ //destination node;
									//cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FORERED << "\u25A0" << RESETTEXT;
									continue;
								}
								else if (j < max_y - 1 && color_map[i - 1][j + 1][1] == 1)
								{ //destination node;
									//cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FORERED << "\u25A0" << RESETTEXT;
									continue;
								}
							}
							if (j > 0)
							{
								if (color_map[i][j - 1][1] == 1)
								{ //destinaiton node
									//cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FORERED << "\u25A0" << RESETTEXT;
									continue;
								}
							}
							if (j < max_y - 1)
							{
								if (color_map[i][j + 1][1] == 1)
								{ //destination node;
									//cout<<FORERED<<"\u25A0\u25A0"<<RESETTEXT;
									cout << FORERED << "\u25A0" << RESETTEXT;
									continue;
								}
							}
						}
						//FINAL PATHS
						if (color_map[i][j][2] == 1)
						{ //path node
							//cout<<FOREBLU<<"\u25A0\u25A0"<<RESETTEXT;
							//cout<<FOREBLU<<"\u25A0"<<RESETTEXT;
							cout << FOREBLU << "*" << RESETTEXT;
							continue;
						}
						//cout<<"..";
						cout << ".";
					}
				}
			}
			cout << endl;
		}
		exit(0);
	}
	//Need to make copy of paths before erasing any of them
	PathMatrix PM_orig = *PM;
	int t_max = optimizing_target_goal_dist(PM);
	optimizing_observer(&PM_orig, nodes, t_max);
}
