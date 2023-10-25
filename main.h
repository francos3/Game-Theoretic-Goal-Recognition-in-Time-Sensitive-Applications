#ifndef MAIN
#define MAIN
#include<vector>
#include <iostream>
#include <chrono>
#include <ctime>
//#include <boost/timer/timer.hpp>
extern std::vector<int> Destinations;


struct vect_comp_int{
  bool operator() (std::vector<int> v1,std::vector<int> v2) const{
    return v1<v2;
    /*//shorter is smaller for our ordering
    if(v1.size()!=v2.size())
      return v1<v2;
    //if a vector is empty use lexicografical comparison
    if(v1.size()==0||v2.size()==0)
      return v1<v2;
    
    //WE COMPARE VECTORS OF SAME SIZE FROM END TO BEGINNING
    //MUCH BETTER WHEN COMPARING PATHS WITH THE SAME ORIGIN POINT!
    for(int i=v1.size()-1;i>=0;i--){
      if(v1[i]!=v2[i])
  return v1[i]<v2[i];
    }
    return true;*/
  }
};
struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        std::hash<int> hasher;
        size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

template <typename T> 
std::ostream& operator<<(std::ostream& os, const std::set<T>& v) 
{ 
    os << "["; 
    for (auto it : v) { 
        os << it; 
        if (it != *v.rbegin()) 
            os << ","; 
    } 
    os << "]\n"; 
    return os; 
} 

struct hash_pair { 
    template <class T1, class T2> 
    size_t operator()(const std::pair<T1, T2>& p) const
    { 
        auto hash1 = std::hash<T1>{}(p.first); 
        auto hash2 = std::hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
};

void elapsed_time(std::string method, std::chrono::_V2::system_clock::time_point start_time);
void calculate_min_path_strategy_AG();
void calculate_Gibbs_strategy_AG(double beta,double Budget);
void simulation();
void simulation_observer();
void calculate_observer_zeta();
void calculate_observer_lambda(int current_node,float current_cost);
float calculate_zeta(std::pair<int,float> s);
int main2(int argc, char *argv[]);
float calculate_q(std::pair<int,float> node,int dest);
float calculate_q_avg(std::pair<int,float> node,int dest);
float calculate_term1(std::pair<int,float> s);
float calculate_term2(std::pair<int,float> s);
float round2(float n);
void from_file_random_dest_and_origin_selection();
void calculate_all_subpaths();
float calculate_q_path(std::pair<int,float> node,int dest,float cost);
float calculate_q_recursive(std::pair<int,float> node,int dest);//Eq 25
void random_orig_and_dest_placement_SaT();
void calculate_min_path_strategy_prefixes(bool skip);
float calculate_q_recursive_prefix(unsigned node);
void simulation_observer_prefix();
std::pair<unsigned, float> dest_predictor_prefix(unsigned node);
float best_avg_target_choice();
void fictitious_play();
void fictitious_play2();
void store_prev_strategies();
void calculate_next_target_strategy();
//std::pair<unsigned, float> prev_dest_predictor_prefix(unsigned node, unsigned iter);
std::pair<float,float> calculate_percentage_difference(std::vector<float> values, int starting_iteration);
std::pair<float,float> calculate_statistics(const std::vector<float> &values);
void store_cutset(std::set<unsigned>& cutset_prefix);
void removeDominatedFromCutset();
void exploreDescendants(int node, std::shared_ptr<std::vector<std::set<unsigned int>>> edges,
                        std::set<int> & nodesToRemove);
void print_cutset_endings();
#endif
