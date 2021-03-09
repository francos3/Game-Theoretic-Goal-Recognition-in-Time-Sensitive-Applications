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
            os << ", "; 
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
void calculate_min_path_strategy_AG();
void simulation();

#endif
