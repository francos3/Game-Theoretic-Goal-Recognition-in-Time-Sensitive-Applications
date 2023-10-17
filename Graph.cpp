///////////////////////////////////////////////////////////////////////////////
///  Graph.cpp
///  <TODO: insert file description here>
///
///  @remarks <TODO: insert remarks here>
///
///  @author Yan Qi @date 8/18/2010
/// 
///  $Id: Graph.cpp 65 2010-09-08 06:48:36Z yan.qi.asu $
///////////////////////////////////////////////////////////////////////////////

#include <limits>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "GraphElements.h"
#include "Graph.h"


const double Graph::DISCONNECT = (numeric_limits<double>::max)();


Graph::Graph( const string& file_name )
{
	_import_from_file(file_name);
}

Graph::Graph( const Graph& graph )
{
	m_nVertexNum = graph.m_nVertexNum;
	m_nEdgeNum = graph.m_nEdgeNum;
	m_vtVertices.assign(graph.m_vtVertices.begin(),graph.m_vtVertices.end());
	m_mpFaninVertices.insert(graph.m_mpFaninVertices.begin(),graph.m_mpFaninVertices.end());
	m_mpFanoutVertices.insert(graph.m_mpFanoutVertices.begin(),graph.m_mpFanoutVertices.end());
	m_mpEdgeCodeWeight.insert(graph.m_mpEdgeCodeWeight.begin(),graph.m_mpEdgeCodeWeight.end());
	m_mpVertexIndex.insert(graph.m_mpVertexIndex.begin(),graph.m_mpVertexIndex.end());
}

Graph::~Graph(void)
{
	clear();
}

///////////////////////////////////////////////////////////////////////////////
///  public  _import_from_file
///  Construct the graph by importing the edges from the input file. 
///
///  @param [in]       file_name const std::string &    The input graph file
///
///  This function doesn't return a value
///
///  @remarks The format of the file is as follows:
///   1. The first line has an integer as the number of vertices of the graph
///   2. Each line afterwards contains a directed edge in the graph:
///		     starting point, ending point and the weight of the edge. 
///		 These values are separated by 'white space'.
///
///  @see <TODO: insert text here>
///
///  @author Yan Qi @date 5/29/2010
///////////////////////////////////////////////////////////////////////////////
void Graph::_import_from_file( const string& input_file_name )
{
	const char* file_name = input_file_name.c_str();

	//1. Check the validity of the file
	ifstream ifs(file_name);
	if (!ifs)
	{
		cerr << "The file " << file_name << " can not be opened!" << endl;
		exit(1);
	}

	//2. Reset the members of the class
	clear();

	//3. Start to read information from the input file. 
	/// Note the format of the data in the graph file.
	//3.1 The first line has an integer as the number of vertices of the graph
	ifs >> m_nVertexNum;

	//3.2 In the following lines, each line contains a directed edge in the graph:
	///   the id of starting point, the id of ending point, the weight of the edge. 
	///   These values are separated by 'white space'. 
	int start_vertex, end_vertex;
	double edge_weight;
	

	while(ifs >> start_vertex)
	{
		if (start_vertex == -1)
		{
			break;
		}
		ifs >> end_vertex;
		ifs >> edge_weight;

		///3.2.1 construct the vertices
		BaseVertex* start_vertex_pt = get_vertex(start_vertex);
		BaseVertex* end_vertex_pt = get_vertex(end_vertex);

		///3.2.2 add the edge weight
		//// note that the duplicate edge would overwrite the one occurring before. 
		m_mpEdgeCodeWeight[get_edge_code(start_vertex_pt, end_vertex_pt)] = edge_weight;

		///3.2.3 update the fan-in or fan-out variables
		//// Fan-in
		get_vertex_set_pt(end_vertex_pt, m_mpFaninVertices)->insert(start_vertex_pt);

		//// Fan-out
		get_vertex_set_pt(start_vertex_pt, m_mpFanoutVertices)->insert(end_vertex_pt);	
	}	

	m_nVertexNum = m_vtVertices.size();
	m_nEdgeNum = m_mpEdgeCodeWeight.size();

	ifs.close();	
}

void Graph::setv()
{
	m_nVertexNum = m_vtVertices.size();
	m_nEdgeNum = m_mpEdgeCodeWeight.size();
}
int Graph::getVertexNum(){
  return m_vtVertices.size();
}


void Graph::set_number_vertices( const int& nVertices )
{
	m_nVertexNum = nVertices;
}

void Graph::add_link( const int& s, const int& d, const double& weight )
{
	///3.2.1 construct the vertices
	BaseVertex* start_vertex_pt = get_vertex(s);
	BaseVertex* end_vertex_pt = get_vertex(d);

	///3.2.2 add the edge weight
	//// note that the duplicate edge would overwrite the one occurring before. 
	m_mpEdgeCodeWeight[get_edge_code(start_vertex_pt, end_vertex_pt)] = weight;

	///3.2.3 update the fan-in or fan-out variables
	//// Fan-in
	get_vertex_set_pt(end_vertex_pt, m_mpFaninVertices)->insert(start_vertex_pt);

	//// Fan-out
	get_vertex_set_pt(start_vertex_pt, m_mpFanoutVertices)->insert(end_vertex_pt);	
}
void Graph::dump_edges(){
  for (auto edge : m_mpEdgeCodeWeight){
    cout<<"edge_code:,"<<edge.first<<",weight:,"<<edge.second<<endl;
  }
}

BaseVertex* Graph::get_vertex( int node_id )
{
	if (m_stRemovedVertexIds.find(node_id) != m_stRemovedVertexIds.end())
	{
		return NULL;
	}
	else
	{
		BaseVertex* vertex_pt = NULL;
		const map<int, BaseVertex*>::iterator pos = m_mpVertexIndex.find(node_id);
		
		if (pos == m_mpVertexIndex.end())
		{
			int vertex_id = m_vtVertices.size();
			vertex_pt = new BaseVertex();
			vertex_pt->setID(node_id);
			m_mpVertexIndex[node_id] = vertex_pt;

			m_vtVertices.push_back(vertex_pt);
		}
		else
		{
			vertex_pt = pos->second;
		}

		return vertex_pt;	
	}
}

void Graph::clear()
{
	m_nEdgeNum = 0;
	m_nVertexNum = 0;

	for( map<BaseVertex*, set<BaseVertex*>*>::const_iterator pos= m_mpFaninVertices.begin();
		 pos!=m_mpFaninVertices.end(); 
		 ++pos)
	{
		delete pos->second;
	}

	m_mpFaninVertices.clear();

	for(map<BaseVertex*, set<BaseVertex*>*>::const_iterator pos=m_mpFanoutVertices.begin();
		pos!=m_mpFanoutVertices.end(); ++pos)
	{
		delete pos->second;
	}
	m_mpFanoutVertices.clear();


	m_mpEdgeCodeWeight.clear();

	//clear the list of vertices objects
	for_each(m_vtVertices.begin(), m_vtVertices.end(), DeleteFunc<BaseVertex>());
	m_vtVertices.clear();
	m_mpVertexIndex.clear();

	m_stRemovedVertexIds.clear();
	m_stRemovedEdge.clear();
}

long Graph::get_edge_code( const BaseVertex* start_vertex_pt, const BaseVertex* end_vertex_pt ) const
{
	/// Note that the computation below works only if 
	/// the result is smaller than the maximum of an integer!
	//cout<<"calculating edge,start_vertex_id:,"<<start_vertex_pt->getID()<<",end_vertex_id:,"<<end_vertex_pt->getID()<<",code:,"<< long(start_vertex_pt->getID())*long(m_nVertexNum)+long(end_vertex_pt->getID())<<endl;
	return long(start_vertex_pt->getID())*long(m_nVertexNum)+long(end_vertex_pt->getID());
}


set<BaseVertex*>* Graph::get_vertex_set_pt( BaseVertex* vertex_, map<BaseVertex*, set<BaseVertex*>*>& vertex_container_index )
{
	BaseVertexPt2SetMapIterator pos = vertex_container_index.find(vertex_);

	if(pos == vertex_container_index.end())
	{
		set<BaseVertex*>* vertex_set = new set<BaseVertex*>();
		pair<BaseVertexPt2SetMapIterator,bool> ins_pos = 
			vertex_container_index.insert(make_pair(vertex_, vertex_set));

		pos = ins_pos.first;
	}

	return pos->second;
}


double Graph::get_edge_weight( const BaseVertex* source, const BaseVertex* sink )
{
	int source_id = source->getID();
	int sink_id = sink->getID();
	//cout<<"\t\t\t\t original edge weight:"<<get_original_edge_weight(source, sink)<<endl;

	if (m_stRemovedVertexIds.find(source_id) != m_stRemovedVertexIds.end()
		|| m_stRemovedVertexIds.find(sink_id) != m_stRemovedVertexIds.end()
		|| m_stRemovedEdge.find(make_pair(source_id, sink_id)) != m_stRemovedEdge.end())
	{
		return DISCONNECT;
	}else
	{
		return get_original_edge_weight(source, sink);
	}
}


void Graph::get_adjacent_vertices( BaseVertex* vertex, set<BaseVertex*>& vertex_set )
{
	int starting_vt_id = vertex->getID();

	if (m_stRemovedVertexIds.find(starting_vt_id) == m_stRemovedVertexIds.end())
	{
		set<BaseVertex*>* vertex_pt_set = get_vertex_set_pt(vertex, m_mpFanoutVertices);
		for(set<BaseVertex*>::const_iterator pos=(*vertex_pt_set).begin();
			pos != (*vertex_pt_set).end(); ++pos)
		{
			int ending_vt_id = (*pos)->getID();
			if (m_stRemovedVertexIds.find(ending_vt_id) != m_stRemovedVertexIds.end()
				|| m_stRemovedEdge.find(make_pair(starting_vt_id, ending_vt_id)) != m_stRemovedEdge.end())
			{
				continue;
			}
			//
			vertex_set.insert(*pos);
		}
	}
}

void Graph::get_precedent_vertices( BaseVertex* vertex, set<BaseVertex*>& vertex_set )
{
	if (m_stRemovedVertexIds.find(vertex->getID()) == m_stRemovedVertexIds.end())
	{
		int ending_vt_id = vertex->getID();
		set<BaseVertex*>* pre_vertex_set = get_vertex_set_pt(vertex, m_mpFaninVertices);
		for(set<BaseVertex*>::const_iterator pos=(*pre_vertex_set).begin(); 
			pos != (*pre_vertex_set).end(); ++pos)
		{
			int starting_vt_id = (*pos)->getID();
			if (m_stRemovedVertexIds.find(starting_vt_id) != m_stRemovedVertexIds.end()
				|| m_stRemovedEdge.find(make_pair(starting_vt_id, ending_vt_id)) != m_stRemovedEdge.end())
			{
				continue;
			}
			//
			vertex_set.insert(*pos);
		}
	}
}

double Graph::get_original_edge_weight( const BaseVertex* source, const BaseVertex* sink )
{
  //cout<<"\t\t\t\t\t source:"<<source->getID()<<",sink:"<<sink->getID()<<",code:"<<get_edge_code(source,sink)<<endl;
	map<long, double>::const_iterator pos = 
		m_mpEdgeCodeWeight.find(get_edge_code(source, sink));

	if (pos != m_mpEdgeCodeWeight.end())
	{
		return pos->second;
	}else
	{
		return DISCONNECT;
	}
}

void PathMatrix::set_W(set<int>  *_W){
  W=*_W;
  cout<<"W:";
  for(auto w : W) cout<<w<<",";
  cout<<endl;
}
void PathMatrix::set_C(int _C){
  C=_C;
}

void PathMatrix::add_paths_set(shared_ptr<set<vector<int>,vect_comp_int > > destination_paths,int dest){
      Path_Matrix.push_back(*destination_paths);
      destinations.push_back(dest);
}
void PathMatrix::eliminate_paths_C(){//Max length of paths
  if(C<0){//no limit given
    return;
  }
  vector<set<vector<int>,vect_comp_int> > New_Path_Matrix;
  
  for(int i=0;i<Path_Matrix.size();i++){
    set<vector<int>,vect_comp_int > surviving_paths;
    for (auto next_path : Path_Matrix[i]){
      if(next_path.size()>C){
	if(debug){
	  cout<<"\t\tEliminating path:";
	  for(auto k : next_path) 
	    cout<<k<<",";
	  cout<<"length:"<<next_path.size()<<",C_limit:"<<C<<endl;
	}
      }
      else{
	surviving_paths.insert(next_path);
      }
    }
    New_Path_Matrix.push_back(surviving_paths);
  }
  Path_Matrix=New_Path_Matrix;
}


void PathMatrix::print_paths(){
  for(int i=0;i<Path_Matrix.size();i++){
    cout<<"Dest:"<<destinations[i]<<endl;
    int j=0;
    for (auto next_path : Path_Matrix[i]){
      cout<<"\tnext_path["<<j++<<"]:";
      for(auto k : next_path) 
	cout<<k<<",";
      cout<<"length:"<<next_path.size()<<endl;
    }
  }
}
void PathMatrix::record_initial_paths_gperf(){
  std::ofstream outfile ("InitialPaths.gperf");
  for(int i=0;i<Path_Matrix.size();i++){
    for (auto next_path : Path_Matrix[i]){
      for(int i=0;i<next_path.size()-1;i++){
	outfile<<next_path[i]<<"_";
      }
      outfile<<next_path.back()<<endl;
    }
  }
  outfile.close();
}
void PathMatrix::record_initial_paths(){
  std::ofstream outfile ("InitialPaths.txt");
  InitialPaths=0;
  for(int i=0;i<Path_Matrix.size();i++){
    outfile<<"Dest:"<<destinations[i]<<endl;
    int j=0;
    InitialPaths+=Path_Matrix[i].size();
    for (auto next_path : Path_Matrix[i]){
      outfile<<"\tnext_path["<<j++<<"]:";
      for(auto k : next_path) 
	outfile<<k<<",";
      outfile<<"length:"<<next_path.size()<<endl;
    }
  }
  outfile.close();
}
void PathMatrix::record_best_solution(){
  //By destinations
  std::ofstream outfile ("BestPathsSet_Directions.txt");
  for(int i=0;i<Path_Matrix.size();i++){
    outfile<<"Dest:"<<destinations[i]<<endl;
    int j=0;
    for (auto next_path : Path_Matrix[i]){
      outfile<<"\tnext_path["<<j++<<"]:";
      for(auto k : next_path) 
	outfile<<k<<",";
      outfile<<"length:"<<next_path.size()<<endl;
    }
  }

  outfile.close();
  //By t values
  std::ofstream outfile2 ("BestPathsSet_T.txt");
  for(int k_index=0;k_index<path_T.size();k_index++){
    if(path_T[k_index].size()==0)
      continue;
    for(int path_index=0;path_index<path_T[k_index].size();path_index++){
      if(path_T[k_index][path_index].size()==0)
	continue;
      outfile2<<"P_"<<iteration<<"_[,"<<path_T[k_index][path_index]<<",";
      outfile2<<"],t_index="<<k_index<<endl;
    }
  }
  outfile2.close();
  //By lambda values
  std::ofstream outfile3 ("BestPathsSet_L.txt");
  final_paths=0;
  for(int k_index=0;k_index<path_L.size();k_index++){
    if(path_L[k_index].size()==0)
      continue;
    final_paths+=path_L[k_index].size();
    for(int path_index=0;path_index<path_L[k_index].size();path_index++){
      if(path_L[k_index][path_index].size()==0)
	continue;
      outfile3<<"P_"<<iteration<<"[,"<<path_L[k_index][path_index]<<",";
      outfile3<<"],Lambda_index="<<k_index<<endl;
    }
  }
  outfile3.close();
}

void PathMatrix::get_safe_nodes(int k){//linear scan to determine which nodes maintain indisclosability at kth index
  if(debug)
    cout<<"Calculating which nodes at index "<<k<<" will mantain undisclosability for all destinations"<<endl;
  safe_nodes.clear();
  set<int> temp_list;

  for(int i=0;i<Path_Matrix.size();i++){
    for (auto next_path : Path_Matrix[i]){
	if(k>next_path.size()){//if k index bigger than path length, skip path
	  continue;
	}
	if(i==0){//first destination, we add all numbers in this path
	  safe_nodes.insert(next_path[k]);
	}
	else{//first check if number exist in safe nodes
	  if(safe_nodes.find(next_path[k])==safe_nodes.end()){//skip number
	    continue;
	  }
	  temp_list.insert(next_path[k]);
	}
    }
    if(i>0){
      safe_nodes=temp_list;//so if any node in safe nodes was not in the current destination, it gets eliminated
    }
  }
  
  //Any nodes in ignore list by the observer will be added to the safe list 
  for (auto node : W){
    if(debug)
      cout<<"\t node:"<<node<<" is safe as it is not in W(nodes with sensors)"<<endl;
    safe_nodes.insert(node);
  }
  
  if(debug){
    cout<<"Safe_nodes for index "<<k<<":";
    for (auto node : safe_nodes)
	cout<<node<<",";
      cout<<endl;
  }
}
void PathMatrix::get_safe_nodes2(int k,int min_dests){//linear scan to determine which nodes maintain indisclosability at kth index
  min_dests=1;
  if(debug)
    cout<<"Calculating which nodes at index "<<k<<" will mantain undisclosability for more than "<< min_dests<<" destinations"<<endl;
  safe_nodes2.clear();
  safe_nodes.clear();
  set<int> temp_list;

  //First get destination count per pat
  for(int i=0;i<Path_Matrix.size();i++){
    set<int> already_found;
    for (auto next_path : Path_Matrix[i]){
      if(k>next_path.size()){//if k index bigger than path length, skip path
	continue;
      }
      auto it = safe_nodes2.find(next_path[k]);
      if(i==0){
	safe_nodes2[next_path[k]]=1;
      }
      else if(it!=safe_nodes2.end()){//add one if found first time
	if(already_found.find(next_path[k])==already_found.end()){//first_time
	  it->second++;
	  already_found.insert(next_path[k]);
	  //cout<<"\t\t\tfound multi-d "<<next_path[k]<<" for first time in d#"<<i<<endl;
	}
	else{
	  //cout<<"\t\t\tfound multi-d"<<next_path[k]<<" for second+ time in d#"<<i<<endl;
	  continue;//skip
	}
      }
      else{
	//cout<<"\t\t\tfound first-d "<<next_path[k]<<" for first time in d#"<<i<<endl;
	safe_nodes2[next_path[k]]=1;
	already_found.insert(next_path[k]);
      }
    }
  }
  //Now remove from safe nodes those nodes whose destination count is lower than min_dests
  map<int,int> temp_map;
  for (auto it : safe_nodes2){
    if(it.second>min_dests){//we keep it
      safe_nodes.insert(it.first);
      if(debug) cout<<"\t\tadding "<<it.first<<" from safe_nodes,reaches dests#:"<<it.second<<",threshold:"<<min_dests<<endl;
    }
    else{
      if(debug) cout<<"\t\tremoving "<<it.first<<" from safe_nodes,reaches dests#:"<<it.second<<",threshold:"<<min_dests<<endl;
    }
  }
  
  //Any nodes in ignore list by the observer will be added to the safe list 
  for (auto node : W){
    if(debug)
      cout<<"\t node:"<<node<<" is safe as it is not in W(nodes with sensors)"<<endl;
    safe_nodes.insert(node);
  }
  
  if(debug){
    cout<<"Potentially safe_nodes for index "<<k<<":";
    for (auto node : safe_nodes)
	cout<<node<<",";
      cout<<endl;
  }
}
//lambda is a meassure of how fat is the target from its actual destination
//at the step on which the final destination is revealed.
//CURRENTLY ASSUMING WEIGHTS ARE ONE, IN THE FUTURE NEED TO ITERATE THROUGH PATH
//TO ADD WEIGHTS
int PathMatrix::calculate_lambda(int t, vector<int> path){
  if(t>path.size()){
    cout<<"t_index cannot be bigger than the actual path_size!"<<endl;
    exit(1);
  }
  int result=path.size()-t-2;
  if(debug)
    cout<<"path:,"<<path<<",lambda:,"<<result<<",t:,"<<t+1<<",length:,"<<path.size()<<endl;
  return result;
}
bool PathMatrix::erase_disclosing_paths(int k,int &disclosing_paths){
  disclosing_paths=0;
  vector<set<vector<int>,vect_comp_int > > New_Path_Matrix;
  bool finished=false;

  if(safe_nodes.size()==0){//we are finished,cannot remove any more paths
    cout<<"safe_nodes="<<safe_nodes.size()<<",finished:"<<finished<<endl;
    finished=true;
    //Add all remaining paths at current k level
    if(k>path_T.size()){
      path_T.resize(k+1);
    }
    for(int i=0;i<Path_Matrix.size();i++){
      for (auto next_path : Path_Matrix[i]){
	path_T[k-1].push_back(next_path);
	int lambda=calculate_lambda(k-1,next_path);
	max_lambda=max(lambda,max_lambda);
	min_lambda=min(lambda,min_lambda);
	avg_lambda+=lambda;
	lambda_points++;
	if(debug)
	  cout<<"\tcalculate_lambda returns:"<<lambda<<",Path_L.size:"<<path_L.size()<<endl;
	if(lambda>=path_L.size()){
	  if(debug)
	    cout<<"\tresizing path_L to new higher lambda of:"<<lambda<<endl;
	  path_L.resize(lambda+1);
	}
	path_L[lambda].push_back(next_path);
      }
    }
    return finished;//we are finished
  }


  //SO there are safe nodes (potentially, lets deal with them)
  for(int i=0;i<Path_Matrix.size();i++){
    set<vector<int>,vect_comp_int > surviving_paths ;
    for (auto next_path : Path_Matrix[i]){
      if((k+1)>=int(next_path.size())){//if k index is length of path then delete path
	if(k>path_T.size()){
	  path_T.resize(k+1);
	}
	path_T[k-1].push_back(next_path);
	avg_lambda+=0;
	lambda_points++;
	if(path_L.size()==0){
	  path_L.resize(1);
	}
	int lambda=calculate_lambda(k-1,next_path);
	path_L[0].push_back(next_path);
        if(debug)
	    cout<<"iteration:"<<iteration<<"path:"<<next_path<<",passed path length:"<<next_path.size()<<",k+2:"<<k+1<<",lambda:"<<lambda<<flush<<endl;
	continue;
      }
      if(safe_nodes.find(next_path[k])!=safe_nodes.end()){//add to surviving paths
	surviving_paths.insert(next_path);
      }
      else{
	if(debug){
	  cout<<"\tEliminating disclosing path:";for(auto j : next_path) cout<<j<<",";
	  cout<<endl;
	}
	if(k>path_T.size()){
	  path_T.resize(k+1);
	}
	path_T[k-1].push_back(next_path);
	int lambda=calculate_lambda(k-1,next_path);
	max_lambda=max(lambda,max_lambda);
	min_lambda=min(lambda,min_lambda);
	avg_lambda+=lambda;
	lambda_points++;

	if(debug)
	  cout<<"\tcalculate_lambda returns:"<<lambda<<",Path_L.size:"<<path_L.size()<<endl;
	if(lambda>=path_L.size()){
	  if(debug)
	    cout<<"\tresizing path_L to new higher lambda of:"<<lambda<<endl;
	  path_L.resize(lambda+1);
	}
	path_L[lambda].push_back(next_path);
	if(debug)
	  cout<<"adding to path_L at"<<lambda<<endl;
	disclosing_paths++;
      }
    }
    New_Path_Matrix.push_back(surviving_paths);
  }

  if(get_remaining_dests()<2){//so at best only one remaining dest
    finished=true;
  }
  

  if(!finished){//re-set paths with backups
    Path_Matrix=New_Path_Matrix;
  }
  erase_false_undisclosing_paths(k);
  return finished;
}
//Just finding "safe nodes" is not enough
//we need to ensure that prefix exist in at least one more path
//going to an alternative destination
void PathMatrix::erase_false_undisclosing_paths(int k){
  vector<set<vector<int>,vect_comp_int > > New_Path_Matrix;
  set<vector<int>,vect_comp_int > surviving_paths;
  for(int i=0;i<Path_Matrix.size();i++){
    surviving_paths.clear();
    for (auto next_path : Path_Matrix[i]){
      if(debug){
	cout<<"\tworking on path "<<next_path<<",prefix:";
	for(int j=0;j<k+1;j++)
	  cout<<next_path[j]<<",";
	cout<<endl;
      }
      bool alternative_exists=false;
      for(int j=0;j<Path_Matrix.size();j++){
	if(i==j)//same destination, skip sub-matrix
	  continue;
	for (auto alt_path : Path_Matrix[j]){//Check for identical sub-path
	  bool identical_prefix=true;
	  for(int l=1;l<k+1;l++){
	    if(next_path[l]!=alt_path[l]){
	      identical_prefix=false;
	      break;
	    }
	  }
	  if(identical_prefix){
	    alternative_exists=true;
	    if(debug){
	      cout<<"\t path:"<<next_path<<" and alt_path:"<<alt_path<<" have same prefix, k:"<<k<<endl;
	      //cout<<"before insert"<<flush<<",next_path.size:"<<next_path.size()<<flush<<endl;
	      }
	    surviving_paths.insert(next_path);
	    break;
	  }
	}
	if(alternative_exists){
	  if(debug)
	    cout<<"alternative_exists"<<flush<<endl;
	  break;
	}
      }
      if(!alternative_exists){
	if(debug)
	  cout<<"\tpath:"<<next_path<<"is actually disclosing, no common prefixes to any alternative destination(s)"<<flush<<endl;
	//cout<<"hola1,path_T.size:"<<path_T.size()<<",k-1:"<<k-1<<flush<<endl;
	if(k>path_T.size()){
	  path_T.resize(k);
	}
	path_T[k-1].push_back(next_path);
	//cout<<"hola2"<<flush<<endl;
	int lambda=calculate_lambda(k-1,next_path);
	max_lambda=max(lambda,max_lambda);
	min_lambda=min(lambda,max_lambda);
	avg_lambda+=lambda;
	lambda_points++;
	if(lambda>=path_L.size()){
	  path_L.resize(lambda+1);
	}
	//cout<<"hola3,lambda:"<<lambda<<",path_L.size:"<<path_L.size()<<flush<<endl;
	path_L[lambda].push_back(next_path);
	//cout<<"hola4"<<flush<<endl;
      }
    }
    //cout<<"hola5"<<flush<<endl;
    New_Path_Matrix.push_back(surviving_paths);
    //cout<<"hola6"<<flush<<endl;
  }
  Path_Matrix=New_Path_Matrix;
}
void PathMatrix::print_t_indexes(int i){
  int total_t_paths=0;
  cout<<"NOTE:ASSUMING UNIT COST FOR EDGES"<<endl;
  cout<<"path_T.size:"<<path_T.size()<<endl;
  for(int k_index=0;k_index<path_T.size();k_index++){
    //set<vector<int>,vect_comp_int > surviving_paths;
    if(path_T[k_index].size()==0)
      continue;
    total_t_paths+=path_T[k_index].size();
    if(!debug){
      cout<<"\tpaths_T_"<<i<<"_["<<k_index<<"]:"<<path_T[k_index].size()<<endl;
      continue;
    }
    for(int path_index=0;path_index<path_T[k_index].size();path_index++){
      if(path_T[k_index][path_index].size()==0)
	continue;
      cout<<"P_"<<i<<"_[,"<<path_T[k_index][path_index]<<",";
      cout<<"],t_index="<<k_index<<endl;
    }
  }
  cout<<"total_t_paths_"<<i<<"_:"<<total_t_paths<<endl;
  print_l_indexes(iteration);
}
void PathMatrix::print_l_indexes(int i,bool full_print){
  cout<<"path_L.size:"<<path_L.size()<<endl;
  if(debug)
    full_print=true;
  int total_l_paths=0;
  for(int k_index=0;k_index<path_L.size();k_index++){
    if(path_L[k_index].size()==0)
      continue;
    total_l_paths+=path_L[k_index].size();
    if(!full_print){
      cout<<"\tpaths_L_"<<i<<"_["<<k_index<<"]:"<<path_L[k_index].size()<<endl;
      continue;
    }
    for(int path_index=0;path_index<path_L[k_index].size();path_index++){
      if(path_L[k_index][path_index].size()==0)
	continue;
      cout<<"P_"<<i<<"[,"<<path_L[k_index][path_index]<<",";
      cout<<"],Lambda_index="<<k_index<<endl;
    }
  }
  cout<<"total_L_paths_"<<i<<"_:"<<total_l_paths<<endl;
}
//PRECALCULATES "U-",the set of paths cover_path protects
//by sharing the same prefix up to k-level.  Paths t and lambda 
//values do not need to be recalculated if at least one path remains
//at k-max level.  This is a one-time operation.  From here on, if a path is removed
//we just have to check on all other paths if it exists at k-max value.  If it does, we remove them
//and check if at least one path remains.  If at least one path remains, no recalculation needed.
void PathMatrix::get_dependent_set_for_path(vector<int> cover_path){
  //record_initial_paths_gperf();exit(1);
  int k=1;
  vector<set<vector<int>,vect_comp_int > > temp_element;

  
  for(int i=0;i<Path_Matrix.size();i++){
    if(destinations[i]==cover_path.back())
      continue;//only covering paths on alternate destinations
    cout<<"\t\tDest:"<<destinations[i]<<",cover_path:"<<cover_path<<endl;
    for (auto next_path : Path_Matrix[i]){
      //cout<<"\t\t\tWorking on candidate_path:"<<next_path;
      //cout<<",cover_path:"<<cover_path<<endl;
      //We check both paths up to the last common vertix
      //from the first node, all paths have origin in common
      //no reason to record this.
      for(k=1;k<min(next_path.size(),cover_path.size());k++){
	if(next_path[k]!=cover_path[k])
	  break;
      }
      if(k==1)
	continue;//only common node is origin
      //cout<<"\t\t\t\t prefix common up to:"<<k-1<<endl;
      //Find if next_path already in cover_map
      auto ret=cover_map.insert(make_pair(next_path,temp_element));
      auto it=ret.first;
      //cout<<"current_size:"<<it->second.size()<<endl;
      //Not add entry to cover set, resizing if necessary
      if(k>it->second.size()){
	it->second.resize(k);
      }
      it->second[k-1].insert(cover_path);
      cout<<"\t\t\tpath:"<<next_path<<" covered by:"<<cover_path<<" at k_level:"<<k-1<<endl;
	//insert(make_pair(cover_path,temp_set);
    }
  }
}
void PathMatrix::populate_cover_matrix(){
  cout<<"Populating cover_matrix"<<endl;
  //auto get_hash = cover_matrix.hash_function();
  int total_paths=get_number_paths_PM();
  indiv_path_values_L.resize(total_paths);
  //initialize cover_matrix
  cover_matrix.reserve(total_paths);
  vector<int> temp_row(total_paths,INT_MAX);
  cover_matrix.assign(total_paths,temp_row);
  //assign an id to each counter, for now using map
  //hash (of vector<ints>) should be more efficient, need to test 
  //how much of a difference it makes
  start_time = std::chrono::high_resolution_clock::now();
  int path_counter=0;
  for(int i=0;i<Path_Matrix.size();i++){
    for (auto next_path : Path_Matrix[i]){
      auto ret=path_ids.emplace(next_path,path_counter);
      if(debug)
	cout<<"path_id["<<next_path<<"]:"<<path_counter<<endl;
      unordered_map<vector<int>,size_t,VectorHash>::iterator it=ret.first;
      active_columns.emplace(path_counter,it);
      active_rows.emplace(path_counter++,it);
    }
  }
  auto end_time = std::chrono::high_resolution_clock::now(); std::chrono::duration<double> diff = end_time-start_time;auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
  cout<<"Time(ms) to populate path_ids:,"<<ms<<endl;
  start_time=end_time;
  //NOW POPULATE ALL LAMBDA VALUES IN COVER MATRIX
  //ROW IS COVERED PATH
  //COLUM IS WHAT COVER VALUE EACH OTHER PATH PROVIDE
  //DEFAULT NON COVER IS INT_MAX
  for (auto cover_path : path_ids){
    if(debug)
      cout<<"working on cover_path:"<<cover_path.first<<flush<<endl;
    for (auto cand_path : path_ids){
      if(cover_path.first.back()==cand_path.first.back())//skipping same destination paths
	continue;
      int k=1;
      for(k=1;k<min(cand_path.first.size(),cover_path.first.size());k++){
	if(cand_path.first[k]!=cover_path.first[k])
	  break;
      }
      if(k==1)
	continue;//only common node is origin
      if(debug)
	cout<<"cover_matrix["<<cover_path.second<<"]["<<cand_path.second<<"]"<<flush<<endl;
      cover_matrix[cover_path.second][cand_path.second]=calculate_lambda(k-1,cover_path.first);
      if(debug)
	cout<<"cover_matrix["<<cover_path.second<<"]["<<cand_path.second<<"]=,"<<cover_matrix[cover_path.second][cand_path.second]<<endl;
    }
  }
  end_time = std::chrono::high_resolution_clock::now(); diff = end_time-start_time; 
  ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
  cout<<"Time(ms) to populate cover_matrix:,"<<ms<<endl;
  MatrixTime=ms;
  start_time=end_time;
}
void PathMatrix::calculate_path_min_lambdas_from_cover_matrix(){
  static bool first_call=true;
  bool skip_recalculation=false;
  start_time = std::chrono::high_resolution_clock::now();
  int path_counter=0;
  int path_min_lambda=INT_MAX;
  min_lambda=INT_MAX;
  max_lambda=0;
  //Now populate min lambda values per path:
  path_L.clear();
  for (auto current_path : active_rows){
    //cout<<"\tgetting path_min_lambda for path["<<path_counter++<<"]"<<endl;
    if(first_call){
      //cout<<"path:"<<current_path.second->first<<endl;exit(1);
      path_min_lambda=current_path.second->first.size()-2;
    }
    else{
      path_min_lambda=INT_MAX;
    }
    if(min_cal&&!first_call){
      skip_recalculation=false;
      //cout<<"removed_paths:"<<removed_paths<<endl;
      int min_lambda_removed=INT_MAX;
      for(auto removed_id : removed_paths){
	min_lambda_removed=min(min_lambda_removed,cover_matrix[current_path.second->second][removed_id]);
	//cout<<"\t\tcurrent_path_id:"<<current_path.second;
	//cout<<",cover_val:"<<cover_matrix[current_path.second->second][removed_id]<<flush<<endl;
      }
      //cout<<"\tmin_lambda_removed:"<<min_lambda_removed<<",current_val:,"<<indiv_path_values_L[current_path.second->second]<<flush<<endl;
      if(min_lambda_removed>indiv_path_values_L[current_path.second->second]){
	skip_recalculation=true;
	path_min_lambda=indiv_path_values_L[current_path.second->second];
      }
    }
    //So skip if all paths removed had higher lambdas than minimal value for the path
    if( (!skip_recalculation) ||first_call){
      for (auto cover_path : active_columns){
	//cout<<"\t\tcurrent_path_id:"<<current_path.second->first<<",cover_path_id:"<<cover_path.second->first;
	path_min_lambda=min(path_min_lambda,cover_matrix[current_path.second->second][cover_path.second->second]);
	//cout<<",path_min_lambda:"<<path_min_lambda<<flush<<",path_L.size:"<<path_L.size()<<endl;
      }
    }
    
    if(path_min_lambda==INT_MAX){//if no cover_paths, lambda is maximal for path
      path_min_lambda=max(0,int(current_path.second->first.size())-2);
      //path_min_lambda=indiv_path_values_L[current_path.first];
    }

    if(path_min_lambda+1>path_L.size()){
      path_L.resize(path_min_lambda+1);
    }
    //populate indiv path values and max/min lambdas
    path_L[path_min_lambda].push_back(current_path.second->first);
    indiv_path_values_L[current_path.first]=path_min_lambda;//usefult to quickly check if path removal changes path best value
    min_lambda=min(min_lambda,path_min_lambda);
    max_lambda=max(max_lambda,path_min_lambda);
  }
  print_l_indexes(iteration,false);
  auto end_time = std::chrono::high_resolution_clock::now(); auto diff = end_time-start_time;
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
  first_call=false;
  cout<<"Time(ms) to populate path_L from cover_matrix:,"<<ms<<endl;
}
//REMOVE ALL PATHS FROM MAXIMAL:
//1)Eliminate form path_L all maximal paths
//2)check path_L still connected,if not we are finished, succesful exit
//3)eliminate removed paths from corresponding active rows and components
//4)return true if overall maximal lambda has been reduced, false otherwise
bool PathMatrix::remove_maximal_lambda_paths(){
  removed_paths.clear();
  iteration++;
  if(path_L.size()==1){
    cout<<"Finished,iteration:"<<iteration<<",max_lambda=0"<<endl;
    print_exit_statement();
    exit(0);
  }
  while(old_max_lambda+1<=path_L.size()){
    if(path_L.size()==0){
      cerr<<"Cannot call elminate_max_lambda_paths_U with no lambda paths!!!"<<endl;
      exit(1);
    }
    cout<<"Eliminating all paths with lambda="<<path_L.size()-1<<endl;
    for( auto next_path : path_L[path_L.size()-1]){
      int id=path_ids[next_path];
      active_rows.erase(id);
      active_columns.erase(id);
      removed_paths.push_back(id);
    }
    path_L.pop_back();
  }
  //Check lambdas still connected
  if(is_lambda_connected()!=destinations.size()){
    auto end_time = std::chrono::high_resolution_clock::now(); std::chrono::duration<double> diff = end_time-orig_start_time;auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
    cout<<"iteration:,"<<iteration<<",time:,"<<ms<<",finished, remaining paths are not fully connected,destinations:"<<destinations<<",reached:"<<is_lambda_connected()<<endl;
    print_exit_statement();
    exit(0);
  }
  //Recalculate affected lambdas in the future
  //for now just recalculate them
  calculate_path_min_lambdas_from_cover_matrix();
  if(path_L.size()<old_max_lambda){
    return true;
  }
  return false;
}

void PathMatrix::get_dependent_set_for_all_paths(){
  for(int i=0;i<Path_Matrix.size();i++){
    cout<<"Dest:"<<destinations[i]<<endl;
    int j=0;
    for (auto next_path : Path_Matrix[i]){
      cout<<"\tnext_path["<<j++<<"]:"<<next_path<<endl;
      get_dependent_set_for_path(next_path);
    }
  }
}
void PathMatrix::print_dependent_set(){
  for(auto el1 : cover_map){
    cout<<"____COVER FOR PATH:"<<el1.first<<endl;
    for(int k=0;k<el1.second.size();k++){
      for(auto el2 : el1.second[k]){
	cout<<"\t  k:"<<k<<",Path:"<<el2<<endl;
      }
    }
  }
}
void PathMatrix::elminate_max_lambda_paths_U(){
  //We are at the next iteration
  iteration++;
  //We always get rid of any path whose lambda value is better than the best ever value
  auto end_time = std::chrono::high_resolution_clock::now(); std::chrono::duration<double> diff = end_time-orig_start_time; auto current_time=diff.count();
  cout<<"starting eliminate_max_lambda_paths,g_time:"<<current_time;
  cout<<"old_max_lambda:"<<old_max_lambda<<"max_lambda="<<max_lambda<<endl;
  //path_L.size is lambda+1 due to lambdas with 0 value
  //but we never substract one from size_t just in case is 0 somehow
  while(old_max_lambda+1<=path_L.size()){
    if(path_L.size()==0){
      cerr<<"Cannot call elminate_max_lambda_paths_U with no lambda paths!!!"<<endl;
      exit(1);
    }
    cout<<"Eliminating all paths with lambda="<<path_L.size()-1<<endl;
    for( auto next_path : path_L[path_L.size()-1]){
      eliminate_max_lambda_path_U(next_path);
    }
    set_path_values_from_cover_paths();
  cout<<"finished eliminating paths with old_max_lambda=,"<<old_max_lambda;
  cout<<",new path_L.size():,"<<path_L.size()<<endl;
  }
}
    
void PathMatrix::eliminate_max_lambda_path_U(vector<int> elim_path){
  vector<int> test_path1{0, 1, 2, 3, 7, 6, 5, 4, 8, 9, 13, 14, 10, 11, 15};
  vector<int> test_path2{0, 1, 2, 3, 7, 6, 5, 4, 8, 9, 13, 12};
  //first delete all alternative paths covering the eliminated path
  cover_map.erase(elim_path);
  //set<vector<int>,vect_comp_int > fully_eliminated_paths;
  //secondly eliminate the path from covering all alternative paths
  bool path_found=false;
  for(auto el1 : cover_map){
    path_found=false;
    for(int k=0;k<el1.second.size();k++){
      for(auto el2 : el1.second[k]){
	if(vect_reverse_comparator(elim_path,el2)){
	  path_found=true;
	  el1.second[k].erase(elim_path);
	  //if(el1.second[k].size()==0)
	    //partial_elim.insert(make_pair(el1.first,k));
	  if(debug)
	    cout<<"\t Eliminated path:"<<elim_path<<" from covert_set"<<el1.first<<endl;
	  if(el1.first==test_path2){
	    //&&elim_path==test_path1)
	    cout<<"UPDATED COVER LIST FOR P:"<<el1.first<<endl;
	    //for(int k_group=0;k_group<el1.second.size();k_group++){
	    int k_group=el1.second.size()-1;k_group=max(0,k_group);
	      for(auto x : el1.second[k_group]){
		cout<<"\t  k:"<<k_group<<",Path:"<<x<<endl;
	      }
	    //}
	  }
	  break;
	  //if(el1.second.size()==0){//last covering path
	  //  fully_eliminated_paths.insert(el1.first);
	  }
	}
      if(path_found){//elim_path can only be present once
	break;
      }
    }
  }
  //for (auto el1 :fully_eliminated_paths){
  //  cover_map.erase(el1);
  //}
}

void PathMatrix::set_path_values_from_cover_paths(){
  cout<<"PRINTING DEPENDENT SET BEFORE RESETTING LAMBDAS,ITERATION:"<<iteration<<endl;
  print_dependent_set();
  path_T.clear();
  path_L.clear();
  clear_lambda_values();
  int number_paths=0;
  for(auto el1 : cover_map){
    int k=el1.second.size();
    //cout<<"path:"<<el1.first<<",t_index:"<<k-1<<endl;
    if(k>path_T.size()){
      path_T.resize(k);
    }
    path_T[k-1].push_back(el1.first);
    number_paths++;
    int lambda=calculate_lambda(k-1,el1.first);
    min_lambda=min(min_lambda,lambda);
    lambda_points++;
    avg_lambda+=lambda;
    //cout<<"lambda:"<<lambda<<flush<<endl;
    if(lambda+1>path_L.size()){
      path_L.resize(lambda+1);
    }
    path_L[lambda].push_back(el1.first);
  }
  max_lambda=path_L.size()-1;


  //print paths
  //bool old_debug=debug;
  //debug=true;
  print_final_results();
  //debug=old_debug;
  
  if(lambda_points==0){
    avg_lambda=0;
  }
  avg_lambda=avg_lambda/double(lambda_points);
  //timings
}

//add to candidate removal list all paths which make
//target distance to goal greatest.
//We get all paths but the ones that make lambda value maximal as the
//new candidate paths
void PathMatrix::get_maximal_lambda_paths(){
  populated_lambda_values=0;
  Candidate_Path_Matrix.clear();
  Candidate_Path_Matrix.resize(destinations.size());

  for(int k_index=0;k_index<path_L.size();k_index++){
    if(path_L[k_index].size()==0){
      cout<<"no paths at Lambda:"<<k_index<<endl;
      continue;
    }
    //removing top values
    if(k_index+1!=path_L.size()){
      if(debug)
	cout<<"\t adding to surviving paths at k_index:"<<k_index<<endl;
      for (auto next_path : path_L[k_index]){
	for(int i=0;i<destinations.size();i++){
	  if(next_path.back()==destinations[i]){
	    populated_lambda_values++;
	    Candidate_Path_Matrix[i].insert(next_path);
	    break;
	  }
	}
      }
    }
  }
  if(Candidate_Path_Matrix.size()<destinations.size()){//Still connected?
    cout<<"Finished,Candidate_Path_Matrix_size is:"<<Candidate_Path_Matrix.size()<<",destinations:"<<destinations.size()<<endl;
    print_exit_statement();
    exit(1);
  }
  if(debug)
    print_Candidate_Path_Matrix();
  cout<<"total surviving candidate paths:"<<get_number_paths(&Candidate_Path_Matrix)<<",populated_lambda_values:"<<populated_lambda_values<<endl;
}
void PathMatrix::clear_lambda_values(){
  min_lambda=INT_MAX;
  max_lambda=0;
  lambda_points=0;
  avg_lambda=0;
}
int PathMatrix::is_lambda_connected(){
  set<int> reached_dests;
  for(auto k : path_L){
    for(auto next_path : k){
      reached_dests.insert(next_path.back());
    }
    if(reached_dests.size()==destinations.size()){
      if(debug)
	cout<<"all destinations:";for(auto i : reached_dests) cout<<","<<i;cout<<" reached"<<endl;
      return destinations.size();
    }
  }
  return reached_dests.size();
}
int PathMatrix::remove_paths_higher_lambda(){
  cout<<"Removing paths whose lambda value>max_lambda:"<<old_max_lambda<<endl;
  path_L.resize(old_max_lambda);
  if(debug)
    print_l_indexes(iteration);
  //CHECK IF GRAPH STILL CONNECTED
  //IF GRAPH IS CONNECTED WE R FINISHED
  int connected_destinations=is_lambda_connected();
  if(connected_destinations!=destinations.size()){
  auto end_time = std::chrono::high_resolution_clock::now(); std::chrono::duration<double> diff = end_time-orig_start_time;auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
    cout<<"iteration:,"<<iteration<<",time:,"<<ms<<",Finished, remaining paths are disconnected,only "<<connected_destinations<<" reached out of "<<destinations.size()<<endl;
    print_exit_statement();
    exit(0);
  }
  else if(debug){
    for (auto dest : destinations) cout<<dest<<",";
    cout<<endl;
  }
  //After removal, we need to recalculate lambdas!
  recalculate_Path_Matrix_from_lambdas();
  return 0;
}
void PathMatrix::recalculate_Path_Matrix_from_lambdas(){
  populated_lambda_values=0;
  Candidate_Path_Matrix.clear();
  Candidate_Path_Matrix.resize(destinations.size());
  for(int k_index=0;k_index<path_L.size();k_index++){
    if(path_L[k_index].size()==0){
      cout<<"no paths at Lambda:"<<k_index<<endl;
      continue;
    }
      
    for (auto next_path : path_L[k_index]){
      for(int i=0;i<destinations.size();i++){
	if(next_path.back()==destinations[i]){
	  populated_lambda_values++;
	  Candidate_Path_Matrix[i].insert(next_path);
	  break;
	}
      }
    }
  }
}

void PathMatrix::check_for_improvement_on_lambda_values(){
  bool disconnected_graph=false;
  bool recalculate_lambdas_now=false;
  for(iteration=1;iteration<iter_limit;iteration++){
    clear_lambda_values();

    if(!recalculate_lambdas_now){
      get_maximal_lambda_paths();
    }
    else{
      recalculate_lambdas_now=false;
    }

    int candidate_paths=0;
    for(auto dest_list : Candidate_Path_Matrix){
      for(auto next_path : dest_list){
	candidate_paths++;
      }
    }
    if(candidate_paths==0){
      cout<<"Finshed after "<<iteration<<" iterations, all surviving candidate paths have same lambda now"<<endl;
      print_exit_statement();
      exit(0);
    }
    Path_Matrix=Candidate_Path_Matrix;
    //clear all t and l values for the surviving paths
    path_L.clear();path_T.clear();
    //Recalculate l and t
    int disclosing_paths=0;
    for(int k=1;k<1000;k++){
      get_safe_nodes2(k,2);
      if(debug)
	cout<<"after get_safe_nodes2 with k "<<k<<flush<<endl;
      if(erase_disclosing_paths(k,disclosing_paths)){
	cout<<"Lambda_Improv,Finished, erasing path would leave graph disconnected,final T score:"<<k<<endl;
	disconnected_graph=true;
	break;
      }
      else if(debug){
	cout<<"check_for_improvment_on_lambda_values,AFTER ELIMINATING DISCLOSING PATHS("<<k<<"), SURVIVING PATHS:"<<endl;
	cout<<"____________________________"<<endl;
	print_paths();
	cout<<"____________________________"<<endl;
      }
    }
      
    if(is_lambda_connected()!=destinations.size()){
      auto end_time = std::chrono::high_resolution_clock::now(); std::chrono::duration<double> diff = end_time-orig_start_time;auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
      cout<<"iteration:,"<<iteration<<",time:,"<<ms<<",finished, remaining paths are not fully connected,destinations:"<<destinations<<",reached:"<<is_lambda_connected()<<endl;
      print_exit_statement();
      exit(0);
    }

    print_t_indexes(iteration);
    if(lambda_points==0){
      avg_lambda=0;
    }
    avg_lambda=avg_lambda/double(lambda_points);

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time-start_time;
    start_time = end_time;//restart
    //OPTIMIZING LAMBDA
    cout<<"iteration:,"<<iteration<<",time(ms):,"<<diff.count()<<",new_avg_lambda:,"<<avg_lambda<<",new_min_lambda:,"<<min_lambda<<",new max_lambda:,"<<max_lambda<<",old_max_lambda:,"<<old_max_lambda<<",remaining_paths:,"<<get_number_paths(&Candidate_Path_Matrix)<<endl;
    if(max_lambda>old_max_lambda){
      cout<<"\titeration:,"<<iteration<<",Removing paths with max_lambda of "<<old_max_lambda<<" has resulted in not reducing lambda "<<max_lambda<<",applying 3c) step, removing all paths whose lambda>"<<old_max_lambda<<endl;
      remove_paths_higher_lambda();
      recalculate_lambdas_now=true;
    }
    else if(max_lambda<old_max_lambda){//we record biggest set of paths with the best value
      old_max_lambda=max_lambda;

      if(debug){
	cout<<"Lambda_Improved, SURVIVING PATHS:"<<endl;
	cout<<"____________________________"<<endl;
	print_paths();
	cout<<"____________________________"<<endl;
      }
      record_best_solution();
    }
  }
}

void PathMatrix::print_Candidate_Path_Matrix(){
  cout<<"___Candidate_Path_Matrix___"<<endl;
  for(int i=0;i<Candidate_Path_Matrix.size();i++){
    cout<<"Dest:"<<destinations[i]<<endl;
    int j=0;
    for (auto next_path : Candidate_Path_Matrix[i]){
      cout<<"\tnext_path["<<j++<<"]:";
      for(auto k : next_path) 
	cout<<k<<",";
      cout<<"length:"<<next_path.size()<<endl;
    }
  }
}
//for easy access from main
int PathMatrix::get_number_paths_PM(){
  int total_paths=0;
  for(auto paths_set : Path_Matrix){
    total_paths+=paths_set.size();
  }
  return total_paths;
}
int PathMatrix::get_remaining_dests(){
  int remaining_dests=0;
  for(auto paths_set : Path_Matrix){
    if(paths_set.size()>0){
      remaining_dests++;
    }
  }
  return remaining_dests;
}
int PathMatrix::get_number_paths(vector<set<vector<int>,vect_comp_int > > *PM){
  int total_paths=0;
  for(auto paths_set : *PM){
    total_paths+=paths_set.size();
  }
  return total_paths;
}
int PathMatrix::get_number_paths_lambda(){
  populated_lambda_values=0;
  for(int k_index=0;k_index<path_L.size();k_index++){
    populated_lambda_values+=path_L[k_index].size();
  }
  return populated_lambda_values;
}
void PathMatrix::print_W(){
  cout<<"All nodes observed excluding:";
  for(auto w : W) cout<<w<<",";
  cout<<endl;
}
void PathMatrix::print_final_results(){
  print_t_indexes(iteration);
  print_W();
  auto end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = end_time-start_time;
  start_time=end_time;//restart
  cout<<"iteration:,"<<iteration<<",time(ms):,"<<diff.count();
  cout<<",min_lambda=,"<<min_lambda<<",max_lambda=,"<<max_lambda;
  avg_lambda=avg_lambda/double(lambda_points);
  cout<<",avg_lambda=,"<<avg_lambda<<",old_max_lambda:"<<old_max_lambda;
  cout<<",remaining_paths:,"<<get_number_paths_lambda()<<endl;
  //old_max_lambda=max_lambda;
}
bool PathMatrix::vect_reverse_comparator(std::vector<int> v1,std::vector<int> v2) {
  //we only care if they are equal or not
  //return v1<v2;
  //shorter is smaller for our ordering
  if(v1.size()!=v2.size())
    return false;
  //if both vectors are empty
  if(v1.size()==0&&v2.size()==0)
    return true;
  
  //WE COMPARE VECTORS OF SAME SIZE FROM END TO BEGINNING
  //MUCH BETTER WHEN COMPARING PATHS WITH THE SAME ORIGIN POINT!
  for(int i=v1.size()-1;i>=0;i--){
    if(v1[i]!=v2[i])
      return v1[i]<v2[i];
  }
  return true;
}
int PathMatrix::get_max_lambda(){return max_lambda;};
int PathMatrix::get_old_max_lambda(){return old_max_lambda;};
void PathMatrix::set_max_lambda(int input){max_lambda=input;};
void PathMatrix::set_old_max_lambda(int input){old_max_lambda=input;};
void PathMatrix::set_InitialLambda(int input){InitialLambda=input;};
void PathMatrix::set_min_cal(bool input){min_cal=input;};
void PathMatrix::print_exit_statement(){
  auto end_time = std::chrono::high_resolution_clock::now(); std::chrono::duration<double> diff = end_time-orig_start_time;auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
  if(InitialLambda==0){//Initiallambda does not get populate if first maximal removal is disconnected
    InitialLambda=old_max_lambda;
  }
  cout<<"Exit statement,N:,"<<N<<",P:,"<<InitialPaths<<",R:,"<<Destinations.size()<<",S:,"<<random_seed<<",Y:,"<<YenTime<<",M:,"<<MatrixTime<<",iterations:,"<<iteration<<",IML:,"<<InitialLambda<<",FL:,"<<old_max_lambda<<",LD:,"<<InitialLambda-old_max_lambda<<",SP:,"<<final_paths<<",OT:,"<<ms;
  cout<<",avg_dest_intradist:,"<<avg_dest_intradist<<",avg_dest_orig_dist:,"<<avg_dest_orig_dest<<endl;//cout<<",LH:,"<<limit_hits<<endl;

    cout<<",Dists_to_orig:,";
  for (auto dest : OrigDests) cout<<dest<<",";
}
