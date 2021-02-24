///////////////////////////////////////////////////////////////////////////////
///  DijkstraShortestPathAlg.cpp
///  The implementation of Dijkstra algorithm to get the shortest path of 
///  a pair of vertices in a graph. 
///
///  @remarks <TODO: insert remarks here>
///
///  @author Yan Qi @date 5/30/2010
/// 
/// $Id: DijkstraShortestPathAlg.cpp 65 2010-09-08 06:48:36Z yan.qi.asu $
///
///////////////////////////////////////////////////////////////////////////////

#include <set>
#include <map>
#include <vector>
#include "GraphElements.h"
#include "Graph.h"
#include "DijkstraShortestPathAlg.h"
BasePath* DijkstraShortestPathAlg::get_shortest_path( BaseVertex* source, BaseVertex* sink )
{

	determine_shortest_paths(source, sink, true);
	//cout<<"source:"<<source->getID()<<",sink:"<<sink->getID()<<endl;

	std::vector<BaseVertex*> vertex_list;
	std::map<BaseVertex*, double>::const_iterator pos = 
		m_mpStartDistanceIndex.find(sink);
	double weight = pos != m_mpStartDistanceIndex.end() ? pos->second : Graph::DISCONNECT;

	if (weight < Graph::DISCONNECT)
	{
		BaseVertex* cur_vertex_pt = sink;
		do 
		{
			vertex_list.insert(vertex_list.begin(), cur_vertex_pt);

			std::map<BaseVertex*, BaseVertex*>::const_iterator pre_pos = 
				m_mpPredecessorVertex.find(cur_vertex_pt);

			if (pre_pos == m_mpPredecessorVertex.end()) break;

			cur_vertex_pt = pre_pos->second;

		} while (cur_vertex_pt != source);

		vertex_list.insert(vertex_list.begin(), source);
	}
	return new BasePath(vertex_list, weight);
}
shared_ptr<BasePath> DijkstraShortestPathAlg::get_shortest_path2( BaseVertex* source, BaseVertex* sink )
{
	determine_shortest_paths(source, sink, true);

	std::vector<BaseVertex*> vertex_list;
	std::map<BaseVertex*, double>::const_iterator pos = 
		m_mpStartDistanceIndex.find(sink);
	double weight = pos != m_mpStartDistanceIndex.end() ? pos->second : Graph::DISCONNECT;

	if (weight < Graph::DISCONNECT)
	{
		BaseVertex* cur_vertex_pt = sink;
		do 
		{
			vertex_list.insert(vertex_list.begin(), cur_vertex_pt);

			std::map<BaseVertex*, BaseVertex*>::const_iterator pre_pos = 
				m_mpPredecessorVertex.find(cur_vertex_pt);

			if (pre_pos == m_mpPredecessorVertex.end()) break;

			cur_vertex_pt = pre_pos->second;

		} while (cur_vertex_pt != source);

		vertex_list.insert(vertex_list.begin(), source);
	}
	return make_shared<BasePath>(vertex_list, weight);
}

void DijkstraShortestPathAlg::determine_all_shortest_paths( BaseVertex* source, float lambda_limit)
{
			      
  current_source=source;
  vector<BaseVertex*> temp_list;temp_list.push_back(source);
  m_mpPerimPaths[source->getID()]=make_shared<BasePath>(temp_list,0);

  cout<<"determine_all_shortest_paths,souce:"<<source->getID()<<",lambda_limit:"<<lambda_limit<<flush<<endl;
	//1. clear the intermediate variables
	clear();

	//2. initiate the local variables
	BaseVertex* start_vertex = source;
	m_mpStartDistanceIndex[start_vertex] = 0;
	//auto best_distances= m_mpStartDistanceIndex;
	start_vertex->Weight(0);
	m_quCandidateVertices.insert(start_vertex);
	//cout<<"before while"<<flush<<endl;

	//3. start searching for the shortest path
	while (!m_quCandidateVertices.empty())
	{
		multiset<BaseVertex*, WeightLess<BaseVertex> >::const_iterator pos = m_quCandidateVertices.begin();

		BaseVertex* cur_vertex_pt = *pos; //m_quCandidateVertices.top();
		//cout<<"\tcur_vertex_pt:"<<cur_vertex_pt->getID()<<",source,:"<<source->getID()<<endl;
		m_quCandidateVertices.erase(pos);
	

		if(cur_vertex_pt->getID() == 31)
		{
			int i = 1;
			i = 100;
		}
		/*if(cur_vertex_pt->getID()==current_source->getID()&&cur_vertex_pt->Weight()>0){//No revisits of origin node
		  cout<<"skipping visiting origin nodes"<<endl;
		  continue;
		}*/
		m_stDeterminedVertices.insert(cur_vertex_pt->getID());
		//cout<<"\t\t\tcur_vertex:"<<cur_vertex_pt->getID()<<flush<<endl;

		improve2vertexKeepPaths(cur_vertex_pt, true,lambda_limit);
		//cout<<"m_mpStartDistanceIndex.size:,"<<m_mpStartDistanceIndex.size()<<endl;
		//if(best_distances.size()<m_mpStartDistanceIndex.size()){
		//  best_distances= m_mpStartDistanceIndex;
		  //cout<<"best_distances.size:,"<<best_distances.size()<<endl;
		//}
	}

	double perim_dist=0;
		for(auto const& x : m_mpStartDistanceIndex ){
		  perim_dist=max(perim_dist,x.second);
		}
		//cout<<"\t Source:,"<<source->getID()<<",Perimeter nodes @ dist:,"<<perim_dist<<",->,"<<endl;;
		for(auto const& x : m_mpStartDistanceIndex ){
		  /*//if(x.second==perim_dist){
		    cout<<"main_graph,"<<x.first->getID()<<",";
		    //Perim_nodes.insert(x.first->getID());
		    //cout<<","<<x.second<<",";
		    cout<<","<<x.second<<endl;
		  //}*/
		}
		//cout<<"All paths:"<<endl;
		//print_paths();
		//for(auto it : m_mpPerimPaths){
		//  cout<<"\t\tPerimNode:"<<it.first<<",OP->";it.second->PrintOut(cout);
		//}
}

void DijkstraShortestPathAlg::determine_shortest_paths( BaseVertex* source, BaseVertex* sink, bool is_source2sink)
{
	//1. clear the intermediate variables
	clear();

	//2. initiate the local variables
	BaseVertex* end_vertex = is_source2sink ? sink : source;
	BaseVertex* start_vertex = is_source2sink ? source : sink;
	m_mpStartDistanceIndex[start_vertex] = 0;
	start_vertex->Weight(0);
	m_quCandidateVertices.insert(start_vertex);
  
	//cout<<"\tdetermine_shortest_paths,source:,"<<source->getID()<<",sink:"<<sink->getID()<<",candidates:"<<m_quCandidateVertices.size()<<endl;

	//bool solution_found=false;
	//3. start searching for the shortest path
	while (!m_quCandidateVertices.empty())
	{
		multiset<BaseVertex*, WeightLess<BaseVertex> >::const_iterator pos = m_quCandidateVertices.begin();

		BaseVertex* cur_vertex_pt = *pos; //m_quCandidateVertices.top();
		m_quCandidateVertices.erase(pos);
		//cout<<"\t\tcur_vertex_pt->getID():"<<cur_vertex_pt->getID()<<endl;
	
		if (cur_vertex_pt == end_vertex){
		 //solution_found=true;
		 break;
		}

		if(cur_vertex_pt->getID() == 31)
		{
			int i = 1;
			i = 100;
		}
		m_stDeterminedVertices.insert(cur_vertex_pt->getID());

		improve2vertex(cur_vertex_pt, is_source2sink);
	}
	/*if(solution_found){
	  //cout<<"\tSolution_found,determine_shortest_paths,sink:"<<sink->getID()<<endl;
	}
	else{
	  //cout<<"\tNo_Solution_found,determine_shortest_paths,sink:"<<sink->getID()<<endl;
	}*/
}

void DijkstraShortestPathAlg::improve2vertexKeepPaths( BaseVertex* cur_vertex_pt, bool is_source2sink, float lambda_limit )
{
	// 1. get the neighboring vertices 
	set<BaseVertex*>* neighbor_vertex_list_pt = new set<BaseVertex*>();
		
	//cout<<"getting neighours"<<endl;
	if(is_source2sink)
	{
		m_pDirectGraph->get_adjacent_vertices(cur_vertex_pt, *neighbor_vertex_list_pt);
	}else
	{
		m_pDirectGraph->get_precedent_vertices(cur_vertex_pt, *neighbor_vertex_list_pt);
	}
	
	/*set<BaseVertex*>* neighbor_adjacent_vertex_list_pt = new set<BaseVertex*>();
	m_pDirectGraph->get_adjacent_vertices(cur_vertex_pt, *neighbor_adjacent_vertex_list_pt);
	cout<<"\t\t\tneighours_adjacent:"<<neighbor_adjacent_vertex_list_pt->size()<<endl;//exit(1);
	set<BaseVertex*>* neighbor_precedent_vertex_list_pt = new set<BaseVertex*>();
	m_pDirectGraph->get_precedent_vertices(cur_vertex_pt, *neighbor_precedent_vertex_list_pt);
	cout<<"\t\t\tneighours_precedent:"<<neighbor_precedent_vertex_list_pt->size()<<endl;//exit(1);*/

	// 2. update the distance passing on the current vertex
	for(set<BaseVertex*>::iterator cur_neighbor_pos=neighbor_vertex_list_pt->begin(); 
		cur_neighbor_pos!=neighbor_vertex_list_pt->end(); ++cur_neighbor_pos)
	{
		//2.1 skip if it has been visited before
		if (m_stDeterminedVertices.find((*cur_neighbor_pos)->getID())!=m_stDeterminedVertices.end())
		{
		  //cout<<"\t\t skipping visited vertex:"<<(*cur_neighbor_pos)->getID()<<endl;
			continue;
		}

		//2.2 calculate the distance
		map<BaseVertex*, double>::const_iterator cur_pos = m_mpStartDistanceIndex.find(cur_vertex_pt);
		double distance =  cur_pos != m_mpStartDistanceIndex.end() ? cur_pos->second : Graph::DISCONNECT;
		//if(cur_vertex_pt->getID()==84){
		 // cout<<"cur_vertex:,"<<cur_vertex_pt->getID()<<",parent distance is:"<<distance<<",lambda_limit:"<<lambda_limit<<endl;
		//}
		//cout<<"cur_vertex:,"<<cur_vertex_pt->getID()<<",distance is:"<<distance<<",lambda_limit:"<<lambda_limit<<endl;
		if(distance>=lambda_limit&&lambda_limit!=-1){
		  //cout<<"cur_vertex:,"<<cur_vertex_pt->getID()<<",distance is:"<<distance<<",lambda_limit:"<<lambda_limit<<",continue"<<endl;
		  continue;
		}

		//cout<<"\t\tcur_pos:"<<(cur_pos->first)->getID()<<",cur_distance:"<<distance<<endl;
		distance += is_source2sink ? m_pDirectGraph->get_edge_weight(cur_vertex_pt, *cur_neighbor_pos) : 
			m_pDirectGraph->get_edge_weight(*cur_neighbor_pos, cur_vertex_pt);
////MASSIVE HACK: SOMETHING GOING FUNNY WITH EDGES WEIGHTS, ALL VALUES ARE 1! CURRENTLY ALL OF OUR EXPERIMENTS ARE WITH UNIT-COST, SO FORCING IT INSTEAD
		//cout<<"\t\tcur_pos:"<<(cur_pos->first)->getID()<<"edge_weight0:"<< m_pDirectGraph->get_edge_weight(cur_vertex_pt, *cur_neighbor_pos)<<",edge_weight2:"<<m_pDirectGraph->get_edge_weight(*cur_neighbor_pos, cur_vertex_pt)<<endl;
		//distance +=1;
		/*if((*cur_neighbor_pos)->getID()==84){
		  cout<<"\t\t\toption 1:"<<m_pDirectGraph->get_edge_weight(cur_vertex_pt, *cur_neighbor_pos)<<",from:,"<<cur_vertex_pt->getID()<<",to:,"<<(*cur_neighbor_pos)->getID()<<",code:,"<<m_pDirectGraph->get_edge_code(cur_vertex_pt,*cur_neighbor_pos )<<endl;
		  cout<<"\t\t\toption 2:"<<m_pDirectGraph->get_edge_weight(*cur_neighbor_pos, cur_vertex_pt)<<",from:,"<<(*cur_neighbor_pos)->getID()<<",to:,"<<cur_vertex_pt->getID()<<endl;
		  cout<<"\t\t\tupdated distance:"<<distance<<endl;
		}*/



		//2.3 update the distance if necessary
		cur_pos = m_mpStartDistanceIndex.find(*cur_neighbor_pos);
		if (cur_pos == m_mpStartDistanceIndex.end() || cur_pos->second > distance)
		{
			m_mpStartDistanceIndex[*cur_neighbor_pos] = distance;
			//cout<<"\tcur_neighours:,"<<(*cur_neighbor_pos)->getID()<<",distance:,"<<distance<<",m_mpStartDistanceIndex.size:,"<<m_mpStartDistanceIndex.size()<<endl;
			m_mpPredecessorVertex[*cur_neighbor_pos] = cur_vertex_pt;

			//if(distance>=lambda_limit||lambda_limit!=-1){
			//KEEP PATHS
			  BaseVertex* cur_vertex_pt2 = *cur_neighbor_pos; //m_quCandidateVertices.top();
			  std::vector<BaseVertex*> vertex_list;
			  do 
			  {
			
			    if(lambda_limit==-1){//regular order, origin to perimeter node
			      vertex_list.insert(vertex_list.begin(), cur_vertex_pt2);
			    }
			    else{//reverse order, perimeter to destination
				  vertex_list.push_back(cur_vertex_pt2);
			    }

				  std::map<BaseVertex*, BaseVertex*>::const_iterator pre_pos = 
					  m_mpPredecessorVertex.find(cur_vertex_pt2);

				  if (pre_pos == m_mpPredecessorVertex.end()) break;

				  cur_vertex_pt2 = pre_pos->second;

			  } while (cur_vertex_pt2 != current_source);
			    if(lambda_limit==-1){//regular order, origin to perimeter node
			      vertex_list.insert(vertex_list.begin(), current_source);
			      //cout<<"\t\tvertex:"<<vertex_list.back()->getID()<<",weight1:"<<distance;
			      m_mpPerimPaths[vertex_list.back()->getID()]=make_shared<BasePath>(vertex_list,distance);
			      //cout<<",distance2="<<m_mpPerimPaths[vertex_list.back()->getID()]->Weight()<<endl;
			      /*if(vertex_list.back()->getID()==84){
				cout<<"84,\t\tvertex:"<<vertex_list.back()->getID()<<",weight1:"<<distance<<endl;
				cout<<"84,,distance2="<<m_mpPerimPaths[vertex_list.back()->getID()]->Weight()<<endl;
			      }*/
			    }
			    else{
			      vertex_list.push_back(current_source);
			      //cout<<"\t\tvertex:"<<vertex_list.front()->getID()<<",weight1:"<<distance;
			      m_mpPerimPaths[vertex_list.front()->getID()]=make_shared<BasePath>(vertex_list,distance);
			      //cout<<",distance2="<<m_mpPerimPaths[vertex_list.front()->getID()]->Weight()<<endl;
			      //cout<<"\t\t\t added "<<vertex_list.front()->getID()<<"to the list of perim nodes"<<endl;
			    }
			  //BasePath current_path(vertex_list,distance);


			  //for(auto vertex : vertex_list) cout<<vertex->getID()<<flush<<",";cout<<endl;//exit(0);

			//}
			
			(*cur_neighbor_pos)->Weight(distance);

			multiset<BaseVertex*, WeightLess<BaseVertex> >::const_iterator pos = m_quCandidateVertices.begin();
			for(; pos != m_quCandidateVertices.end(); ++pos)
			{
				if ((*pos)->getID() == (*cur_neighbor_pos)->getID())
				{
					break;
				}
			}
			if(pos != m_quCandidateVertices.end())
			{
				m_quCandidateVertices.erase(pos);
			}
			m_quCandidateVertices.insert(*cur_neighbor_pos);
		}
	}
}
void DijkstraShortestPathAlg::improve2vertex( BaseVertex* cur_vertex_pt, bool is_source2sink, float lambda_limit )
{
	// 1. get the neighboring vertices 
	set<BaseVertex*>* neighbor_vertex_list_pt = new set<BaseVertex*>();
		
	if(is_source2sink)
	{
		m_pDirectGraph->get_adjacent_vertices(cur_vertex_pt, *neighbor_vertex_list_pt);
	}else
	{
		m_pDirectGraph->get_precedent_vertices(cur_vertex_pt, *neighbor_vertex_list_pt);
	}
	//cout<<"\t\t source:"<<cur_vertex_pt->getID()<<",neighour_size:"<<neighbor_vertex_list_pt->size()<<endl;
	//cout<<"\t\t\tneighours:"<<neighbor_vertex_list_pt->size()<<endl;

	// 2. update the distance passing on the current vertex
	for(set<BaseVertex*>::iterator cur_neighbor_pos=neighbor_vertex_list_pt->begin(); 
		cur_neighbor_pos!=neighbor_vertex_list_pt->end(); ++cur_neighbor_pos)
	{
	  //cout<<"\t\t\tneihgbor:"<<(*cur_neighbor_pos)->getID()<<endl;

		//2.1 skip if it has been visited before
		if (m_stDeterminedVertices.find((*cur_neighbor_pos)->getID())!=m_stDeterminedVertices.end())
		{
			continue;
		}

		//2.2 calculate the distance
		map<BaseVertex*, double>::const_iterator cur_pos = m_mpStartDistanceIndex.find(cur_vertex_pt);
		double distance =  cur_pos != m_mpStartDistanceIndex.end() ? cur_pos->second : Graph::DISCONNECT;
		//cout<<"\t\t\tparent distance:"<<distance<<endl;
		if(distance>=lambda_limit){
		  continue;
		}

		distance += is_source2sink ? m_pDirectGraph->get_edge_weight(cur_vertex_pt, *cur_neighbor_pos) : 
			m_pDirectGraph->get_edge_weight(*cur_neighbor_pos, cur_vertex_pt);
		//cout<<"\t\t\toption 1:"<<m_pDirectGraph->get_edge_weight(cur_vertex_pt, *cur_neighbor_pos)<<endl;
		//cout<<"\t\t\toption 2:"<<m_pDirectGraph->get_edge_weight(*cur_neighbor_pos, cur_vertex_pt)<<endl;
		//cout<<"\t\t\tupdated distance:"<<distance<<endl;
		//MASSIVE HACK! SOMETHING WRONG WITH WEIGHTS, fording edge weight to one because we r doing unit cost, FIX ME!!!
		//distance +=1;
		//cout<<"cur_neighour:"<<(*cur_neighbor_pos)->getID()<<",distance:"<<distance<<endl;

		//2.3 update the distance if necessary
		cur_pos = m_mpStartDistanceIndex.find(*cur_neighbor_pos);
		if (cur_pos == m_mpStartDistanceIndex.end() || cur_pos->second > distance)
		{
			m_mpStartDistanceIndex[*cur_neighbor_pos] = distance;
			//cout<<"cur_vertex:,"<<cur_vertex_pt->getID()<<",distance:,"<<distance<<",m_mpStartDistanceIndex.size:,"<<m_mpStartDistanceIndex.size()<<endl;
			m_mpPredecessorVertex[*cur_neighbor_pos] = cur_vertex_pt;

			
			(*cur_neighbor_pos)->Weight(distance);

			multiset<BaseVertex*, WeightLess<BaseVertex> >::const_iterator pos = m_quCandidateVertices.begin();
			for(; pos != m_quCandidateVertices.end(); ++pos)
			{
				if ((*pos)->getID() == (*cur_neighbor_pos)->getID())
				{
					break;
				}
			}
			if(pos != m_quCandidateVertices.end())
			{
				m_quCandidateVertices.erase(pos);
			}
			m_quCandidateVertices.insert(*cur_neighbor_pos);
		}
	}
}

void DijkstraShortestPathAlg::clear()
{
	m_stDeterminedVertices.clear();
	m_mpPredecessorVertex.clear();
	m_mpStartDistanceIndex.clear();
	m_quCandidateVertices.clear();
	m_mpPerimPaths.clear();
	m_mpFullPaths.clear();
}
DijkstraShortestPathAlg::DijkstraShortestPathAlg( void)
{
  clear();
}

BasePath* DijkstraShortestPathAlg::update_cost_forward( BaseVertex* vertex )
{
	double cost = Graph::DISCONNECT;

 	// 1. get the set of successors of the input vertex
	set<BaseVertex*>* adj_vertex_set = new set<BaseVertex*>();
	m_pDirectGraph->get_adjacent_vertices(vertex, *adj_vertex_set);
 
 	// 2. make sure the input vertex exists in the index
	map<BaseVertex*, double>::iterator pos4vertexInStartDistIndex = m_mpStartDistanceIndex.find(vertex);
	if(pos4vertexInStartDistIndex == m_mpStartDistanceIndex.end())
 	{
		pos4vertexInStartDistIndex = 
			(m_mpStartDistanceIndex.insert(make_pair(vertex, Graph::DISCONNECT))).first;
 	}

 	// 3. update the distance from the root to the input vertex if necessary
 	for(set<BaseVertex*>::const_iterator pos=adj_vertex_set->begin(); pos!=adj_vertex_set->end();++pos)
 	{
 		// 3.1 get the distance from the root to one successor of the input vertex
		map<BaseVertex*, double>::const_iterator cur_vertex_pos = m_mpStartDistanceIndex.find(*pos);
		double distance = cur_vertex_pos == m_mpStartDistanceIndex.end() ?
			Graph::DISCONNECT : cur_vertex_pos->second;
 
 		// 3.2 calculate the distance from the root to the input vertex
		distance += m_pDirectGraph->get_edge_weight(vertex, *pos);
	
 		// 3.3 update the distance if necessary 
		double cost_of_vertex = pos4vertexInStartDistIndex->second;
 		if(cost_of_vertex > distance)
 		{
			m_mpStartDistanceIndex[vertex] = distance;
			m_mpPredecessorVertex[vertex] = cur_vertex_pos->first;
 			cost = distance;
 		}
 	}

 	// 4. create the sub_path if exists
	BasePath* sub_path = NULL;
	if(cost < Graph::DISCONNECT) 
 	{
		vector<BaseVertex*> vertex_list;
		vertex_list.push_back(vertex);

		map<BaseVertex*, BaseVertex*>::const_iterator pos4PredVertexMap =
			m_mpPredecessorVertex.find(vertex);
		
		while(pos4PredVertexMap != m_mpPredecessorVertex.end())
		{
			BaseVertex* pred_vertex_pt = pos4PredVertexMap->second;
			vertex_list.push_back(pred_vertex_pt);
			pos4PredVertexMap = m_mpPredecessorVertex.find(pred_vertex_pt);
		}

		sub_path = new BasePath(vertex_list, cost);
 	}
 	return sub_path;
}

void DijkstraShortestPathAlg::correct_cost_backward( BaseVertex* vertex )
{
 	// 1. initialize the list of vertex to be updated
	vector<BaseVertex*> vertex_pt_list;
	vertex_pt_list.push_back(vertex);

	// 2. update the cost of relevant precedents of the input vertex
	while(!vertex_pt_list.empty())
 	{
		BaseVertex* cur_vertex_pt = *(vertex_pt_list.begin());
		vertex_pt_list.erase(vertex_pt_list.begin());

 		double cost_of_cur_vertex = m_mpStartDistanceIndex[cur_vertex_pt];

		set<BaseVertex*> pre_vertex_set;
		m_pDirectGraph->get_precedent_vertices(cur_vertex_pt, pre_vertex_set);
		for(set<BaseVertex*>::const_iterator pos=pre_vertex_set.begin(); pos!=pre_vertex_set.end();++pos)
		{
			map<BaseVertex*,double>::const_iterator pos4StartDistIndexMap = 
				m_mpStartDistanceIndex.find(*pos);
			double cost_of_pre_vertex = m_mpStartDistanceIndex.end() == pos4StartDistIndexMap ?
				Graph::DISCONNECT : pos4StartDistIndexMap->second;

			double fresh_cost = cost_of_cur_vertex + m_pDirectGraph->get_edge_weight(*pos, cur_vertex_pt);
			if(cost_of_pre_vertex > fresh_cost)
			{
				m_mpStartDistanceIndex[*pos] = fresh_cost;
				m_mpPredecessorVertex[*pos] = cur_vertex_pt;
				vertex_pt_list.push_back(*pos);
			}
		}
	}
}
	
void DijkstraShortestPathAlg::make_full_paths(std::map<int,shared_ptr<BasePath> > *FullOriginPaths){
  for (auto it :  m_mpPerimPaths){
    if(FullOriginPaths->find(it.first)==FullOriginPaths->end()){
      //cerr<<"node "<<it.first<<" has no path from origin"<<endl;
      BasePath FullPath=*it.second;
      m_mpFullPaths[it.first]=make_shared<BasePath>(FullPath);
      continue;
    }
    //cout<<"\tnode:"<<it.first<<", perim_path:";it.second->PrintOut(cout);
    BasePath FullPath=*((*FullOriginPaths)[it.first]);
    FullPath.append_to_path(it.second);

    //std::map<int,shared_ptr<BasePath> > m_mpFullPaths;
    m_mpFullPaths[it.first]=make_shared<BasePath>(FullPath);
    //*(it.second)=FullPath;
    //cout<<"\t\t full_path:";it.second->PrintOut(cout);
    //cout<<"\t\t full_path:"<<m_mpFullPaths[it.first]->get_path_string()<<endl;
  }
}
void DijkstraShortestPathAlg::print_paths(){
  //cout<<"hola print_paths"<<endl;
  for (auto it :  m_mpPerimPaths){
    cout<<"\tnode:"<<it.first<<", perim_path:";it.second->PrintOut(cout);
  }
}
	
void DijkstraShortestPathAlg::get_perim_nodes(set<int> *perim_nodes,int lambda_limit){
  for (auto it :  m_mpPerimPaths){
    if(it.second->Weight()>=lambda_limit){
      perim_nodes->insert(it.first);
    }
  }
}
//check if we can improve distances
//void DijkstraShortestPathAlg::improve_distances(unordered_map<pair<int,int>,float,hash_pair > *current_distances){
void DijkstraShortestPathAlg::improve_distances(unordered_set<int> *FO_nodes, vector<Link> *FO_edges,unordered_map<pair<int,int>,float,hash_pair> *best_d){
  //Note: We only ever get one distance for the subgraphs because we only call this function once per visible node
  //which is the origin node for this call.
  for (auto it :  m_mpPerimPaths){
    //auto current_pair=make_pair<int,int>(current_source->getID(),int(it.first));
    //if(current_distances->find(current_pair)==current_distances->end()){//populate if no value
      //(*current_distances)[current_pair]=it.second->Weight();
      if(FO_nodes->find(it.first)==FO_nodes->end()){//node is not visible so we do not need to update the distance
	continue;
      }
      //Ensure that we do not have already a better distance!
      auto temp_pair=make_pair<int,int>(current_source->getID(),int(it.first));
      auto got=best_d->find(temp_pair);
      if(got!=best_d->end()){
	if(got->second>=it.second->Weight())
	  continue;
      }
      (*best_d)[temp_pair]=it.second->Weight();
	
      Link temp_link{current_source->getID(),it.first,it.second->Weight()};
      FO_edges->push_back(temp_link);
      //cout<<"source:,"<<current_source->getID()<<","<<it.first<<","<<it.second->Weight()<<endl;
    //}
    //else if((*current_distances)[current_pair]>it.second->Weight()){//else update if found smaller path
    //  cout<<"\t found improvement for origin:,"<<
    //  (*current_distances)[current_pair]=it.second->Weight();
     //}
  }
}
int DijkstraShortestPathAlg::get_perim_size(){
  return m_mpPerimPaths.size();
}
shared_ptr<BasePath> DijkstraShortestPathAlg::recover_shortest_perim_path(int sink){
  if(m_mpPerimPaths.find(sink)!=m_mpPerimPaths.end()){
      return m_mpPerimPaths[sink];
      }
  else{ 
    cerr<<"Cannot recover shortest perimeter path for destination:"<<sink<<","<<",source:"<<current_source->getID();
    cerr<<",size:"<<m_mpPerimPaths.size()<<flush<<",stored_nodes:";
    for(auto it : m_mpPerimPaths)
      cerr<<it.first<<",";
    cerr<<endl;
  }
  exit(1);
}
bool DijkstraShortestPathAlg::is_node_reachable(int sink){
  if(m_mpPerimPaths.find(sink)!=m_mpPerimPaths.end()){
    return true;
  }
  return false;
}
//Cost from start to node
double DijkstraShortestPathAlg::get_cost(int sink){
  if(m_mpPerimPaths.find(sink)!=m_mpPerimPaths.end()){
    return m_mpPerimPaths[sink]->Weight();
  }
  else{
    //cout<<"node:"<<sink<<" has no associated cost from start, returning 0"<<endl;
    return 0;
  }
}
void DijkstraShortestPathAlg::print_perim_paths(){
  for (auto path : m_mpPerimPaths){
	std::map<int,shared_ptr<BasePath> > m_mpPerimPaths;
    cout<<"\t\tperim_path,sink:"<<path.first<<",source:"<<path.second->get_node(0)<<endl;
  }
}

shared_ptr<BasePath> DijkstraShortestPathAlg::recover_shortest_full_path(int sink){
  if(m_mpFullPaths.find(sink)!=m_mpFullPaths.end()){
      return m_mpPerimPaths[sink];
      }
  else{ 
    cerr<<"Cannot recover shortest path for destination:"<<sink<<",source:"<<current_source<<endl;
  }
  exit(1);
}
	
void DijkstraShortestPathAlg::populate_final_paths(set<BasePath> *Final_paths){
  for(auto it : m_mpFullPaths){
    Final_paths->insert(*(it.second));
  }
}
void DijkstraShortestPathAlg::dump_edges(){
  m_pDirectGraph->dump_edges();
}
//def stuck(x)
//   if x == t
//     return False
//   for each neighbor y of x
//     if y not in seen
//       insert y in seen
//       if !stuck(y)
//         return False
//   return True

//def search(x)
//  if x == t
//    print path
//  seen = set(path)
//  if stuck(x)
//    return
//  for each neighbor y of x
//    push y on the path
//    search(y)
//    pop y from the path