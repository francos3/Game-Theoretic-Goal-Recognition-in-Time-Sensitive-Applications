///////////////////////////////////////////////////////////////////////////////
///  GraphElements.h
///  <TODO: insert file description here>
///
///  @remarks <TODO: insert remarks here>
///
///  @author Yan Qi @date 5/28/2010
///  modified by Santiago Franco @date 17/06/2019
///
///  $Id: GraphElements.h 65 2010-09-08 06:48:36Z yan.qi.asu $
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <sstream>
#include <memory>
#include <iostream>
#include "main.h"



template<class T>
class WeightGreater
{
public:
	// Determine priority.
	bool operator()(const T& a, const T& b) const
	{
		return a.Weight() > b.Weight();
	}

	bool operator()(const T* a, const T* b) const
	{
		return a->Weight() > b->Weight();
	}
};

template<class T>
class WeightLess
{
public:
	// Determine priority.
	bool operator()(const T& a, const T& b) const
	{
		return a.Weight() < b.Weight();
	}

	bool operator()(const T* a, const T* b) const
	{
		return a->Weight() < b->Weight();
	}
};

//////////////////////////////////////////////////////////////////////////
// A class for the object deletion
//////////////////////////////////////////////////////////////////////////
template<class T>
class DeleteFunc
{
public:
	void operator()(const T* it) const
	{
		delete it;
	}
};



/**************************************************************************
*  BaseVertex
*  <TODO: insert class description here>
*
*
*  @remarks <TODO: insert remarks here>
*
*  @author Yan Qi @date 6/6/2010
**************************************************************************/
class BaseVertex
{
	int m_nID;
	double m_dWeight;
//	BaseVertex* parent;	

public:

	int getID() const { return m_nID; }
	void setID(int ID_) { m_nID = ID_; }

//	BaseVertex* getParent() const { return parent; }
//	void setParent(BaseVertex* parent_) { parent = parent_; }

	double Weight() const { return m_dWeight; }
	void Weight(double val) { m_dWeight = val; }

	void PrintOut(std::ostream& out_stream)
	{
	  std::string test;
	  std::stringstream ss;
	  ss << m_nID;
	  test = ss.str();
		//out_stream << m_nID;
	  out_stream << test;
	}
};


/**************************************************************************
*  BasePath
*  <TODO: insert class description here>
*
*
*  @remarks <TODO: insert remarks here>
*
*  @author Yan Qi @date 6/6/2010
**************************************************************************/
class BasePath
{
protected:

	int m_nLength; 
	double m_dWeight;
	std::vector<BaseVertex*> m_vtVertexList;
public:
	BasePath(const std::vector<BaseVertex*>& vertex_list, double weight)
		:m_dWeight(weight)
	{
		m_vtVertexList.assign(vertex_list.begin(), vertex_list.end());
		m_nLength = m_vtVertexList.size();
		//m_dWeight=m_nLength-1;
		//std::cout<<"BasePath constructed, size:"<<vertex_list.size()<<std::endl;
	}
	
	~BasePath(void){}

	double Weight() const { return m_dWeight; }
	void Weight(double val) { m_dWeight = val; }

	int length() { return m_nLength; }

	BaseVertex* GetVertex(int i)
	{
		return m_vtVertexList.at(i);
	}

	bool SubPath(std::vector<BaseVertex*>& sub_path, BaseVertex* ending_vertex_pt)
	{
		for (std::vector<BaseVertex*>::const_iterator pos = m_vtVertexList.begin(); 
			pos != m_vtVertexList.end(); ++pos)
		{
			if (*pos != ending_vertex_pt)
			{
				sub_path.push_back(*pos);
			}else
			{
				//break;
				return true;
			}
		}

		return false;
	}

	// display the content
	int get_node(int lambda_limit){
	  //std::cout<<"hi get_node"<<std::flush<<std::endl;
	  auto temp_node=m_vtVertexList.begin()+ m_vtVertexList.size()-lambda_limit-1;
	  return (*temp_node)->getID();
	}
	void PrintOut(std::ostream& out_stream) const
	{
		//out_stream << "Cost: " <<std::flush<< m_dWeight <<std::flush<< ", Length: " << std::flush<<m_vtVertexList.size() << std::endl;
		
		int size = m_vtVertexList.size(); 
		int i = 0;
		for(std::vector<BaseVertex*>::const_iterator pos=m_vtVertexList.begin(); 
			pos!=m_vtVertexList.end();
			++pos)
		{
			i++;
			(*pos)->PrintOut(out_stream);

			if ( i < size ) out_stream << "-";
		}
		out_stream << ",Cost:," << m_dWeight << ",Length:," << m_vtVertexList.size() << std::endl;
	}
	//Santiago: we pass each paths to PathMatrix corresponding destination set
	void add_paths_set(std::shared_ptr<std::set<std::vector<int>,vect_comp_int > > current_paths){
	  std::vector<int> path;
	  for(std::vector<BaseVertex*>::const_iterator pos=m_vtVertexList.begin(); pos!=m_vtVertexList.end();++pos)
	  {
	    path.push_back((*pos)->getID());
	  }
	  current_paths->insert(path);
	  //std::cout<<"current_paths size:"<<current_paths->size()<<std::endl;
	}
	void append_to_path(std::shared_ptr<BasePath> tail){
	  //for(std::vector<BaseVertex*>::const_iterator pos=tail->m_vtVertexList.begin();  pos!=tail->m_vtVertexList.end();pos++)
	  for(auto pos=tail->m_vtVertexList.begin()+1;  pos!=tail->m_vtVertexList.end();++pos){
	    m_vtVertexList.push_back(*pos);
	  }
	  m_nLength=m_vtVertexList.size();
	  m_dWeight+=tail->Weight();
	}
	std::string get_path_string() const
	{
	  std::string output_string;
		//out_stream << "Cost: " <<std::flush<< m_dWeight <<std::flush<< ", Length: " << std::flush<<m_vtVertexList.size() << std::endl;
		
		int size = m_vtVertexList.size(); 
		int i = 0;
		for(std::vector<BaseVertex*>::const_iterator pos=m_vtVertexList.begin(); 
			pos!=m_vtVertexList.end();
			++pos)
		{
			i++;
			output_string+=",";
			output_string+=std::to_string((*pos)->getID());
			//(*pos)->PrintOut(out_stream);

			//if ( i < size ) output_string.length() << "-";
		}
		//out_stream << ",Cost:," << m_dWeight << ",Length:," << m_vtVertexList.size() << std::endl;
		return output_string;
	}


	bool operator<(const BasePath& other) const {
	  if(m_nLength!=other.m_nLength){
	    return m_nLength<other.m_nLength;
	  }
	  else{
	    for (size_t i=0;i<m_nLength;i++){
	      if((m_vtVertexList[i]->getID())==(other.m_vtVertexList[i]->getID())){
		continue;
	      }
	      else return (m_vtVertexList[i]->getID()<other.m_vtVertexList[i]->getID());
	    }
	  }
	  return false;
	}

};
