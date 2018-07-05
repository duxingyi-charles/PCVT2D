#ifndef __WEIGHTED_TREE_H__
#define __WEIGHTED_TREE_H__

#include <map>

namespace Geex {

class TreeNode {
public:
	TreeNode() {
		v_index = -1 ;
		v_weight = 0 ;
		total_weight = 0 ;
		lchild = NULL ;
		rchild = NULL ;
		parent = NULL ;
	}
	~TreeNode() {
		if(lchild!=NULL) delete lchild ;
		if(rchild!=NULL) delete rchild ;
	}

	int    v_index ; // corresponding vertex index in the DT
	double v_weight ; 
	double total_weight ;
	TreeNode *lchild ;
	TreeNode *rchild ;
	TreeNode *parent ;
} ;

class WeightedTree {
public:
	WeightedTree(int v=-1, double area=0) ;
	~WeightedTree() ;
//gpointer new_tree(GtsVertex *v, double area);
	int select_vertex_from_tree();
//GtsVertex *select_vertex_from_tree(weighted_tree *WT);
	void add_to_tree(int v, double area);
//void add_to_tree(weighted_tree *WT, GtsVertex *v, double area);
	void update_tree(int v, double new_area);
//void update_tree(weighted_tree *WT, GtsVertex *v, double new_area);
	int in_tree(int v);
//int in_tree(weighted_tree *WT, GtsVertex *v);
	void print_verts();
//void print_verts(weighted_tree *WT);
	int no_area_left();
//int no_area_left(weighted_tree *WT);
	void check_verts(double exclude_radius);
//void check_verts(weighted_tree *WT, double exclude_radius);
	double area_left();
//double area_left(weighted_tree *WT);

protected:
	TreeNode *root_ ;
	std::map<int, TreeNode*> table_ ;
} ;

} 

#endif 
