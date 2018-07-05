#include "weighted_tree.h"
#include "geometry.h"

#define printf

namespace Geex {

WeightedTree::WeightedTree(int v, double area) {
	root_ = new TreeNode() ;
	root_->v_index = v ;
	root_->v_weight = area ;
	root_->total_weight = area ;
	table_[v] = root_ ;
}

WeightedTree::~WeightedTree() {
	if(root_!=NULL)
		delete root_ ;
}

double rand48() {
	double d = Numeric::random_float64() ;
	while (d==0.0) {
		d = Numeric::random_float64() ;
	} ;
	return d ;
}
static int select_vertex_from_tree_w(TreeNode *tn, double w)
{
  //treenode *tn = TREENODE(n);
  //GNode *left_child = g_node_first_child(n);
	TreeNode *left_child = tn->lchild ;

  printf("selecting from tree with total weight %f\n", tn->total_weight);

//  gx_assert(w > 0.0) ;

  if ((w <= tn->v_weight) || (left_child == NULL)) {
	  if(tn->v_weight == 0) {
		  std::cerr << "tree node " << tn->v_index 
			        << " total weight " << tn->total_weight 
					<< " w "  << w << std::endl ;
	  }
//	  gx_assert(tn->v_weight > 0) ;
    return (tn->v_index);
  }

  w -= tn->v_weight;

  if (w <= (left_child)->total_weight)
    return (select_vertex_from_tree_w(left_child,
                                     (left_child)->total_weight * rand48()));

  if (tn->rchild!=NULL)
    return (select_vertex_from_tree_w(tn->rchild,
                                     (tn->rchild)->total_weight *  rand48()));

  return (select_vertex_from_tree_w(left_child, 
                                     (left_child)->total_weight *  rand48()));
}

int WeightedTree::select_vertex_from_tree() {
	int v = select_vertex_from_tree_w(root_, root_->total_weight*rand48());
	if(table_[v]->v_weight == 0) {
		std::cerr << "vertex " << v << " with 0 area is selected...total area is " 
			      << root_->total_weight << std::endl ;
	}
	while(table_[v]->v_weight == 0 && root_->total_weight!=0) {
		v = select_vertex_from_tree_w(root_, root_->total_weight*rand48());
	}
	return v ;
}

static double weight_of_node(TreeNode *n)
{
  double w = (n)->v_weight;
  if (n->lchild != NULL) {
	  w += (n)->lchild->total_weight;
    if ((n)->rchild != NULL) 
      w += (n)->rchild->total_weight;
  }

  return (w);
}

void add_to_tree_t(TreeNode *t, TreeNode *new_node)
{
	if ((t)->lchild == NULL) {
    // make this the left (first) child of this node
		t->lchild = new_node;
		new_node->parent = t ;
	} else if (t->rchild == NULL) {
    // make this the right (last) child of this node
		t->rchild = new_node;
		new_node->parent = t ;
	} else {
		if (Numeric::random_float64() > 0.5) {
			add_to_tree_t(t->lchild, new_node);
		} else {
			add_to_tree_t(t->rchild, new_node);
	    }
  }

  (t)->total_weight = weight_of_node(t);
}

void WeightedTree::add_to_tree(int v, double area) {
//treenode *tn = g_malloc(sizeof(treenode));
//GNode *vn;
	TreeNode *tn = new TreeNode() ;
	tn->v_index = v;
	tn->v_weight = area ;
	tn->total_weight = area ;

	//vn = g_node_new(tn);
	add_to_tree_t(root_, tn);

	//g_hash_table_insert(WT->table, v, vn);
	table_[v] = tn ;
}

void WeightedTree::update_tree(int v, double new_area){
//  GNode *n = g_hash_table_lookup(WT->table, v);
	TreeNode *n = table_[v] ;
  
  // boundary vertices
  if (n == NULL) return;

  (n)->v_weight = new_area;
  
  while (n) {
    (n)->total_weight = weight_of_node(n);
    n = n->parent;
  }
}

#define IND for (i = 0; i < indent; i++) printf(" ");

void tree_area(TreeNode *n, int indent)
{
  int i;
  IND; printf("total  area: %lf\n", (n)->total_weight);
  IND; printf("vertex area: %lf\n", (n)->v_weight);
  if (n->lchild) {
    IND; printf("left:");
	tree_area((n)->lchild, indent + 2);
	if (n->rchild) {
      IND; printf("right:");
	  tree_area(n->rchild, indent + 2);
    }
  }
}

int WeightedTree::in_tree(int v){
	return table_.find(v) != table_.end() ;
//	return (g_hash_table_lookup(WT->table, v) != NULL);
}

void WeightedTree::print_verts(){

}

int WeightedTree::no_area_left(){
	return (root_->total_weight <= 0.0);
}

void WeightedTree::check_verts(double exclude_radius){
}

double WeightedTree::area_left(){
	return  root_->total_weight;
}

} //