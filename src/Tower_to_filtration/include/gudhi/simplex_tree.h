/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2018  TU Graz (Austria)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SIMPLEX_TREE_H
#define SIMPLEX_TREE_H

#include <vector>
#include <unordered_map>
#include <queue>
#include <stack>
#include <algorithm>

namespace Gudhi {
namespace tower_to_filtration {

class Simplex_tree
{
public:
    class Node;
    using vertex = long long;
    using index = long long;
    using simplex_handle = Node*;
    using size_type = long long;
    using simplex_vertex_range = std::vector<vertex>;
    using label_dictionary = std::unordered_map<vertex, Node*>;

    Simplex_tree();
    ~Simplex_tree();

    class Node
    {
    public:
        Node(vertex label, index insertionIndex, int dim, Node *parent);
        ~Node();

        vertex get_label() const;
        index get_insertion_index() const;

        Node* get_parent() const;
        void set_parent(Node *value);
        label_dictionary *get_children() const;

        Node* get_next() const;
        void set_next(Node *value);
        Node* get_prev() const;
        void set_prev(Node *value);

        int get_dim() const;

    private:
        vertex label_;
        index insertionIndex_;
        int dim_;
        Node *parent_;
        label_dictionary *children_;
        Node *next_;
        Node *prev_;
    };

    /* Important for tower_converter.h */

    simplex_handle insert_simplex(simplex_vertex_range &simplex, simplex_handle *handle = nullptr);
    bool insert_simplex_and_faces(simplex_vertex_range &simplex, std::vector<simplex_handle> *addedSimplices = nullptr);
    bool insert_edge_and_expand(vertex u, vertex v, int maxDim = -1, std::vector<simplex_handle> *addedSimplices = nullptr);	// u < v !!! ; maxDim == -1 -> no limit
    void remove_simplex(simplex_vertex_range &simplex, std::vector<simplex_handle> *removedIndices = nullptr);
    void remove_simplex(simplex_handle &simplex, std::vector<simplex_handle> *removedIndices = nullptr);
    simplex_handle get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_handle> *closedStar);    //returns vertex with smallest closed star; closedStar is ordered
    void get_boundary(simplex_handle &simplex, std::vector<simplex_handle> *boundary);   //ordered by increasing insertion numbers
    simplex_vertex_range get_vertices(simplex_handle &simplex);

    /* Other */

    size_type get_size() const;
    int get_dimension(simplex_handle handle);
    simplex_handle get_simplex_handle(simplex_vertex_range &simplex);
    size_type get_max_size() const;
    int get_max_dimension() const;
    void get_cofaces(simplex_vertex_range &simplex, std::vector<simplex_vertex_range*> *cofaces);
    void get_cofacets(simplex_vertex_range &simplex, std::vector<simplex_vertex_range*> *cofacets);
    bool contains(simplex_vertex_range &simplex);
    void print();

private:
    std::vector<label_dictionary*> *dictionaries_; // circular list for each label at each height
    std::unordered_map<vertex, vertex> *verticesTranslation_;
    size_type numberOfSimplices_;
    index maxIndex_;
    size_type maxSize_;
    int maxDim_;

    Node* insert_simplex_in_tree(simplex_vertex_range &simplex);
    Node* insert_union(Node *simplex, vertex v);
    void insert_node_in_dictionary(Node *simplex);
    void delete_node_from_dictionary(Node *simplex);
    void expand_node(Node *simplex, int maxDim, std::vector<simplex_handle> *addedSimplices);
    Node* find(simplex_vertex_range &simplex);
    void get_cofaces(Node *simplex, std::vector<Node*> *cofaces);
    void get_cofacets(Node *simplex, std::vector<Node*> *cofacets);
    bool is_coface(Node *node, Node *simplex);
    bool are_all_facets_of_union_inserted(Node *facet, vertex endLabel);
    void get_simplices_in_subtree(Node *node, std::vector<Node*> *simplices);
    vertex get_smallest_star(vertex v, vertex u, std::queue<Node*> *qv, std::queue<Node*> *qu);
    bool get_smallest_star_find_next_node(vertex v, vertex toAvoid, std::vector<label_dictionary*>::size_type &currentHeight,
                                          Node *&currentRoot, Node *&startRoot, Node *&currentNode, bool &currentNodeIsRoot,
                                          label_dictionary *&children, label_dictionary::iterator &it, std::queue<Node*> *&tail);
    bool get_smallest_star_find_next_root(vertex v, vertex toAvoid, std::vector<label_dictionary*>::size_type &currentHeight,
                                          Node *&currentRoot, Node *&startRoot, Node *&currentNode, bool &currentNodeIsRoot);
    Node* get_opposite_facet(Node *simplex, vertex v);
    void print_from_node(Node *node, std::string prefix);
};

inline Simplex_tree::Simplex_tree() : numberOfSimplices_(0), maxIndex_(-1), maxSize_(0), maxDim_(0)
{
    verticesTranslation_ = new std::unordered_map<vertex, vertex>();
    dictionaries_ = new std::vector<label_dictionary*>();
}

inline Simplex_tree::~Simplex_tree()
{
    delete verticesTranslation_;
    for (std::vector<label_dictionary*>::size_type i = 0; i < dictionaries_->size(); i++){
        for (label_dictionary::iterator it = dictionaries_->at(i)->begin(); it != dictionaries_->at(i)->end(); it++) delete it->second;
        delete dictionaries_->at(i);
    }
    delete dictionaries_;
}

inline Simplex_tree::simplex_handle Simplex_tree::insert_simplex(simplex_vertex_range &simplex, simplex_handle *handle)
{
    Node *simplexNode = insert_simplex_in_tree(simplex);
    if (simplexNode == nullptr) return simplexNode;

    if (handle != nullptr) *handle = simplexNode;

    int dim = simplex.size() - 1;

    numberOfSimplices_++;
    maxIndex_++;
    if (maxSize_ < numberOfSimplices_) maxSize_ = numberOfSimplices_;
    if (maxDim_ < dim) maxDim_ = dim;

    return simplexNode;
}

inline bool Simplex_tree::insert_simplex_and_faces(simplex_vertex_range &simplex, std::vector<simplex_handle> *addedSimplices)
{
    int dim;
    std::vector<Node*> nodes;
    std::vector<Node*> tmpNodes;
    Node *vertexNode;

    if (find(simplex) != nullptr) return false;

    if (numberOfSimplices_ == 0) dictionaries_->push_back(new label_dictionary());
    label_dictionary *vertices = dictionaries_->front();

    for (simplex_vertex_range::size_type i = 0; i < simplex.size(); ++i){
	vertex v = simplex.at(i);
	tmpNodes.clear();

	if (vertices->find(v) != vertices->end()){
	    vertexNode = vertices->at(v);
	} else {
	    vertexNode = new Node(v, maxIndex_ + 1, 0, nullptr);
	    vertices->emplace(v, vertexNode);
	    ++numberOfSimplices_;
	    ++maxIndex_;
	    if (addedSimplices != nullptr) addedSimplices->push_back(vertexNode);
	}
	tmpNodes.push_back(vertexNode);

	for (Node *node : nodes){
	    label_dictionary* children = node->get_children();
	    dim = node->get_dim();
	    if (children->find(v) == children->end()){
		if ((int)dictionaries_->size() == dim + 1) dictionaries_->push_back(new label_dictionary());
		Node *vNode = new Node(v, maxIndex_ + 1, dim + 1, node);
		children->emplace(v, vNode);
		insert_node_in_dictionary(vNode);
		++numberOfSimplices_;
		++maxIndex_;
		if (addedSimplices != nullptr) addedSimplices->push_back(vNode);
	    }
	    tmpNodes.push_back(children->at(v));
	}
	for (Node *node : tmpNodes){
	    nodes.push_back(node);
	}
    }

    dim = simplex.size() - 1;
    if (maxSize_ < numberOfSimplices_) maxSize_ = numberOfSimplices_;
    if (maxDim_ < dim) maxDim_ = dim;

    return true;
}

inline bool Simplex_tree::insert_edge_and_expand(vertex u, vertex v, int maxDim, std::vector<simplex_handle> *addedSimplices)
{
    simplex_vertex_range vect;
    vect.push_back(u);
    Node *vertexNode = insert_simplex(vect);
    if (vertexNode && addedSimplices != nullptr) addedSimplices->push_back(vertexNode);
    vect.at(0) = v;
    vertexNode = insert_simplex(vect);
    if (vertexNode && addedSimplices != nullptr) addedSimplices->push_back(vertexNode);
    vect.at(0) = u;
    vect.push_back(v);
    Node *edgeNode = insert_simplex(vect);

    if (edgeNode) {
	if (addedSimplices != nullptr) addedSimplices->push_back(edgeNode);
	if (maxDim > 1 || maxDim == -1) expand_node(edgeNode, maxDim, addedSimplices);
	return true;
    }

    return false;
}

inline void Simplex_tree::remove_simplex(simplex_vertex_range &simplex, std::vector<simplex_handle> *removedIndices)
{
    Node *simplexNode = find(simplex);
    remove_simplex(simplexNode, removedIndices);
}

inline void Simplex_tree::remove_simplex(simplex_handle &simplex, std::vector<simplex_handle> *removedIndices)
{
    std::vector<Node*> cofaces;
    get_cofaces(simplex, &cofaces);

    std::sort(cofaces.begin(), cofaces.end(), [](Node *n1, Node *n2){return n1->get_dim() > n2->get_dim();});

    for (std::vector<Node*>::size_type i = 0; i < cofaces.size(); i++){
	Node *toDelete = cofaces.at(i);
	if (removedIndices != nullptr) removedIndices->push_back(toDelete);

	if (toDelete->get_dim() > 0) toDelete->get_parent()->get_children()->erase(toDelete->get_label());
	delete_node_from_dictionary(toDelete);

	numberOfSimplices_--;

	delete toDelete;
    }
}

inline Simplex_tree::simplex_handle Simplex_tree::get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_handle> *closedStar)
{
    std::queue<Node*> qv;
    std::queue<Node*> qu;
    Node *s;

    if (get_smallest_star(v, u, &qv, &qu) == v){
        while (!qv.empty()){
            s = qv.front();
            qv.pop();
	    if (s->get_parent() != nullptr){
		closedStar->push_back(get_opposite_facet(s, v));
            }
	    closedStar->push_back(s);
        }

	return dictionaries_->front()->at(v);
    } else {
        while (!qu.empty()){
            s = qu.front();
            qu.pop();
	    if (s->get_parent() != nullptr){
		closedStar->push_back(get_opposite_facet(s, u));
            }
	    closedStar->push_back(s);
        }

	return dictionaries_->front()->at(u);
    }
}

inline Simplex_tree::simplex_handle Simplex_tree::get_simplex_handle(simplex_vertex_range &simplex)
{
    Node *simplexNode = find(simplex);
    return simplexNode;
}

inline void Simplex_tree::get_boundary(simplex_handle &simplexNode, std::vector<simplex_handle> *boundary)
{
    if (simplexNode->get_parent() == nullptr) return;

    simplex_vertex_range tail;
    tail.push_back(simplexNode->get_label());
    Node *it = simplexNode->get_parent();
    simplex_vertex_range::reverse_iterator tailIt;
    Node *itb;

    while (it != nullptr){
	itb = it;
	tailIt = tail.rbegin();
	tailIt++;
	while (tailIt != tail.rend()){
	    itb = itb->get_children()->at(*tailIt);
	    tailIt++;
	}
	boundary->push_back(itb);
	tail.push_back(it->get_label());
	it = it->get_parent();
    }
    itb = dictionaries_->front()->at(tail.at(tail.size() - 2));
    tailIt = tail.rbegin();
    tailIt += 2;
    while (tailIt != tail.rend()){
	itb = itb->get_children()->at(*tailIt);
	tailIt++;
    }
    boundary->push_back(itb);

    std::sort(boundary->begin(), boundary->end(), [](const Node *n1, const Node *n2) {
	return n1->get_insertion_index() < n2->get_insertion_index();
    });
}

inline Simplex_tree::simplex_vertex_range Simplex_tree::get_vertices(simplex_handle &simplex)
{
    simplex_vertex_range vector;
    Node *nodeIt = simplex;
    while (nodeIt != nullptr){
	vector.push_back(nodeIt->get_label());
	nodeIt = nodeIt->get_parent();
    }
    std::reverse(vector.begin(), vector.end());
    return vector;
}

inline Simplex_tree::size_type Simplex_tree::get_size() const
{
    return numberOfSimplices_;
}

inline int Simplex_tree::get_dimension(simplex_handle handle)
{
    return handle->get_dim();
}

inline Simplex_tree::size_type Simplex_tree::get_max_size() const
{
    return maxSize_;
}

inline int Simplex_tree::get_max_dimension() const
{
    return maxDim_;
}

inline void Simplex_tree::get_cofaces(simplex_vertex_range &simplex, std::vector<simplex_vertex_range*> *cofaces)
{
    /*Node *simplexNode = find(simplex);
    std::vector<Node*> cofaceNodes;
    get_cofaces(simplexNode, &cofaceNodes);

    for (std::vector<Node*>::size_type i = 0; i < cofaceNodes.size(); i++){
	cofaces->push_back(node_to_vector(cofaceNodes.at(i)));
    }*/
}

inline void Simplex_tree::get_cofacets(simplex_vertex_range &simplex, std::vector<simplex_vertex_range*> *cofacets)
{
    /*Node *simplexNode = find(simplex);
    std::vector<Node*> cofacetNodes;
    get_cofacets(simplexNode, &cofacetNodes);

    for (std::vector<Node*>::size_type i = 0; i < cofacetNodes.size(); i++){
        cofacets->push_back(node_to_vector(cofacetNodes.at(i)));
    }*/
}

inline bool Simplex_tree::contains(simplex_vertex_range &simplex)
{
    return find(simplex) != nullptr;
}

inline void Simplex_tree::print()
{
    Node *vertex;
    for (auto it = dictionaries_->front()->begin(); it != dictionaries_->front()->end(); ++it){
	vertex = it->second;
	print_from_node(vertex, "");
    }
}

inline Simplex_tree::Node* Simplex_tree::insert_simplex_in_tree(simplex_vertex_range &simplex)
{
    int dim = simplex.size() - 1;

    if (dim == 0){
        if (numberOfSimplices_ == 0) dictionaries_->push_back(new label_dictionary());
        label_dictionary *vertices = dictionaries_->front();
	if (vertices->find(simplex.at(0)) != vertices->end()) return nullptr;
	Node *simplexNode = new Node(simplex.at(0), maxIndex_ + 1, 0, nullptr);
	vertices->emplace(simplex.at(0), simplexNode);
        return simplexNode;
    }

    Node *it = dictionaries_->front()->at(simplex.at(0));
    for (int i = 1; i < dim; i++){
	it = it->get_children()->at(simplex.at(i));
    }

    if (it->get_children()->find(simplex.back()) != it->get_children()->end()) return nullptr;

    if ((int)dictionaries_->size() == dim) dictionaries_->push_back(new label_dictionary());
    Node *simplexNode = new Node(simplex.back(), maxIndex_ + 1, dim, it);
    it->get_children()->emplace(simplex.back(), simplexNode);
    insert_node_in_dictionary(simplexNode);

    return simplexNode;
}

inline Simplex_tree::Node* Simplex_tree::insert_union(Node *simplex, vertex v)
{
    Node *it = simplex;
    std::stack<vertex> tail;
    vertex endLabel = v;

    while (it != nullptr && it->get_label() > v){
	tail.push(it->get_label());
	it = it->get_parent();
    }
    if (it != nullptr && it->get_label() == v) return nullptr;
    if (it == nullptr){
	it = dictionaries_->front()->at(v);
    } else if (it->get_children()->find(v) != it->get_children()->end()){
	it = it->get_children()->at(v);
    }

    if (it != simplex){
	while (!tail.empty() && it->get_children()->find(tail.top()) != it->get_children()->end()){
	    it = it->get_children()->at(tail.top());
	    tail.pop();
	}
	if (tail.empty() || tail.size() != 1) return nullptr;
	endLabel = tail.top();
    }

    if (!are_all_facets_of_union_inserted(it, endLabel)) return nullptr;

    if ((int)dictionaries_->size() == simplex->get_dim() + 1) dictionaries_->push_back(new label_dictionary());
    Node *simplexNode = new Node(endLabel, maxIndex_ + 1, simplex->get_dim() + 1, it);
    it->get_children()->emplace(endLabel, simplexNode);
    insert_node_in_dictionary(simplexNode);

    numberOfSimplices_++;
    maxIndex_++;
    if (maxSize_ < numberOfSimplices_) maxSize_ = numberOfSimplices_;
    if (maxDim_ < simplex->get_dim() + 1) maxDim_ = simplex->get_dim() + 1;

    return simplexNode;
}

inline void Simplex_tree::insert_node_in_dictionary(Node *simplex)
{
    label_dictionary *dict = dictionaries_->at(simplex->get_dim());

    if (dict->find(simplex->get_label()) == dict->end()) {
        dict->emplace(simplex->get_label(), simplex);
        simplex->set_next(simplex);
        simplex->set_prev(simplex);
    } else {
        Node *head = dict->at(simplex->get_label());
        simplex->set_next(head->get_next());
        simplex->set_prev(head);
        head->get_next()->set_prev(simplex);
        head->set_next(simplex);
    }
}

inline void Simplex_tree::delete_node_from_dictionary(Node *simplex)
{
    label_dictionary *dict = dictionaries_->at(simplex->get_dim());

    if (simplex->get_next() == simplex) dict->erase(simplex->get_label());
    else {
        simplex->get_next()->set_prev(simplex->get_prev());
        simplex->get_prev()->set_next(simplex->get_next());
        if (dict->at(simplex->get_label()) == simplex) dict->at(simplex->get_label()) = simplex->get_next();
    }
}

inline void Simplex_tree::expand_node(Node *simplex, int maxDim, std::vector<simplex_handle> *addedSimplices)
{
    std::vector<Node*> cofacets;
    get_cofacets(simplex->get_parent(), &cofacets);

    for (Node *cofacet : cofacets){
	if (cofacet == simplex) continue;
	Node *node = insert_union(cofacet, simplex->get_label());
	if (node != nullptr){	//if node could be inserted, i.e. all its facets were there and it-self was not already inserted
	    if (addedSimplices != nullptr) addedSimplices->push_back(node);
	    if (node->get_dim() < maxDim || maxDim == -1) expand_node(node, maxDim, addedSimplices);
	}
    }
}

inline Simplex_tree::Node *Simplex_tree::find(simplex_vertex_range &simplex)
{
    if (simplex.empty() || dictionaries_->empty()) return nullptr;

    vertex label = simplex.front();
    label_dictionary *vertices = dictionaries_->front();

    if (vertices->find(label) == vertices->end()) return nullptr;

    Node *it = vertices->at(label);
    for (simplex_vertex_range::size_type i = 1; i < simplex.size(); ++i) {
	label = simplex.at(i);
        label_dictionary::iterator next = it->get_children()->find(label);

	if (next == it->get_children()->end()) return nullptr;
        else it = next->second;
    }
    return it;
}

inline void Simplex_tree::get_cofaces(Node *simplex, std::vector<Node *> *cofaces)
{
    for (int d = dictionaries_->size() - 1; d >= simplex->get_dim(); d--){
        if (dictionaries_->at(d)->find(simplex->get_label()) != dictionaries_->at(d)->end()){
            Node *it = dictionaries_->at(d)->at(simplex->get_label());
            Node *startingNode = it;
            if (is_coface(it, simplex)) get_simplices_in_subtree(it, cofaces);
            it = it->get_next();
            while (it != startingNode){
                if (is_coface(it, simplex)) get_simplices_in_subtree(it, cofaces);
                it = it->get_next();
            }
        }
    }
}

inline void Simplex_tree::get_cofacets(Node *simplex, std::vector<Node*> *cofacets)
{
    if ((vertex)dictionaries_->size() > simplex->get_dim() &&
	    dictionaries_->at(simplex->get_dim() + 1)->find(simplex->get_label()) != dictionaries_->at(simplex->get_dim() + 1)->end()){
	Node *it = dictionaries_->at(simplex->get_dim() + 1)->at(simplex->get_label());
	Node *startingNode = it;
	if (is_coface(it, simplex)) cofacets->push_back(it);
	it = it->get_next();
	while (it != startingNode){
	    if (is_coface(it, simplex)) cofacets->push_back(it);
	    it = it->get_next();
	}
    }

    for (auto childIt = simplex->get_children()->begin(); childIt != simplex->get_children()->end(); childIt++){
        cofacets->push_back(childIt->second);
    }
}

inline bool Simplex_tree::is_coface(Node *node, Node *simplex)
{
    Node *nodeIt = node;
    Node *simplexIt = simplex;

    while (simplexIt != nullptr && nodeIt != nullptr){
        if (simplexIt->get_label() == nodeIt->get_label()) simplexIt = simplexIt->get_parent();
        nodeIt = nodeIt->get_parent();
    }

    if (simplexIt == nullptr) return true;
    else return false;
}

inline bool Simplex_tree::are_all_facets_of_union_inserted(Node *facet, vertex endLabel)
{
    if (facet->get_parent() == nullptr) return dictionaries_->front()->find(endLabel) != dictionaries_->front()->end();

    simplex_vertex_range tail;
    tail.push_back(endLabel);
    Node *it = facet;
    simplex_vertex_range::reverse_iterator tailIt;
    Node *itb;

    while (it != nullptr){
	itb = it;
	tailIt = tail.rbegin();
	tailIt++;
	while (tailIt != tail.rend()){
	    if (itb->get_children()->find(*tailIt) == itb->get_children()->end()) return false;
	    itb = itb->get_children()->at(*tailIt);
	    tailIt++;
	}
	tail.push_back(it->get_label());
	it = it->get_parent();
    }
    if (dictionaries_->front()->find(tail.at(tail.size() - 2)) == dictionaries_->front()->end()) return false;
    itb = dictionaries_->front()->at(tail.at(tail.size() - 2));
    tailIt = tail.rbegin();
    tailIt += 2;
    while (tailIt != tail.rend()){
	if (itb->get_children()->find(*tailIt) == itb->get_children()->end()) return false;
	itb = itb->get_children()->at(*tailIt);
	tailIt++;
    }

    return true;
}

inline void Simplex_tree::get_simplices_in_subtree(Node *node, std::vector<Node *> *simplices)
{
    simplices->push_back(node);
    for (auto it = node->get_children()->begin(); it != node->get_children()->end(); it++){
        get_simplices_in_subtree(it->second, simplices);
    }
}

inline Simplex_tree::vertex Simplex_tree::get_smallest_star(vertex v, vertex u, std::queue<Node*> *qv, std::queue<Node*> *qu)
{
    label_dictionary *childrenV;
    label_dictionary *childrenU;
    label_dictionary::iterator vit;
    label_dictionary::iterator uit;
    std::queue<Node*> *tv = new std::queue<Node*>();
    std::queue<Node*> *tu = new std::queue<Node*>();
    std::vector<label_dictionary*>::size_type currentHeightV = 0;
    std::vector<label_dictionary*>::size_type currentHeightU = 0;
    Node *startRootV = dictionaries_->at(0)->at(v);
    Node *startRootU = dictionaries_->at(0)->at(u);
    Node *currentRootV = startRootV;
    Node *currentRootU = startRootU;
    Node *currentNodeV = startRootV;
    Node *currentNodeU = startRootU;
    bool currentNodeIsRootV = true;
    bool currentNodeIsRootU = true;
    bool continueV = true;
    bool continueU = true;

    while (continueV && continueU){
        qv->push(currentNodeV);
        if (!currentNodeIsRootV) tv->push(currentNodeV);
        qu->push(currentNodeU);
        if (!currentNodeIsRootU) tu->push(currentNodeU);

        continueV = get_smallest_star_find_next_node(v, u, currentHeightV, currentRootV, startRootV, currentNodeV, currentNodeIsRootV, childrenV, vit, tv);
        continueU = get_smallest_star_find_next_node(u, v, currentHeightU, currentRootU, startRootU, currentNodeU, currentNodeIsRootU, childrenU, uit, tu);
    }

    delete tv;
    delete tu;

    if (!continueU) return u;
    else return v;
}

inline bool Simplex_tree::get_smallest_star_find_next_node(vertex v, vertex toAvoid, std::vector<label_dictionary*>::size_type &currentHeight, Node *&currentRoot,
                                                    Node *&startRoot, Node *&currentNode, bool &currentNodeIsRoot,
                                                    label_dictionary *&children, label_dictionary::iterator &it, std::queue<Node*> *&tail)
{
    if (currentNodeIsRoot){
        children = currentRoot->get_children();
        it = children->begin();
        currentNodeIsRoot = false;
    } else {
        it++;
    }

    if (it != children->end() && it->first == toAvoid) it++;
    while (it == children->end()){
        if (!tail->empty()){
            children = tail->front()->get_children();
            it = children->begin();
            tail->pop();
        } else {
            return get_smallest_star_find_next_root(v, toAvoid, currentHeight, currentRoot, startRoot, currentNode, currentNodeIsRoot);
        }
    }

    if (it->first == toAvoid){
        return get_smallest_star_find_next_node(v, toAvoid, currentHeight, currentRoot, startRoot, currentNode, currentNodeIsRoot, children, it, tail);
    }

    currentNode = it->second;
    return true;
}

inline bool Simplex_tree::get_smallest_star_find_next_root(vertex v, vertex toAvoid, std::vector<label_dictionary*>::size_type &currentHeight, Node *&currentRoot,
							   Node *&startRoot, Node *&currentNode, bool &currentNodeIsRoot)
{
    if (currentRoot->get_next() == startRoot){
        currentHeight++;
        while (currentHeight < dictionaries_->size() && dictionaries_->at(currentHeight)->find(v) == dictionaries_->at(currentHeight)->end()) currentHeight++;
        if (currentHeight == dictionaries_->size()) return false;
        startRoot = dictionaries_->at(currentHeight)->at(v);
        currentRoot = startRoot;
    } else {
        currentRoot = currentRoot->get_next();
    }

    if (is_coface(currentRoot, dictionaries_->at(0)->at(toAvoid)))
        return get_smallest_star_find_next_root(v, toAvoid, currentHeight, currentRoot, startRoot, currentNode, currentNodeIsRoot);
    currentNode = currentRoot;
    currentNodeIsRoot = true;
    return true;
}

inline Simplex_tree::Node *Simplex_tree::get_opposite_facet(Node *simplex, vertex v)
{
    if (simplex->get_parent() == nullptr) return nullptr;

    Node *trav = simplex;
    std::stack<vertex> tail;

    while (trav != nullptr && trav->get_label() > v) {
        tail.push(trav->get_label());
        trav = trav->get_parent();
    }

    if (trav == nullptr || trav->get_label() != v) return nullptr;

    if (trav->get_parent() != nullptr){
        trav = trav->get_parent();
    } else {
        trav = dictionaries_->at(0)->at(tail.top());
        tail.pop();
    }
    while (!tail.empty()){
        trav = trav->get_children()->at(tail.top());
        tail.pop();
    }

    return trav;
}

inline void Simplex_tree::print_from_node(Node *node, std::string prefix)
{
    prefix += std::to_string(node->get_label());

    std::cout << prefix << "\n";

    for (auto it = node->get_children()->begin(); it != node->get_children()->end(); ++it){
	print_from_node(it->second, prefix);
    }
}

inline Simplex_tree::Node::Node(vertex label, index insertionIndex, int dim, Node *parent) :
    label_(label), insertionIndex_(insertionIndex), dim_(dim), parent_(parent), next_(this), prev_(this)
{
    children_ = new label_dictionary();
}

inline Simplex_tree::Node::~Node()
{
    delete children_;
}

inline Simplex_tree::vertex Simplex_tree::Node::get_label() const
{
    return label_;
}

inline Simplex_tree::index Simplex_tree::Node::get_insertion_index() const
{
    return insertionIndex_;
}

inline Simplex_tree::Node *Simplex_tree::Node::get_parent() const
{
    return parent_;
}

inline void Simplex_tree::Node::set_parent(Node *value)
{
    parent_ = value;
}

inline Simplex_tree::label_dictionary *Simplex_tree::Node::get_children() const
{
    return children_;
}

inline Simplex_tree::Node *Simplex_tree::Node::get_next() const
{
    return next_;
}

inline void Simplex_tree::Node::set_next(Node *value)
{
    next_ = value;
}

inline Simplex_tree::Node *Simplex_tree::Node::get_prev() const
{
    return prev_;
}

inline void Simplex_tree::Node::set_prev(Node *value)
{
    prev_ = value;
}

inline int Simplex_tree::Node::get_dim() const
{
    return dim_;
}

}
}

#endif // SIMPLEX_TREE_H
