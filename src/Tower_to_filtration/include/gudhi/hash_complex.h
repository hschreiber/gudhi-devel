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

#ifndef COMPLEX_H
#define COMPLEX_H

#include <unordered_map>
#include <vector>
#include <queue>
#include <algorithm>

namespace Gudhi {
namespace tower_to_filtration {

class Hash_complex
{
public:
    using vertex = long long;
    using index = long long;
    using simplex_handle = index;
    using simplex_vertex_range = std::vector<vertex>;
    using size_type = long long;

    Hash_complex();
    ~Hash_complex();

    class Simplex
    {
    public:
	Simplex(index num, simplex_vertex_range *vertices);
	~Simplex();

	index get_insertion_num() const;
	void set_insertion_num(const index &insertionNum);

	void add_cofacet(Simplex *coface, vertex v);
	std::unordered_map<vertex, Simplex*>* get_cofacets() const;

	simplex_vertex_range *get_vertices() const;

    private:
	index insertionNum_;
	std::unordered_map<vertex, Simplex*> *cofacets_;
	simplex_vertex_range* vertices_;
    };

    struct Key_hasher {
	std::size_t operator()(const std::pair<simplex_vertex_range*, int> *k) const
	{
	    std::size_t seed;
	    if (k->second < 0) seed = k->first->size();
	    else seed = k->first->size() - 1;

	    for (int i = 0; i < (int)k->first->size(); i++) {
		if (i != k->second) seed ^= (std::size_t)(k->first->at(i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	    }
	    return seed;
	}
    };

    struct Simplices_equals : std::binary_function<const std::pair<simplex_vertex_range*, int>*, const std::pair<simplex_vertex_range*, int>*, bool> {
	bool operator()(const std::pair<simplex_vertex_range*, int> *s1, const std::pair<simplex_vertex_range*, int> *s2) const
	{
	    const std::pair<simplex_vertex_range*, int> *key;
	    const std::pair<simplex_vertex_range*, int> *inMap;
	    simplex_vertex_range::size_type size;

	    if (s1->second > -1) {
		key = s1;
		inMap = s2;
	    } else {
		key = s2;
		inMap = s1;
	    }

	    if (key->second < 0) size = key->first->size();
	    else size = key->first->size() - 1;
	    if (size != inMap->first->size()) return false;
	    int j = 0;
	    for (simplex_vertex_range::size_type i = 0; i < size; i++){
		if (j == key->second) j++;
		if (key->first->at(j) != inMap->first->at(i)) {
		    return false;
		}
		j++;
	    }
	    return true;
	}
    };

    /* Important for tower_converter.h */

    bool insert_simplex(simplex_vertex_range &simplex, simplex_handle *handle = nullptr);
    bool insert_simplex_and_faces(simplex_vertex_range &simplex, std::vector<simplex_handle> *addedSimplices = nullptr);
    bool insert_edge_and_expand(vertex u, vertex v, int maxDim = -1, std::vector<simplex_handle> *addedSimplices = nullptr);	// u < v !!! ; maxDim == -1 -> no limit
    bool remove_simplex(simplex_vertex_range &simplex, std::vector<simplex_handle> *removedIndices = nullptr);
    bool remove_simplex(simplex_handle &simplex, std::vector<simplex_handle> *removedIndices = nullptr);
    simplex_handle get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_handle> *closedStar);    //returns vertex with smallest closed star; closedStar is ordered
    void get_boundary(simplex_handle &handle, std::vector<simplex_handle> *boundary);	//ordered by increasing insertion numbers
    simplex_vertex_range& get_vertices(simplex_handle &handle) const;

    /* Other */

    simplex_handle get_simplex_handle(simplex_vertex_range &simplex);
    int get_dimension(simplex_handle handle);
    size_type get_size() const;
    size_type get_max_size() const;
    int get_max_dimension() const;
    void get_cofaces(simplex_vertex_range &simplex, std::vector<simplex_vertex_range*> *cofaces);
    void get_cofacets(simplex_vertex_range &simplex, std::vector<simplex_vertex_range *> *cofacets);
    bool contains(simplex_vertex_range &simplex);
    void print();

private:
    index maxIndex_;
    size_type maxSize_;
    int maxDim_;
    std::unordered_map<
	std::pair<simplex_vertex_range*, int>*,
	Simplex*,
	Key_hasher,
	Simplices_equals
    > *simplices_;
    std::unordered_map<simplex_handle, Simplex*> *handleToSimplex_;

    vertex get_smallest_star(vertex v, vertex u, std::queue<Simplex*> *qv, std::queue<Simplex*> *qu);
    int get_vertex_index(simplex_vertex_range *simplex, vertex v);
    void expand_simplex(simplex_vertex_range *vectSimplex, int maxDim, std::vector<simplex_handle> *addedSimplices);
    Simplex* insert_union(Simplex *simplex, vertex v);
    bool remove_simplex(Simplex *simplex, std::vector<simplex_handle> *removedIndices);
};

inline Hash_complex::Hash_complex() : maxIndex_(-1), maxSize_(0), maxDim_(0)
{
    simplices_ = new std::unordered_map<
		    std::pair<simplex_vertex_range*, int>*,
		    Simplex*,
		    Key_hasher,
		    Simplices_equals
		>();
    handleToSimplex_ = new std::unordered_map<simplex_handle, Simplex*>();
}

inline Hash_complex::~Hash_complex()
{
    for (auto it = simplices_->begin(); it != simplices_->end(); ++it){
        delete it->first;
        delete it->second;
    }
    delete simplices_;
    delete handleToSimplex_;
}

inline bool Hash_complex::insert_simplex(simplex_vertex_range &simplex, simplex_handle *handle)
{
    simplex_vertex_range *vs = new simplex_vertex_range(simplex);
    std::pair<simplex_vertex_range*,int> *p = new std::pair<simplex_vertex_range*,int>(vs, -1);
    Simplex *splx = new Simplex(maxIndex_ + 1, vs);

    if (simplices_->emplace(p, splx).second == false) {
	p->first = nullptr;
	delete splx;
	delete p;
	return false;
    }

    handleToSimplex_->emplace(++maxIndex_, splx);
    if ((int)vs->size() - 1 > maxDim_) maxDim_ = vs->size() - 1;
    if (maxSize_ < (size_type)simplices_->size()) maxSize_ = simplices_->size();

    if (handle != nullptr) *handle = maxIndex_;

    if (vs->size() <= 1) return true;

    for (simplex_vertex_range::size_type i = 0; i < vs->size(); i++){
	p->second = i;
	simplices_->at(p)->add_cofacet(splx, vs->at(i));
    }
    p->second = -1;

    return true;
}

inline bool Hash_complex::insert_simplex_and_faces(simplex_vertex_range &simplex, std::vector<simplex_handle> *addedSimplices)
{
    std::pair<simplex_vertex_range*,int> p(&simplex, -1);
    if (simplices_->find(&p) != simplices_->end()) return false;

    std::vector<simplex_vertex_range> simplices;
    std::vector<simplex_vertex_range> tmpSimplices;
    simplex_handle currentIndex;
    simplices.push_back(simplex_vertex_range());

    for (simplex_vertex_range::size_type i = 0; i < simplex.size(); ++i){
	vertex v = simplex.at(i);
	tmpSimplices.clear();
	for (simplex_vertex_range &prefix : simplices){
	    tmpSimplices.push_back(simplex_vertex_range(prefix));
	    tmpSimplices.back().push_back(v);
	    if (insert_simplex(tmpSimplices.back(), &currentIndex) && addedSimplices != nullptr) addedSimplices->push_back(currentIndex);
	}
	for (auto it = tmpSimplices.begin(); it != tmpSimplices.end(); ++it){
	    simplices.push_back(*it);
	}
    }

    return true;
}

inline bool Hash_complex::insert_edge_and_expand(vertex u, vertex v, int maxDim, std::vector<simplex_handle> *addedSimplices)
{
    simplex_vertex_range vect;
    simplex_handle handle;
    std::pair<simplex_vertex_range*,int> p(&vect, -1);
    vect.push_back(u);
    if (insert_simplex(vect, &handle) && addedSimplices != nullptr) addedSimplices->push_back(handle);
    vect.at(0) = v;
    if (insert_simplex(vect, &handle) && addedSimplices != nullptr) addedSimplices->push_back(handle);
    vect.at(0) = u;
    vect.push_back(v);

    if (insert_simplex(vect, &handle)) {
	if (addedSimplices != nullptr) addedSimplices->push_back(handle);
	if (maxDim > 1 || maxDim == -1) expand_simplex(&vect, maxDim, addedSimplices);
	return true;
    }

    return false;
}

inline bool Hash_complex::remove_simplex(simplex_vertex_range &simplex, std::vector<simplex_handle> *removedIndices)
{
    std::pair<simplex_vertex_range*,int> p(&simplex, -1);
    if (simplices_->find(&p) == simplices_->end()) return false;
    return remove_simplex(simplices_->at(&p), removedIndices);
}

inline bool Hash_complex::remove_simplex(simplex_handle &simplex, std::vector<simplex_handle> *removedIndices)
{
    if (handleToSimplex_->find(simplex) == handleToSimplex_->end()) return false;
    return remove_simplex(handleToSimplex_->at(simplex), removedIndices);
}

inline bool Hash_complex::remove_simplex(Simplex *simplex, std::vector<simplex_handle> *removedIndices)
{
    std::unordered_map<vertex, Simplex*> *cofacets = simplex->get_cofacets();
    simplex_vertex_range *vertices = simplex->get_vertices();
    std::pair<simplex_vertex_range*,int> p(vertices, -1);

    if (vertices->size() > 1){
	for (int i = 0; i < (int)vertices->size(); i++){
	    p.second = i;
	    simplices_->at(&p)->get_cofacets()->erase(vertices->at(i));
	}
    }
    p.second = -1;

    while (!cofacets->empty()) remove_simplex(cofacets->begin()->second, removedIndices);

    std::pair<simplex_vertex_range*,int> *tmp = simplices_->find(&p)->first;
    if (removedIndices != nullptr) removedIndices->push_back(simplex->get_insertion_num());
    simplices_->erase(&p);
    handleToSimplex_->erase(simplex->get_insertion_num());
    delete tmp;
    delete simplex;
    return true;
}

inline Hash_complex::simplex_handle Hash_complex::get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_handle> *closedStar)
{
    std::queue<Simplex*> qv;
    std::queue<Simplex*> qu;
    Simplex *s;
    std::pair<simplex_vertex_range*,int> p;
    simplex_vertex_range vv(1, v);
    simplex_vertex_range vu(1, u);

    if (get_smallest_star(v, u, &qv, &qu) == v){
        while (!qv.empty()){
            s = qv.front();
            qv.pop();
	    p.first = s->get_vertices();
	    p.second = get_vertex_index(s->get_vertices(), v);
	    if (s->get_vertices()->size() > 1){
		closedStar->push_back(simplices_->at(&p)->get_insertion_num());
            }
            p.second = -1;
	    closedStar->push_back(s->get_insertion_num());
        }

	p.first = &vv;
    } else {
        while (!qu.empty()){
            s = qu.front();
            qu.pop();
	    p.first = s->get_vertices();
	    p.second = get_vertex_index(s->get_vertices(), u);
	    if (s->get_vertices()->size() > 1){
		closedStar->push_back(simplices_->at(&p)->get_insertion_num());
            }
            p.second = -1;
	    closedStar->push_back(s->get_insertion_num());
        }

	p.first = &vu;
    }

    p.second = -1;
    return simplices_->at(&p)->get_insertion_num();
}

inline void Hash_complex::get_boundary(simplex_handle &handle, std::vector<simplex_handle> *boundary)
{
    Simplex *simplex = handleToSimplex_->at(handle);
    std::pair<simplex_vertex_range*,int> p(simplex->get_vertices(), -1);

    if (simplex->get_vertices()->size() == 1) return;

    for (simplex_vertex_range::size_type i = 0; i < simplex->get_vertices()->size(); i++){
	p.second = i;
	boundary->push_back(simplices_->at(&p)->get_insertion_num());
    }
    std::sort(boundary->begin(), boundary->end());
}

inline Hash_complex::simplex_vertex_range& Hash_complex::get_vertices(simplex_handle &handle) const
{
    return *(handleToSimplex_->at(handle)->get_vertices());
}

inline Hash_complex::simplex_handle Hash_complex::get_simplex_handle(simplex_vertex_range &simplex)
{
    std::pair<simplex_vertex_range*,int> p(&simplex, -1);
    return simplices_->at(&p)->get_insertion_num();
}

inline int Hash_complex::get_dimension(simplex_handle handle)
{
    return handleToSimplex_->at(handle)->get_vertices()->size() - 1;
}

inline Hash_complex::size_type Hash_complex::get_size() const
{
    return simplices_->size();
}

inline Hash_complex::size_type Hash_complex::get_max_size() const
{
    return maxSize_;
}

inline int Hash_complex::get_max_dimension() const
{
    return maxDim_;
}

inline void Hash_complex::get_cofaces(simplex_vertex_range &simplex, std::vector<simplex_vertex_range*> *cofaces)
{
    auto hash = [](simplex_vertex_range* const& k) {
        std::size_t seed = k->size();
	for (simplex_vertex_range::size_type i = 0; i < k->size(); i++) {
            seed ^= (std::size_t)(k->at(i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    };
    auto comp = [](const simplex_vertex_range* s1, const simplex_vertex_range* s2) {
	simplex_vertex_range::size_type size = s1->size();
        if (size != s2->size()) return false;
	for (simplex_vertex_range::size_type i = 0; i < size; i++){
            if (s1->at(i) != s2->at(i)) return false;
        }
        return true;
    };

    std::unordered_map<simplex_vertex_range*, bool, decltype(hash), decltype(comp)> visited(100, hash, comp);
    std::pair<simplex_vertex_range*,int> p(&simplex, -1);
    std::unordered_map<vertex, Simplex*> *cofacets;
    std::unordered_map<vertex, Simplex*>::iterator it;
    std::queue<simplex_vertex_range*> q;

    q.push(&simplex);

    while (!q.empty()){
        p.first = q.front();
        q.pop();
        cofacets = simplices_->at(&p)->get_cofacets();
        for (it = cofacets->begin(); it != cofacets->end(); it++){
            p.first = it->second->get_vertices();
            if (visited.find(p.first) == visited.end()){
                visited.emplace(p.first, true);
		cofaces->push_back(new simplex_vertex_range(*(p.first)));
                q.push(p.first);
            }
        }
    }
}

inline void Hash_complex::get_cofacets(simplex_vertex_range &simplex, std::vector<simplex_vertex_range *> *cofacets)
{
    std::pair<simplex_vertex_range*,int> p(&simplex, -1);
    std::unordered_map<vertex, Simplex*> *mapedCofacets = simplices_->at(&p)->get_cofacets();
    for (auto it = mapedCofacets->begin(); it != mapedCofacets->end(); it++){
	cofacets->push_back(new simplex_vertex_range(*(it->second->get_vertices())));
    }
}

inline bool Hash_complex::contains(simplex_vertex_range &simplex)
{
    std::pair<simplex_vertex_range*,int> p(&simplex, -1);
    return (simplices_->find(&p) != simplices_->end());
}

inline void Hash_complex::print()
{
    for (auto it = simplices_->begin(); it != simplices_->end(); ++it){
	simplex_vertex_range *current = it->second->get_vertices();
	for (simplex_vertex_range::size_type i = 0; i < current->size(); ++i){
	    std::cout << current->at(i);
	}
	std::cout << "\n";
    }
}

inline Hash_complex::vertex Hash_complex::get_smallest_star(vertex v, vertex u, std::queue<Simplex*> *qv, std::queue<Simplex*> *qu)
{
    auto hash = [](simplex_vertex_range* const& k) {
        std::size_t seed = k->size();
	for (simplex_vertex_range::size_type i = 0; i < k->size(); i++) {
            seed ^= (std::size_t)(k->at(i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    };
    auto comp = [](const simplex_vertex_range* s1, const simplex_vertex_range* s2) {
	simplex_vertex_range::size_type size = s1->size();
        if (size != s2->size()) return false;
	for (simplex_vertex_range::size_type i = 0; i < size; i++){
            if (s1->at(i) != s2->at(i)) return false;
        }
        return true;
    };

    std::unordered_map<simplex_vertex_range*, bool, decltype(hash), decltype(comp)> visitedV(100, hash, comp);
    std::unordered_map<simplex_vertex_range*, bool, decltype(hash), decltype(comp)> visitedU(100, hash, comp);
    std::pair<simplex_vertex_range*,int> p;
    std::unordered_map<vertex, Simplex*> *cofacetsV;
    std::unordered_map<vertex, Simplex*> *cofacetsU;
    std::unordered_map<vertex, Simplex*>::iterator vit;
    std::unordered_map<vertex, Simplex*>::iterator uit;
    std::queue<simplex_vertex_range*> tv;
    std::queue<simplex_vertex_range*> tu;
    simplex_vertex_range sv(1, v);
    simplex_vertex_range su(1, u);

    p.first = &sv;
    p.second = -1;
    qv->push(simplices_->at(&p));
    cofacetsV = simplices_->at(&p)->get_cofacets();
    vit = cofacetsV->begin();

    p.first = &su;
    qu->push(simplices_->at(&p));
    cofacetsU = simplices_->at(&p)->get_cofacets();
    uit = cofacetsU->begin();

    if (vit != cofacetsV->end() && vit->first == u) vit++;
    if (uit != cofacetsU->end() && uit->first == v) uit++;

    while (vit != cofacetsV->end() && uit != cofacetsU->end()) {
        p.first = vit->second->get_vertices();
        visitedV.emplace(p.first, true);
	qv->push(vit->second);
        tv.push(p.first);

        p.first = uit->second->get_vertices();
        visitedU.emplace(p.first, true);
	qu->push(uit->second);
        tu.push(p.first);

        vit++;
        while ((vit == cofacetsV->end() && !tv.empty()) || (vit != cofacetsV->end() && (vit->first == u || visitedV.find(vit->second->get_vertices()) != visitedV.end()))){
            if (vit == cofacetsV->end() && !tv.empty()){
                p.first = tv.front();
                tv.pop();
                cofacetsV = simplices_->at(&p)->get_cofacets();
                vit = cofacetsV->begin();
            }
            if (vit != cofacetsV->end() && (vit->first == u || visitedV.find(vit->second->get_vertices()) != visitedV.end())) vit++;
        }

        uit++;
        while ((uit == cofacetsU->end() && !tu.empty()) || (uit != cofacetsU->end() && (uit->first == v || visitedU.find(uit->second->get_vertices()) != visitedU.end()))){
            if (uit == cofacetsU->end() && !tu.empty()){
                p.first = tu.front();
                tu.pop();
                cofacetsU = simplices_->at(&p)->get_cofacets();
                uit = cofacetsU->begin();
            }
            if (uit != cofacetsU->end() && (uit->first == v || visitedU.find(uit->second->get_vertices()) != visitedU.end())) uit++;
        }
    }

    if (uit == cofacetsU->end()) return u;
    else return v;
}

inline int Hash_complex::get_vertex_index(simplex_vertex_range *simplex, vertex v)
{
    int i = 0;
    while (i < (int)simplex->size() && simplex->at(i) != v){
        i++;
    }
    if (i == (int)simplex->size()) return -1;
    return i;
}

inline void Hash_complex::expand_simplex(simplex_vertex_range *vectSimplex, int maxDim, std::vector<simplex_handle> *addedSimplices)
{
    std::pair<simplex_vertex_range*,int> p(vectSimplex, -1);
    Simplex *simplex = simplices_->at(&p);
    p.second = vectSimplex->size() - 1;
    Simplex *facet = simplices_->at(&p);
    std::unordered_map<vertex, Simplex*> *cofacets = facet->get_cofacets();

    for (auto it = cofacets->begin(); it != cofacets->end(); ++it){
	Simplex *current = it->second;
	if (current == simplex) continue;
	Simplex *union_simplex = insert_union(current, vectSimplex->back());
	if (union_simplex != nullptr){	//if union could be inserted, i.e. all its facets were there and it-self was not already inserted
	    if (addedSimplices != nullptr) addedSimplices->push_back(union_simplex->get_insertion_num());
	    if ((int)union_simplex->get_vertices()->size() - 1 < maxDim  || maxDim == -1) expand_simplex(union_simplex->get_vertices(), maxDim, addedSimplices);
	}
    }
}

inline Hash_complex::Simplex *Hash_complex::insert_union(Simplex *simplex, vertex v)
{
    simplex_vertex_range unionVect;
    simplex_vertex_range *simplexVect = simplex->get_vertices();
    simplex_vertex_range::size_type i = 0;
    while (i < simplexVect->size() && simplexVect->at(i) < v){
	unionVect.push_back(simplexVect->at(i));
	++i;
    }
    if (i < simplexVect->size() && simplexVect->at(i) == v) return nullptr;
    unionVect.push_back(v);
    while (i < simplexVect->size()) {
	unionVect.push_back(simplexVect->at(i));
	++i;
    }

    i = 0;
    std::pair<simplex_vertex_range*,int> p(&unionVect, -1);
    while (i < unionVect.size()){
	p.second = i;
	if (simplices_->find(&p) == simplices_->end()) return nullptr;
	++i;
    }

    insert_simplex(unionVect);
    p.second = -1;
    return simplices_->at(&p);
}

inline Hash_complex::Simplex::Simplex(index num, simplex_vertex_range *vertices) : insertionNum_(num), vertices_(vertices)
{
    cofacets_ = new std::unordered_map<vertex, Simplex*>();
}

inline Hash_complex::Simplex::~Simplex()
{
    delete cofacets_;
    delete vertices_;
}

inline Hash_complex::index Hash_complex::Simplex::get_insertion_num() const
{
    return insertionNum_;
}

inline void Hash_complex::Simplex::set_insertion_num(const index &insertionNum)
{
    insertionNum_ = insertionNum;
}

inline void Hash_complex::Simplex::add_cofacet(Simplex *coface, vertex v)
{
    cofacets_->emplace(v, coface);
}

inline std::unordered_map<Hash_complex::vertex, Hash_complex::Simplex*>* Hash_complex::Simplex::get_cofacets() const
{
    return cofacets_;
}

inline Hash_complex::simplex_vertex_range *Hash_complex::Simplex::get_vertices() const
{
    return vertices_;
}

}
}

#endif // COMPLEX_H
