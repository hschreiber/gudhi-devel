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
#include <iostream>
#include <functional>

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
	using simplex_key = std::pair<simplex_vertex_range*, int>;

	Hash_complex();

	class Simplex
	{
	public:
		Simplex(index num, simplex_vertex_range &vertices);
		Simplex(Simplex&& toMove) noexcept;

		index get_insertion_num() const;
		void set_insertion_num(const index &insertionNum);

		void add_cofacet(Simplex &coface, vertex v);
		std::unordered_map<vertex, Simplex*>& get_cofacets();

		simplex_vertex_range &get_vertices();

	private:
		index insertionNum_;
		std::unordered_map<vertex, Simplex*> cofacets_;
		simplex_vertex_range vertices_;
	};

	struct Key_hasher {
		std::size_t operator()(const simplex_key &k) const
		{
			std::size_t seed;
			if (k.second < 0) seed = k.first->size();
			else seed = k.first->size() - 1;

			for (int i = 0; i < static_cast<int>(k.first->size()); i++) {
				if (i != k.second) seed ^= static_cast<std::size_t>(k.first->at(i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};

	struct Simplices_equals : std::binary_function<const simplex_key&, const simplex_key&, bool> {
		bool operator()(const simplex_key &s1, const simplex_key &s2) const
		{
			const simplex_key &key = s1.second > -1 ? s1 : s2;
			const simplex_key &inMap = s1.second > -1 ? s2 : s1;
			simplex_vertex_range::size_type size;

			if (key.second < 0) size = key.first->size();
			else size = key.first->size() - 1;
			if (size != inMap.first->size()) return false;
			int j = 0;
			for (simplex_vertex_range::size_type i = 0; i < size; i++){
				if (j == key.second) j++;
				if (key.first->at(j) != inMap.first->at(i)) {
					return false;
				}
				j++;
			}
			return true;
		}
	};

	/* Important for tower_converter.h */

	std::vector<simplex_handle> insert_simplex(simplex_vertex_range &simplex);
	std::vector<simplex_handle> insert_simplex_and_faces(simplex_vertex_range &simplex);
	std::vector<simplex_handle> insert_edge_and_expand(vertex u, vertex v, int maxDim = -1);	// u < v !!! ; maxDim == -1 -> no limit
	std::vector<index> remove_simplex(simplex_vertex_range &simplex);
	std::vector<index> remove_simplex(simplex_handle simplex);
	simplex_handle get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_handle> &closedStar);    //returns vertex with smallest closed star; closedStar is ordered
	std::vector<simplex_handle> get_boundary(simplex_handle simplex);   //ordered by increasing insertion numbers
	const simplex_vertex_range& get_vertices(simplex_handle simplex);
	index get_key(simplex_handle simplex);

	/* Other */

	size_type get_size() const;
	int get_dimension(simplex_handle handle);
	simplex_handle get_simplex_handle(simplex_vertex_range &simplex);
	size_type get_max_size() const;
	int get_max_dimension() const;
	std::vector<simplex_handle> get_cofaces(simplex_handle simplex);
	std::vector<simplex_handle> get_cofacets(simplex_handle simplex);
	bool contains(simplex_vertex_range &simplex);
	void print();

private:
	index maxIndex_;
	size_type maxSize_;
	int maxDim_;
	std::unordered_map<
			  simplex_key,
			  Simplex*,
			  Key_hasher,
			  Simplices_equals
			 > simplices_;
	std::unordered_map<simplex_handle, Simplex> handleToSimplex_;

	vertex _get_smallest_star(vertex v, vertex u, std::queue<Simplex*> &qv, std::queue<Simplex*> &qu);
	int _get_vertex_index(simplex_vertex_range &simplex, vertex v);
	void _expand_simplex(simplex_vertex_range &vectSimplex, int maxDim, std::vector<simplex_handle> &addedSimplices);
	Simplex* _insert_union(Simplex *simplex, vertex v);
	void _remove_simplex(Simplex &simplex, std::vector<index>& removedIndices);
};

inline Hash_complex::Hash_complex() : maxIndex_(-1), maxSize_(0), maxDim_(0)
{}

inline std::vector<Hash_complex::simplex_handle> Hash_complex::insert_simplex(simplex_vertex_range &simplex)
{
	if (simplices_.find(simplex_key(&simplex, -1)) != simplices_.end()) {
		return std::vector<Hash_complex::simplex_handle>();
	}

	++maxIndex_;
	auto res = handleToSimplex_.emplace(maxIndex_, Simplex(maxIndex_, simplex));
	Simplex &splx = res.first->second;
	simplex_vertex_range &vertices = splx.get_vertices();
	simplex_key key(&vertices, -1);
	simplices_.emplace(key, &splx);

	if (static_cast<int>(vertices.size()) - 1 > maxDim_)
		maxDim_ = vertices.size() - 1;

	if (maxSize_ < static_cast<size_type>(simplices_.size()))
		maxSize_ = simplices_.size();

	if (vertices.size() <= 1)
		return std::vector<Hash_complex::simplex_handle>(1, maxIndex_);

	for (simplex_vertex_range::size_type i = 0; i < vertices.size(); i++){
		key.second = i;
		simplices_.at(key)->add_cofacet(splx, vertices.at(i));
	}

	return std::vector<Hash_complex::simplex_handle>(1, maxIndex_);
}

inline std::vector<Hash_complex::simplex_handle> Hash_complex::insert_simplex_and_faces(simplex_vertex_range &simplex)
{
	std::vector<simplex_handle> addedSimplices;

	if (simplices_.find(simplex_key(&simplex, -1)) != simplices_.end())
		return addedSimplices;

	std::vector<simplex_vertex_range> simplices(1);

	for (simplex_vertex_range::size_type i = 0; i < simplex.size(); ++i){
		vertex v = simplex.at(i);
		std::vector<simplex_vertex_range> tmpSimplices(simplices);

		for (unsigned int i = 0; i < simplices.size(); ++i){
			tmpSimplices[i].push_back(v);
			std::vector<Hash_complex::simplex_handle> currentIndex = insert_simplex(tmpSimplices[i]);
			if (!currentIndex.empty()) addedSimplices.push_back(currentIndex.front());
		}

		simplices.insert(simplices.end(), std::make_move_iterator(tmpSimplices.begin()), std::make_move_iterator(tmpSimplices.end()));
	}

	return addedSimplices;
}

inline std::vector<Hash_complex::simplex_handle> Hash_complex::insert_edge_and_expand(vertex u, vertex v, int maxDim)
{
	std::vector<simplex_handle> addedSimplices;
	simplex_vertex_range vect;

	vect.push_back(u);
	std::vector<Hash_complex::simplex_handle> handles = insert_simplex(vect);
	if (!handles.empty())
		addedSimplices.push_back(handles.front());

	vect.at(0) = v;
	handles = insert_simplex(vect);
	if (!handles.empty())
		addedSimplices.push_back(handles.front());

	vect.at(0) = u;
	vect.push_back(v);
	handles = insert_simplex(vect);
	if (!handles.empty()) {
		addedSimplices.push_back(handles.front());
		if (maxDim > 1 || maxDim == -1) _expand_simplex(vect, maxDim, addedSimplices);
		return addedSimplices;
	}

	return addedSimplices;
}

inline std::vector<Hash_complex::simplex_handle> Hash_complex::remove_simplex(simplex_vertex_range &simplex)
{
	std::vector<index> removedIndices;
	simplex_key key(&simplex, -1);
	if (simplices_.find(key) == simplices_.end()) return std::vector<simplex_handle>();
	_remove_simplex(*simplices_.at(key), removedIndices);
	return removedIndices;
}

inline std::vector<Hash_complex::simplex_handle> Hash_complex::remove_simplex(simplex_handle simplex)
{
	std::vector<index> removedIndices;
	if (handleToSimplex_.find(simplex) == handleToSimplex_.end()) return std::vector<Hash_complex::simplex_handle>();
	_remove_simplex(handleToSimplex_.at(simplex), removedIndices);
	return removedIndices;
}

inline void Hash_complex::_remove_simplex(Simplex &simplex, std::vector<index>& removedIndices)
{
	std::unordered_map<vertex, Simplex*> &cofacets = simplex.get_cofacets();
	simplex_vertex_range &vertices = simplex.get_vertices();
	simplex_key key(&vertices, -1);

	if (vertices.size() > 1){
		for (int i = 0; i < (int)vertices.size(); i++){
			key.second = i;
			simplices_.at(key)->get_cofacets().erase(vertices.at(i));
		}
	}
	key.second = -1;

	while (!cofacets.empty())
		_remove_simplex(*(cofacets.begin()->second), removedIndices);

	removedIndices.push_back(simplex.get_insertion_num());
	simplices_.erase(key);
	handleToSimplex_.erase(simplex.get_insertion_num());
}

inline Hash_complex::simplex_handle Hash_complex::get_smallest_closed_star(vertex v, vertex u, std::vector<simplex_handle> &closedStar)
{
	std::queue<Simplex*> qv;
	std::queue<Simplex*> qu;
	Simplex *s;
	simplex_key key;
	simplex_vertex_range vAsRange(1, v);
	simplex_vertex_range uAsRange(1, u);

	if (_get_smallest_star(v, u, qv, qu) == v){
		while (!qv.empty()){
			s = qv.front();
			qv.pop();
			key.first = &s->get_vertices();
			key.second = _get_vertex_index(s->get_vertices(), v);
			if (s->get_vertices().size() > 1){
				closedStar.push_back(simplices_.at(key)->get_insertion_num());
			}
			key.second = -1;
			closedStar.push_back(s->get_insertion_num());
		}

		key.first = &vAsRange;
	} else {
		while (!qu.empty()){
			s = qu.front();
			qu.pop();
			key.first = &s->get_vertices();
			key.second = _get_vertex_index(s->get_vertices(), u);
			if (s->get_vertices().size() > 1){
				closedStar.push_back(simplices_.at(key)->get_insertion_num());
			}
			key.second = -1;
			closedStar.push_back(s->get_insertion_num());
		}

		key.first = &uAsRange;
	}

	key.second = -1;
	return simplices_.at(key)->get_insertion_num();
}

inline std::vector<Hash_complex::simplex_handle> Hash_complex::get_boundary(simplex_handle handle)
{
	std::vector<simplex_handle> boundary;
	Simplex &simplex = handleToSimplex_.at(handle);
	simplex_key key(&simplex.get_vertices(), -1);

	if (simplex.get_vertices().size() == 1) return boundary;

	for (simplex_vertex_range::size_type i = 0; i < simplex.get_vertices().size(); i++){
		key.second = i;
		boundary.push_back(simplices_.at(key)->get_insertion_num());
	}
	std::sort(boundary.begin(), boundary.end());

	return boundary;
}

inline const Hash_complex::simplex_vertex_range& Hash_complex::get_vertices(simplex_handle handle)
{
	return handleToSimplex_.at(handle).get_vertices();
}

inline Hash_complex::index Hash_complex::get_key(Hash_complex::simplex_handle simplex)
{
	return simplex;
}

inline Hash_complex::simplex_handle Hash_complex::get_simplex_handle(simplex_vertex_range &simplex)
{
	return simplices_.at(simplex_key(&simplex, -1))->get_insertion_num();
}

inline int Hash_complex::get_dimension(simplex_handle handle)
{
	return handleToSimplex_.at(handle).get_vertices().size() - 1;
}

inline Hash_complex::size_type Hash_complex::get_size() const
{
	return simplices_.size();
}

inline Hash_complex::size_type Hash_complex::get_max_size() const
{
	return maxSize_;
}

inline int Hash_complex::get_max_dimension() const
{
	return maxDim_;
}

inline std::vector<Hash_complex::simplex_handle> Hash_complex::get_cofaces(simplex_handle handle)
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

	std::vector<simplex_handle> cofaces;
	Simplex& simplex = handleToSimplex_.at(handle);
	std::unordered_map<simplex_vertex_range*, bool, decltype(hash), decltype(comp)> visited(100, hash, comp);
	simplex_key key(&simplex.get_vertices(), -1);
	std::unordered_map<vertex, Simplex*> cofacets;
	std::queue<simplex_vertex_range*> q;

	q.push(&simplex.get_vertices());

	while (!q.empty()){
		key.first = q.front();
		q.pop();
		cofacets = simplices_.at(key)->get_cofacets();
		for (auto it : cofacets){
			Simplex* s = it.second;
			key.first = &s->get_vertices();
			if (visited.find(key.first) == visited.end()){
				visited.emplace(key.first, true);
				cofaces.push_back(s->get_insertion_num());
				q.push(key.first);
			}
		}
	}

	return cofaces;
}

inline std::vector<Hash_complex::simplex_handle> Hash_complex::get_cofacets(simplex_handle handle)
{
	std::unordered_map<vertex, Simplex*> &mapedCofacets = handleToSimplex_.at(handle).get_cofacets();
	std::vector<simplex_handle> cofacets(mapedCofacets.size());
	unsigned int i = 0;
	for (auto it : mapedCofacets){
		cofacets[i++] = it.second->get_insertion_num();
	}
	return cofacets;
}

inline bool Hash_complex::contains(simplex_vertex_range &simplex)
{
	return simplices_.find(simplex_key(&simplex, -1)) != simplices_.end();
}

inline void Hash_complex::print()
{
	for (auto it = simplices_.begin(); it != simplices_.end(); ++it){
		simplex_vertex_range &current = it->second->get_vertices();
		for (simplex_vertex_range::size_type i = 0; i < current.size(); ++i){
			std::cout << current.at(i);
		}
		std::cout << "\n";
	}
}

inline Hash_complex::vertex Hash_complex::_get_smallest_star(vertex v, vertex u, std::queue<Simplex*> &qv, std::queue<Simplex*> &qu)
{
	auto hash = [](simplex_vertex_range* const& k) {
		std::size_t seed = k->size();
		for (simplex_vertex_range::size_type i = 0; i < k->size(); i++) {
			seed ^= static_cast<std::size_t>(k->at(i)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
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
	simplex_key key;
	std::unordered_map<vertex, Simplex*>::iterator vit;
	std::unordered_map<vertex, Simplex*>::iterator uit;
	std::queue<simplex_vertex_range*> tv;
	std::queue<simplex_vertex_range*> tu;
	simplex_vertex_range sv(1, v);
	simplex_vertex_range su(1, u);

	key.first = &sv;
	key.second = -1;
	qv.push(simplices_.at(key));
	std::unordered_map<vertex, Simplex*>* cofacetsV = &simplices_.at(key)->get_cofacets();
	vit = cofacetsV->begin();

	key.first = &su;
	qu.push(simplices_.at(key));
	std::unordered_map<vertex, Simplex*>* cofacetsU = &simplices_.at(key)->get_cofacets();
	uit = cofacetsU->begin();

	if (vit != cofacetsV->end() && vit->first == u) vit++;
	if (uit != cofacetsU->end() && uit->first == v) uit++;

	while (vit != cofacetsV->end() && uit != cofacetsU->end()) {
		key.first = &vit->second->get_vertices();
		visitedV.emplace(key.first, true);
		qv.push(vit->second);
		tv.push(key.first);

		key.first = &uit->second->get_vertices();
		visitedU.emplace(key.first, true);
		qu.push(uit->second);
		tu.push(key.first);

		vit++;
		while ((vit == cofacetsV->end() && !tv.empty()) ||
			   (vit != cofacetsV->end() && (vit->first == u || visitedV.find(&vit->second->get_vertices()) != visitedV.end())))
		{
			if (vit == cofacetsV->end() && !tv.empty()){
				key.first = tv.front();
				tv.pop();
				cofacetsV = &simplices_.at(key)->get_cofacets();
				vit = cofacetsV->begin();
			}
			if (vit != cofacetsV->end() &&
					(vit->first == u || visitedV.find(&vit->second->get_vertices()) != visitedV.end()))
				vit++;
		}

		uit++;
		while ((uit == cofacetsU->end() && !tu.empty()) ||
			   (uit != cofacetsU->end() && (uit->first == v || visitedU.find(&uit->second->get_vertices()) != visitedU.end())))
		{
			if (uit == cofacetsU->end() && !tu.empty()){
				key.first = tu.front();
				tu.pop();
				cofacetsU = &simplices_.at(key)->get_cofacets();
				uit = cofacetsU->begin();
			}
			if (uit != cofacetsU->end() &&
					(uit->first == v || visitedU.find(&uit->second->get_vertices()) != visitedU.end()))
				uit++;
		}
	}

	if (uit == cofacetsU->end()) return u;
	else return v;
}

inline int Hash_complex::_get_vertex_index(simplex_vertex_range &simplex, vertex v)
{
	int i = 0;
	while (i < static_cast<int>(simplex.size()) && simplex.at(i) != v){
		i++;
	}
	if (i == static_cast<int>(simplex.size())) return -1;
	return i;
}

inline void Hash_complex::_expand_simplex(simplex_vertex_range &vectSimplex, int maxDim, std::vector<simplex_handle> &addedSimplices)
{
	simplex_key p(&vectSimplex, -1);
	Simplex *simplex = simplices_.at(p);
	p.second = vectSimplex.size() - 1;
	Simplex *facet = simplices_.at(p);
	std::unordered_map<vertex, Simplex*> &cofacets = facet->get_cofacets();

	for (auto it = cofacets.begin(); it != cofacets.end(); ++it){
		Simplex *current = it->second;
		if (current == simplex) continue;
		Simplex *union_simplex = _insert_union(current, vectSimplex.back());
		if (union_simplex != nullptr){	//if union could be inserted, i.e. all its facets were there and it-self was not already inserted
			addedSimplices.push_back(union_simplex->get_insertion_num());
			if (static_cast<int>(union_simplex->get_vertices().size()) - 1 < maxDim  || maxDim == -1)
				_expand_simplex(union_simplex->get_vertices(), maxDim, addedSimplices);
		}
	}
}

inline Hash_complex::Simplex* Hash_complex::_insert_union(Simplex *simplex, vertex v)
{
	simplex_vertex_range unionVect;
	simplex_vertex_range &simplexVect = simplex->get_vertices();
	simplex_vertex_range::size_type i = 0;
	while (i < simplexVect.size() && simplexVect.at(i) < v){
		unionVect.push_back(simplexVect.at(i));
		++i;
	}
	if (i < simplexVect.size() && simplexVect.at(i) == v) return nullptr;
	unionVect.push_back(v);
	while (i < simplexVect.size()) {
		unionVect.push_back(simplexVect.at(i));
		++i;
	}

	i = 0;
	std::pair<simplex_vertex_range*,int> p(&unionVect, -1);
	while (i < unionVect.size()){
		p.second = i;
		if (simplices_.find(p) == simplices_.end()) return nullptr;
		++i;
	}

	insert_simplex(unionVect);
	p.second = -1;
	return simplices_.at(p);
}

inline Hash_complex::Simplex::Simplex(index num, simplex_vertex_range &vertices)
	: insertionNum_(num), vertices_(vertices)
{}

inline Hash_complex::Simplex::Simplex(Simplex &&toMove) noexcept
	: insertionNum_(std::exchange(toMove.insertionNum_, 0)),
	  cofacets_(std::move(toMove.cofacets_)),
	  vertices_(std::move(toMove.vertices_))
{}

inline Hash_complex::index Hash_complex::Simplex::get_insertion_num() const
{
	return insertionNum_;
}

inline void Hash_complex::Simplex::set_insertion_num(const index &insertionNum)
{
	insertionNum_ = insertionNum;
}

inline void Hash_complex::Simplex::add_cofacet(Simplex &coface, vertex v)
{
	cofacets_.emplace(v, &coface);
}

inline std::unordered_map<Hash_complex::vertex, Hash_complex::Simplex*>& Hash_complex::Simplex::get_cofacets()
{
	return cofacets_;
}

inline Hash_complex::simplex_vertex_range& Hash_complex::Simplex::get_vertices()
{
	return vertices_;
}

}
}

#endif // COMPLEX_H
