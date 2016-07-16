#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <ctime>

#include <list>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <utility>

namespace misc
{

typedef const char *cstr_t;

template <class T>
struct identity : std::unary_function<T, T>
{
	const T &operator () (const T &v) const {return v;}
}; // struct identity

} // namespace misc

//#endif // COMMON_H

#ifndef LINEAR_EQN_H
#define LINEAR_EQN_H

namespace linear
{

typedef std::vector<double> v_double;

bool Gauss2(size_t N, v_double eqn[], double rhs0[], double rhs1[]);

} // namespace linear

#endif // LINEAR_EQN_H


#ifndef GEOMETRY_H
#define GEOMETRY_H

namespace misc
{


struct dbl_pair
{
	double first, second;

	dbl_pair(double f = 0, double s = 0) : first(f), second(s) {}

	template <class T, class U>
	dbl_pair(const std::pair<T, U> &p) : first(p.first), second(p.second) {}

	// operators
	dbl_pair &operator += (const dbl_pair &r)
		{first += r.first; second += r.second; return *this;}

	dbl_pair &operator -= (const dbl_pair &r)
		{first -= r.first; second -= r.second; return *this;}

	dbl_pair &operator *= (const dbl_pair &r)
		{first *= r.first; second *= r.second; return *this;}

	dbl_pair &operator /= (const dbl_pair &r)
		{first /= r.first; second /= r.second; return *this;}

	dbl_pair &operator += (double r)
		{first += r; second += r; return *this;}

	dbl_pair &operator -= (double r)
		{first -= r; second -= r; return *this;}

	dbl_pair &operator *= (double r)
		{first *= r; second *= r; return *this;}

	dbl_pair &operator /= (double r)
		{first /= r; second /= r; return *this;}

}; // struct dbl_pair

inline dbl_pair operator + (dbl_pair l, const dbl_pair &r) {l += r; return l;}
inline dbl_pair operator - (dbl_pair l, const dbl_pair &r) {l -= r; return l;}
inline dbl_pair operator * (dbl_pair l, const dbl_pair &r) {l *= r; return l;}
inline dbl_pair operator / (dbl_pair l, const dbl_pair &r) {l /= r; return l;}
inline dbl_pair operator + (dbl_pair l, double r) {l += r; return l;}
inline dbl_pair operator - (dbl_pair l, double r) {l -= r; return l;}
inline dbl_pair operator * (dbl_pair l, double r) {l *= r; return l;}
inline dbl_pair operator / (dbl_pair l, double r) {l /= r; return l;}
inline dbl_pair operator + (double l, const dbl_pair &r) {dbl_pair ll(l, l); ll += r; return ll;}
inline dbl_pair operator - (double l, const dbl_pair &r) {dbl_pair ll(l, l); ll -= r; return ll;}
inline dbl_pair operator * (double l, const dbl_pair &r) {dbl_pair ll(l, l); ll *= r; return ll;}
inline dbl_pair operator / (double l, const dbl_pair &r) {dbl_pair ll(l, l); ll /= r; return ll;}

inline dbl_pair max(const dbl_pair &l, const dbl_pair &r)
{
	return dbl_pair(
		(l.first < r.first)? r.first: l.first,
		(l.second < r.second)? r.second: l.second);
}

inline dbl_pair min(const dbl_pair &l, const dbl_pair &r)
{
	return dbl_pair(
		(l.first > r.first)? r.first: l.first,
		(l.second > r.second)? r.second: l.second);
}

} // namespace misc

#endif // GEOMETRY_H


#ifndef TOKEN_PARSER_H
#define TOKEN_PARSER_H

namespace misc
{
	
// delimit chars are spaces
char *skip_space(char *buf, cstr_t delimit);
char *goto_space(char *buf, cstr_t delimit);

// in-place split by delimit chars
void split_tokens(char *buf, std::vector<cstr_t> *ptokens, cstr_t delimit);

// split to string vector
void split_tokens(char *buf, std::vector<std::string> *ptokens, cstr_t delimit);

// string table
size_t string_table_get_id(cstr_t str);
cstr_t string_table_get_str(size_t id);

} // namespace misc

#endif // TOKEN_PARSER_H

#ifndef bm_H
#define bm_H


namespace bm
{

struct signal_t;

struct module_t
{
	std::string name;
	std::vector<signal_t *> nets;
	double w, h; // weight/height for blocks and ratio for terminals
	double x, y; // center position
	int id;
}; // struct module_t

struct signal_t
{
	std::vector<module_t *> pins;
}; // struct signal_t

struct detailed_signal_t
{
	typedef std::pair<module_t *, misc::dbl_pair> pin_t;
	std::vector<pin_t> pins;
}; // struct detailed_signal_t

// read a bm benchmark, return # of blocks
int read_bm(const std::string &bm,
	std::vector<module_t> &modules,
	std::vector<signal_t> &signals,
	std::vector<detailed_signal_t> *detailed = 0,
	misc::dbl_pair *shape = 0);

} // namespace bm

#endif // bm_H

#ifndef NETS_H
#define NETS_H

namespace fplan
{

// nets that contain only modules
class module_nets_t
{
	std::vector<size_t> net_modules_;
	std::vector<size_t> nets_;

public:
	module_nets_t();

	size_t num_nets() const {return nets_.size()-1;}
	size_t num_modules(size_t inet) const
	{
		assert(inet < num_nets());
		return nets_[inet+1]-nets_[inet];
	}
	size_t module(size_t inet, size_t im) const
	{
		assert(inet < num_nets());
		assert(im < num_modules(inet));
		return net_modules_[nets_[inet]+im];
	}

	// add nets one by one
	void append_net();
	void append_module(size_t im);

		// WL function
	double HPWL(const misc::dbl_pair centers[]) const;
	
}; // class module_nets_t

// nets that contain pins
class pin_nets_t
{
public:
	typedef std::pair<size_t, misc::dbl_pair> pin_t;

private:
	std::vector<pin_t> net_pins_;
	std::vector<size_t> nets_;

public:
	pin_nets_t();

	size_t num_nets() const {return nets_.size()-1;}
	size_t num_pins(size_t inet) const
	{
		assert(inet < num_nets());
		return nets_[inet+1]-nets_[inet];
	}
	pin_t pin(size_t inet, size_t ip) const
	{
		assert(inet < num_nets());
		assert(ip < num_pins(inet));
		return net_pins_[nets_[inet]+ip];
	}

	// add nets one by one
	void append_net();
	void append_pin(size_t ip, const misc::dbl_pair &loc);

	// WL functions
	double HPWL(size_t nm, const misc::dbl_pair shapes[],
		const size_t ignions[], const misc::dbl_pair centers[]) const;


}; // class pin_nets_t

} // namespace fplan

#endif // NETS_H


#ifndef bm_INFO_H
#define bm_INFO_H


namespace fplan
{

void bm_info(int n,
	const std::vector<bm::module_t> &modules,
	const std::vector<bm::signal_t> &signals,
	std::vector<misc::dbl_pair> &module_shapes,
	std::vector<misc::cstr_t> &module_names,
	std::vector<misc::dbl_pair> &pin_locations,
	module_nets_t &mnets,
	const std::vector<bm::detailed_signal_t> *detailed = 0,
	pin_nets_t *pnets = 0);

} // namespace fplan

#endif // bm_INFO_H



#ifndef SHARED_PTR_H
#define SHARED_PTR_H


namespace misc
{

namespace detail
{

// policies for copy and deallocate

template <class T>
struct shared_ptr_policy_new_delete
{
	static void release(T *t) {delete t;}
	static T *clone(T *t) {return new T(t);}
}; // struct shared_ptr_policy_new_delete<T>

template <class T>
struct shared_ptr_policy_clone_release
{
	static void release(T *t) {t->release();}
	static T *clone(T *t) {return t->clone();}
}; // struct shared_ptr_policy_clone_release<T>

} // namespace detail

template <class T,
	class Policy = detail::shared_ptr_policy_clone_release<T>,
	class Count = int>
class shared_ptr
{
	T *p_;
	Count *count_;

public:
	typedef shared_ptr<T, Policy, Count> this_type;

	explicit shared_ptr(T *p) : p_(p), count_(new Count(1)) {}

	shared_ptr(const this_type &s) : p_(s.p_), count_(s.count_) {++*count_;}
	
	~shared_ptr()
	{
		if (--*count_ == 0)
		{
			if (p_)
				Policy::release(p_);
			delete count_;
		}
	}

	void swap(this_type &r)
	{
		std::swap(count_, r.count_);
		std::swap(p_, r.p_);
	}

	this_type &operator = (this_type r) {swap(r); return *this;}

	T &operator * () const {return *p_;}
	T *operator -> () const {return p_;}
	T *get() const {return p_;}

	this_type copy() const
	{
		if (p_)
			return this_type(Policy::clone(p_));
		else
			return this_type(0);
	}

}; // shared_ptr<T>

} // namespace misc

#endif // SHARED_PTR_H


#ifndef COST_H
#define COST_H


namespace fplan
{

namespace detail
{
	struct cost_i
	{
		virtual double cost() = 0;
		virtual void release() = 0;

	protected:
		~cost_i() {} // don't delete
	}; // struct cost_i
} // namespace detail

typedef misc::shared_ptr<detail::cost_i> cost_ptr;

} // namespace fplan

namespace misc
{

// unary
fplan::cost_ptr double_cost(double d);


// binary
fplan::cost_ptr pow(double d, const fplan::cost_ptr &c);
fplan::cost_ptr pow(const fplan::cost_ptr &c, double d);
fplan::cost_ptr pow(const fplan::cost_ptr &c, const fplan::cost_ptr &cc);

fplan::cost_ptr operator + (double d, const fplan::cost_ptr &c);
fplan::cost_ptr operator - (double d, const fplan::cost_ptr &c);
fplan::cost_ptr operator * (double d, const fplan::cost_ptr &c);
fplan::cost_ptr operator / (double d, const fplan::cost_ptr &c);
fplan::cost_ptr operator + (const fplan::cost_ptr &c, double d);
fplan::cost_ptr operator - (const fplan::cost_ptr &c, double d);
fplan::cost_ptr operator * (const fplan::cost_ptr &c, double d);
fplan::cost_ptr operator / (const fplan::cost_ptr &c, double d);
fplan::cost_ptr operator + (const fplan::cost_ptr &c, const fplan::cost_ptr &cc);
fplan::cost_ptr operator - (const fplan::cost_ptr &c, const fplan::cost_ptr &cc);
fplan::cost_ptr operator * (const fplan::cost_ptr &c, const fplan::cost_ptr &cc);
fplan::cost_ptr operator / (const fplan::cost_ptr &c, const fplan::cost_ptr &cc);

} // namespace misc

#endif // COST_H

#ifndef NETS_H
#define NETS_H



namespace fplan
{

// nets that contain only modules
class module_nets_t
{
	std::vector<size_t> net_modules_;
	std::vector<size_t> nets_;

public:
	module_nets_t();

	size_t num_nets() const {return nets_.size()-1;}
	size_t num_modules(size_t inet) const
	{
		assert(inet < num_nets());
		return nets_[inet+1]-nets_[inet];
	}
	size_t module(size_t inet, size_t im) const
	{
		assert(inet < num_nets());
		assert(im < num_modules(inet));
		return net_modules_[nets_[inet]+im];
	}

	// add nets one by one
	void append_net();
	void append_module(size_t im);
		
	// WL functions
	double HPWL(const misc::dbl_pair centers[]) const;
	
}; // class module_nets_t

// nets that contain pins
class pin_nets_t
{
public:
	typedef std::pair<size_t, misc::dbl_pair> pin_t;

private:
	std::vector<pin_t> net_pins_;
	std::vector<size_t> nets_;

public:
	pin_nets_t();

	size_t num_nets() const {return nets_.size()-1;}
	size_t num_pins(size_t inet) const
	{
		assert(inet < num_nets());
		return nets_[inet+1]-nets_[inet];
	}
	pin_t pin(size_t inet, size_t ip) const
	{
		assert(inet < num_nets());
		assert(ip < num_pins(inet));
		return net_pins_[nets_[inet]+ip];
	}

	// add nets one by one
	void append_net();
	void append_pin(size_t ip, const misc::dbl_pair &loc);

	// WL functions
	double HPWL(size_t nm, const misc::dbl_pair shapes[],
		const size_t ignions[], const misc::dbl_pair centers[]) const;

}; // class pin_nets_t

} // namespace fplan

#endif // NETS_H

#ifndef FPLAN_H
#define FPLAN_H

namespace fplan
{

// floorplan interface
struct fplan_i
{
	virtual size_t get_num_modules() = 0;
	virtual misc::dbl_pair get_shape() = 0;
	virtual misc::dbl_pair get_centers(size_t ids[], misc::dbl_pair centers[]) = 0;
	virtual misc::dbl_pair get_physical_info(
		size_t ids[], misc::dbl_pair centers[],
		misc::dbl_pair lls[], misc::dbl_pair urs[],
		size_t ign[], misc::dbl_pair ign_shapes[]) = 0;

	cost_ptr area_cost();
	cost_ptr fixed_cost(double w, double h, bool rotate = true);
	cost_ptr square_cost();

	// assume ids are 0, 1, ...
	cost_ptr HPWL_cost(const module_nets_t &mnets, size_t np, const misc::dbl_pair pins[]);
	
protected:
	~fplan_i() {} // no delete
}; // struct fplan_i

} // namespace fplan

#endif // FPLAN_H


#ifndef ANNEALING_H
#define ANNEALING_H


namespace fplan
{

namespace detail
{


// interface for perturbations
struct perturb_i
{
	// allow temperature aware perturbations
	virtual void perturbation(double T) = 0;
	virtual void reject() = 0;

	virtual void record_best() = 0;
	virtual void swap_best() = 0;

	virtual void release() = 0;

protected:
	~perturb_i() {} // no delete

}; // struct perturb_i

} // namespace detail

typedef misc::shared_ptr<detail::perturb_i> perturb_ptr;

int anneal(cost_ptr c, perturb_ptr pt, int N, FILE *fp = 0, int k = 10,
	double r = 0.90, double term_T = 0.1, double max_reject_rate = 0.90);


} // namespace fplan

#endif // ANNEALING_H

#ifndef ACG_BASE_H
#define ACG_BASE_H

//#include "../misc/geometry.h"

namespace acg
{

// fanins/fanouts, from/to
enum pos_t {POS_IN = 0, POS_OUT = 1};
inline pos_t other(pos_t pos) {return pos_t(1-pos);}

// horizental/vertical
enum color_t {COLOR_RED = 0, COLOR_BLACK = 1};
inline color_t other(color_t c) {return color_t(1-c);}

// forwards
class node_t;

class edge_t
{
	friend class edge_list_t;
	friend class acg_dev_t;

	// no copy, no assign
	edge_t(const edge_t &);
	edge_t &operator = (const edge_t &);

	edge_t *next_[2], *prev_[2]; // 0: same fanin, 1: same fanout
	
	node_t *node_[2]; // 0: from, 1: to

	edge_t() : m_data(0)
	{
		node_[0] = node_[1] = 0;
		next_[0] = next_[1] = prev_[0] = prev_[1] = this;
	}

	~edge_t() {}

public:

	mutable void *m_data; // scratching data for edges

	edge_t *next(pos_t pos) {return next_[pos];}
	edge_t *prev(pos_t pos) {return prev_[pos];}
	node_t *node(pos_t pos) {return node_[pos];}

	const edge_t *next(pos_t pos) const {return next_[pos];}
	const edge_t *prev(pos_t pos) const {return prev_[pos];}
	const node_t *node(pos_t pos) const {return node_[pos];}

	bool valid() const;

}; // class edge_t

class edge_list_t
{
	friend class node_t;

	// no copy, no assign
	edge_list_t(const edge_list_t &);
	edge_list_t &operator = (const edge_list_t &);

	edge_t head_; // two double linked circular list, cross linked with other lists

	// need this to init array elements...
	void init(node_t *n)
	{
		head_.node_[0] = head_.node_[1] = n;
	}

public:
	edge_list_t() {}

	node_t *node() {return head_.node(POS_IN);}
	const node_t *node() const {return head_.node(POS_IN);}

	edge_t *begin(pos_t pos) {return head_.next(pos);}
	edge_t *end(pos_t pos) {return &head_;}
	const edge_t *begin(pos_t pos) const {return head_.next(pos);}
	const edge_t *end(pos_t pos) const {return &head_;}

	bool empty(pos_t pos) const {return begin(pos) == end(pos);}

	typedef std::pair<edge_t *, edge_t *> range_t;
	typedef std::pair<const edge_t *, const edge_t *> const_range_t;
	range_t all(pos_t pos) {return std::make_pair(begin(pos), end(pos));}
	const_range_t all(pos_t pos) const {return std::make_pair(begin(pos), end(pos));}

	edge_t *front(pos_t pos) {assert(!empty(pos)); return head_.next(pos);}
	edge_t *back(pos_t pos) {assert(!empty(pos)); return head_.prev(pos);}
	const edge_t *front(pos_t pos) const {assert(!empty(pos)); return head_.next(pos);}
	const edge_t *back(pos_t pos) const {assert(!empty(pos)); return head_.next(pos);}

	edge_t *pop_front(pos_t pos);
	edge_t *pop_back(pos_t pos);
	void push_front(edge_t *edge, pos_t pos);
	void push_back(edge_t *edge, pos_t pos);

	// The following functions are all static because we don't keep the size of the edge list.

	// Insert empty 'edge' before 'after' on 'pos' side; return 'edge'.
	static void insert(edge_t *after, edge_t *edge, pos_t pos);

	// Erase 'edge' from 'pos' side.
	static void erase(edge_t *edge, pos_t pos);

	// Splice ['first', 'last') before 'after' at 'pos' side.
	static void splice(edge_t *after, edge_t *first, edge_t *last, pos_t pos);

	bool valid() const;

}; // class edge_list_t

class node_t
{
	friend class acg_dev_t;

	// no copy, no assign
	node_t(const node_t &);
	node_t &operator = (const node_t &);

	edge_list_t edges_[2]; // COLOR_BLACK, COLOR_RED

	double width_[2];
	double bb_[2][2]; // color/pos

	node_t();
	node_t(double w, double h, size_t id);
	~node_t() {}

public:
	size_t m_id; // user managed ids

	mutable void *m_data; // scratching data for nodes
	// more scratching data from edge list
	void *&edge_data(color_t c) const {return edges_[c].head_.m_data;}

	edge_list_t *edges(color_t c) {return &edges_[c];}
	const edge_list_t *edges(color_t c) const {return &edges_[c];}

	// several shortcuts to get boundary node directly
	// return self if no such node exists
	
	node_t *nearest_node(color_t c, pos_t pos)
		{return edges(c)->begin(pos)->node(other(pos));}
	const node_t *nearest_node(color_t c, pos_t pos) const
		{return edges(c)->begin(pos)->node(other(pos));}

	node_t *farest_node(color_t c, pos_t pos)
		{return edges(c)->end(pos)->prev(pos)->node(other(pos));}
	const node_t *farest_node(color_t c, pos_t pos) const
		{return edges(c)->end(pos)->prev(pos)->node(other(pos));}

	bool valid() const;

	// shape and packing
	// use functions so that you don't need to know how to index

	double &width(color_t c) {return width_[c];}
	double width(color_t c) const {return width_[c];}

	double &bb(color_t c, pos_t pos) {return bb_[c][pos];}
	double bb(color_t c, pos_t pos) const {return bb_[c][pos];}

	double boundary(color_t c, pos_t pos) const {return bb(c, pos)+width(c);}

	// tuples
	misc::dbl_pair width() const
		{return misc::dbl_pair(width_[COLOR_RED], width_[COLOR_BLACK]);}
	misc::dbl_pair bb(pos_t pos) const
		{return misc::dbl_pair(bb_[COLOR_RED][pos], bb_[COLOR_BLACK][pos]);}
	misc::dbl_pair boundary(pos_t pos) const {return bb(pos)+width();}
	misc::dbl_pair center(pos_t pos) const {return bb(pos)+width()/2;}

	// copy/swap shape, bb, and id
	void copy_content(const node_t &src);
	void swap_content(node_t &src);

}; // class node_t

class acg_base_t
{
	friend class acg_dev_t;

	size_t n_edges_;
	std::vector<node_t *> nodes_;

public:
	// create an empty one
	acg_base_t();
	// initialized as empty for direct building
	acg_base_t(size_t n, misc::dbl_pair shapes[], size_t ids[] = 0);

	~acg_base_t();
	void clear();

	// swap, copy, and assign
	void swap(acg_base_t &r);
	// node_t::m_data and edge_t::m_data will be overwritten
	acg_base_t(const acg_base_t &r);
	acg_base_t &operator = (acg_base_t r) {swap(r); return *this;}

	// statistics
	size_t num_nodes() const {return nodes_.size();}
	size_t num_edges() const {return n_edges_;}
	bool empty() const {return nodes_.empty();}

	// individual node visiting
	node_t *at(size_t i) {assert(i < num_nodes()); return nodes_[i];}
	const node_t *at(size_t i) const {assert(i < num_nodes()); return nodes_[i];}
	node_t *operator[](size_t i) {return at(i);}
	const node_t *operator[](size_t i) const {return at(i);}

	// write modules and edges using pointer as names
	// type default to ACG-BASE
	void xml_write(FILE *fp, const char *type = "ACG-BASE");
	// However, edge read in progress...
	// bool? xml_read(FILE *fp, ???);

	// node perturbations
	void node_swap(size_t i, size_t j);
	void node_rotate(size_t where);
	void node_resize(size_t where, double w, double h);

}; // class acg_base_t

} // namespace acg

namespace std
{

template<> inline void swap(acg::acg_base_t &a, acg::acg_base_t &b) {a.swap(b);}

} // namespace std

#endif // ACG_BASE_H

#ifndef ACG_DEV_H
#define ACG_DEV_H

namespace acg
{

class acg_dev_t
{
	size_t &n_edges_;
	std::vector<node_t *> &nodes_;

	// helper functions

	void dup_edges_POS_IN(color_t c, const node_t *src);
	void dup_edges_POS_OUT(color_t c, const node_t *src);

	void dump_node_edges(const node_t *node, color_t c, pos_t pos, FILE *fp) const;

public:
	explicit acg_dev_t(acg_base_t &base)
		: n_edges_(base.n_edges_), nodes_(base.nodes_)
		{}

	// statistics
	size_t num_nodes() const {return nodes_.size();}
	size_t num_edges() const {return n_edges_;}
	bool empty() const {return nodes_.empty();}

	// individual node visiting
	
	node_t *&at(size_t i) {assert(i < num_nodes()); return nodes_[i];}
	const node_t *at(size_t i) const {assert(i < num_nodes()); return nodes_[i];}
	node_t *&operator[](size_t i) {return at(i);}
	const node_t *operator[](size_t i) const {return at(i);}

	// count # of edges here
	edge_t *alloc_edge();
	void free_edge(edge_t *edge);

	// free a group of edges

	// free [first, last) at pos side
	void free_edges(edge_t *first, edge_t *last, pos_t pos);

	// free edges associated with a node
	void free_edges(node_t *node, pos_t pos);
	void free_edges(node_t *node, color_t c);
	void free_edges(node_t *node, color_t c, pos_t pos);

	void free_all_edges();

	// allocator for nodes
	node_t *alloc_node();
	node_t *alloc_node(double w, double h, size_t id);
	void free_node(node_t *node);

	// insert/erase nodes

	// POS_IN: front, POS_OUT: back
	node_t *front(pos_t pos);
	void push(node_t *node, pos_t pos);
	void pop(pos_t pos);

	void insert(size_t where, node_t *node);
	void erase(size_t where);

	// update nodes order
	// order[0] is the old pos of the new 0
	void update_order(const size_t order[]);

	// dumping for cg (acg/racg/rcg) debugging use
	// node_t::m_data will be overwritten
	// see acg::xml for better dumping
	void dump_cg_edges(FILE *fp) const;

	// duplicate edges
	// node_t::m_data and edge_t::m_data will be overwritten
	// assume dest has the same number of nodes as src and 0 edges
	void duplicate_edges(const acg_base_t &src);

}; // class acg_dev_t

} // namespace acg

#endif // ACG_DEV_H

#ifndef MALLOC_POOL_H
#define MALLOC_POOL_H

namespace misc
{

template <size_t n>
class malloc_pool
{
	void *pool_;

	static void *&next(void *p)
	{
		return *reinterpret_cast<void **>(p);
	}

public:
	malloc_pool() : pool_(0)
	{
		// compile time check
		// make sure the pool works
		char n_is_not_big_enough[(n < sizeof(void *))? 0: 1];
		(void)n_is_not_big_enough;
	}
	
	~malloc_pool()
	{
		for (; pool_ != 0;)
		{
			void *tmp = next(pool_);
			::free(pool_);
			pool_ = tmp;
		}

		// most of time a pool is a singleton
		// clear pool_ in case a singleton is resurrected
		pool_ = 0;
	}

	void free(void *p)
	{
		next(p) = pool_;
		pool_ = p;
	}

	void *alloc()
	{
		if (pool_ == 0)
			return malloc(n);
		else
		{
			void *tmp = pool_;
			pool_ = next(pool_);
			return tmp;
		}
	}
}; // class malloc_pool<n>

} // namespace misc

#endif // MALLOC_POOL_H


#ifndef ACG_RACG_H
#define ACG_RACG_H


namespace acg
{

// direct build
void racg_build(acg_base_t &base);

// head/tail insert/erase
void racg_push(acg_base_t &base,
	double w, double h, size_t id,
	color_t c, pos_t pos);
void racg_pop(acg_base_t &base, pos_t pos);

// general insert/erase
void racg_insert(acg_base_t &base, size_t where,
	double w, double h, size_t id);
void racg_erase(acg_base_t &base, size_t where);

// all the perturbations here
// switch edge type between red/black
// 1 adjacent edge/class 1
// 2 last non-alter edge/last class 2
// 3 first alter edge/class 3
void racg_switch_adjacent(acg_base_t &base, size_t where, pos_t pos);
void racg_switch_last(acg_base_t &base, size_t where, pos_t pos);
void racg_switch_alter(acg_base_t &base, size_t where, pos_t pos);

// Reduced ACG <=> ACG
typedef std::vector<edge_t *> racg_acg_context_t;
void racg_to_acg(acg_base_t &base, racg_acg_context_t &ctx);
void acg_to_racg(acg_base_t &base, racg_acg_context_t &ctx);

} // namespace acg

#endif // ACG_RACG_H



#ifndef CG_PACK_H
#define CG_PACK_H

//#include "acg_base.h"

namespace acg
{

// longest path packing
// assume partial order in params
// increase count by the number of the edges used
misc::dbl_pair packing_longest_path(acg_base_t &base, size_t &count);
misc::dbl_pair packing_longest_path(acg_base_t &base);

// compute longest path bounding box, which gives slack
void longest_path_bb(acg_base_t &base, misc::dbl_pair shape);

// pack at c direction, don't assume partial order
double packing_lp_unsort(acg_base_t &base, color_t c);

} // namespace acg

#endif // CG_PACK_H


#ifndef FPLAN_RCG_H
#define FPLAN_RCG_H

namespace fplan
{

// Use the get_fplan_interface function
// instead of a public inheritance here
// since probably we want to inherit the
// implementation of rcg_fp_t.
class rcg_fp_t : protected fplan_i
{
protected:
	friend class rcg_perturb_t;

	acg::acg_base_t acg_base_;
	std::vector<size_t> ign_;

	misc::dbl_pair shape_;
	bool mapped_;

	virtual void do_packing();
	virtual void compute_bb();

	// by pass the build process so derived classes could work
	rcg_fp_t();

public:
	rcg_fp_t(size_t n, misc::dbl_pair shapes[], size_t ids[] = 0,
		size_t ign[] = 0, bool build_rcg = true); // last param is for derived class
	~rcg_fp_t();

	void swap(rcg_fp_t &r);

	// default ones work here
	// rcg_fp_t(const rcg_fp_t &r);
	// rcg_fp_t &operator = (rcg_fp_t r);

	virtual size_t get_num_modules();
	virtual misc::dbl_pair get_shape();
	virtual misc::dbl_pair get_centers(size_t ids[], misc::dbl_pair centers[]);
	virtual misc::dbl_pair get_physical_info(
		size_t ids[], misc::dbl_pair centers[],
		misc::dbl_pair lls[], misc::dbl_pair urs[],
		size_t ign[], misc::dbl_pair ign_shapes[]);

	fplan_i *get_fplan_interface() {return this;}

	bool dump(const std::string &name, misc::cstr_t module_names[] = 0);

	perturb_ptr create_perturbation(double ign, double exchangeP);

}; // class rcg_fp_t

} // namespace fplan

namespace std
{

template<> inline void swap(fplan::rcg_fp_t &a, fplan::rcg_fp_t &b) {a.swap(b);}

} // namespace std

#endif // FPLAN_RCG_H

#ifndef FPLAN_ACG_H
#define FPLAN_ACG_H


namespace fplan
{

class acg_fp_t : protected rcg_fp_t
{
protected:
	friend class acg_perturb_t;

	acg::racg_acg_context_t ctx_;
	bool reduced_;

	void require_acg();
	void require_racg();

	virtual void do_packing();
	virtual void compute_bb();

public:
	acg_fp_t(size_t n, misc::dbl_pair shapes[], size_t ids[] = 0, size_t ign[] = 0);
	~acg_fp_t();

	void swap(acg_fp_t &r);

	acg_fp_t(acg_fp_t &r);
	acg_fp_t &operator = (acg_fp_t r);

	using rcg_fp_t::get_num_modules;
	using rcg_fp_t::get_shape;
	using rcg_fp_t::get_centers;
	using rcg_fp_t::get_physical_info;

	using rcg_fp_t::get_fplan_interface;

	using rcg_fp_t::dump;

	perturb_ptr create_perturbation(double ign, double exchangeP);

}; // class acg_fp_t

} // namespace fplan

namespace std
{

template<> inline void swap(fplan::acg_fp_t &a, fplan::acg_fp_t &b) {a.swap(b);}

} // namespace std

#endif // FPLAN_ACG_H
#endif // COMMON_H
