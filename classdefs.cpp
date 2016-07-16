#include "hfile.h"
#include <iostream>
using namespace std;
extern time_t now_g;
namespace misc
{

char *skip_space(char *buf, cstr_t delimit)
{
	for (; *buf != 0; ++buf)
	{
		size_t i;
		for (i = 0; delimit[i] != 0; ++i)
			if (*buf == delimit[i])
				break;
		if (delimit[i] == 0)
			return buf;
	}
	return buf;
}

char *goto_space(char *buf, cstr_t delimit)
{
	for (; *buf != 0; ++buf)
		for (size_t i = 0; delimit[i] != 0; ++i)
			if (*buf == delimit[i])
				return buf;
	return buf;
}

void split_tokens(char *buf, std::vector<cstr_t> *ptokens, cstr_t delimit)
{
	for (char *begin = buf; *begin != 0;)
	{
		char *end = goto_space(begin, delimit);
		ptokens->push_back(begin);
		if (*end == 0)
			return;
		*end = 0;
		begin = skip_space(end+1, delimit);
	}
}

// split to string vector
void split_tokens(char *buf, std::vector<std::string> *ptokens, cstr_t delimit)
{
	for (char *begin = buf; *begin != 0;)
	{
		char *end = goto_space(begin, delimit);
		ptokens->push_back(std::string(begin, end));
		if (*end == 0)
			return;
		begin = skip_space(end+1, delimit);
	}
}

// string table
typedef std::map<std::string, size_t> map_str2id;
typedef std::vector<const char *> str_v;

static map_str2id st_str_tbl;
static str_v st_strs;

size_t string_table_get_id(const char *str)
{
	map_str2id::iterator i = st_str_tbl.find(str);
	if (i == st_str_tbl.end())
	{
		i = st_str_tbl.insert(std::make_pair(str, st_strs.size())).first;
		st_strs.push_back(i->first.c_str());
	}
	return i->second;
}

const char *string_table_get_str(size_t id)
{
	assert(id < st_strs.size());

	return st_strs[id];
}

} // namespace misc

namespace linear
{

typedef std::vector<size_t> v_size_t;

bool Gauss2(size_t N, v_double eqn[], double rhs0[], double rhs1[])
{
	v_size_t col_order(N);
	const double eps = 1e-10;

	for (size_t i = 0; i < N; ++i)
		col_order[i] = i;

	// deduce
	for (size_t row = 0; row < N; ++row)
	{
		// select row max
		size_t cdmax = row;
		size_t colmax = col_order[cdmax];
		double valmax = fabs(eqn[row][colmax]);
		for (size_t i = row+1; i < N; ++i)
		{
			size_t col = col_order[i];
			double val = fabs(eqn[row][col]);
			if (valmax < val)
				valmax = val, colmax = col, cdmax = i;
		}
		if (valmax < eps)
		{
			fprintf(stderr, "Zero at row %u\n", row);
			return false;
		}

		// booking row max
		col_order[cdmax] = col_order[row];
		col_order[row] = colmax;

		// solve
		for (size_t j = row+1; j < N; ++j)
		{
			if (fabs(eqn[j][colmax]) < eps)
				continue;

			double mul = eqn[j][colmax]/eqn[row][colmax];
			for (size_t i = row+1; i < N; ++i)
			{
				size_t col = col_order[i];
				eqn[j][col] -= mul*eqn[row][col];
			}
			rhs0[j] -= mul*rhs0[row];
			rhs1[j] -= mul*rhs1[row];
		}
	}

	// solve for rhs
	v_double r0(rhs0, rhs0+N), r1(rhs1, rhs1+N);
	for (int j = int(N-1); j >= 0; --j)
	{
		size_t col = col_order[j];
		double rt0 = r0[j];
		double rt1 = r1[j];
		for (size_t i = j+1; i < N; ++i)
		{
			size_t coll = col_order[i];
			rt0 -= rhs0[coll]*eqn[j][coll];
			rt1 -= rhs1[coll]*eqn[j][coll];
		}
		rhs0[col] = rt0/eqn[j][col];
		rhs1[col] = rt1/eqn[j][col];
	}
	
	return true;
}

} // namespace linear


namespace fplan
{

void bm_info(int n,
	const std::vector<bm::module_t> &modules,
	const std::vector<bm::signal_t> &signals,
	std::vector<misc::dbl_pair> &module_shapes,
	std::vector<misc::cstr_t> &module_names,
	std::vector<misc::dbl_pair> &pin_locations,
	module_nets_t &mnets,
	const std::vector<bm::detailed_signal_t> *detailed,
	pin_nets_t *pnets)
{
	const size_t N = size_t(n);

	assert(modules.size() >= N);

	module_shapes.reserve(N);
	module_names.reserve(N);
	pin_locations.reserve(modules.size()-N);

	for (size_t i = 0; i < N; ++i)
	{
		module_shapes.push_back(misc::dbl_pair(modules[i].w, modules[i].h));
		module_names.push_back(modules[i].name.c_str());
	}
	for (size_t j = N; j < modules.size(); ++j)
		pin_locations.push_back(misc::dbl_pair(modules[j].w, modules[j].h));
	for (size_t k = 0; k < signals.size(); ++k)
	{
		mnets.append_net();
		for (size_t i = 0; i < signals[k].pins.size(); ++i)
			mnets.append_module(signals[k].pins[i]-&modules[0]);
	}
	if (detailed != 0)
	{
		for (size_t i = 0 ; i < detailed->size(); ++i)
		{
			pnets->append_net();
			for (size_t j = 0; j < detailed->at(i).pins.size(); ++j)
				pnets->append_pin(
					detailed->at(i).pins[j].first-&modules[0], 
					detailed->at(i).pins[j].second);
		}
	}
}


} // namespace fplan

namespace bm
{

using namespace misc;

typedef std::map<std::string, module_t *> map_modules;

int bm_blocks(const std::string &blocks,
	std::vector<module_t> &modules, map_modules &names);

bool bm_nets(const std::string &nets,
	std::vector<signal_t> &signals, map_modules &names,
	std::vector<detailed_signal_t> *detailed);
bool bm_pl(const std::string &pl,
	std::vector<module_t> &modules, map_modules &names);

int read_bm(const std::string &bm,
	std::vector<module_t> &modules,
	std::vector<signal_t> &signals,
	std::vector<detailed_signal_t> *detailed,
	misc::dbl_pair *shape)
{
	assert(modules.size() == 0);
	assert(signals.size() == 0);

	std::string blocks = bm+".blocks";
	std::string nets = bm+".nets";
	
	map_modules names;

	int nblocks = bm_blocks(blocks, modules, names);

	if (nblocks == 0)
	{
		modules.clear();
		signals.clear();
		return 0;
	}

	if (!bm_nets(nets, signals, names, detailed))
		signals.clear();

	// calculating terminal positions
	double maxx = -1e9, minx = 1e9, maxy = -1e9, miny = 1e9;
	for (size_t i = nblocks; i < modules.size(); ++i)
	{
		if (modules[i].x > maxx) maxx = modules[i].x;
		if (modules[i].y > maxy) maxy = modules[i].y;
		if (modules[i].x < minx) minx = modules[i].x;
		if (modules[i].y < miny) miny = modules[i].y;
	}
	if (maxx > minx)
		for (size_t i = nblocks; i < modules.size(); ++i)
			modules[i].w = (modules[i].x-minx)/(maxx-minx);
	if (maxy > miny)
		for (size_t i = nblocks; i < modules.size(); ++i)
			modules[i].h = (modules[i].y-miny)/(maxy-miny);

	if (shape != 0)
	{
		shape->first = maxx-minx;
		shape->second = maxy-miny;
	}

	return nblocks;
}

const cstr_t bm_delimit = " \r\n\t,()=:";

int bm_blocks(const std::string &blocks,
	std::vector<module_t> &modules, map_modules &names)
{
	FILE *fp = fopen(blocks.c_str(), "r");
	if (fp == 0)
	{
		fprintf(stderr, "%s: no blocks file\n", blocks.c_str());
		return 0;
	}

	int nblocks = -1, nterminals = -1, iblock = 0, iterminal = 0;
	bool init = false;

	char linebuf[1000];
	size_t line = 0;
	while (fgets(linebuf, sizeof(linebuf), fp) != 0)
	{
		line++;
		char *temp = skip_space(linebuf, bm_delimit);
		if ((*temp == 0) || (*temp == '#'))
			continue; // empty line and comments

		std::vector<cstr_t> tokens;
		split_tokens(temp, &tokens, bm_delimit);

		if (!init && (nblocks > 0) && (nterminals >= 0))
		{
			init = true;
			modules.resize(nblocks+nterminals);
			iblock = 0;
			iterminal = nblocks;
		}

		if (!init && (tokens.size() == 2)
			&& (strcmp(tokens[0], "NumHardRectilinearBlocks") == 0))
			nblocks = atoi(tokens[1]);
		else if (!init && (tokens.size() == 2)
			&& (strcmp(tokens[0], "NumTerminals") == 0))
			nterminals = atoi(tokens[1]);
		else if (init && (tokens.size() == 11)
			&& (strcmp(tokens[1], "hardrectilinear") == 0)
			&& (strcmp(tokens[2], "4") == 0))
		{
			if (iblock == nblocks)
			{
				fprintf(stderr, "%s:%d: error: extra hard block\n", blocks.c_str(), line);
				fclose(fp);
				return 0;
			}
			
			module_t &module = modules[iblock];
			module.name = tokens[0];
			module.w = atof(tokens[7]);
			module.h = atof(tokens[8]);
			module.x = module.y = 0;
			module.id = 0;

			names[tokens[0]] = &module;
			iblock++;
		}
		else if (init && (tokens.size() == 2)
			&& (strcmp(tokens[1], "terminal") == 0))
		{
			if (iterminal == (int)modules.size())
			{
				fprintf(stderr, "%s:%d: error: extra terminal\n", blocks.c_str(), line);
				fclose(fp);
				return 0;
			}

			module_t &module = modules[iterminal];
			module.name = tokens[0];
			module.w = module.h = 0;
			module.x = module.y = 0;
			module.id = -1;

			names[tokens[0]] = &module;
			iterminal++;
		}
	}

	if (iblock != nblocks)
	{
		fprintf(stderr,
			"%s: error: missing %d hard blocks\n",
			blocks.c_str(), nblocks-iblock);
		fclose(fp);
		return 0;
	}
			
	if (iterminal != (int)modules.size())
	{
		fprintf(stderr,
			"%s: error: missing %d terminals\n",
			blocks.c_str(), modules.size()-iterminal);
		fclose(fp);
		return 0;
	}

	fclose(fp);
	return nblocks;
}

bool bm_nets(const std::string &nets,
	std::vector<signal_t> &signals, map_modules &names,
	std::vector<detailed_signal_t> *detailed)
{
	FILE *fp = fopen(nets.c_str(), "r");
	if (fp == 0)
	{
		fprintf(stderr, "%s: no nets file\n", nets.c_str());
		return false;
	}

	int nnets = -1, npins = -1, inet = 0, ipin = 0, degree = 0;
	bool init = false;

	char linebuf[1000];
	size_t line = 0;
	while (fgets(linebuf, sizeof(linebuf), fp) != 0)
	{
		line++;
		char *temp = skip_space(linebuf, bm_delimit);
		if ((*temp == 0) || (*temp == '#'))
			continue; // empty line and comments

		std::vector<cstr_t> tokens;
		split_tokens(temp, &tokens, bm_delimit);

		if (!init && (nnets > 0) && (npins > 0))
		{
			init = true;
			signals.resize(nnets);
			if (detailed != 0)
				detailed->resize(nnets);
			inet = -1;
			ipin = 0;
		}

		if (!init && (tokens.size() == 2)
			&& (strcmp(tokens[0], "NumNets") == 0))
			nnets = atoi(tokens[1]);
		else if (!init && (tokens.size() == 2)
			&& (strcmp(tokens[0], "NumPins") == 0))
			npins = atoi(tokens[1]);
		else if (init && (tokens.size() == 2)
			&& (strcmp(tokens[0], "NetDegree") == 0))
		{
			if (degree != 0)
			{
				fprintf(stderr,
					"%s:%d: error: missing %d pins",
					nets.c_str(), line, degree);
				fclose(fp);
				return false;
			}
			degree = atoi(tokens[1]);
			++inet;
			if (inet == nnets)
			{
				fprintf(stderr, "%s:%d: error: extra nets\n", nets.c_str(), line);
				fclose(fp);
				return false;
			}
		}
		else if (init && (degree > 0)
			&& (tokens.size() >= 2))
		{
			map_modules::iterator it = names.find(tokens[0]);
			if (it == names.end())
			{
				fprintf(stderr, "%s:%d: error: unknown pin %s\n",
					nets.c_str(), line, tokens[0]);
				fclose(fp);
				return false;
			}
			signals[inet].pins.push_back(it->second);
			it->second->nets.push_back(&signals[inet]);

			if (detailed != 0) // add detailed pin info
			{
				detailed_signal_t::pin_t pin;
				pin.first = it->second;

				if (tokens.size() == 4)
				{
					if ((tokens[2][0] != '%') || (tokens[3][0] != '%'))
					{
						fprintf(stderr, "%s:%d: error: pin position %s\n",
							nets.c_str(), line, tokens[0]);
						fclose(fp);
						return false;
					}
					pin.second.first = atof(tokens[2]+1)/100;
					pin.second.second = atof(tokens[3]+1)/100;
//					printf("%s: %.2f, %.2f\n", tokens[0], pin.second.first, pin.second.second);
				}
				else if (tokens.size() != 2)
				{
					fprintf(stderr, "%s:%d: error: pin %s\n",
						nets.c_str(), line, tokens[0]);
					fclose(fp);
					return false;
				}

				detailed->at(inet).pins.push_back(pin);
			}
			
			--degree;
			++ipin;
		}
	}

	++inet;
	if (inet != nnets)
	{
		fprintf(stderr,
			"%s: error: missing %d nets\n",
			nets.c_str(), nnets-inet);
		fclose(fp);
		return false;
	}
			
	if (ipin != npins)
	{
		fprintf(stderr,
			"%s: error: missing %d pins\n",
			nets.c_str(), npins-ipin);
		fclose(fp);
		return false;
	}

	fclose(fp);
	return true;
}

} // namespace bm



namespace fplan
{

namespace detail
{

class anneal_t
{
protected:
	cost_ptr &cost_;
	perturb_ptr &pt_;
	FILE * const fp_;

	const int N_, k_;

	double best_cost_, prev_cost_;
	int counter_;

	double initial_temperature();

	// return reject rate
	double round(double T, int n, double &avg);

public:
	anneal_t(cost_ptr &cost, perturb_ptr &pt, int N, FILE *fp, int k);

	int go(double r, double term_T, double max_rr);

}; // class anneal_t

} // namespace detail

int anneal(cost_ptr c, perturb_ptr pt, int N, FILE *fp, int k,
	double r, double term_T, double max_reject_rate)
{
	detail::anneal_t an(c, pt, N, fp, k);

	return an.go(r, term_T, max_reject_rate);
}

namespace detail
{

anneal_t::anneal_t(cost_ptr &cost, perturb_ptr &pt, int N, FILE *fp, int k)
	: cost_(cost), pt_(pt), fp_(fp), N_(N), k_(k),
	best_cost_(0), prev_cost_(0), counter_(0)
{
}

int anneal_t::go(double r, double term_T, double max_rr)
{
	FILE *fplog = fopen("log.txt", "w");
	// accept initial solution
	prev_cost_ = best_cost_ = cost_->cost();
	pt_->record_best();
	
	// initial temperature
	counter_ = 0;
	double T = initial_temperature();

	// annealing
	for (int iter = 1;; ++iter)
	{
		
		double avg = 0;
		double rr = round(T, N_*k_, avg);
		
		time_t then1 = time(0);
		
		time_t diff1 = then1 - now_g;
		
		//cout << "Annealing in Progress ... Running since: " << diff1 << "seconds" << endl;
		fprintf(fplog, "iteration: %d T: %f scaled_area: %f rejectrate: %f exec_time: %d \n", iter, T, cost_->cost(),rr, diff1);
	

		if (rr >= max_rr)
		{
			if (rr >= 0.95)
				if (fp_ != 0) fprintf(fp_, "convergent\n");
			break;
		}
		else if (T <= term_T)
		{
			if (fp_ != 0) fprintf(fp_, "cool enough, T=%f\n", T);
			break;
		}

		T *= r;
	}

	pt_->swap_best();

	fclose(fplog);
	return counter_;
}

double anneal_t::initial_temperature()
{
	double avg = 0;
	int uphill = 0;
	for (int j = 0; j < N_; ++j, ++counter_)
	{
		pt_->perturbation(-1);
		
		double cost = cost_->cost();

		if (cost > prev_cost_)
		{
			avg += prev_cost_-cost;
			++uphill;
		}
		else if (cost < best_cost_)
		{
			best_cost_ = cost;
			pt_->record_best();
			/*if (fp_ != 0) fprintf(fp_, "best cost %.3f\n", best_cost_);*/
		}
		prev_cost_ = cost;
	}
	avg /= uphill;

	// initial temperature
	double T = avg/log(0.85);
	/*if (fp_ != 0) fprintf(fp_, "uphill = %d, avg = %e, T = %f\n", uphill, avg, T);*/

	return T;
}

double anneal_t::round(double T, int n, double &avg)
{
	int MT = 0, uphill = 0, reject = 0;
	for (; (uphill < n) && (MT < 2*n); ++MT, ++counter_)
	{
		pt_->perturbation(T);
		double cost = cost_->cost();

		if ((cost < prev_cost_) || double(rand())/RAND_MAX < exp((prev_cost_-cost)/T))
		{
			if (prev_cost_ < cost)
			{
				avg += cost-prev_cost_;
				uphill++;
			}
			if (cost < best_cost_)
			{
				best_cost_ = cost;
				pt_->record_best();
				/*if (fp_ != 0) fprintf(fp_, "best cost %.3f\n", best_cost_);*/
			}
			prev_cost_ = cost;
		}
		else
		{
			reject++;
			pt_->reject();
		}
	}

	double rr = double(reject)/MT;
	/*if (fp_ != 0) fprintf(fp_,
		"uphill = %d, reject = %d, reject_rate = %f\n",
		uphill, reject, rr);*/

	return rr;
}

} // namespace detail

} // namespace fplan

namespace misc
{

using fplan::cost_ptr;

namespace detail
{

typedef double (*func_t)(double, double);

class cost_binary_t : public fplan::detail::cost_i
{
	cost_ptr c_;
	cost_ptr cc_;
	func_t func_;

	virtual double cost() {return func_(c_->cost(), cc_->cost());}
	virtual void release() {delete this;}

	~cost_binary_t() {}

public:
	cost_binary_t(const cost_ptr &c, const cost_ptr &cc, func_t func)
		: c_(c), cc_(cc), func_(func)
		{}
}; // class cost_binary_t

double add(double a, double b) {return a+b;}
double sub(double a, double b) {return a-b;}
double mul(double a, double b) {return a*b;}
double div(double a, double b) {return a/b;}

class double_cost_t : public fplan::detail::cost_i
{
	double d_;
	
	virtual double cost() {return d_;}
	virtual void release() {delete this;}

	~double_cost_t() {}

public:
	double_cost_t(double d) : d_(d) {}
}; // class double_cost_t

} // namespace detail

cost_ptr double_cost(double d)
{
	return cost_ptr(new detail::double_cost_t(d));
}


cost_ptr operator + (const cost_ptr &c, const cost_ptr &cc)
{
	return cost_ptr(new detail::cost_binary_t(c, cc, detail::add));
}

cost_ptr operator - (const cost_ptr &c, const cost_ptr &cc)
{
	return cost_ptr(new detail::cost_binary_t(c, cc, detail::sub));
}

cost_ptr operator * (const cost_ptr &c, const cost_ptr &cc)
{
	return cost_ptr(new detail::cost_binary_t(c, cc, detail::mul));
}

cost_ptr operator / (const cost_ptr &c, const cost_ptr &cc)
{
	return cost_ptr(new detail::cost_binary_t(c, cc, detail::div));
}

extern "C" double pow(double, double);
cost_ptr pow(const cost_ptr &c, const cost_ptr &cc)
{
	return cost_ptr(new detail::cost_binary_t(c, cc, pow));
}

cost_ptr operator + (double d, const cost_ptr &c)
{
	return double_cost(d)+c;
}

cost_ptr operator + (const cost_ptr &c, double d)
{
	return c+double_cost(d);
}

cost_ptr operator - (double d, const cost_ptr &c)
{
	return double_cost(d)-c;
}

cost_ptr operator - (const cost_ptr &c, double d)
{
	return c-double_cost(d);
}

cost_ptr operator * (double d, const cost_ptr &c)
{
	return double_cost(d)*c;
}

cost_ptr operator * (const cost_ptr &c, double d)
{
	return c*double_cost(d);
}

cost_ptr operator / (double d, const cost_ptr &c)
{
	return double_cost(d)/c;
}

cost_ptr operator / (const cost_ptr &c, double d)
{
	return c/double_cost(d);
}

cost_ptr pow(double d, const cost_ptr &c)
{
	return pow(double_cost(d), c);
}

cost_ptr pow(const cost_ptr &c, double d)
{
	return pow(c, double_cost(d));
}

} // namespace misc


namespace fplan
{

namespace detail
{

class area_cost_t : public cost_i
{
	fplan_i *fp_;

	virtual double cost()
	{
		misc::dbl_pair shape = fp_->get_shape(); // eval
		return shape.first*shape.second;
	}

	virtual void release() {delete this;}

public:
	area_cost_t(fplan_i *fp) : fp_(fp) {}
}; // class area_cost_t

class fixed_cost_t : public cost_i
{
	fplan_i *fp_;
	double w_, h_;
	bool rotate_;

	virtual double cost()
	{
		misc::dbl_pair shape = fp_->get_shape(); // eval
			
		double w = shape.first, h = shape.second;

		if (rotate_ && (w < h))
			std::swap(w, h);

		w = (w > w_)? w-w_: 0;
		h = (h > h_)? h-h_: 0;

		return exp(h/h_+w/w_);
	}

	virtual void release() {delete this;}

public:
	fixed_cost_t(fplan_i *fp, double w, double h, bool rotate)
		: fp_(fp), w_(w), h_(h), rotate_(rotate)
	{
		if (rotate_ && (w_ < h_))
			std::swap(w_, h_);
	}
}; // class fixed_cost_t


class mnets_cost_t : public cost_i
{
	fplan_i *fp_;

	const module_nets_t &mnets_;

	typedef double (module_nets_t::*WL_func_t)(const misc::dbl_pair centers[]) const;
	WL_func_t wl_func_;

	const size_t nm_, np_;
	const misc::dbl_pair *pins_;

	std::vector<misc::dbl_pair> centers_;

	virtual double cost()
	{
		// you should create a new cost_ptr if the floorplan changed
		assert(nm_ == fp_->get_num_modules());

		// get shape and centers
		misc::dbl_pair shape = fp_->get_centers(0, &centers_[0]);

		// compute pin positions
		for (size_t i = 0; i < np_; ++i)
			centers_[i+nm_] = pins_[i]*shape;
			

		return (mnets_.*wl_func_)(&centers_[0]);
	}

	virtual void release() {delete this;}

public:
	mnets_cost_t(fplan_i *fp, const module_nets_t &mnets,
		WL_func_t wl_func, size_t np, const misc::dbl_pair pins[])
		: fp_(fp), mnets_(mnets), wl_func_(wl_func),
		nm_(fp->get_num_modules()), np_(np), pins_(pins),
		centers_(fp->get_num_modules()+np)
		{}
}; // class mnets_cost_t

} // namespace detail

cost_ptr fplan_i::area_cost()
{
	return cost_ptr(new detail::area_cost_t(this));
}


cost_ptr fplan_i::fixed_cost(double w, double h, bool rotate)
{
	return cost_ptr(new detail::fixed_cost_t(this, w, h, rotate));
}

cost_ptr fplan_i::HPWL_cost(const module_nets_t &mnets,
	size_t np, const misc::dbl_pair pins[])
{
	return cost_ptr(new detail::mnets_cost_t(this, mnets, &module_nets_t::HPWL, np, pins));
}


} // namespace fplan


namespace fplan
{

class acg_perturb_t : public detail::perturb_i
{
protected:
	acg_fp_t *fp_;
	double ign_, exchangeP_;

	acg::acg_base_t best_, temp_;
	std::vector<size_t> best_ign_;

	int choice_;
	size_t n1_, n2_;

	virtual void perturbation(double T);
	virtual void reject();
	virtual void record_best();
	virtual void swap_best();

	// you should use the shared_ptr interface
	virtual void release() {delete this;}
	~acg_perturb_t() {}

public:
	acg_perturb_t(acg_fp_t *fp, double ign, double exchangeP);

}; // acg_perturb_t

acg_fp_t::acg_fp_t(size_t n, misc::dbl_pair shapes[], size_t ids[], size_t ign[])
	: rcg_fp_t(n, shapes, ids, ign, false), reduced_(true)
{
	acg::racg_build(acg_base_);
}

acg_fp_t::acg_fp_t(acg_fp_t &r)
	: reduced_(true)
{
	r.require_racg();

	acg_base_ = r.acg_base_;
	ign_ = r.ign_;
	shape_ = r.shape_;
	mapped_ = r.mapped_;
}

acg_fp_t::~acg_fp_t()
{
}

void acg_fp_t::swap(acg_fp_t &r)
{
	acg_base_.swap(r.acg_base_);
	ign_.swap(r.ign_);
	std::swap(shape_, r.shape_);
	std::swap(mapped_, r.mapped_);

	ctx_.swap(r.ctx_);
	std::swap(reduced_, r.reduced_);
}

acg_fp_t &acg_fp_t::operator = (acg_fp_t r)
{
	swap(r);

	return *this;
}

void acg_fp_t::require_acg()
{
	if (reduced_)
	{
		acg::racg_to_acg(acg_base_, ctx_);
		reduced_ = false;
	}
}

void acg_fp_t::require_racg()
{
	if (!reduced_)
	{
		acg::acg_to_racg(acg_base_, ctx_);
		reduced_ = true;
	}
}

void acg_fp_t::do_packing()
{
	if (!mapped_)
	{
		require_acg();

		shape_ = acg::packing_longest_path(acg_base_);

		mapped_ = true;
	}
}

void acg_fp_t::compute_bb()
{
	do_packing();

	require_acg();

	acg::longest_path_bb(acg_base_, shape_);
}

perturb_ptr acg_fp_t::create_perturbation(double ign, double exchangeP)
{
	return perturb_ptr(new acg_perturb_t(this, ign, exchangeP));
}

acg_perturb_t::acg_perturb_t(acg_fp_t *fp, double ign, double exchangeP)
	: fp_(fp),
	ign_(ign), exchangeP_(exchangeP),
	choice_(-1), n1_(0), n2_(0)
{
}

void acg_perturb_t::perturbation(double T)
{
	const size_t n = fp_->acg_base_.num_nodes();

	assert(n != 0);

	double P = double(rand())/RAND_MAX;

	fp_->mapped_ = false;

	if ((P < ign_) || (n == 1))
	{
		choice_ = 0;
		
		n1_ = rand()%n;

		fp_->acg_base_.node_rotate(n1_);
		fp_->ign_[n1_] = 1-fp_->ign_[n1_];
	}
	else if (P < ign_+exchangeP_)
	{
		choice_ = 1;
	
		n1_ = rand()%fp_->acg_base_.num_nodes();
	
		do
		{
			n2_ = (int)rand()%n;
		}
		while (n2_ == n1_);

		fp_->acg_base_.node_swap(n1_, n2_);
		std::swap(fp_->ign_[n1_], fp_->ign_[n2_]);
	}
	else
	{
		choice_ = 2;

		P -= ign_+exchangeP_;
		double ref = 1-ign_-exchangeP_;

		fp_->require_racg();
		temp_ = fp_->acg_base_; // force racg and copy
	
		acg::pos_t pos = (rand()%2 == 0)? acg::POS_IN: acg::POS_OUT;
		if (P < 0.34*ref)
			acg::racg_switch_adjacent(fp_->acg_base_, n1_, pos);
		else if (P < 0.67*ref)
			acg::racg_switch_alter(fp_->acg_base_, n1_, pos);
		else
			acg::racg_switch_last(fp_->acg_base_, n1_, pos);
	}
}

void acg_perturb_t::reject()
{
	fp_->mapped_ = false;

	switch (choice_)
	{
	case 0:
		fp_->acg_base_.node_rotate(n1_);
		fp_->ign_[n1_] = 1-fp_->ign_[n1_];
		break;
	case 1:
		fp_->acg_base_.node_swap(n1_, n2_);
		std::swap(fp_->ign_[n1_], fp_->ign_[n2_]);
		break;
	case 2:
		fp_->acg_base_.swap(temp_);
		fp_->reduced_ = true;
		fp_->ctx_.clear();
		break;
	default:
		assert(false);
	}
}

void acg_perturb_t::record_best()
{
	fp_->require_racg();

	best_ = fp_->acg_base_;
	best_ign_ = fp_->ign_;
}

void acg_perturb_t::swap_best()
{
	best_.swap(fp_->acg_base_);
	best_ign_.swap(fp_->ign_);
	fp_->mapped_ = false;

	fp_->ctx_.clear();
	fp_->reduced_ = true;
}

} // namespace fplan


namespace acg
{

namespace detail
{

color_t node_color(const node_t *node, pos_t pos); // from acg_racg.cpp

void rcg_remove_transitive(acg_dev_t &base,
	node_t *nA, node_t *nB, pos_t pos, color_t c);

void rcg_add_original(acg_dev_t &base,
	node_t *nA, node_t *nB, pos_t pos, color_t c, edge_t *eAB);

} // namespace detail

void rcg_build(acg_base_t &base)
{
	racg_build(base);

	racg_acg_context_t ctx;

	racg_to_acg(base, ctx);
}

void rcg_swap(acg_base_t &b, size_t where)
{
	using detail::node_color;
	using detail::rcg_remove_transitive;
	using detail::rcg_add_original;

	if (where == b.num_nodes()-1)
		return;

	acg_dev_t base(b);

	node_t *node0 = base[where];
	node_t *node1 = base[where+1];

	color_t c = node_color(node0, POS_IN);
	assert(c == node_color(node1, POS_OUT));
	assert(node0 == node1->nearest_node(c, POS_OUT));
	assert(node1 == node0->nearest_node(c, POS_IN));

	rcg_remove_transitive(base, node0, node1, POS_OUT, other(c));
	rcg_remove_transitive(base, node1, node0, POS_IN, other(c));

	edge_t *change = node0->edges(c)->front(POS_IN);
	assert(change->node(POS_OUT) == node1);

	rcg_add_original(base, node0, node1, POS_OUT, c, change);
	rcg_add_original(base, node1, node0, POS_IN, c, change);

	// re-assign change
	edge_list_t::erase(change, POS_IN);
	edge_list_t::erase(change, POS_OUT);
	node0->edges(other(c))->push_front(change, POS_OUT);
	node1->edges(other(c))->push_front(change, POS_IN);
	
	// swap the node
	std::swap(base[where], base[where+1]);
}

namespace detail
{

/*
 * A--B      A
 * :  :      :
 * ?  :  =>  B
 * :  :      :
 * CCCC      C
 */
void rcg_remove_transitive(acg_dev_t &base,
	node_t *nA, node_t *nB, pos_t pos, color_t c)
{
	assert(nB->nearest_node(other(c), pos) == nA);
	assert(nA->nearest_node(other(c), other(pos)) == nB);

	for (edge_list_t::range_t rB = nB->edges(c)->all(pos);
		rB.first != rB.second; rB.first = rB.first->next(pos))
	{
		edge_t *eCB = rB.first;
		node_t *nC = eCB->node(other(pos)); (void)nC;

		// there is a c path from C to A
		assert(eCB != nC->edges(c)->front(other(pos)));

		edge_t *eCA = eCB->prev(other(pos));
		if (eCA->node(pos) != nA) // possiblely not eCA
			break;

		// remove eCA
		edge_list_t::erase(eCA, POS_IN);
		edge_list_t::erase(eCA, POS_OUT);
		base.free_edge(eCA);
	}
}

/*
 * C--A--B      C-----A
 * C  :  B  =>  C     :
 * C--D??B      C--D  B
 */
void rcg_add_original(acg_dev_t &base,
	node_t *nA, node_t *nB, pos_t pos, color_t c, edge_t *eAB)
{
	assert(eAB->node(pos) == nB);
	assert(eAB->node(other(pos)) == nA);

	for (edge_list_t::range_t rA = nA->edges(c)->all(pos);
		rA.first != rA.second; rA.first = rA.first->next(pos))
	{
		edge_t *eCA = rA.first;
		node_t *nC = eCA->node(other(pos));

		if (eCA != nC->edges(c)->front(other(pos)))
		{
			node_t *nD = eCA->prev(other(pos))->node(pos);

			// oc path from D to B iff eDB is oc

			/* this won't work since eDA is erased if eDB is there :(
			// looking for eDA
			edge_list_t::range_t rD = nD->edges(other(c))->all(other(pos));
			for (;; rD.first = rD.first->next(other(pos)))
			{
				assert(rD.first != rD.second);

				if (rD.first->node(pos) == nA)
					break;
			}
			edge_t *eDA = rD.first;

			// if eDB is oc, it must be the one after eDA
			// should get out of the loop if eDB is not there
			// since there must be a c path from D to B
			if (eDA->next(other(pos))->node(pos) != nB)
				break;
			*/

			edge_list_t::range_t rD = nD->edges(other(c))->all(other(pos));
			for (; rD.first != rD.second; rD.first = rD.first->next(other(pos)))
				if (rD.first->node(pos) == nB)
					break;
			if (rD.first == rD.second) // eDB is not there
				break;
		}

		// add edge from C to B
		edge_t *eCB = base.alloc_edge();
		edge_list_t::insert(eCA, eCB, other(pos));
		edge_list_t::insert(eAB, eCB, pos);
	}
}

} // namespace detail

} // namespace acg


namespace fplan
{

class rcg_perturb_t : public detail::perturb_i
{
protected:
	rcg_fp_t *fp_;
	double ign_, exchangeP_;

	acg::acg_base_t best_;
	std::vector<size_t> best_ign_;

	int choice_;
	size_t n1_, n2_;

	virtual void perturbation(double T);
	virtual void reject();
	virtual void record_best();
	virtual void swap_best();

	// you should use the shared_ptr interface
	virtual void release() {delete this;}
	~rcg_perturb_t() {}

public:
	rcg_perturb_t(rcg_fp_t *fp, double ign, double exchangeP);

}; // class rcg_perturb_t


rcg_fp_t::rcg_fp_t()
	: mapped_(false)
{
}

rcg_fp_t::rcg_fp_t(size_t n, misc::dbl_pair shapes[], size_t ids[],
	size_t ign[], bool build_rcg)
	: acg_base_(n, shapes, ids), ign_(n, 0), mapped_(false)
{
	if (ign != 0)
		for (size_t i = 0; i < n; ++i)
		{
			ign_[i] = ign[i];
			if (ign[i] != 0)
				acg_base_.node_rotate(i);
		}

	if (build_rcg)
		acg::rcg_build(acg_base_);
}

rcg_fp_t::~rcg_fp_t()
{
}

void rcg_fp_t::swap(rcg_fp_t &r)
{
	acg_base_.swap(r.acg_base_);
	ign_.swap(r.ign_);
	std::swap(shape_, r.shape_);
	std::swap(mapped_, r.mapped_);
}

perturb_ptr rcg_fp_t::create_perturbation(double ign, double exchangeP)
{
	return perturb_ptr(new rcg_perturb_t(this, ign, exchangeP));
}

void rcg_fp_t::do_packing()
{
	if (!mapped_)
	{
		shape_ = acg::packing_longest_path(acg_base_);

		mapped_ = true;
	}
}

void rcg_fp_t::compute_bb()
{
	do_packing();

	acg::longest_path_bb(acg_base_, shape_);
}

size_t rcg_fp_t::get_num_modules()
{
	return acg_base_.num_nodes();
}

misc::dbl_pair rcg_fp_t::get_shape()
{
	do_packing();

	return shape_;
}

misc::dbl_pair rcg_fp_t::get_centers(size_t ids[], misc::dbl_pair centers[])
{
	do_packing();
	
	if (ids == 0)
		
		for (size_t i = 0; i < acg_base_.num_nodes(); ++i)
			centers[acg_base_[i]->m_id] = acg_base_[i]->center(acg::POS_IN);
	else
		for (size_t i = 0; i < acg_base_.num_nodes(); ++i)
		{
			ids[i] = acg_base_[i]->m_id;
			centers[i] = acg_base_[i]->center(acg::POS_IN);
		}

	return shape_;
}

misc::dbl_pair rcg_fp_t::get_physical_info(
	size_t ids[], misc::dbl_pair centers[],
	misc::dbl_pair lls[], misc::dbl_pair urs[],
	size_t ign[], misc::dbl_pair ign_shapes[])
{
	compute_bb();

	for (size_t i = 0; i < acg_base_.num_nodes(); ++i)
	{
		size_t index = (ids == 0)? acg_base_[i]->m_id: i;

		acg::node_t *node = acg_base_[i];
		if (centers != 0)
			centers[index] = node->center(acg::POS_IN);
		if (lls != 0)
			lls[index] = node->bb(acg::POS_IN);
		if (urs != 0)
			urs[index] = node->bb(acg::POS_OUT);
		if (ign_shapes != 0)
			ign_shapes[index] = node->width();

		if (ids != 0)
			ids[index] = acg_base_[i]->m_id;

		if (ign != 0)
			ign[index] = ign_[i];
	}

	return shape_;
}


bool rcg_fp_t::dump(const std::string &name, misc::cstr_t module_names[])
{
	compute_bb();

	FILE *fp = fopen(name.c_str(), "w");
	if (fp == 0)
		return false;

	acg_base_.xml_write(fp);

	fclose(fp);

	return true;
}

rcg_perturb_t::rcg_perturb_t(rcg_fp_t *fp, double ign, double exchangeP)
	: fp_(fp),
	ign_(ign), exchangeP_(exchangeP),
	choice_(-1), n1_(0), n2_(0)
{
}

void rcg_perturb_t::perturbation(double)
{
	const size_t n = fp_->acg_base_.num_nodes();

	assert(n != 0);

	double P = double(rand())/RAND_MAX;

	fp_->mapped_ = false;

	if ((P < ign_) || (n == 1))
	{
		choice_ = 0;
		
		n1_ = rand()%n;

		fp_->acg_base_.node_rotate(n1_);
		fp_->ign_[n1_] = 1-fp_->ign_[n1_];
	}
	else if (P < ign_+exchangeP_)
	{
		choice_ = 1;
	
		n1_ = rand()%fp_->acg_base_.num_nodes();
	
		do
		{
			n2_ = (int)rand()%n;
		}
		while (n2_ == n1_);

		fp_->acg_base_.node_swap(n1_, n2_);
		std::swap(fp_->ign_[n1_], fp_->ign_[n2_]);
	}
	else
	{
		choice_ = 2;

		n1_ = rand()%(n-1);

		acg::rcg_swap(fp_->acg_base_, n1_);
		std::swap(fp_->ign_[n1_], fp_->ign_[n1_+1]);
	}
}

void rcg_perturb_t::reject()
{
	fp_->mapped_ = false;

	switch (choice_)
	{
	case 0:
		fp_->acg_base_.node_rotate(n1_);
		fp_->ign_[n1_] = 1-fp_->ign_[n1_];
		break;
	case 1:
		fp_->acg_base_.node_swap(n1_, n2_);
		std::swap(fp_->ign_[n1_], fp_->ign_[n2_]);
		break;
	case 2:
		acg::rcg_swap(fp_->acg_base_, n1_);
		std::swap(fp_->ign_[n1_], fp_->ign_[n1_+1]);
		break;
	default:
		assert(false);
	}
}

void rcg_perturb_t::record_best()
{
	best_ = fp_->acg_base_;
	best_ign_ = fp_->ign_;
}

void rcg_perturb_t::swap_best()
{
	best_.swap(fp_->acg_base_);
	best_ign_.swap(fp_->ign_);
	fp_->mapped_ = false;
}

} // namespace fplan


namespace acg
{

bool edge_t::valid() const
{
	for (int i = 0; i < 2; ++i)
	{
		pos_t pos = pos_t(i);

		if (next(pos)->prev(pos) != this)
			return false;
		if (prev(pos)->next(pos) != this)
			return false;

		if (node(pos) != next(pos)->node(pos))
			return false;
		if (node(pos) != prev(pos)->node(pos))
			return false;
	}

	return true;
}

edge_t *edge_list_t::pop_front(pos_t pos)
{
	edge_t *tmp = front(pos);
	erase(tmp, pos);
	return tmp;
}

edge_t *edge_list_t::pop_back(pos_t pos)
{
	edge_t *tmp = back(pos);
	erase(tmp, pos);
	return tmp;
}

void edge_list_t::push_front(edge_t *edge, pos_t pos)
{
	insert(head_.next(pos), edge, pos);
}

void edge_list_t::push_back(edge_t *edge, pos_t pos)
{
	insert(&head_, edge, pos);
}

void edge_list_t::insert(edge_t *after, edge_t *edge, pos_t pos)
{
	assert(after->valid());
	assert(edge->valid());

	assert(edge->node_[pos] == 0);

	edge->node_[pos] = after->node_[pos];

	edge->prev_[pos] = after->prev_[pos];
	edge->next_[pos] = after;
	after->prev_[pos]->next_[pos] = edge;
	after->prev_[pos] = edge;
}

void edge_list_t::erase(edge_t *edge, pos_t pos)
{
	assert(edge->valid());

	assert(edge->node_[pos] != 0);
	assert(edge != edge->next_[pos]);

	edge->node_[pos] = 0;

	edge->prev_[pos]->next_[pos] = edge->next_[pos];
	edge->next_[pos]->prev_[pos] = edge->prev_[pos];
	edge->prev_[pos] = edge->next_[pos] = edge;
}

void edge_list_t::splice(edge_t *after, edge_t *first, edge_t *last, pos_t pos)
{
	assert(after->valid());
	assert(first->valid());
	assert(last->valid());

	assert(after->node_[pos] != 0);
	assert(first->node(pos) == last->node(pos));

	if (first == last)
		return;
	
	// worth splice, relink
	edge_t *prev_first = first->prev_[pos];
	edge_t *prev_after = after->prev_[pos];
	edge_t *prev_last = last->prev_[pos];

	// prev_after<->after, prev_first<->first...prev_last<->last
	// prev_after<->first...prev_last<->after, prev_first<->last

	prev_after->next_[pos] = first;
	first->prev_[pos] = prev_after;

	prev_last->next_[pos] = after;
	after->prev_[pos] = prev_last;

	prev_first->next_[pos] = last;
	last->prev_[pos] = prev_first;

	// reset node
	for (; first != after; first = first->next(pos))
		first->node_[pos] = after->node_[pos];
}

bool edge_list_t::valid() const
{
	if (!head_.valid())
		return false;

	for (int i = 0; i < 2; ++i)
	{
		pos_t pos = pos_t(i);

		for (const_range_t r = all(pos);
			r.first != r.second;
			r.first = r.first->next(pos))
			if (!r.first->valid())
				return false;
	}

	if (head_.node(POS_IN) != head_.node(POS_OUT))
		return false;

	return true;
}

node_t::node_t() : m_id(0), m_data(0)
{
	edges_[0].init(this);
	edges_[1].init(this);

	width_[0] = width_[1] = 0;
	bb_[0][0] = bb_[0][1] = bb_[1][0] = bb_[1][1] = 0;
}

node_t::node_t(double w, double h, size_t id)
	: m_id(id), m_data(0)
{
	edges_[0].init(this);
	edges_[1].init(this);

	width_[0] = w; width_[1] = h;
	bb_[0][0] = bb_[0][1] = bb_[1][0] = bb_[1][1] = 0;
}

bool node_t::valid() const
{
	if (!edges(COLOR_RED)->valid())
		return false;
	if (!edges(COLOR_BLACK)->valid())
		return false;

	if (edges(COLOR_RED)->node() != this)
		return false;
	if (edges(COLOR_BLACK)->node() != this)
		return false;

	return true;
}

void node_t::copy_content(const node_t &src)
{
	m_id = src.m_id;
	memcpy(width_, src.width_, sizeof(width_));
	memcpy(bb_, src.bb_, sizeof(bb_));
}

void node_t::swap_content(node_t &src)
{
	std::swap(m_id, src.m_id);
	std::swap(width_[0], src.width_[0]);
	std::swap(width_[1], src.width_[1]);
	std::swap(bb_[0][0], src.bb_[0][0]);
	std::swap(bb_[0][1], src.bb_[0][1]);
	std::swap(bb_[1][0], src.bb_[1][0]);
	std::swap(bb_[1][1], src.bb_[1][1]);
}

acg_base_t::acg_base_t()
	: n_edges_(0)
{
}

acg_base_t::acg_base_t(size_t n, misc::dbl_pair shapes[], size_t ids[])
	: n_edges_(0), nodes_(n)
{
	acg_dev_t dev(*this);

	for (size_t i = 0; i < n; ++i)
		nodes_[i] = dev.alloc_node(
			shapes[i].first, shapes[i].second,
			(ids == 0)? i: ids[i]);
}

acg_base_t::acg_base_t(const acg_base_t &r)
	: n_edges_(0), nodes_(r.num_nodes())
{
	acg_dev_t dev(*this);

	for (size_t i = 0; i < r.num_nodes(); ++i)
	{
		nodes_[i] = dev.alloc_node();
		nodes_[i]->copy_content(*r[i]);
	}

	dev.duplicate_edges(r);
}

acg_base_t::~acg_base_t()
{
	clear();
}

void acg_base_t::clear()
{
	acg_dev_t dev(*this);

	dev.free_all_edges();

	for (size_t i = 0; i < num_nodes(); ++i)
		dev.free_node(nodes_[i]);
	nodes_.clear();
}

void acg_base_t::swap(acg_base_t &r)
{
	nodes_.swap(r.nodes_);
	std::swap(n_edges_, r.n_edges_);
}

void acg_base_t::node_swap(size_t i, size_t j)
{
	at(i)->swap_content(*at(j));
}

void acg_base_t::node_rotate(size_t where)
{
	node_t *node = at(where);
	std::swap(node->width(COLOR_RED), node->width(COLOR_BLACK));
}

void acg_base_t::node_resize(size_t where, double w, double h)
{
	// could be useful for soft modules
	node_t *node = at(where);
	node->width(COLOR_RED) = w;
	node->width(COLOR_BLACK) = h;
}

namespace detail
{

void xml_write_node_edges(FILE *fp, const node_t *node, color_t c)
{
	const char *color = (c == COLOR_RED)? "red": "black";

	for (edge_list_t::const_range_t r = node->edges(c)->all(POS_IN);
		r.first != r.second; r.first = r.first->next(POS_IN))
		fprintf(fp,
			"    <edge color='%s' from='%p' to='%p' />\n",
			color, node, r.first->node(POS_OUT));
}

} // namespace detail

void acg_base_t::xml_write(FILE *fp, const char *type)
{
	using detail::xml_write_node_edges;

	const size_t n = num_nodes();
	
	FILE *fpmoduleinfo = fopen("moduleinfo.txt", "w");
	
	for (size_t i = 0; i < n; ++i)
	{
		fprintf(fp,
			"sb%d %f %f",
			at(i)->m_id, 
			at(i)->bb(acg::COLOR_RED, acg::POS_IN),
			at(i)->bb(acg::COLOR_BLACK, acg::POS_IN));
		fprintf(fp,"\n");
		
		fprintf(fpmoduleinfo, "%d\t%f\t%f\t%f\t%f\n",
			at(i)->m_id, 
			at(i)->width(acg::COLOR_RED),
			at(i)->width(acg::COLOR_BLACK),
			at(i)->bb(acg::COLOR_RED, acg::POS_OUT),
			at(i)->bb(acg::COLOR_BLACK, acg::POS_OUT));
		
			
	}
	fclose(fpmoduleinfo);
	
}

} // namespace acg


namespace acg
{

namespace detail
{

typedef misc::malloc_pool<sizeof(edge_t)> edge_pool_t;
typedef misc::malloc_pool<sizeof(node_t)> node_pool_t;

// singletons for pools
edge_pool_t *get_edge_pool()
{
	static edge_pool_t pool;
	return &pool;
}

node_pool_t *get_node_pool()
{
	static node_pool_t pool;
	return &pool;
}

} // namespace detail

//Diagnostic Purposes
edge_t *acg_dev_t::alloc_edge()
{
	++n_edges_;

	// there will be a leakage for ctor exceptions
	// however, edge_t is so simple that it will never
	// generate an exception... which applies to node_t too.
	return new (detail::get_edge_pool()->alloc()) edge_t;
}

void acg_dev_t::free_edge(edge_t *edge)
{
	--n_edges_;

	edge->~edge_t();

	detail::get_edge_pool()->free(edge);
}

node_t *acg_dev_t::alloc_node()
{
	return new (detail::get_node_pool()->alloc()) node_t;
}

node_t *acg_dev_t::alloc_node(double w, double h, size_t id)
{
	return new (detail::get_node_pool()->alloc()) node_t(w, h, id);
}

void acg_dev_t::free_node(node_t *node)
{
	node->~node_t();

	detail::get_node_pool()->free(node);
}

void acg_dev_t::insert(size_t where, node_t *node)
{
	assert(where <= nodes_.size());

	nodes_.insert(nodes_.begin()+where, node);
}

void acg_dev_t::erase(size_t where)
{
	assert(where < nodes_.size());

	assert(nodes_[where]->edges(COLOR_RED)->empty(POS_IN));
	assert(nodes_[where]->edges(COLOR_RED)->empty(POS_OUT));
	assert(nodes_[where]->edges(COLOR_BLACK)->empty(POS_IN));
	assert(nodes_[where]->edges(COLOR_BLACK)->empty(POS_OUT));

	nodes_.erase(nodes_.begin()+where);
}

void acg_dev_t::push(node_t *node, pos_t pos)
{
	if (pos == POS_IN)
		insert(0, node);
	else
		nodes_.push_back(node);
}

void acg_dev_t::pop(pos_t pos)
{
	if (pos == POS_IN)
		erase(0);
	else
	{
		assert(!nodes_.empty());

		assert(nodes_.back()->edges(COLOR_RED)->empty(POS_IN));
		assert(nodes_.back()->edges(COLOR_RED)->empty(POS_OUT));
		assert(nodes_.back()->edges(COLOR_BLACK)->empty(POS_IN));
		assert(nodes_.back()->edges(COLOR_BLACK)->empty(POS_OUT));
		
		nodes_.pop_back();
	}
}

node_t *acg_dev_t::front(pos_t pos)
{
	assert(num_nodes() > 0);

	if (pos == POS_IN)
		return nodes_.front();
	else
		return nodes_.back();
}

void acg_dev_t::free_edges(edge_t *first, edge_t *last, pos_t pos)
{
	assert(first->valid());
	assert(last->valid());

	assert(first->node(pos) == last->node(pos));

	for (; first != last;)
	{
		edge_t *tmp = first; first = first->next(pos);

		edge_list_t::erase(tmp, POS_IN);
		edge_list_t::erase(tmp, POS_OUT);

		free_edge(tmp);
	}
}

void acg_dev_t::free_edges(node_t *node, pos_t pos)
{
	free_edges(node, COLOR_RED, pos);
	free_edges(node, COLOR_BLACK, pos);
}

void acg_dev_t::free_edges(node_t *node, color_t c)
{
	free_edges(node, c, POS_IN);
	free_edges(node, c, POS_OUT);
}

void acg_dev_t::free_edges(node_t *node, color_t c, pos_t pos)
{
	edge_list_t *el = node->edges(c);

	free_edges(el->begin(pos), el->end(pos), pos);
}

void acg_dev_t::free_all_edges()
{
	if (n_edges_ != 0)
		for (size_t i = 0; i < num_nodes(); ++i)
			free_edges(nodes_[i], POS_IN);
	assert(n_edges_ == 0);
}

void acg_dev_t::duplicate_edges(const acg_base_t &src)
{
	assert(num_nodes() == src.num_nodes());

	const size_t n = num_nodes();

	// clear all the dest nodes first
	for (size_t i = 0; i < n; ++i)
	{
		assert(at(i)->valid());
		assert(src[i]->valid());

		assert(at(i)->edges(COLOR_RED)->empty(POS_IN));
		assert(at(i)->edges(COLOR_RED)->empty(POS_OUT));
		assert(at(i)->edges(COLOR_BLACK)->empty(POS_IN));
		assert(at(i)->edges(COLOR_BLACK)->empty(POS_OUT));

		src[i]->m_data = at(i);
	}

	// create edges and link them at POS_IN side
	for (size_t from = 0; from < n; ++from)
	{
		dup_edges_POS_IN(COLOR_RED, src[from]);
		dup_edges_POS_IN(COLOR_BLACK, src[from]);
	}

	// link edges at POS_OUT side
	for (size_t to = 0; to < n; ++to)
	{
		dup_edges_POS_OUT(COLOR_RED, src[to]);
		dup_edges_POS_OUT(COLOR_BLACK, src[to]);
	}

	assert(num_edges() == src.num_edges());
}


void acg_dev_t::update_order(const size_t order[])
{
	const size_t n = num_nodes();
	std::vector<node_t *> nodes(n);

	for (size_t i = 0; i < n; ++i)
		nodes[i] = at(order[i]);

	nodes_.swap(nodes);
}

void acg_dev_t::dup_edges_POS_IN(color_t c, const node_t *src)
{
	for (edge_list_t::const_range_t r = src->edges(c)->all(POS_IN);
		r.first != r.second; r.first = r.first->next(POS_IN))
	{
		edge_t *new_edge = alloc_edge();

		r.first->m_data = new_edge;

		reinterpret_cast<node_t *>(src->m_data)
			->edges(c)->push_back(new_edge, POS_IN);
	}
}

void acg_dev_t::dup_edges_POS_OUT(color_t c, const node_t *src)
{
	for (edge_list_t::const_range_t r = src->edges(c)->all(POS_OUT);
		r.first != r.second; r.first = r.first->next(POS_OUT))
	{
		edge_t *new_edge = reinterpret_cast<edge_t *>(r.first->m_data);

		reinterpret_cast<node_t *>(src->m_data)
			->edges(c)->push_back(new_edge, POS_OUT);
	}
}


} // namespace acg



namespace acg
{

namespace detail
{

void racg_create_edges(acg_dev_t &base, node_t *node, node_t *next, color_t c, pos_t pos);

edge_t *racg_relink_remove(acg_dev_t &base, edge_list_t *edges, edge_t *AB, pos_t pos);

color_t node_color(const node_t *node, pos_t pos);

void racg_erase_case1(acg_dev_t &base, node_t *node, color_t c);
void racg_erase_case2(acg_dev_t &base, node_t *node, color_t out);

} // namespace detail

void racg_build(acg_base_t &b)
{
	using detail::racg_create_edges;

	acg_dev_t base(b);

	for (size_t i = 1; i < base.num_nodes(); ++i)
		racg_create_edges(base, base[i], base[i-1], color_t(rand()%2), POS_OUT);
}

void racg_push(acg_base_t &b,
	double w, double h, size_t id,
	color_t c, pos_t pos)
{
	using detail::racg_create_edges;

	acg_dev_t base(b);

	node_t *node = base.alloc_node(w, h, id);

	if (!base.empty())
		racg_create_edges(base, node, base.front(pos), c, pos);

	base.push(node, pos);
}

void racg_pop(acg_base_t &b, pos_t pos)
{
	acg_dev_t base(b);

	node_t *tmp = base.front(pos);

	base.free_edges(tmp, pos);

	base.pop(pos);

	base.free_node(tmp);
}

/*
 * refer to erase case 1
 *
 * X--A--B--C
 * X  :  :  :
 * X--DDDDDDD--Y <= X--A--B--CCCCCCC
 *    :  :  :  Y       :  :  :  :  :
 *    E--F--G--Y       EEEEEEE--F--G--Y
 */
void racg_insert(acg_base_t &b, size_t where,
	double w, double h, size_t id)
{
	using detail::node_color;

	if (where == 0)
	{
		racg_push(b, w, h, id, color_t(rand()%2), POS_IN);
		return;
	}
	if (where == b.num_nodes())
	{
		racg_push(b, w, h, id, color_t(rand()%2), POS_OUT);
		return;
	}

	acg_dev_t base(b);

	node_t *node = base.alloc_node(w, h, id);

	node_t *nD = node;
	node_t *nE = base[where];

	color_t cCE = node_color(nE, POS_OUT);
	edge_t *eCE = nE->edges(cCE)->front(POS_OUT);
	
	for (int i = 0; i < 2; ++i)
	{
		pos_t pos = pos_t(i);

		// following variables named according to pos == POS_IN

		node_t *nE = eCE->node(other(pos));

		// create D->E
		edge_t *eDE = base.alloc_edge();
		nD->edges(cCE)->push_front(eDE, pos);
		nE->edges(cCE)->push_front(eDE, other(pos));

		for (edge_t *eCF = eCE->next(pos);
			eCF != eCE->node(pos)->edges(cCE)->end(pos);
			eCF = eCE->next(pos)) // C->F, C->G, etc
		{
			if (rand()%2 == 0)
				break;

			// C->F becomes D->F
			edge_list_t::erase(eCF, pos);
			nD->edges(cCE)->push_back(eCF, pos);
		}

		// create D->Y if needed
		node_t *nG = nD->farest_node(cCE, pos);
		assert(nG != nD);
		if (!nG->edges(other(cCE))->empty(pos))
		{
			edge_t *eGY = nG->edges(other(cCE))->front(pos);

			edge_t *eDY = base.alloc_edge();
			nD->edges(other(cCE))->push_front(eDY, pos);
			edge_list_t::insert(eGY->next(other(pos)), eDY, other(pos));
		}
	}
	
	// remove C->E
	edge_list_t::erase(eCE, POS_IN);
	edge_list_t::erase(eCE, POS_OUT);
	base.free_edge(eCE);
	
	// insert to base
	base.insert(where, node);
}

void racg_erase(acg_base_t &b, size_t where)
{
	using detail::node_color;
	using detail::racg_erase_case1;
	using detail::racg_erase_case2;

	if (where == 0)
	{
		racg_pop(b, POS_IN);
		return;
	}
	if (where == b.num_nodes()-1)
	{
		racg_pop(b, POS_OUT);
		return;
	}

	acg_dev_t base(b);

	node_t *node = base[where];

	color_t in = node_color(node, POS_IN);
	color_t out = node_color(node, POS_OUT);

	if (in == out)
		racg_erase_case1(base, node, in);
	else
		racg_erase_case2(base, node, out);

	base.erase(where);

	base.free_node(node);
}

/*
 *          S
 *          :                      SSSSSSSSSS
 * X--A--B--CCCCCCC--Y             :  :  :  :
 * X  :  :  :  :  :  Y => X--A--B--C--D--E--F--Y
 * X--DDDDDDD--E--F--Y       :  :  :  :
 *          :                TTTTTTTTTT
 *          T
 */
void racg_switch_adjacent(acg_base_t &b, size_t where, pos_t pos)
{
	using detail::node_color;
	using detail::racg_relink_remove;

	acg_dev_t base(b);

	node_t *node = base[where];
	color_t cCD = node_color(node, pos);

	if (node->edges(cCD)->empty(pos))
		return;

	edge_t *eCD = node->edges(cCD)->front(pos);
	node_t *nC = eCD->node(POS_IN);
	node_t *nD = eCD->node(POS_OUT);

	edge_t *eSE = racg_relink_remove(base, nC->edges(cCD), eCD, POS_IN);
	edge_t *eBT = racg_relink_remove(base, nD->edges(cCD), eCD, POS_OUT);

/*
 *          SSSSSSS
 *          :  :  :
 * X--A--B--C--------Y
 * X  :  :  :  :  :  Y
 * X--------D--E--F--Y
 *    :  :  :
 *    TTTTTTT
 */

	edge_list_t *otherC = nC->edges(other(cCD));
	edge_list_t *otherD = nD->edges(other(cCD));

	// remove alter edges
	if (!otherC->empty(POS_IN))
	{
		edge_t *eCY = otherC->pop_front(POS_IN);
		edge_list_t::erase(eCY, POS_OUT);
		base.free_edge(eCY);
	}
	if (!otherD->empty(POS_OUT))
	{
		edge_t *eXD = otherD->pop_front(POS_OUT);
		edge_list_t::erase(eXD, POS_IN);
		base.free_edge(eXD);
	}

	// remove adjacent edge
	edge_list_t::erase(eCD, POS_IN);
	edge_list_t::erase(eCD, POS_OUT);

/*
 *          SSSSSSSSSS
 *          :     :  :
 * X--A--B--C  D--E--F--Y
 *    :  :     :
 *    TTTTTTTTTT
 */

	assert(nC->edges(cCD)->empty(POS_IN));
	assert(otherC->empty(POS_IN));
	assert(nD->edges(cCD)->empty(POS_OUT));
	assert(otherD->empty(POS_OUT));

	// new adjacent edge, using the other color
	otherC->push_back(eCD, POS_IN); // to an empty list
	otherD->push_back(eCD, POS_OUT); // to an empty list

	// add two more edges
	if (eSE != 0)
	{
		edge_t *eSD = base.alloc_edge();
		edge_list_t::insert(eSE, eSD, POS_IN);
		nD->edges(cCD)->push_back(eSD, POS_OUT); // to an empty list
	}
	if (eBT != 0)
	{
		edge_t *eCT = base.alloc_edge();
		edge_list_t::insert(eBT, eCT, POS_OUT);
		nC->edges(cCD)->push_back(eCT, POS_IN); // to an empty list
	}
}

/*
 * S             SSSSSSS
 * :             :     :
 * AAAAAAA--E => AAAA--D
 * :  :  :  E    :  :  D
 * B--C--D--E    B--C--D--E
 */
void racg_switch_last(acg_base_t &b, size_t where, pos_t pos)
{
	using detail::node_color;

	acg_dev_t base(b);

	node_t *nA = base[where];
	color_t cAB = node_color(nA, pos);
	edge_list_t *edgeA = nA->edges(cAB);

	if (edgeA->empty(pos)) // no such edge, it is the last one
		return;
	
	edge_t *eAD = edgeA->back(pos);
	if (edgeA->front(pos) == eAD) // should call swap_adjacent
		return;

	node_t *nD = eAD->node(other(pos));

	edge_list_t::erase(eAD, POS_IN);
	edge_list_t::erase(eAD, POS_OUT);

	assert(nD->edges(cAB)->empty(other(pos)));

	edge_list_t *otherA = nA->edges(other(cAB));
	if (!otherA->empty(pos))
	{
		edge_t *eAE = otherA->pop_front(pos);
		edge_list_t::erase(eAE, other(pos));
		base.free_edge(eAE);
	}
	assert(otherA->empty(pos));

	otherA->push_back(eAD, pos);
	nD->edges(other(cAB))->push_back(eAD, other(pos));

/*
 * S
 * :
 * AAAA--D
 * :  :  D
 * B--C--D--E
 */

	if (!edgeA->empty(other(pos)))
	{
		edge_t *eSA = edgeA->front(other(pos));

		edge_t *eSD = base.alloc_edge();
		edge_list_t::insert(eSA->next(pos), eSD, pos);
		nD->edges(cAB)->push_back(eSD, other(pos)); // to an empty list
	}
}

/*
 * C--D--E--F--T => C--D--E--F--T
 * :  :  :  F       :  :  :  :  T
 * BBBBBBB--F       BBBBBBBBBB--T
 * :        F       :           T
 * A--------F       A-----------T
 * :        :       :
 * XXXXXXXXXX       X
 */
void racg_switch_alter(acg_base_t &b, size_t where, pos_t pos)
{
	using detail::node_color;
	using detail::racg_relink_remove;

	acg_dev_t base(b);

	node_t *nB = base[where];
	color_t cBC = node_color(nB, pos);
	edge_list_t *otherB = nB->edges(other(cBC));

	if (otherB->empty(pos)) // no such edge
		return;

	edge_t *eBF = otherB->front(pos);
	node_t *nF = eBF->node(other(pos));
	
	edge_t *eAT = racg_relink_remove(base, nF->edges(other(cBC)), eBF, other(pos));

/*
 * C--D--E--F--T
 * :  :  :  F  T
 * BBBBBBB--F  T
 * :        :  T
 * A-----------T
 * :        :
 * XXXXXXXXXX
 */

	edge_list_t *edgeF = nF->edges(cBC);

	// remove alter edges
	if (!edgeF->empty(other(pos)))
	{
		edge_t *eXF = edgeF->pop_front(other(pos));
		edge_list_t::erase(eXF, pos);
		base.free_edge(eXF);
	}

	edge_list_t::erase(eBF, POS_IN);
	edge_list_t::erase(eBF, POS_OUT);

/*
 * C--D--E--F--T
 * :  :  :     T
 * BBBBBBB     T
 * :           T
 * A-----------T
 * :
 * X
 */

	assert(edgeF->empty(other(pos)));
	assert(otherB->empty(pos));

	// relink eBF using the other color
	nB->edges(cBC)->push_back(eBF, pos); // must be the last one
	edgeF->push_back(eBF, other(pos)); // to an empty list

	// add one more edge
	if (eAT != 0)
	{
		edge_t *eBT = base.alloc_edge();
		otherB->push_back(eBT, pos); // to an empty list
		edge_list_t::insert(eAT, eBT, other(pos));
	}
}

void racg_to_acg(acg_base_t &b, racg_acg_context_t &ctx)
{
	using detail::node_color;

	assert(ctx.empty());

	acg_dev_t base(b);

	if (base.num_nodes() < 2)
		return;

	ctx.reserve(2*base.num_nodes());

	for (size_t i = 1; i < base.num_nodes(); ++i)
	{
		node_t *node = base[i];
		color_t c = node_color(node, POS_OUT);
		edge_list_t *edges = node->edges(c);
		edge_list_t *others = node->edges(other(c));

		if (others->empty(POS_OUT))
			continue;

		assert(others->front(POS_OUT) == others->back(POS_OUT)); // iacg

		edge_t *bridge = others->front(POS_OUT)->prev(POS_IN);

		for (;;)
		{
			assert(bridge->node(POS_OUT) == edges->back(POS_OUT)->node(POS_IN));

			bridge = bridge->next(POS_OUT);
			if (bridge == bridge->node(POS_OUT)->edges(other(c))->end(POS_OUT))
				break;

			// add acg edge
			edge_t *acg_edge = base.alloc_edge();
			bridge->node(POS_IN)->edges(c)->push_back(acg_edge, POS_IN);
			edges->push_back(acg_edge, POS_OUT);

			ctx.push_back(acg_edge);

			// go next
			c = other(c);
			std::swap(edges, others);
			bridge = acg_edge->prev(POS_IN);
		}
	}
}

void acg_to_racg(acg_base_t &b, racg_acg_context_t &ctx)
{
	acg_dev_t base(b);

	for (size_t i = 0; i < ctx.size(); ++i)
	{
		edge_list_t::erase(ctx[i], POS_IN);
		edge_list_t::erase(ctx[i], POS_OUT);
		base.free_edge(ctx[i]);
	}
	ctx.clear();
}

namespace detail
{

void racg_create_edges(acg_dev_t &base, node_t *node, node_t *next, color_t c, pos_t pos)
{
	node_t *prev = 0;
	for (color_t cur = c; prev != next;
		prev = next, next = next->nearest_node(other(c), pos))
	{
		assert(next->edges(c)->empty(other(pos))); // Reduced ACG

		edge_t *edge = base.alloc_edge();
		node->edges(cur)->push_back(edge, pos);
		next->edges(cur)->push_back(edge, other(pos));
		
		if (cur != c)
			return; // alter edge
		else
			cur = color_t(rand()%2);
	}
}

/*
 * S                SSSSSSSSSSSSS
 * :                :  :  :  :  :
 * AAAAAAAAAAAAA => A  :  :  :  :
 * :  :  :  :  :    :  :  :  :  :
 * B--C--D--E--F    B--C--D--E--F
 */
edge_t *racg_relink_remove(acg_dev_t &base, edge_list_t *edges, edge_t *AB, pos_t pos)
{
	edge_t *AC = AB->next(pos);

	if (edges->empty(other(pos)))
		// remove
	{
		base.free_edges(AC, edges->end(pos), pos);

		return 0;
	}
	else
		// relink
	{
		edge_t *SA = edges->front(other(pos));
		edge_list_t::splice(SA->next(pos), AC, edges->end(pos), pos);

		return SA->next(pos);
	}
}

color_t node_color(const node_t *node, pos_t pos)
{
	const node_t *red = node->nearest_node(COLOR_RED, pos);

	if (red == node)
		return COLOR_BLACK;

	if (red->nearest_node(COLOR_RED, other(pos)) == node)
		return COLOR_RED;
	else
		return COLOR_BLACK;
}

/*
 * case 1
 *
 * X--A--B--C
 * X  :  :  :
 * X--DDDDDDD--Y => X--A--B--CCCCCCC
 *    :  :  :  Y       :  :  :  :  :
 *    E--F--G--Y       EEEEEEE--F--G--Y
 */
void racg_erase_case1(acg_dev_t &base, node_t *node, color_t c)
{
	node_t *nD = node;

	edge_list_t *edgeD = nD->edges(c);
	edge_list_t *otherD = nD->edges(other(c));

	// remove X->D and D->Y
	if (!otherD->empty(POS_IN))
	{
		edge_t *eDY = otherD->pop_front(POS_IN);
		edge_list_t::erase(eDY, POS_OUT);
		base.free_edge(eDY);
	}
	if (!otherD->empty(POS_OUT))
	{
		edge_t *eXD = otherD->pop_front(POS_OUT);
		edge_list_t::erase(eXD, POS_IN);
		base.free_edge(eXD);
	}

/*
 * X--A--B--C
 *    :  :  :
 *    DDDDDDD
 *    :  :  :
 *    E--F--G--Y
 */

	assert(otherD->empty(POS_IN) && otherD->empty(POS_OUT));

	edge_t *eCD = edgeD->front(POS_OUT);
	edge_t *eDE = edgeD->front(POS_IN);

	// D->F, D->G, etc become C->F, C->G, etc
	edge_list_t::splice(eCD->next(POS_IN),
		eDE->next(POS_IN), edgeD->end(POS_IN), POS_IN);
	// B->D, A->D, etc become B->E, A->E, etc
	edge_list_t::splice(eDE->next(POS_OUT),
		eCD->next(POS_OUT), edgeD->end(POS_OUT), POS_OUT);

/*
 * X--A--B--CCCCCCC
 *    :  :  :  :  :
 *    :  :  D  :  :
 *    :  :  :  :  :
 *    EEEEEEE--F--G--Y
 */

	// reuse C->D to be C->E
	edge_list_t::erase(eCD, POS_OUT);
	edge_list_t::insert(eDE, eCD, POS_OUT);

	// free D->E
	edge_list_t::erase(eDE, POS_IN);
	edge_list_t::erase(eDE, POS_OUT);
	base.free_edge(eDE);

	assert(edgeD->empty(POS_IN));
	assert(edgeD->empty(POS_OUT));
}

/*
 * case 2
 *
 * X--A--B--CCCC    X--A--B--C
 * X  :  :  :  :    X  :  :  :
 * X--DDDDDDD--E => X--EEEEEEE
 *    DDDDDDD  :    X        :
 *    DDDDDDD--F    X--------F
 *    DDDDDDD  :    X        :
 *    DDDDDDD--G    X--------G
 *    :        :             :
 *    YYYYYYYYYY             Y
 */
void racg_erase_case2(acg_dev_t &base, node_t *node, color_t out)
{
	color_t cCD = out;
	pos_t pos = POS_OUT;
	node_t *nD = node;
	node_t *nC = nD->nearest_node(cCD, pos);

	if (nC->nearest_node(cCD, other(pos)) == nC->farest_node(cCD, other(pos)))
	{
		cCD = other(cCD);
		pos = POS_IN;
	}

	edge_list_t *edges = nD->edges(cCD);
	edge_list_t *other_edges = nD->edges(other(cCD));

	// free D->Y
	if (!edges->empty(other(pos)))
	{
		edge_t *eDY = edges->pop_front(other(pos));
		edge_list_t::erase(eDY, pos);
		base.free_edge(eDY);
	}

	edge_t *eCD = edges->front(pos);
	edge_t *eCE = eCD->next(other(pos));

	// B->D, A->D, etc become B->E, A->E, etc
	edge_list_t::splice(eCE->next(pos), eCD->next(pos), edges->end(pos), pos);

	// free C->D
	edge_list_t::erase(eCD, POS_IN);
	edge_list_t::erase(eCD, POS_OUT);
	base.free_edge(eCD);

	assert(edges->empty(POS_IN));
	assert(edges->empty(POS_OUT));

	edge_list_t::range_t r = other_edges->all(other(pos));
	if (other_edges->empty(pos))
		// no X, free D->E, D->F, D->G, etc
		base.free_edges(r.first, r.second, other(pos));
	else
	{
		edge_t *eXD = other_edges->front(pos);

		// has X, D->E, D->F, D->G, etc become X->E, X->F, X->G, etc
		edge_list_t::splice(eXD->next(other(pos)),
			r.first, r.second, other(pos));

		// free X->D
		edge_list_t::erase(eXD, POS_IN);
		edge_list_t::erase(eXD, POS_OUT);
		base.free_edge(eXD);
	}
	assert(other_edges->empty(POS_IN));
	assert(other_edges->empty(POS_OUT));
}

} // namespace detail

} // namespace acg



namespace acg
{

namespace detail
{

void nodes_bb_reset(acg_base_t &base, pos_t pos);

misc::dbl_pair packing_longest_path_impl(acg_base_t &base, pos_t pos, size_t &count);

size_t longest_path_propagate(node_t *node, pos_t pos, color_t c, double &bound);

void nodes_bb_reset_unsort(acg_base_t &base, color_t c);

double calc_boundary_unsort(node_t *node, color_t c);

} // namespace detail

misc::dbl_pair packing_longest_path(acg_base_t &base, size_t &count)
{
	return detail::packing_longest_path_impl(base, POS_IN, count);
}

misc::dbl_pair packing_longest_path(acg_base_t &base)
{
	size_t count = 0;
	return detail::packing_longest_path_impl(base, POS_IN, count);
}

void longest_path_bb(acg_base_t &base, misc::dbl_pair shape)
{
	size_t count = 0;

	detail::packing_longest_path_impl(base, POS_OUT, count);

	const size_t n = base.num_nodes();
	for (size_t i = 0; i < n; ++i)
	{
		base[i]->bb(COLOR_RED, POS_OUT)
			= shape.first-base[i]->bb(COLOR_RED, POS_OUT);
		base[i]->bb(COLOR_BLACK, POS_OUT)
			= shape.second-base[i]->bb(COLOR_BLACK, POS_OUT);
	}
}

double packing_lp_unsort(acg_base_t &base, color_t c)
{
	detail::nodes_bb_reset_unsort(base, c);

	const size_t n = base.num_nodes();
	double width = 0;

	for (size_t i = 0; i < n; ++i)
	{
		double t = detail::calc_boundary_unsort(base[i], c);
		if (width < t)
			width = t;
	}
	return width;
}

namespace detail
{

void nodes_bb_reset(acg_base_t &base, pos_t pos)
{
	const size_t n = base.num_nodes();

	for (size_t i = 0; i < n; ++i)
		base[i]->bb(COLOR_RED, pos) = base[i]->bb(COLOR_BLACK, pos) = 0;
}

size_t longest_path_propagate(node_t *node, pos_t pos, color_t c, double &bound)
{
	size_t count = 0;

	double t = node->bb(c, pos)+node->width(c);

	if (node->edges(c)->empty(pos))
	{
		if (bound < t)
			bound = t;
	}
	else
		for (edge_list_t::range_t r = node->edges(c)->all(pos);
			r.first != r.second; r.first = r.first->next(pos))
		{
			node_t *to = r.first->node(other(pos));
            
			if (to->bb(c, pos) < t)
				to->bb(c, pos) = t;

			++count;
		}

	return count;
}

misc::dbl_pair packing_longest_path_impl(acg_base_t &base, pos_t pos, size_t &count)
{
	const size_t n = base.num_nodes();

	nodes_bb_reset(base, pos);

	double width = 0, height = 0;

	if (pos == POS_IN)
		for (size_t i = 0; i < n; ++i)
		{
			count += longest_path_propagate(base[i], pos, COLOR_RED, width);
			count += longest_path_propagate(base[i], pos, COLOR_BLACK, height);
		}
	else
		for (size_t i = n; i > 0; --i)
		{
			count += longest_path_propagate(base[i-1], pos, COLOR_RED, width);
			count += longest_path_propagate(base[i-1], pos, COLOR_BLACK, height);
		}

	return std::make_pair(width, height);
}

void nodes_bb_reset_unsort(acg_base_t &base, color_t c)
{
	const size_t n = base.num_nodes();

	for (size_t i = 0; i < n; ++i)
		base[i]->bb(c, POS_IN) = base[i]->bb(c, POS_OUT) = -1;
}

double calc_boundary_unsort(node_t *node, color_t c)
{
	if (node->bb(c, POS_IN) < -0.5)
		// -1 is always less than and 0 is always greater than -0.5
	{
		node->bb(c, POS_IN) = 0;

		for (edge_list_t::range_t r = node->edges(c)->all(POS_OUT);
			r.first != r.second; r.first = r.first->next(POS_OUT))
		{
			double t = calc_boundary_unsort(r.first->node(POS_IN), c);
			
			if (node->bb(c, POS_IN) < t)
				node->bb(c, POS_IN) = t;
		}
	}
	return node->boundary(c, POS_IN);
}

} // namespace detail

} // namespace acg


namespace fplan
{

module_nets_t::module_nets_t()
{
	nets_.push_back(0);
}

void module_nets_t::append_net()
{
	nets_.push_back(net_modules_.size());
}

void module_nets_t::append_module(size_t im)
{
	net_modules_.push_back(im);
	++nets_.back();
}

double module_nets_t::HPWL(const misc::dbl_pair centers[]) const
{
	double wl = 0;

	for (size_t i = 0; i < nets_.size()-1; ++i)
	{
		size_t from = nets_[i];
		size_t to = nets_[i+1];

		if (from+2 == to) // two pin nets
		{
			const misc::dbl_pair pt0 = centers[net_modules_[from]];
			const misc::dbl_pair pt1 = centers[net_modules_[from+1]];
			wl += fabs(pt0.first-pt1.first)+fabs(pt0.second-pt1.second);
			continue;
		}
		else if (from+3 == to) // three pin nets
		{
			const misc::dbl_pair pt0 = centers[net_modules_[from]];
			const misc::dbl_pair pt1 = centers[net_modules_[from+1]];
			const misc::dbl_pair pt2 = centers[net_modules_[from+2]];
			const misc::dbl_pair ll = misc::min(pt0, misc::min(pt1, pt2));
			const misc::dbl_pair ur = misc::max(pt0, misc::max(pt1, pt2));
			wl += ur.first-ll.first+ur.second-ll.second;
			continue;
		}
		else if (from+4 == to) // four pin nets
		{
			const misc::dbl_pair pt0 = centers[net_modules_[from]];
			const misc::dbl_pair pt1 = centers[net_modules_[from+1]];
			const misc::dbl_pair pt2 = centers[net_modules_[from+2]];
			const misc::dbl_pair pt3 = centers[net_modules_[from+3]];
			const misc::dbl_pair ll = misc::min(misc::min(pt0, pt1), misc::min(pt2, pt3));
			const misc::dbl_pair ur = misc::max(misc::max(pt0, pt1), misc::max(pt2, pt3));
			wl += ur.first-ll.first+ur.second-ll.second;
			continue;
		}
		else if ((from == to) || (from+1 == to))
			continue;

		misc::dbl_pair ll = centers[net_modules_[from]];
		misc::dbl_pair ur = ll;

		for (++from; from != to; ++from)
		{
			const misc::dbl_pair pt = centers[net_modules_[from]];
			
			ll = misc::min(pt, ll);
			ur = misc::max(pt, ur);
		}

		wl += ur.first-ll.first+ur.second-ll.second;
	}

	return wl;
}

pin_nets_t::pin_nets_t()
{
	nets_.push_back(0);
}

void pin_nets_t::append_net()
{
	nets_.push_back(net_pins_.size());
}

void pin_nets_t::append_pin(size_t ip, const misc::dbl_pair &loc)
{
	net_pins_.push_back(std::make_pair(ip, loc));
	++nets_.back();
}

double pin_nets_t::HPWL(size_t nm,
	const misc::dbl_pair shapes[],
	const size_t ignions[],
	const misc::dbl_pair centers[]) const
{
	double wl = 0;

	for (size_t i = 0; i < nets_.size()-1; ++i)
	{
		size_t from = nets_[i];
		size_t to = nets_[i+1];

		if ((from == to) || (from+1 == to))
			continue;

		size_t old_from = from;

		misc::dbl_pair ll, ur;

		for (; from != to; ++from)
		{
			const pin_t &pin = net_pins_[from];
			misc::dbl_pair pt = centers[pin.first];
			const misc::dbl_pair &s = shapes[pin.first];
			const misc::dbl_pair &p = pin.second;
			
			if (pin.first < nm)
				switch (ignions[pin.first])
				{
				case 0:
					pt.first += s.first*p.first;
					pt.second += s.second*p.second;
					break;
				case 1:
					pt.first += s.second*p.second;
					pt.second -= s.first*p.first;
					break;
				case 2:
					pt.first -= s.first*p.first;
					pt.second -= s.second*p.second;
					break;
				case 3:
					pt.first -= s.second*p.second;
					pt.second += s.first*p.first;
					break;
				case 4:
					pt.first += s.first*p.first;
					pt.second -= s.second*p.second;
					break;
				case 5:
					pt.first += s.second*p.second;
					pt.second += s.first*p.first;
					break;
				case 6:
					pt.first -= s.first*p.first;
					pt.second += s.second*p.second;
					break;
				case 7:
					pt.first -= s.second*p.second;
					pt.second -= s.first*p.first;
					break;
				default:
					assert(false);
				}

			if (from == old_from)
				ll = ur = pt;
			else
			{
				ll = misc::min(pt, ll);
				ur = misc::max(pt, ur);
			}
		}

		wl += ur.first-ll.first+ur.second-ll.second;
	}

	return wl;
}

} // namespace fplan


