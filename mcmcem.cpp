// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <numeric>
#include <string>
#include <math.h>
#include <random>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>


using namespace Rcpp;

boost::random::mt19937 gen;

enum Trans_type {SPLIT, MERGE, SWITCH, SWAP, 
                 UPDATE_IN, UPDATE_OUT, SWITCH1, SWITCH2, 
                 EXTRACT, ABSORB,
                 INVALID};

struct Data {
  int n_obs;
  int n_groups;
  IntegerVector group;
  IntegerVector group_sizes;
  IntegerVector group_start_inds;
  NumericVector time;
  IntegerVector loc;
  NumericVector loc_prior;
  IntegerVector scoreSample_first;
  IntegerVector scoreSample_second;
  
  boost::random::discrete_distribution<> dist_prop;
  boost::random::discrete_distribution<> dist_flat;
  
  Data (IntegerVector group_,
        IntegerVector group_sizes_,
        IntegerVector group_start_inds_,
        NumericVector time_,
        IntegerVector loc_,
        NumericVector loc_prior_,
        IntegerVector scoreSample_first_,
        IntegerVector scoreSample_second_):
    n_obs(group_.size()), n_groups(group_sizes_.size()),
    group(group_), group_sizes(group_sizes_),
    group_start_inds(group_start_inds_),
    time(time_), loc(loc_), loc_prior(loc_prior_),
    scoreSample_first(scoreSample_first_), scoreSample_second(scoreSample_second_)
  {
    std::vector<double> group_weights(n_groups);
    std::vector<double> group_pos(n_groups);
    for (int g = 0; g < n_groups; g++) {
      if (group_sizes[g] == 1) {
        group_weights[g] = 0;
        group_pos[g] = 0;
      } else {
        group_weights[g] = (double)group_sizes[g] * ((double)group_sizes[g] - 1);
        group_pos[g] = (double)group_sizes[g] * ((double)group_sizes[g] - 1);
      }
    }
    dist_prop = boost::random::discrete_distribution<>(group_weights.begin(), group_weights.end());
    dist_flat = boost::random::discrete_distribution<>(group_pos.begin(), group_pos.end());
  }
  
};

class MH_chain; // pre declare so Model class knows it exists

class Model {
public:
  double rho;
  double nb_shape1;
  double nb_shape2;
  double A;
  double discount;
  double zeta;
  int S;
  
  Model (const Data& data, List params) : nb_shape1(0.25), nb_shape2(.01), 
  A(-.64 + .46), discount(0.64), zeta(0.32) {
    S = data.n_obs;
    rho = S / 20.0;
    if (params["nb_shape1"] != R_NilValue)
      nb_shape1 = params["nb_shape1"];
    if (params["nb_shape2"] != R_NilValue)
      nb_shape2 = params["nb_shape2"];
    if (params["crp_A"] != R_NilValue)
      A = params["crp_A"];
    if (params["crp_discount"] != R_NilValue)
      discount = params["crp_discount"];
    if (params["dircat_discount"] != R_NilValue)
      zeta = params["dircat_discount"];
    if (params["rho"] != R_NilValue)
      rho = params["rho"];
  }
  void max_exp_loglik(const MH_chain& chain, const Data& data);
  
  double get_rho() {
    return(rho);
  }
};

class Correspondence {
private:
  std::vector<int> L1;
  std::vector<int> L2;
public:
  Correspondence(const Data& data, IntegerVector init_L1) {
    L1.resize(data.n_obs);
    L2.resize(data.n_obs);
    
    if (init_L1.size() == data.n_obs) {
      Rcout << "initializing correspondence with init_L1" << std::endl;
      for (int id = 1; id <= data.n_obs; id++) {
        if (init_L1[id-1] != 0) {
          L1[id-1] = init_L1[id-1];
          L2[init_L1[id-1]-1] = id;
        }
      }
    } else {
      // default initialization
      for (int g = 0; g < data.group_sizes.size(); g++) {
        std::map<int, std::set<int>> ids_per_loc;
        for (int id =
             data.group_start_inds[g];
             id < data.group_start_inds[g] + data.group_sizes[g] - 1;
             id++) {
          ids_per_loc[data.loc[id-1]].insert(id);
        }
        
        for (auto it = ids_per_loc.begin(); it != ids_per_loc.end(); ++it) {
          std::set<int>::iterator it2 = (it->second).begin();
          int id_prev, id_cur;
          id_prev = *it2;
          id_cur = *it2;
          ++it2;
          while (it2 != (it->second).end()) {
            id_cur = *it2;
            L1[id_prev-1] = id_cur;
            L2[id_cur-1] = id_prev;
            ++it2;
            id_prev = id_cur;
          }
          L1[id_cur-1] = 0;
        }
      }
    }
    
  }
  
  std::vector<int> get_L1() {
    return(L1);
  }
  
  bool connected(int id1, int id2) {
    if (id1 > id2)
      std::swap(id1, id2);
    
    int id = id1;
    while (id < id2 && id != 0) {
      id = L1[id-1];
      if (id == id2)
        return(true);
    }
    return(false);
  }
  
  int get_prev(int id) { return(L2[id-1]); }
  int get_next(int id) { return(L1[id-1]); }
  
  bool check_homogeneous(int id1, const Data& data) {
    int loc = data.loc[id1-1];
    int id = L2[id1-1];
    while (id != 0) {
      if (data.loc[id-1] != loc)
        return(false);
      id = L2[id-1];
    }
    id = L1[id1-1];
    while (id != 0) {
      if (data.loc[id-1] != loc)
        return(false);
      id = L1[id-1];
    }
    return(true);
  }
  
  
  int count_loc(int loc, int id1, const Data& data) {
    int c = 0;
    int id = id1;
    while (id != 0) {
      if (data.loc[id-1] == loc)
        c++;
      id = L1[id-1];
    }
    id = L2[id1-1];
    while (id != 0) {
      if (data.loc[id-1] == loc)
        c++;
      id = L2[id-1];
    }
    return(c);
  }
  
  int n_locs(int id1, const Data& data) {
    std::unordered_set<int> locs;
    int id = id1;
    while (id != 0) {
      locs.insert(data.loc[id-1]);
      id = L1[id-1];
    }
    id = L2[id1-1];
    while (id != 0) {
      locs.insert(data.loc[id-1]);
      id = L2[id-1];
    }
    return(locs.size());
  }
  
  int upper_bound(int id1, double t_split, const Data& data) {
    int id = L1[id1-1];
    int old_id = id1;
    while (id != 0) {
      if (data.time[id-1] > t_split) {
        return(old_id);
      }
      old_id = id;
      id = L1[id-1];
    }
    return(old_id);
  }
  
  int lower_bound(int id2, double t_split, const Data& data) {
    int id = L2[id2-1];
    int old_id = id2;
    while (id != 0) {
      if (data.time[id-1] < t_split) {
        return(old_id);
      }
      old_id = id;
      id = L2[id-1];
    }
    return(old_id);
  }
  
  int n_before(int id2) {
    int n = 0;
    int id = id2;
    while (id != 0) {
      n++;
      id = L2[id-1];
    }
    return(n);
  }
  
  int n_after(int id1) {
    int n = 0;
    int id = id1;
    while (id != 0) {
      n++;
      id = L1[id-1];
    }
    return(n);
  }
  
  void insert(int id1, int id2, Trans_type trans_type, const Data& data) {
    switch (trans_type) {
    case SPLIT:
      L1[id1-1] = 0;
      L2[id2-1] = 0;
      break;
    case MERGE:
      L1[id1-1] = id2;
      L2[id2-1] = id1;
      break;
    case SWITCH:
      if (L1[id1-1] != 0)
        L2[L1[id1-1]-1] = L2[id2-1];
      if (L2[id2-1] != 0)
        L1[L2[id2-1]-1] = L1[id1-1];
      L1[id1-1] = id2;
      L2[id2-1] = id1;
      break;
    case SWAP:
      if (L2[id1-1] != 0)
        L1[L2[id1-1]-1] = L1[id1-1];
      if (L1[id1-1] != 0)
        L2[L1[id1-1]-1] = L2[id1-1];
      
      if (id2 < id1) {
        L2[id1-1] = id2;
        L1[id1-1] = L1[id2-1];
        if (L1[id2-1] != 0) {
          L2[L1[id2-1]-1] = id1;
        }
        L1[id2-1] = id1;
      } else if (id2 > id1) {
        L1[id1-1] = id2;
        L2[id1-1] = L2[id2-1];
        if (L2[id2-1] != 0) {
          L1[L2[id2-1]-1] = id1;
        }
        L2[id2-1] = id1;
      }
      break;
    case UPDATE_IN:
      if (id2 < id1) {
        if (L1[id2-1] != 0)
          L2[L1[id2-1]-1] = id1;
        L2[id1-1] = id2;
        L1[id1-1] = L1[id2-1];
        L1[id2-1] = id1;
      } else if (id2 > id1) {
        if (L2[id2-1] != 0)
          L1[L2[id2-1]-1] = id1;
        L1[id1-1] = id2;
        L2[id1-1] = L2[id2-1];
        L2[id2-1] = id1;
      }
      break;
    case UPDATE_OUT:
      if (L1[id1-1] != 0)
        L2[L1[id1-1]-1] = L2[id1-1];
      if (L2[id1-1] != 0)
        L1[L2[id1-1]-1] = L1[id1-1];
      L1[id1-1] = 0;
      L2[id1-1] = 0;
      break;
    case SWITCH1:
      L2[L1[id1-1]-1] = 0;
      L1[id1-1] = id2;
      L2[id2-1] = id1;
      break;
    case SWITCH2:
      L1[L2[id2-1]-1] = 0;
      L1[id1-1] = id2;
      L2[id2-1] = id1;
      break;
    case EXTRACT:
    {
      std::set<int> id_list1;
      std::set<int> id_list2;
      
      int id = id1;
      while (id != 0) {
        if (data.loc[id-1] == data.loc[id1-1]) {
          id_list2.insert(id);
        } else {
          id_list1.insert(id);
        }
        id = L1[id-1];
      }
      id = L2[id1-1];
      while (id != 0) {
        if (data.loc[id-1] == data.loc[id1-1]) {
          id_list2.insert(id);
        } else {
          id_list1.insert(id);
        }
        id = L2[id-1];
      }
      
      if (id_list1.size() == 0 || id_list2.size() == 0) {
        Rcout << "id_list1.size() = " << id_list1.size() << std::endl;
        Rcout << "id_list2.size() = " << id_list2.size() << std::endl;
        Rf_error("extract error");
      }
      
      auto it = id_list1.begin();
      id = *it;
      L2[id-1] = 0;
      ++it;
      int id_next;
      while (it != id_list1.end()) {
        id_next = *it;
        L2[id_next-1] = id;
        L1[id-1] = id_next;
        id = id_next;
        ++it;
      }
      L1[id-1] = 0;
      
      auto it2 = id_list2.begin();
      id = *it2;
      L2[id-1] = 0;
      ++it2;
      while (it2 != id_list2.end()) {
        id_next = *it2;
        L2[id_next-1] = id;
        L1[id-1] = id_next;
        id = id_next;
        ++it2;
      }
      L1[id-1] = 0;
      
    }
      
      break;
    case ABSORB:
    {
      std::set<int> id_list_absorbed;
      
      int id = id1;
      while (id != 0) {
        id_list_absorbed.insert(id);
        id = L1[id-1];
      }
      id = L2[id1-1];
      while (id != 0) {
        id_list_absorbed.insert(id);
        id = L2[id-1];
      }
      id = id2;
      while (id != 0) {
        id_list_absorbed.insert(id);
        id = L1[id-1];
      }
      id = L2[id2-1];
      while (id != 0) {
        id_list_absorbed.insert(id);
        id = L2[id-1];
      }
      
      auto it = id_list_absorbed.begin();
      id = *it;
      L2[id-1] = 0;
      ++it;
      int id_next;
      while (it != id_list_absorbed.end()) {
        id_next = *it;
        L2[id_next-1] = id;
        L1[id-1] = id_next;
        id = id_next;
        ++it;
      }
      L1[id-1] = 0;
      
    }
      break;
    case INVALID:
      Rf_error("tried to insert invalid trans type");
      break;
    }
  }
};

class MH_chain {
private:
  Correspondence correspondence;
  unsigned int N_avail;
public:
  int n_samples;
  unsigned long long int N_avail_sum;
  std::vector<int> N_group;
  
  std::vector<std::pair<std::pair<int, int>, unsigned long int>> score_sample;
  
  MH_chain(const Data& data, const Model& model, IntegerVector init_L1) : 
    correspondence(data, init_L1), N_avail(data.n_obs) {
    N_group.resize(data.n_groups);
    N_group.assign(data.n_groups, 0);
    
    N_avail = 0;
    for (int g = 0; g < data.n_groups; g++) {
      for (int id = data.group_start_inds[g];
           id <= data.group_start_inds[g] + data.group_sizes[g] - 1;
           id++) {
        if (correspondence.get_next(id) == 0) {
          N_group[g]++;
          N_avail++;
        }
      }
    }
  }
  
  void reset(int n_samples_, const Data& data, const Model& model) {
    n_samples = n_samples_;
    N_avail_sum = 0;
  }
  
  std::vector<int> get_L1() {
    return(correspondence.get_L1());
  }
  
  std::pair<int, int> sample_pair(const Data& data) {
    int id1, id2, group_id, group_size, offset, tmp;
    
    group_id = data.dist_flat(gen) + 1;
    group_size = data.group_sizes[group_id-1];
    offset = data.group_start_inds[group_id-1];
    
    boost::random::uniform_int_distribution<> dist2(1 + offset, offset + group_size);
    id1 = dist2(gen);
    
    tmp = 1 + rand() % (group_size - 1);
    if (tmp < (id1-offset)) {
      id2 = tmp + offset;
    } else {
      id2 = tmp + offset + 1;
    }
    
    return(std::make_pair(id1, id2));
  }
  
  std::pair<int, int> sample_pair_prop(const Data& data) {
    // sample pair weighted by group size
    int id1, id2, group_id, group_size, offset, tmp;
    
    group_id = data.dist_prop(gen) + 1;
    group_size = data.group_sizes[group_id-1];
    offset = data.group_start_inds[group_id-1];
    
    boost::random::uniform_int_distribution<> dist2(1 + offset, offset + group_size);
    id1 = dist2(gen);
    
    tmp = 1 + rand() % (group_size - 1);
    if (tmp < (id1-offset)) {
      id2 = tmp + offset;
    } else {
      id2 = tmp + offset + 1;
    }
    
    return(std::make_pair(id1, id2));
  }
  
  NumericVector get_score_counts() {
    CharacterVector pair_names(score_sample.size());
    NumericVector probs(score_sample.size());
    std::stringstream ss;
    
    for (size_t i = 0; i < score_sample.size(); i++) {
      ss.str("");
      ss.clear();
      ss << score_sample[i].first.first << "-" << score_sample[i].first.second;
      pair_names[i] = ss.str();
      probs[i] = 1000.0*((double)score_sample[i].second / n_samples);
    }
    probs.names() = pair_names;
    return(probs);
  }
  
  void insert(int id1, int id2, Trans_type trans_type,
              int group_id, const Data& data, const Model& model) {
    switch (trans_type) {
    case SPLIT:
      N_avail++;
      N_group[group_id-1]++;
      break;
    case MERGE:
      N_avail--;
      N_group[group_id-1]--;
      break;
    case SWITCH:
      break;
    case SWAP:
      break;
    case UPDATE_IN:
      N_avail--;
      N_group[group_id-1]--;
      break;
    case UPDATE_OUT:
      N_avail++;
      N_group[group_id-1]++;
      break;
    case SWITCH1:
      break;
    case SWITCH2:
      break;
    case EXTRACT:
      N_avail++;
      N_group[group_id-1]++;
      break;
    case ABSORB:
      N_avail--;
      N_group[group_id-1]--;
      break;
    case INVALID:
      Rf_error("tried to insert invalid trans type");
    }
  correspondence.insert(id1, id2, trans_type, data);
  }
  
  double log_hastings_ratio(const Data& data, int id1, int id2, Trans_type trans_type, bool search_model_);
  double log_prob_ratio(const Model& model, const Data& data, int id1, int id2, Trans_type trans_type, 
                        bool estimate_rho_, bool search_model_);
  
  void sample(const Data& data, const Model& model, 
              bool init = false, bool score = false, 
              bool estimate_rho = false, bool search_model = false) {
    int acceptance_rate = 0;
    int split_count = 0;
    int split_accept = 0;
    int merge_count = 0;
    int merge_accept = 0;
    int switch_count = 0;
    int switch_accept = 0;
    int swap_count = 0;
    int swap_accept = 0;
    int update_in_count = 0;
    int update_in_accept = 0;
    int update_out_count = 0;
    int update_out_accept = 0;
    int extract_count = 0;
    int extract_accept = 0;
    int absorb_count = 0;
    int absorb_accept = 0;
    int iter = 0;
    
    if (score) {
      for (int i = 0; i < data.scoreSample_first.size(); i++) {
        score_sample.push_back(std::make_pair(std::make_pair(
            data.scoreSample_first[i], data.scoreSample_second[i]), 0));
      }
    }
    
    while(iter < n_samples) {
      std::pair<int, int> id_pair;
      int id1, id2, id3, group_id;
      Trans_type trans_type = INVALID;
      double accept_prob_log = -99;
      
      if (init) {
        id_pair = sample_pair_prop(data);
      } else {
        id_pair = sample_pair(data);
      }
      id1 = id_pair.first;
      id2 = id_pair.second;
      group_id = data.group[id1-1];
      
      double c = ((double) rand() / (RAND_MAX));
      if (c <= 0.75) {
        // SWAP or UPDATE
        double tmp;
        if (id2 < id1) std::swap(id1, id2);
        if (((double) rand() / (RAND_MAX)) <= .5) {
          std::swap(id1, id2);
          tmp = id2;
          id2 = correspondence.upper_bound(id2, data.time[id1-1], data);
        } else {
          tmp = id2;
          id2 = correspondence.lower_bound(id2, data.time[id1-1], data);
        }
        if (id1==id2 || 
            correspondence.get_prev(id2) == id1 || 
            correspondence.get_next(id2) == id1) {
          trans_type = UPDATE_OUT;
          if (search_model) {
            int c_loc = correspondence.count_loc(data.loc[id1-1], id1, data);
            int nn = correspondence.n_before(id1) + correspondence.n_after(id1) - 1;
            
            if (c_loc == 1) {
              trans_type = UPDATE_OUT;
            } else {
              if (((double) rand() / (RAND_MAX)) <= ((double)c_loc / nn)) {
                trans_type = UPDATE_OUT;
              } else {
                trans_type = EXTRACT;
              }
            }
          }
        } else if (correspondence.get_prev(id1) == 0 && 
          correspondence.get_next(id1) == 0) {
          trans_type = UPDATE_IN;
        } else {
          trans_type = SWAP;
          if (search_model) {
            int c_loc = correspondence.count_loc(data.loc[id1-1], id1, data);
            int c_loc2 = correspondence.count_loc(data.loc[id1-1], id2, data);
            
            if (correspondence.check_homogeneous(id1, data) &&
                (c_loc2 == 0)) {
              // absorb or swap
              if (((double) rand() / (RAND_MAX)) <= (1.0 / 10)) {
                trans_type = SWAP;
              } else {
                trans_type = ABSORB;
              }
            } else {
              trans_type = SWAP;
              if (c_loc == 1 &&
                  data.loc[id1-1] == data.loc[id2-1] &&
                  correspondence.check_homogeneous(id2, data)) {
                trans_type = SWAP;
              }
            }
          }
        }
        
      } else {
        // SWITCH or SPLIT/MERGE
        if (id2 < id1) {
          std::swap(id1, id2);
          if (correspondence.get_prev(id2) != 0) {
            trans_type = SPLIT;
            id1 = correspondence.get_prev(id2);
          } else { 
            id3 = correspondence.get_next(id1);
            if (id3 != 0) {
              trans_type = SWITCH1;
            } else {
              trans_type = MERGE;
            }
          }
        } else {
          if (correspondence.get_next(id1) != 0) {
            trans_type = SPLIT;
            id2 = correspondence.get_next(id1);
          } else {
            id3 = correspondence.get_prev(id2);
            if (id3 != 0) {
              trans_type = SWITCH2;
            } else {
              trans_type = MERGE;
            }
          }
        }
        
      }
      
      if (trans_type == SPLIT) {
        split_count++;
      } else if (trans_type == MERGE) {
        merge_count++;
      } else if (trans_type == SWAP) {
        swap_count++;
      } else if (trans_type == UPDATE_IN) {
        update_in_count++;
      } else if (trans_type == UPDATE_OUT) {
        update_out_count++;
      } else if (trans_type == SWITCH1) {
        switch_count++;
      } else if (trans_type == SWITCH2) {
        switch_count++;
      } else if (trans_type == EXTRACT) {
        extract_count++;
      } else if (trans_type == ABSORB) {
        absorb_count++;
      } else if (trans_type == INVALID) {
        Rf_error("invalid trans_type");
      }
      
      accept_prob_log = log_hastings_ratio(data, id1, id2, trans_type, search_model);
      if (!init) {
        accept_prob_log += log_prob_ratio(model, data, id1, id2, trans_type, estimate_rho, search_model);
      }
      
      if ((rand() / ((float)RAND_MAX)) <= exp(accept_prob_log)) {
        acceptance_rate++;
        if (trans_type == SPLIT) {
          split_accept++;
        } else if (trans_type == MERGE) {
          merge_accept++;
        } else if (trans_type == SWITCH1 || trans_type == SWITCH2) {
          switch_accept++;
        } else if (trans_type == SWAP) {
          swap_accept++;
        } else if (trans_type == UPDATE_IN) {
          update_in_accept++;
        } else if (trans_type == UPDATE_OUT) {
          update_out_accept++;
        } else if (trans_type == EXTRACT) {
          extract_accept++;
        } else if (trans_type == ABSORB) {
          absorb_accept++;
        } else if (trans_type == INVALID) {
          Rf_error("invalid trans_type");
        }
        insert(id1, id2, trans_type, group_id, data, model);
      }
      
      // update sufficient statistics
      N_avail_sum += N_avail;
      
      if (score && (iter % 1000) == 0) {
        if (iter % 1000000 == 0) {
          Rcout << "Score iter " << iter << std::endl;
        }
        for (size_t i = 0; i < score_sample.size(); i++) {
          if (correspondence.connected(score_sample[i].first.first, 
                                       score_sample[i].first.second)) {
            score_sample[i].second++;
          }
        }
      }
      
      iter++;
    }
    
    Rcout << "acceptance rate = " << acceptance_rate / (double)n_samples << std::endl;
    Rcout << "split_count = " << split_count << "\t accepted = " << split_accept << std::endl;
    Rcout << "merge_count = " << merge_count << "\t accepted = " << merge_accept << std::endl;
    Rcout << "switch_count = " << switch_count << "\t accepted = " << switch_accept << std::endl;
    Rcout << "swap_count = " << swap_count << "\t accepted = " << swap_accept << std::endl;
    Rcout << "update_in_count = " << update_in_count << "\t accepted = " << update_in_accept << std::endl;
    Rcout << "update_out_count = " << update_out_count << "\t accepted = " << update_out_accept << std::endl;
    Rcout << "extract_count = " << extract_count << "\t accepted = " << extract_accept << std::endl;
    Rcout << "absorb_count = " << absorb_count << "\t accepted = " << absorb_accept << std::endl;
    Rcout << "N = " << N_avail << std::endl;
    
  }
  
};

void Model::max_exp_loglik(const MH_chain& chain, const Data& data) {
  rho = chain.N_avail_sum / (S * (double)chain.n_samples);
  Rcout << "rho = " << rho << std::endl;
}

double MH_chain::log_hastings_ratio(const Data& data,
                                    int id1, int id2,
                                    Trans_type trans_type,
                                    bool search_model_) {
  
  double alpha = 0.0;
  switch (trans_type) {
  case SPLIT:
  {
    int n1 = correspondence.n_before(id1);
    int n2 = correspondence.n_after(id2);
    alpha += boost::math::lgamma(n1 + 1) + boost::math::lgamma(n2 + 1) - boost::math::lgamma(n1 + n2 + 1);
    alpha += log(N_avail + 1);
    
    int offset = data.group_start_inds[data.group[id1-1]-1];
    int group_size = data.group_sizes[data.group[id1-1]-1];
    int nb2 = id2 - offset - 1;
    int na1 = group_size + offset - id1;
    
    alpha += log(2);
    alpha -= log(na1 + nb2);
  }
    break;
  case MERGE:
  {
    int n1 = correspondence.n_before(id1);
    int n2 = correspondence.n_after(id2);
    alpha += -boost::math::lgamma(n1 + 1) - boost::math::lgamma(n2 + 1) + boost::math::lgamma(n1 + n2 + 1);
    alpha += -log(N_avail);
    
    int offset = data.group_start_inds[data.group[id1-1]-1];
    int group_size = data.group_sizes[data.group[id1-1]-1];
    int nb2 = id2 - offset - 1;
    int na1 = group_size + offset - id1;
    
    alpha -= log(2);
    alpha += log(na1 + nb2);
  }
    break;
  case SWITCH:
  {
    int n1b = correspondence.n_before(id1);
    int n1a = correspondence.n_after(id1) - 1;
    int n2a = correspondence.n_after(id2);
    int n2b = correspondence.n_before(id2) - 1;
    alpha += boost::math::lgamma(n1b + n2a + 1) + boost::math::lgamma(n2b + n1a + 1) -
      boost::math::lgamma(n1b + n1a + 1) - boost::math::lgamma(n2b + n2a + 1);
  }
    break;
  case SWAP:
  {
    int n1 = correspondence.n_before(id1) + correspondence.n_after(id1) - 1;
    int n2 = correspondence.n_before(id2) + correspondence.n_after(id2) - 1;
    
    alpha += log(n1 - 1) - log(n2); // proposal selection ratio
    alpha += log(n2 + 1) - log(n1); // order within chain
    
    if (search_model_) {
      int c_loc = correspondence.count_loc(data.loc[id1-1], id1, data);
      int c_loc2 = correspondence.count_loc(data.loc[id1-1], id2, data);
      if (correspondence.check_homogeneous(id1, data) &&
          c_loc2 == 0) {
        alpha -= log(1.0/10);
      } else if (data.loc[id2-1] == data.loc[id1-1] &&
        correspondence.check_homogeneous(id2, data) &&
        c_loc == 1) {
        alpha += log(1.0/10);
      }
    }
  }
    break;
  case UPDATE_IN:
  {
    int nn = correspondence.n_before(id2) + correspondence.n_after(id2) - 1;
    alpha += log(nn + 1);
    alpha += -log(N_avail);
    
    if (search_model_) {
      int n_loc = correspondence.count_loc(data.loc[id1-1], id2, data);
      if (n_loc > 0) {
        alpha += (log(n_loc + 1) - log(nn + 1));
      }
    }
    
  }
    break;
  case UPDATE_OUT:
  {
    int nn = correspondence.n_before(id1) + correspondence.n_after(id1) - 1;
    alpha += -log(nn);
    alpha += log(N_avail + 1);
    
    if (search_model_) {
      int n_loc = correspondence.count_loc(data.loc[id1-1], id1, data);
      if (n_loc > 1) {
        alpha -= (log(n_loc) - log(nn));
      }
    }
    
  }
    break;
  case SWITCH1:
  {
    int n1 = correspondence.n_before(id1);
    int n2 = correspondence.n_after(id2);
    int n3 = correspondence.n_after(id1) - 1;
    
    alpha += boost::math::lgamma(n1+n2+1) + boost::math::lgamma(n3+1) - 
      boost::math::lgamma(n1+n3+1) - boost::math::lgamma(n2+1);
    
  }
    break;
  case SWITCH2:
  {
    int n1 = correspondence.n_before(id1);
    int n2 = correspondence.n_after(id2);
    int n3 = correspondence.n_before(id2) - 1;
    
    alpha += boost::math::lgamma(n1+n2+1) + boost::math::lgamma(n3+1) - 
      boost::math::lgamma(n2+n3+1) - boost::math::lgamma(n1+1);
  }
    break;
  case EXTRACT:
  {
    alpha += log(N_avail + 1);
    
    int nn = correspondence.n_before(id1) + correspondence.n_after(id1) - 1;
    int n_loc = correspondence.count_loc(data.loc[id1-1], id1, data);
    
    alpha += boost::math::lgamma(n_loc+1) + boost::math::lgamma(nn - n_loc + 1) -
      boost::math::lgamma(nn+1);
    
    alpha -= (log(nn - 1) + log(nn - n_loc) - log(nn));
    alpha += log(nn - n_loc);
    
    alpha += log(9.0/10);
  }
    break;
  case ABSORB:
  {
    alpha += -log(N_avail);
    int n2 = correspondence.n_before(id2) + correspondence.n_after(id2) - 1;
    int n_loc = correspondence.n_before(id1) + correspondence.n_after(id1) - 1;
    
    alpha += boost::math::lgamma(n2 + n_loc + 1) -
      boost::math::lgamma(n2 + 1) - boost::math::lgamma(n_loc + 1);
    
    alpha -= log(n2);
    alpha += (log(n2 + n_loc - 1) + log(n2) - log(n2 + n_loc));
    
    alpha -= log(9.0/10);
  }
    break;
  case INVALID:
    Rf_error("tried to insert invalid trans type");
  }
  
  return(alpha);
}


double MH_chain::log_prob_ratio(const Model& model, const Data& data, 
                                int id1, int id2,
                                Trans_type trans_type,
                                bool estimate_rho,
                                bool search_model_) {
  double alpha = 0.0;
  
  switch (trans_type) {
  case SPLIT:
    if (estimate_rho) {
      alpha += log(model.rho) - log(1.0 - model.rho);
      alpha += log(model.S - N_avail) - log(N_avail + 1);
    }
    
    alpha += log(N_group[data.group[id1-1]-1] - model.discount) - log(N_avail + model.A);
    {
      int n_before1 = correspondence.n_before(id1);
      int n_after2 = correspondence.n_after(id2);
      double p = 1.0 / (1 + model.nb_shape2);
      
      alpha += model.nb_shape1 * log(1 - p) - log(p) +
        boost::math::lgamma(n_before1-1 + model.nb_shape1) -
        boost::math::lgamma(n_before1) +
        boost::math::lgamma(n_after2-1 + model.nb_shape1) -
        boost::math::lgamma(n_after2) -
        boost::math::lgamma(n_before1 + n_after2 - 1 + model.nb_shape1) +
        boost::math::lgamma(n_before1 + n_after2) -
        boost::math::lgamma(model.nb_shape1);
      
      if (search_model_) {
        std::unordered_map<int, int> nz1;
        std::unordered_map<int, int> nz2;
        std::unordered_map<int, int> nz12;
        int id;
        id = id1;
        while (id != 0) {
          nz1[data.loc[id-1]]++;
          id = correspondence.get_prev(id);
        }
        nz12 = nz1;
        id = id2;
        while (id != 0) {
          nz2[data.loc[id-1]]++;
          nz12[data.loc[id-1]]++;
          id = correspondence.get_next(id);
        }
        
        alpha += (boost::math::lgamma(model.zeta) - boost::math::lgamma(model.zeta + n_before1)) +
          (boost::math::lgamma(model.zeta) - boost::math::lgamma(model.zeta + n_after2)) -
          (boost::math::lgamma(model.zeta) - boost::math::lgamma(model.zeta + n_before1 + n_after2));
        
        for (auto it = nz1.begin(); it != nz1.end(); ++it) {
          alpha += boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
          alpha -= boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
        }
        for (auto it = nz2.begin(); it != nz2.end(); ++it) {
          alpha += boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
          alpha -= boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
        }
        for (auto it = nz12.begin(); it != nz12.end(); ++it) {
          alpha -= boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
          alpha += boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
        }
        
      }
    }
    
    break;
  case MERGE:
    if (estimate_rho) {
      alpha += log(1 - model.rho) - log(model.rho);
      alpha += log(N_avail) - log(model.S - N_avail + 1);
    }
    
    alpha += log(N_avail + model.A - 1) - log(N_group[data.group[id1-1]-1] - 1 - model.discount);
    {
      int n_before1 = correspondence.n_before(id1);
      int n_after2 = correspondence.n_after(id2);
      double p = 1 / (1 + model.nb_shape2);
      
      alpha -= (model.nb_shape1 * log(1 - p) - log(p) +
        boost::math::lgamma(n_before1-1 + model.nb_shape1) -
        boost::math::lgamma(n_before1) +
        boost::math::lgamma(n_after2-1 + model.nb_shape1) -
        boost::math::lgamma(n_after2) -
        boost::math::lgamma(n_before1 + n_after2 - 1 + model.nb_shape1) +
        boost::math::lgamma(n_before1 + n_after2) -
        boost::math::lgamma(model.nb_shape1));
      
      if (search_model_) {
        std::unordered_map<int, int> nz1;
        std::unordered_map<int, int> nz2;
        std::unordered_map<int, int> nz12;
        int id;
        id = id1;
        while (id != 0) {
          nz1[data.loc[id-1]]++;
          id = correspondence.get_prev(id);
        }
        nz12 = nz1;
        id = id2;
        while (id != 0) {
          nz2[data.loc[id-1]]++;
          nz12[data.loc[id-1]]++;
          id = correspondence.get_next(id);
        }
        
        alpha += -(boost::math::lgamma(model.zeta) -
          boost::math::lgamma(model.zeta + n_before1) -
          boost::math::lgamma(model.zeta + n_after2) +
          boost::math::lgamma(model.zeta + n_before1 + n_after2));
        
        for (auto it = nz1.begin(); it != nz1.end(); ++it) {
          alpha -= boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
          alpha += boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
        }
        for (auto it = nz2.begin(); it != nz2.end(); ++it) {
          alpha -= boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
          alpha += boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
        }
        for (auto it = nz12.begin(); it != nz12.end(); ++it) {
          alpha += boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
          alpha -= boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
        }
        
      }
      
    }
    
    break;
  case SWITCH:
  {
    int n_before1 = correspondence.n_before(id1);
    int n_after1 = correspondence.n_after(id1) - 1;
    int n_before2 = correspondence.n_before(id2) - 1;
    int n_after2 = correspondence.n_after(id2);
    
    alpha += boost::math::lgamma(n_before1 + n_after2 - 1 + model.nb_shape1) -
      boost::math::lgamma(n_before1 + n_after2) +
      boost::math::lgamma(n_before2 + n_after1 - 1 + model.nb_shape1) -
      boost::math::lgamma(n_before2 + n_after1) -
      boost::math::lgamma(n_before2 + n_after2 - 1 + model.nb_shape1) +
      boost::math::lgamma(n_before2 + n_after2) -
      boost::math::lgamma(n_before1 + n_after1 - 1 + model.nb_shape1) +
      boost::math::lgamma(n_before1 + n_after1);
    
    if (search_model_) {
      Rf_error("not implemented");
    }
    
  }
    break;
  case SWITCH1:
  {
    int n1 = correspondence.n_before(id1);
    int n2 = correspondence.n_after(id2);
    int n3 = correspondence.n_after(id1) - 1;
    
    alpha += boost::math::lgamma(n1 + n2 - 1 + model.nb_shape1) +
      boost::math::lgamma(n3 - 1 + model.nb_shape1) +
      boost::math::lgamma(n1 + n3) +
      boost::math::lgamma(n2) -
      boost::math::lgamma(n1 + n3 - 1 + model.nb_shape1) -
      boost::math::lgamma(n2 - 1 + model.nb_shape1) -
      boost::math::lgamma(n1 + n2) -
      boost::math::lgamma(n3);
    
    if (search_model_) {
      int id3 = correspondence.get_next(id1);
      std::unordered_map<int, int> nz12;
      std::unordered_map<int, int> nz13;
      std::unordered_map<int, int> nz2;
      std::unordered_map<int, int> nz3;
      int id;
      id = id1;
      while (id != 0) {
        nz12[data.loc[id-1]]++;
        id = correspondence.get_prev(id);
      }
      nz13 = nz12;
      id = id3;
      while (id != 0) {
        nz13[data.loc[id-1]]++;
        nz3[data.loc[id-1]]++;
        id = correspondence.get_next(id);
      }
      id = id2;
      while (id != 0) {
        nz12[data.loc[id-1]]++;
        nz2[data.loc[id-1]]++;
        id = correspondence.get_next(id);
      }
      
      alpha += -boost::math::lgamma(model.zeta + n1 + n2) -
        boost::math::lgamma(model.zeta + n3) +
        boost::math::lgamma(model.zeta + n1 + n3) +
        boost::math::lgamma(model.zeta + n2);
      
      for (auto it = nz13.begin(); it != nz13.end(); ++it) {
        alpha -= boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
        alpha += boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
      }
      for (auto it = nz2.begin(); it != nz2.end(); ++it) {
        alpha -= boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
        alpha += boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
      }
      for (auto it = nz3.begin(); it != nz3.end(); ++it) {
        alpha += boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
        alpha -= boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
      }
      for (auto it = nz12.begin(); it != nz12.end(); ++it) {
        alpha += boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
        alpha -= boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
      }
    }
  } 
    break;
  case SWITCH2:
  {
    int n1 = correspondence.n_before(id1);
    int n2 = correspondence.n_after(id2);
    int n3 = correspondence.n_before(id2) - 1;
    
    alpha += boost::math::lgamma(n1 + n2 - 1 + model.nb_shape1) +
      boost::math::lgamma(n3 - 1 + model.nb_shape1) +
      boost::math::lgamma(n2 + n3) +
      boost::math::lgamma(n1) -
      boost::math::lgamma(n2 + n3 - 1 + model.nb_shape1) -
      boost::math::lgamma(n1 - 1 + model.nb_shape1) -
      boost::math::lgamma(n1 + n2) -
      boost::math::lgamma(n3);
    
    if (search_model_) {
      int id3 = correspondence.get_prev(id2);
      std::unordered_map<int, int> nz12;
      std::unordered_map<int, int> nz32;
      std::unordered_map<int, int> nz3;
      std::unordered_map<int, int> nz1;
      int id;
      id = id3;
      while (id != 0) {
        nz32[data.loc[id-1]]++;
        id = correspondence.get_prev(id);
      }
      nz3 = nz32;
      id = id1;
      while (id != 0) {
        nz1[data.loc[id-1]]++;
        id = correspondence.get_prev(id);
      }
      nz12 = nz1;
      id = id2;
      while (id != 0) {
        nz12[data.loc[id-1]]++;
        nz32[data.loc[id-1]]++;
        id = correspondence.get_next(id);
      }
      
      alpha += -boost::math::lgamma(model.zeta + n1 + n2) -
        boost::math::lgamma(model.zeta + n3) +
        boost::math::lgamma(model.zeta + n2 + n3) +
        boost::math::lgamma(model.zeta + n1);
      
      for (auto it = nz32.begin(); it != nz32.end(); ++it) {
        alpha -= boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
        alpha += boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
      }
      for (auto it = nz1.begin(); it != nz1.end(); ++it) {
        alpha -= boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
        alpha += boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
      }
      for (auto it = nz3.begin(); it != nz3.end(); ++it) {
        alpha += boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
        alpha -= boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
      }
      for (auto it = nz12.begin(); it != nz12.end(); ++it) {
        alpha += boost::math::lgamma(it->second + model.zeta*data.loc_prior[it->first-1]);
        alpha -= boost::math::lgamma(model.zeta * data.loc_prior[it->first-1]);
      }
    }
  } 
    break;
    
  case SWAP:
  {
    int nn1 = correspondence.n_before(id1) + correspondence.n_after(id1) - 1;
    int nn2 = correspondence.n_before(id2) + correspondence.n_after(id2) - 1;
    
    alpha += -log(nn1 + model.nb_shape1 - 2) + log(nn1 - 1) - log(nn2) + log(nn2 + model.nb_shape1 - 1);
    
    if (search_model_) {
      int count_in_first = 0;
      int count_in_second = 0;
      
      int id;
      id = correspondence.get_prev(id1);
      while (id != 0) {
        if (data.loc[id-1] == data.loc[id1-1])
          count_in_first++;
        id = correspondence.get_prev(id);
      }
      id = correspondence.get_next(id1);
      while (id != 0) {
        if (data.loc[id-1] == data.loc[id1-1])
          count_in_first++;
        id = correspondence.get_next(id);
      }
      id = id2;
      while (id != 0) {
        if (data.loc[id-1] == data.loc[id1-1])
          count_in_second++;
        id = correspondence.get_prev(id);
      }
      id = correspondence.get_next(id2);
      while (id != 0) {
        if (data.loc[id-1] == data.loc[id1-1])
          count_in_second++;
        id = correspondence.get_next(id);
      }
      
      alpha += log(count_in_second + model.zeta * data.loc_prior[data.loc[id1-1]-1]) - 
        log(model.zeta + nn2);
      
      alpha -= (log(count_in_first + model.zeta * data.loc_prior[data.loc[id1-1]-1]) - 
        log(model.zeta + nn1 - 1));
    }
  }
    break;
  case UPDATE_IN:
    if (estimate_rho) {
      alpha += log(1 - model.rho) - log(model.rho);
      alpha += log(N_avail) - log(model.S - N_avail + 1);
    }
    
    alpha += log(N_avail + model.A - 1) - log(N_group[data.group[id1-1]-1] - 1 - model.discount);
    {
      int nn = correspondence.n_before(id2) + correspondence.n_after(id2) - 1;
      double p = 1 / (1 + model.nb_shape2);
      alpha -= (model.nb_shape1 * log(1 - p) - log(p)
                  + log(nn) - log(nn - 1 + model.nb_shape1));
      
      if (search_model_) {
        int count_in_second = 0;
        int id;
        id = id2;
        while (id != 0) {
          if (data.loc[id-1] == data.loc[id1-1])
            count_in_second++;
          id = correspondence.get_prev(id);
        }
        id = correspondence.get_next(id2);
        while (id != 0) {
          if (data.loc[id-1] == data.loc[id1-1])
            count_in_second++;
          id = correspondence.get_next(id);
        }
        
        alpha += log(count_in_second + model.zeta * data.loc_prior[data.loc[id1-1]-1]) - 
          log(model.zeta + nn);
        alpha -= log(data.loc_prior[data.loc[id1-1]-1]);
      }
    }
    break;
  case UPDATE_OUT:
    if (estimate_rho) {
      alpha += log(model.rho) - log(1 - model.rho);
      alpha += log(model.S - N_avail) - log(N_avail + 1);
    }
    
    alpha += log(N_group[data.group[id1-1]-1] - model.discount) - log(N_avail + model.A);
    {
      int nn = correspondence.n_before(id1) + correspondence.n_after(id1) - 1;
      double p = 1 / (1 + model.nb_shape2);
      alpha += model.nb_shape1 * log(1 - p) - log(p) + 
        log(nn - 1) - log(nn + model.nb_shape1 - 2);
      
      if (search_model_) {
        int count_in_first = 0;
        int id;
        id = correspondence.get_prev(id1);
        while (id != 0) {
          if (data.loc[id-1] == data.loc[id1-1])
            count_in_first++;
          id = correspondence.get_prev(id);
        }
        id = correspondence.get_next(id1);
        while (id != 0) {
          if (data.loc[id-1] == data.loc[id1-1])
            count_in_first++;
          id = correspondence.get_next(id);
        }
        
        alpha += -log(count_in_first + model.zeta * data.loc_prior[data.loc[id1-1]-1]) +
          log(model.zeta + nn - 1);
        alpha += log(data.loc_prior[data.loc[id1-1]-1]);
      }
    }
    break;
  case EXTRACT:
  {
    alpha += log(N_group[data.group[id1-1]-1] - model.discount) - log(N_avail + model.A);
    
    int n2 = correspondence.count_loc(data.loc[id1-1], id1, data);
    int n1 = correspondence.n_before(id1) + correspondence.n_after(id1) - 1 - n2;
    double p = 1.0 / (1 + model.nb_shape2);
    
    alpha += model.nb_shape1 * log(1 - p) - log(p) +
      boost::math::lgamma(n1-1 + model.nb_shape1) -
      boost::math::lgamma(n1) +
      boost::math::lgamma(n2-1 + model.nb_shape1) -
      boost::math::lgamma(n2) -
      boost::math::lgamma(n1 + n2 - 1 + model.nb_shape1) +
      boost::math::lgamma(n1 + n2) -
      boost::math::lgamma(model.nb_shape1);
    
    if (search_model_) {
      alpha -= (boost::math::lgamma(model.zeta + n2) - boost::math::lgamma(model.zeta));
      alpha += (boost::math::lgamma(model.zeta + n1 + n2) - boost::math::lgamma(model.zeta + n1));
    }
  }
    break;
  case ABSORB:
  {
    alpha += log(N_avail + model.A - 1) - log(N_group[data.group[id1-1]-1] - 1 - model.discount);
    
    int n1 = correspondence.n_before(id1) + correspondence.n_after(id1) - 1;
    int n2 = correspondence.n_before(id2) + correspondence.n_after(id2) - 1;
    double p = 1 / (1 + model.nb_shape2);
    
    alpha -= (model.nb_shape1 * log(1 - p) - log(p) +
      boost::math::lgamma(n1-1 + model.nb_shape1) -
      boost::math::lgamma(n1) +
      boost::math::lgamma(n2-1 + model.nb_shape1) -
      boost::math::lgamma(n2) -
      boost::math::lgamma(n1 + n2 - 1 + model.nb_shape1) +
      boost::math::lgamma(n1 + n2) -
      boost::math::lgamma(model.nb_shape1));
    
    if (search_model_) {
      alpha -= (boost::math::lgamma(model.zeta + n1 + n2) - boost::math::lgamma(model.zeta + n2));
      alpha += (boost::math::lgamma(model.zeta + n1) - boost::math::lgamma(model.zeta));
    }
  }
    break;
  case INVALID:
    Rf_error("tried to insert invalid trans type");
  }
  return(alpha);
}


// [[Rcpp::export]]
List mcmcem(int max_iter,
            int n_samples,
            bool search_model,
            bool estimate_rho,
            List params,
            IntegerVector group,
            IntegerVector group_sizes,
            IntegerVector group_start_inds,
            NumericVector time,
            IntegerVector loc,
            NumericVector loc_prior,
            IntegerVector init_L1,
            IntegerVector scoreSample_first,
            IntegerVector scoreSample_second) {
  
  Rcout << "initialize data structures" << std::endl;
  
  Rcout << "Data object" << std::endl;
  Data data(group, group_sizes, group_start_inds, time, 
            loc, loc_prior, scoreSample_first, scoreSample_second);
  
  Rcout << "Model object" << std::endl;
  Model model(data, params);
  
  Rcout << "MH_chain" << std::endl;
  MH_chain chain(data, model, init_L1);
  
  // randomly intialize correspondence
  if (init_L1.size() != data.n_obs) {
    Rcout << "start random initialization..." << std::endl;
    chain.reset(5e4, data, model);
    chain.sample(data, model, true, false, false, false);
    if (estimate_rho)
      model.max_exp_loglik(chain, data);
  }
  
  bool score = false;
  if (scoreSample_first.size() > 1)
    score = true;
  
  std::vector<double> rho_vec(max_iter);
  
  for (int i=0; i<max_iter; i++) {
    Rcout << "iter " << i << std::endl;
    chain.reset(n_samples, data, model);
    chain.sample(data, model, false, score, estimate_rho,search_model);
    if (estimate_rho && !score) {
      if (i % 2 == 0)
        model.max_exp_loglik(chain, data);
      rho_vec[i] = model.get_rho();
    }
  }
  Rcout << "Finished sampling" << std::endl;
  
  // return some stuff to R
  List out;
  out["N_group"] = chain.N_group;
  out["ScoreSample"] = chain.get_score_counts();
  out["L1"] = chain.get_L1();
  if (estimate_rho) {
    out["rho_vec"] = rho_vec;
  }
  
  return(out);
}
