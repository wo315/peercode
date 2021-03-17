/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <thread>

#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"

//--functionality_2
//--simulation failed, let me know if it worked for you
//--END


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  Point initial_pos; //< Node's initial position
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double len;  //< Edge length
  double k;       //< Spring constant
  EdgeData() : len(1.0), k(100.0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

class Force {

  public:
    Force() {}
    virtual Point operator()(Node n_, double d_) {
      (void) n_;
      (void) d_;
      return Point(0,0,0);
    }

};

class GravityForce : public Force {

  public:
    GravityForce() {}
    Point operator() (Node n_, double d_) {
      (void) d_;
      return n_.value().mass*Point(0,0,-grav);
   }

};

class MassSpringForce : public Force {

  public:
    MassSpringForce() {}
    Point operator() (Node n_, double d_) {
      (void) d_;
      Point cur_node_pos = n_.position();
      double k;
      double length;
      Point frc = Point(0,0,0);

      for(auto it = n_.edge_begin(); it != n_.edge_end(); ++it) {
        length = (*it).value().len;
        k = (*it).value().k;
        Point adj_node_pos = (*it).node2().position();
        Point cur_frc = cur_node_pos - adj_node_pos;
        cur_frc = -k*(cur_frc/norm(cur_frc))*(norm(cur_frc) - length);
        frc += cur_frc;
      }
      return frc;
    }

};

class DampingForce : public Force {

  public:
    DampingForce(double c) : const_(c) {}
    Point operator() (Node n_, double d_) {
      (void) d_;
      Point vel = n_.value().vel;
      Point frc = -const_*vel;
      return frc;
    }

  private:
    double const_;

};

class CombinedForce {

  public:
    CombinedForce(std::vector<Force*> vf) : force_vec(vf) {}
    Point operator() (Node n_, double d_) {
      Point frc = Point(0,0,0);
      for(unsigned int i = 0; i < force_vec.size(); i++) {
        frc += (*force_vec[i])(n_, d_);
      }
      return frc;
    }

  private:
    std::vector<Force*> force_vec;

};

//--design_1
//--argument objects passed by value are destroyed after function finishes
//--START
template<typename force1, typename force2>
CombinedForce make_combined_force(force1 f1, force2 f2) {
  std::vector<Force*> vf{&f1, &f2};
  return CombinedForce(vf);
}
//--END

template<typename force1, typename force2, typename force3>
CombinedForce make_combined_force(force1 f1, force2 f2, force3 f3) {
  std::vector<Force*> vf{&f1, &f2, &f3};
  return CombinedForce(vf);
}




/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {

    (void) t;

    CombinedForce cf = make_combined_force(MassSpringForce(), GravityForce());
    return cf(n, t);

  }
};


class Constraint {

  public:
    Constraint() {}
    virtual void operator()(Graph<NodeData, EdgeData> &g, double t) {
      (void) g;
      (void) t;
    }

};


class PinConstraint : public Constraint {

  public:
    PinConstraint(Point p1, Point p2) : p1_(p1), p2_(p2) {}
    void operator()(Graph<NodeData, EdgeData> &g, double t) {
      for(auto it = g.node_begin(); it != g.node_end(); ++it) {
        if((*it).value().initial_pos == p1_) {
          (*it).position() = p1_;
        } else if ((*it).value().initial_pos == p2_) {
          (*it).position() = p2_;
        }
      }
    (void) t;
    }

  private:
    Point p1_;
    Point p2_;

};

class PlaneConstraint : public Constraint {

  public:
    PlaneConstraint(double z) : z_(z) {}
    void operator()(Graph<NodeData, EdgeData> &g, double t) {
      for(auto it = g.node_begin(); it != g.node_end(); ++it) {
        Point cur_pos = (*it).position();
        if(inner_prod(cur_pos, Point(0,0,1)) < z_) {
          (*it).value().vel[2] = 0;
          (*it).position()[2] = z_;
        }
      }
    (void) t;
    }

  private:
    double z_;

};

class SphereConstraint : public Constraint {

  public:
    SphereConstraint(Point c, double r) : center(c), radius(r) {}
    void operator()(Graph<NodeData, EdgeData> &g, double t) {
      for(auto it = g.node_begin(); it != g.node_end(); ++it) {
        Point cur_pos = (*it).position();
        double dist = norm(cur_pos - center);
        if(dist < radius) {
          Point r = (cur_pos - center)/dist;
          Point cur_vel = (*it).value().vel;
          (*it).value().vel = cur_vel - inner_prod(cur_vel, r)*r;
          (*it).position() = center + (radius/dist)*(cur_pos - center);
        }
      }
      (void) t;
    }

  private:
    Point center;
    double radius;

};

class TearConstraint : public Constraint {

  public:
    TearConstraint(Point c, double r) : center(c), radius(r) {}
    void operator()(Graph<NodeData, EdgeData> &g, double t) {
      for(auto it = g.node_begin(); it != g.node_end(); ){
        Point cur_pos = (*it).position();
        double dist = norm(cur_pos - center);
        if(dist < radius) {
          it = g.remove_node(it);
        } else {
          ++it;
        }
      }
      (void) t;
    }

  private:
    Point center;
    double radius;

};

class CombinedConstraints {

  public:
    CombinedConstraints(std::vector<Constraint*> vc) : constr_vec(vc) {}
    void operator() (Graph<NodeData, EdgeData> &g, double t) {
      for(unsigned int i = 0; i < constr_vec.size(); i++) {
        (*constr_vec[i])(g, t);
      }
      return;
    }

  private:
    std::vector<Constraint*> constr_vec;

};


template<typename constr1, typename constr2>
CombinedConstraints make_combined_constraints(constr1 c1, constr2 c2) {
  std::vector<Constraint*> vc{&c1, &c2};
  return CombinedConstraints(vc);
}
