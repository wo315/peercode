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


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  Point initial_location;       //< store the initial location 
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1), initial_location(0) {}
};
/** Custom structure of data to store with Edges*/
struct EdgeData {
  double K;       //< Edge K
  double leng;     //< Edge length
  EdgeData() : K(0), leng(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/** Force function object for HW2 #1. */
class Force {
  // * A parent class, for forces (hopefully for good) which can be applied the cloths of the world.
  // * Implement virtual Point operator()(NODE n, double t) it's children should.
  public:
    virtual Point operator()(Node n, double t) = 0;

    virtual void name() {
      std::cout<<"Force"<<std::endl;
    }
};

class GravityForce: public Force {
  public:
    virtual Point operator()(Node n, double t) {
      auto gravity = Point(0,0,grav);
      return -1 * gravity * n.value().mass;
    }

    virtual void name() {
      std::cout<<"GravityForce"<<std::endl;
    }
};

class MassSpringForce: public Force {
  public:
    virtual Point operator()(Node n, double t) {
      auto my_loc = n.position();
      Point net_force = Point(0,0,0);
      for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
        auto edge = (*e);
        auto L = edge.value().leng;
        auto K = edge.value().K;
        auto x_j = edge.node2();
        if (x_j == n) x_j = edge.node1();

        auto diff = my_loc - x_j.position();
        net_force -=  K * (diff/norm(diff))*(norm(diff) - L);
      }
      return net_force;
    }

    virtual void name() {
      std::cout<<"MassSpringForce"<<std::endl;
    }
};

class DampingForce: public Force {
  public:
    virtual Point operator()(Node n, double t) {
      float c = 0.0005;
      return -1 * c *n.value().vel;
    }

    virtual void name() {
      std::cout<<"DampingForce"<<std::endl;
    }
};

class CombinedForce {   
  public: 
    std::vector<Force*> forces; 
    CombinedForce(std::vector<Force*> fs) : forces(fs) {  } 
  
    Point operator () (Node n, double t) const { 
        auto net_force = Point(0);
        for (auto f : forces) {
          auto delta = (*f)(n, t);
          net_force += delta;
        }

        return net_force;
    } 
};

template <typename F1, typename F2>
CombinedForce make_combined_force(F1 f1, F2 f2) {
  std::vector<Force*> vect{&f1, &f2};
  auto combined_force = CombinedForce(vect);
  return combined_force;
}

template <typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 f1, F2 f2, F3 f3) {
  std::vector<Force*> vect{&f1, &f2, &f3};
  auto combined_force = CombinedForce(vect);
  return combined_force;
}


class Constraint {
public:
  virtual void operator() (GraphType& g, double t) {}
};

class PinConstraint: public Constraint {
public:
  virtual void operator() (GraphType& g, double t) {
    for (auto n_ptr = g.node_begin(); n_ptr != g.node_end(); ++n_ptr) {
      auto n = *n_ptr;
      auto p = n.value().initial_location;
      if (p == Point(0,0,0) || p == Point(1,0,0)) {
        n.position() = Point(p.x, p.y, p.z);
        n.value().vel = Point(0);
      }
    }
    
  }
};

class PlaneConstraint: public Constraint {
  float z = -0.75;
public:
  virtual void operator() (GraphType& g, double t) {
    for (auto n_ptr = g.node_begin(); n_ptr != g.node_end(); ++n_ptr) {
      auto n = *n_ptr;
      if (n.position().z < z) {
        n.position().z = z;
        n.value().vel.z = 0;
      }
    }
  }
};


class SphereConstraint: public Constraint {
private:
  Point c = Point(0.5, 0.5, -0.5);
  float r = 0.15;
public:
  virtual void operator() (GraphType& g, double t) {
    for (auto n_ptr = g.node_begin(); n_ptr != g.node_end(); ++n_ptr) {
      auto n = *n_ptr;
      if (norm(n.position() - c) < r) {
        auto unit_direction = (n.position() - c) / norm(n.position() - c);
        n.position() = c + unit_direction * r;
        n.value().vel = n.value().vel - dot(n.value().vel, unit_direction)* unit_direction;
      }
    }
  }
};

class TearSphereConstraint: public Constraint {
private:
  Point c = Point(0.5, 0.5, -0.5);
  float r = 0.15;

public:
  virtual void operator() (GraphType& g, double t) {
    for (auto n_ptr = g.node_begin(); n_ptr != g.node_end(); ++n_ptr) {
      auto n = *n_ptr;
      if (norm(n.position() - c) < r) {
        g.remove_node(n);
        n_ptr = g.node_begin();
      }
    }
  }
};


class CombinedConstraint {   
  public: 
    std::vector<Constraint*> constraints; 
    CombinedConstraint(std::vector<Constraint*> cs) : constraints(cs) {  } 
  
    virtual void operator() (GraphType& g, double t) const {
        for (auto c : constraints) {
          (*c)(g,t);
        }
    } 
};

class EmptyConstraint: public Constraint {
public:
  virtual void operator() (GraphType& g, double t) {}
};


template <typename C1, typename C2>
CombinedConstraint make_combined_constraint(C1 c1, C2 c2) {
  std::vector<Constraint*> vect{&c1, &c2};
  auto combined_constraint = CombinedConstraint(vect);
  return combined_constraint;
}

//--functionality_1
//--I get a segfault when using these, probably because the pointers are pointing at objects that fall out of range
//--START
template <typename C1, typename C2, typename C3>
CombinedConstraint make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  std::vector<Constraint*> vect{&c1, &c2, &c3};
  auto combined_constraint = CombinedConstraint(vect);
  return combined_constraint;
}
//--END

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


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. 
   * just 4 kicks, this is defined in problem 3, 
   * cuz we good docs write,
   * with this the spring-mass force being*/
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE

    if (n.position () ==  Point (0,0,0) || n.position () ==  Point (1,0,0)) return  Point (0,0,0);

    int K = 100;
    auto my_loc = n.position();
    Point net_force = Point(0,0,0);
    for (auto e = n.edge_begin(); e != n.edge_end(); ++e) {
      auto edge = (*e);
      auto L = edge.length();
      auto x_j = edge.node2();
      if (x_j == n) x_j = edge.node1();

      auto diff = my_loc - x_j.position();
      net_force -=  K * (diff/norm(diff))*(norm(diff) - L);
    }
    // std::cout<<net_force<<std::endl;
    auto gravity = Point(0,0,grav);
    net_force -= gravity * n.value().mass;
    // std::cout<<net_force<<std::endl;
    return net_force;

  }
};







