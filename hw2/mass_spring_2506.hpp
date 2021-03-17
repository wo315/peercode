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
  double mass;     //< Node mass
  Point initial_pos;
  NodeData() : vel(0), mass(1), initial_pos(0) {}
};

struct EdgeData {
  double K; //spring constant
  double L; //rest length
  EdgeData() : K(100), L(0) {}
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


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    // (void) n; (void) t; (void) grav;    // silence compiler warnings
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }
    Point f_spring = Point(0);
    for (auto it = n.edge_begin(); it != n.edge_end();++it){
      auto e = *it;
      auto x1 = e.node1().position();
      auto x2 = e.node2().position();
      Point diff = x1 - x2;
      double dist = norm_2(diff);
      f_spring += (-(e.value().K)*(dist - e.value().L)/dist)*diff;
    }
    return f_spring + Point(0,0,-grav*n.value().mass);
  }
};

class Force {
  public:
    // virtual Point operator()(Node n, double t){
    //   return Point(0);
    // }
    virtual Point operator()(Node n, double t) = 0;
};

class GravityForce: public Force {
  public:
    virtual Point operator()(Node n, double t){
      (void) t;
      //std::cout << Point(0,0,-grav*n.value().mass) << std::endl;
      return Point(0,0,-grav*n.value().mass);
    }

    virtual ~GravityForce() {}
};

class MassSpringForce: public Force {
  public:
    virtual Point operator()(Node n, double t){
      Point f_spring = Point(0);
      for (auto it = n.edge_begin(); it != n.edge_end();++it){
        auto e = *it;
        auto x1 = e.node1().position();
        auto x2 = e.node2().position();
        Point diff = x1 - x2;
        double dist = norm_2(diff);
        f_spring += (-(e.value().K)*(dist - e.value().L)/dist)*diff;
      }
      return f_spring;
    }

  virtual ~MassSpringForce() {}

};

class DampingForce: public Force {
  private:
    double c;

  public:
    virtual Point operator()(Node n, double t){
      Point f_damping = Point(0);
      f_damping += (-c)*n.value().vel;
      return f_damping;
    }

    DampingForce() : c(1) {}

  virtual ~DampingForce() {}

};

struct CombinedForce {
  private:
    std::vector<Force*> forces;

  public:
    virtual Point operator()(Node n, double t){
      Point combined_force = Point(0);
      for (auto it = forces.begin(); it != forces.end(); ++it){
        combined_force += (*(*it))(n,t);
      }
      return combined_force;
    }

    CombinedForce(std::vector<Force*> forces_) : forces(forces_) {}

    virtual ~CombinedForce() {}
};

//--design_1
//--argument objects are destroyed after function finishes
//--START
template<typename F1, typename F2>
CombinedForce make_combined_force(F1 f1, F2 f2){
  std::vector<Force*> forces;
  forces.push_back(&f1);
  forces.push_back(&f2);
  //forces.push_back(f3);
  return CombinedForce(forces);
}
//--END

template<typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 f1, F2 f2, F3 f3){
  std::vector<Force*> forces;
  forces.push_back(&f1);
  forces.push_back(&f2);
  forces.push_back(&f3);
  return CombinedForce(forces);
}

class Constraint {
  public:
    // virtual Point operator()(Node n, double t){
    //   return Point(0);
    // }
    virtual void operator()(GraphType& g, double t) = 0;
};

class PinConstraint: public Constraint {
  public:
    virtual void operator()(GraphType& g, double t){
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
          auto n = *it;
          if (( n.value().initial_pos == Point(0,0,0)) || (n.value().initial_pos == Point(1,0,0))){
            n.position() = n.value().initial_pos;
          }
        }
    }
};

class PlaneConstraint: public Constraint {
  public:
    virtual void operator()(GraphType& g, double t){
      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        auto n = *it;
        if (dot(n.position(),Point(0,0,1)) < -0.75){
          n.position().z = -0.75;
          n.value().vel.z = 0;
        }
      }
    }
};

class SphereConstraint: public Constraint {
  public:
    virtual void operator()(GraphType& g, double t){
      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        auto n = *it;
        Point c = Point(0.5,0.5,-0.5);
        if (norm(n.position()-c) < 0.15){
          Point R = (n.position() - c)/norm(n.position() - c);
          n.position() = c + 0.15*R;
          n.value().vel -= dot(n.value().vel, R)*R;
        }
      }
    }
};

class TearConstraint: public Constraint {
  public:
    virtual void operator()(GraphType& g, double t){
      //std::cout  << "start" << std::endl;
      Point c = Point(0.5,0.5,-0.5);
      auto it = g.node_begin();
      while (it != g.node_end()){
        auto n = *it;
        if (norm(n.position()-c) < 0.15){
          //std::cout << "before" << std::endl;
          g.remove_node(it);
          //std::cout << "after" << std::endl;
        }
        else {
          ++it;
        }
      }
      //std::cout << "end" << std::endl;
      // std::cout << nodes_to_del.size() << std::endl;
      // for (auto it2 = nodes_to_del.begin(); it2 != nodes_to_del.end(); ++it2){
      //   if ( ((*(*it2)).value().initial_pos == Point(0,0,0)) || ((*(*it2)).value().initial_pos == Point(1,0,0))){
      //     std::cout << "hi" << std::endl;
      //   }
      //   g.remove_node(*(*it2));
      // }
    }
};

struct CombinedConstraints {
  private:
    std::vector<Constraint*> constraints;

  public:
    virtual void operator()(GraphType& g, double t){
      for (auto it = constraints.begin(); it != constraints.end(); ++it){
        (*(*it))(g,t);
      }
    }

    CombinedConstraints(std::vector<Constraint*> constraints_) : constraints(constraints_) {}

    virtual ~CombinedConstraints() {}
};

template<typename C1, typename C2>
CombinedConstraints make_combined_constraint(C1 c1, C2 c2){
  std::vector<Constraint*> constraints;
  constraints.push_back(&c1);
  constraints.push_back(&c2);
  return CombinedConstraints(constraints);
}

template<typename C1, typename C2, typename C3>
CombinedConstraints make_combined_constraint(C1 c1, C2 c2, C3 c3){
  std::vector<Constraint*> constraints;
  constraints.push_back(&c1);
  constraints.push_back(&c2);
  constraints.push_back(&c3);
  return CombinedConstraints(constraints);
}
