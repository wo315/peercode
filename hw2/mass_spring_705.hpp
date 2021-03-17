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
  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, double>;
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



class Force{
public:
  template <typename Node>
  Force(Node p):force_(p){};
  virtual Point operator()(double t) {
    (void)t;
    return Point(0); }
protected:
  Node force_;
};

class GravityForce: public Force{
public:
  template <typename Node>
  GravityForce(Node n1): Force(n1){  };
  virtual Point operator()(double t){
    (void)t;
    return force_.value().mass * Point(0,0,-grav);
  }
};


class MassSpringForce: public Force{
public:
  template <typename Node>
  MassSpringForce(Node n1): Force(n1){};
  virtual Point operator()(double t){
    (void)t;
    unsigned K=100;
    Point f_spring=Point(0);
    Point n_pos=force_.position();
    for(auto it=force_.edge_begin(); it!=force_.edge_end(); ++it){
      Node j=(*it).node2();
      Point j_pos=j.position();
      double L=(*it).value();
      Point temp=-(K*(n_pos-j_pos)/norm(n_pos-j_pos))*(norm(n_pos-j_pos)-L);
      f_spring=f_spring+temp;
    }
    return f_spring;
  }
};



class DampingForce: public Force{
public:
  template <typename Node>
  DampingForce(Node n1): Force(n1), c(1.0){};
  DampingForce(Node n1, double x): Force(n1), c(x){};
  virtual Point operator()(double t){
    (void)t;
    return -force_.value().vel * c;
  }
private:
  double c;
};


struct CombinedForce{
  CombinedForce(std::vector<Force*> f): forcevec(f){}; 
  Point operator()(Node n, double t){
    (void)t;
    Point total=n.position();
    for (auto it=forcevec.begin(); it!=forcevec.end(); ++it){
      auto cur=*it;
      total=total+cur(1.0);
    }
    return total;
  }
  private:
    std::vector<Force*> forcevec;
};

template<typename force1, typename force2>
CombinedForce make_combined_force(force1 f1, force2 f2){
  std::vector<Force*> f;
  f.push_back(&f1);
  f.push_back(&f2);
  return CombinedForce(f);
}

template<typename force1, typename force2, typename force3>
CombinedForce make_combined_force(force1 f1, force2 f2, force3 f3){
  std::vector<Force*> f;
  f.push_back(&f1);
  f.push_back(&f2);
  f.push_back(&f3);
  return CombinedForce(f);
}







// This does not work. But I cannot figure out why..

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void)t;
    if(n.position()==Point(0,0,0)||n.position()==Point(1,0,0)){
      return Point(0,0,0);
    }
    unsigned K=100;
    Point f_spring=Point(0);
    Point n_pos=n.position();
    for(auto it=n.edge_begin(); it!=n.edge_end(); ++it){
      Node j=(*it).node2();
      Point j_pos=j.position();
      double L=(*it).value();
      Point temp=-(K*(n_pos-j_pos)/norm(n_pos-j_pos))*(norm(n_pos-j_pos)-L);
      f_spring=f_spring+temp;
    }
    Point cur=n.value().mass*Point(0,0,-grav);
    Point f=f_spring+cur;
    return f;
  }
};



