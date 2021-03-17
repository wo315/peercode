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
static constexpr double K=100;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double K; //spring constant
  double L; //length??
  EdgeData(): K(100), L(0) {}
};
// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
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

template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force,C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  //velocity not yet been updated
  // position updated

  constraint(g,t);
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

  //initalize this with k and l
  double K_;
  Problem1Force(const double K) : K_(K) {};
  Problem1Force(): K_(100) {};
  template <typename NODE>
  Point operator()(NODE n, double t) {
      // n =(0,0,0)
      // n=(0.25,0,0)
    // HW2 #1: YOUR CODE HERE
    if (n.position() == Point(0) || n.position()==Point(1,0,0))
      return Point(0);

    Point Fgravity = Point(0,0,-grav)*n.value().mass;
    Point Fspring = Point(0,0,0);
    Point xi = n.position();

    // each edge will add Fspring twice

    for (auto i = n.edge_begin(); i!= n.edge_end(); ++i)
    {
      Point xj = (*i).node2().position();
      Fspring += (-K_ *(xi-xj)/norm(xi-xj))*(norm(xi-xj)-(*i).value().L);
      /*
      std::cout<< "source node = " << xi[0] << " " << xi[1] << " " << xi[2] <<std::endl;
      std::cout<< "destination node " << xj[0] << " " <<xj[1] << " " << xj[2] << std::endl;
      std::cout<< "t = "<< t << " x = " << Fspring[0] << " y= " << Fspring[1] << "z = " << Fspring[2]<< std::endl;
      std::cout<< "k =" << K_ << "L = " <<  (*i).value().L << std::endl;
      std::cout<< "norm = "<< norm(xi-xj) << std::endl; */
    }
    (void) t;  // silence compiler warnings
    return Fgravity+Fspring;
  }
};

//--functionality_1
//--The dynamic of your plot is not what is expected
//--(not symmetric). You might have a force that is computed
//--wrong, or a constraint that enforces some position when
//--it shouldn't
//--END

// parent class for all force children
struct Force{
  Force(){}
  virtual Point operator()(Node n, double t){
    (void) t; (void) n;
    return Point(0);
  }
  virtual ~Force(){}
};

// Gravity force, inherit from Force

struct GravityForce: public Force{
  Point operator()(Node n, double t) { 
    (void) t; //prevents warning
    return Point(0,0,-grav)*n.value().mass;}
};

struct MassSpringForce: public Force{
  double K_;
  MassSpringForce(const double K) : K_(K) {};
  MassSpringForce(): K_(100) {}; 
  Point operator()(Node n, double t) {
    (void) t;
    Point Fspring = Point(0);
    Point xi = n.position();
    for (auto i = n.edge_begin(); i!= n.edge_end(); ++i)
    {
      Point xj = (*i).node2().position();
      Fspring += (-K_ *(xi-xj)/norm(xi-xj))*(norm(xi-xj)-(*i).value().L);
    }
    return Fspring;
  }
};

struct DampingForce: public Force{
  double c;
  DampingForce(): c(1.0) {};
  DampingForce(const double tmp): c(tmp) {};
  Point operator()(Node n, double t) {
      (void) t;
      return -c* n.value().vel;
  }
};

// doesn't work :/
/*
struct make_combined_force: public Force{
  Force f1,f2,f3;
  make_combined_force(Force* one,Force* two): f1(one),f2(two),f3(Force()) {};
  make_combined_force(Force one,Force two, Force three): f1(one),f2(two),f3(three) {};
  Point operator()(Node n, double t) {
      return f1(n,t)+f2(n,t)+f3(n,t);
  }
}; */

struct CombinedForces{
  std::vector<Force*> internal_vec;
  //constructor
  CombinedForces(std::vector<Force*> v):internal_vec(v){};

  Point operator()(Node n, double t) {
    Point p1 = Point(0);
    //iterate through vec and call function to add on to it
    for (auto iter = internal_vec.begin(); iter != internal_vec.end(); ++iter){
      p1 += (*(*iter))(n,t); //first star get us pointer to forces, second star gives us forces
    }
    return p1;
  }
};

//--style_1
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START
template <typename F1, typename F2>
CombinedForces make_combined_force(F1 one, F2 two){
    std::vector<Force*> v = {&one, &two}; //vector of pointer of forces
    return CombinedForces(v);
}
//--END


template <typename F1, typename F2, typename F3>
CombinedForces make_combined_force(F1 one, F2 two, F3 three){
    std::vector<Force*> v = {&one, &two, &three};
    return CombinedForces(v);
}

struct Constraint {
  Constraint(){}
  virtual void operator()(GraphType& g, double t){
    //do nothing
    (void) g; (void) t;
  }
  virtual ~Constraint(){};
};

struct PinConstraint: public Constraint{
  void operator()(GraphType& g, double t){
    (void) t;
    static bool firstTime = true;
    static int n1 = -1 , n2=-1;
    // n.position() += n.value().vel * dt;

    if (firstTime) {
      for (auto i = g.node_begin(); i!= g.node_end(); ++i)
      {
        if ((*i).position() == Point(0)){
          n1 = (*i).index();
        } 
        
        if ((*i).position()==Point(1,0,0))
          n2 = (*i).index();
      }
      firstTime = false;
    }
    else {
      Node(&g,n1).position() = Point(0);
      Node(&g,n2).position() = Point(1,0,0);
    }
  }
};

struct PlaneConstraint: public Constraint{
  void operator()(GraphType& g, double t){
    (void) t;
    for (auto i = g.node_begin(); i!= g.node_end(); ++i)
    {
      Node n = *i;
      if(dot(n.position(),Point(0,0,1))<-0.75){
        n.position() = Point(0,0,n.position()[2]);
        n.value().vel[2]=0;
      }
    }
  }
};

struct SphereConstraint: public Constraint{
  void operator()(GraphType& g, double t){
    (void) t;
    Point c = Point(0.5,0.5,-0.5);
    for (auto i = g.node_begin(); i!= g.node_end(); ++i)
    {
      Node n = *i;
      Point diff = n.position() - c;
      if(norm(diff)< 0.15){
        diff = diff/norm(diff);
        //set position
        n.position() = 0.15*diff+c;
        //set component of the vel that is normal to sphere surface to 0
        n.value().vel -= dot(n.value().vel, diff)*diff;
      }
    }
  }
};

struct HoleyConstraint: public Constraint{
  Point c;
  double r;
  HoleyConstraint() :c(Point(0.5,0.5,-0.5)), r(0.15) {};
  void operator()(GraphType& g, double t){
    (void) t;
//--design_1
//--You're skipping some nodes with this for loop
//--START
    for(auto n = g.node_begin(); n!= g.node_end(); ++n){
      if (norm((*n).position()-c)<r) g.remove_node(*n);
    }
//--END
  }
};

/*
struct CombinedConstraints{
  Constraint c1,c2,c3;
  CombinedConstraints(Constraint one, Constraint two): c1(one),c2(two),c3(Constraint()) {};
  CombinedConstraints(Constraint one,Constraint two,Constraint three): c1(one),c2(two),c3(three) {};
  void operator()(GraphType& g, double t){
    c1(g,t);
    c2(g,t);
    c3(g,t);
  }
}; */


struct CombinedConstraints{
  std::vector<Constraint*> internal_con;
  //constructor
  CombinedConstraints(std::vector<Constraint*> v) : internal_con(v){};

  void operator()(GraphType& g, double t){
    for (unsigned int i = 0; i < internal_con.size(); ++i){
      (*(internal_con[i]))(g, t);
    }
  }
};

template <typename C1, typename C2>
CombinedConstraints make_combined_constraints(C1 c1,C2 c2){
    std::vector<Constraint*> vec;
    vec = {&c1, &c2};
    return CombinedConstraints(vec);
}

template <typename C1, typename C2, typename C3>
CombinedConstraints make_combined_constraints(C1 c1,C2 c2 ,C3 c3){
    std::vector<Constraint*> vec;
    vec = {&c1, &c2, &c3};
    return CombinedConstraints(vec);
}
