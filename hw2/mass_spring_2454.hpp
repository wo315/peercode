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
  Point orig;
  NodeData() : vel(0), mass(1) {}
  NodeData( double v, double m, Point o ) : vel(v), mass(m), orig(o) {}
};

struct EdgeData {
  double K;
  double L;
  EdgeData() : K(100), L(1) {}
  EdgeData( double k, double l) : K(k), L(l) {}
};


// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/**
Base Force Type
*/
struct Force {
  /**
  The () operator for the basic Force type applied to node @a n at time @a t. Meant to be trivial and have no real physical meaning
  @return a zero-valued force
  */
  virtual Point operator()( Node n, double t ) { (void) n, (void) t; return Point(0, 0, 0); }
};

/**
  A force representing gravitational force (z-direction)
*/
struct GravityForce : public Force {
  /**
    @return the gravitational force applied to Node @a n at time @a t (although time has no meaning in this scenario).
  */
  Point operator()(Node n, double t) { (void) t; return n.value().mass * Point(0, 0, -1*grav); };
};

/**
  A force representing Spring force on some node
*/
struct MassSpringForce : public Force {
  /**
  @return the mass spring force applied to Node @a n at time @a t. Rest length and spring constant are assumed to be stored in the edges incidient to @a n. 
  */
  Point operator()(Node n, double t) {
    (void) t;
    Point result = Point(0,0,0);
    for( auto it = n.edge_begin(); it != n.edge_end(); ++it ) {
      double dist = norm( (*it).node2().position() - n.position() );
      Point diff = n.position() - (*it).node2().position();
      
      result += -1*(*it).value().K * diff/dist *(dist - (*it).value().L);
    }
    return result;
  }

};


/**
  A force representing frictional force
*/
struct DampingForce : public Force {
  float c;
  /**
  Constructs a force class with a default friction constant 
  */
  DampingForce(float k = 0.1) : c(k) {}

  /**
  @return a frictional force applied to Node @a n at time @a t with the constructed constant
  */
  Point operator()( Node n, double t ) {
    (void) t;
    return -1*c*n.value().vel;
  }

};

/**
  Class that represents the combination of some forces
*/
struct CombinedForce {
  std::vector<Force*> forces;
  /**
  Construct a combined force with a vector of the forces that will compose it
  */
  CombinedForce( std::vector<Force*> f ) : forces(f) {}

  /**
  @return a force representing the combination of all forces passed in constructor
  */
  Point operator()(Node n, double t ) {
    Point result = Point(0, 0, 0);
    for(size_t i = 0; i < forces.size(); i++){
      result += (*forces[i])(n, t);
    }
    return result;
  }

};

/**
@return the combined force of @a f1 and @a f2
*/
CombinedForce make_combined_force( Force&& f1,  Force&& f2 ) {
  std::vector<Force*> forces;
  forces.push_back(&f1);
  forces.push_back(&f2);
  return CombinedForce( forces );
}

/**
@return the combined force of @a f1, @a f2, and @a f3
*/
CombinedForce make_combined_force( Force&& f1,  Force&& f2, Force&& f3 ) {
  std::vector<Force*> forces;
  forces.push_back(&f1);
  forces.push_back(&f2);
  forces.push_back(&f3);
  return CombinedForce( forces );
}


/**
Class representing the generic constraint
*/

struct Constraint {
  /** 
  The base constraint does nothing when using the () operator
  */
  virtual void operator()( GraphType* g, double t ) { (void)g, (void)t;}
};

/**
  A constraint that pins the points at (0,0,0) and (1,0,0) fixed in place
*/
struct PinConstraint : public Constraint {

  /**
  Checks the graph @a g at time @a t to see if anything violates the constraint and modifies points to ensure that it follows the constraints
  */

  void operator()( GraphType* g, double t ) {
    (void) t;
    for( auto it = g->node_begin(); it != g->node_end(); ++it ) {
      if( (*it).position() != Point(1,0,0) && (*it).value().orig == Point(1,0,0) ) (*it).position() = Point(1,0,0);
      if( (*it).position() != Point(0,0,0) && (*it).value().orig == Point(0,0,0) ) (*it).position() = Point(0,0,0);
    }
  }

};

/**
Checks to see if we are violating the plane constraint at some plane z = k value
*/
struct PlaneConstraint: public Constraint {
  float k;
  /**
  Sets up a plane constraint at the selected value Z = @a c 
  */
  PlaneConstraint( float c = -0.75 ) : k(c) {}

  /**
  Applies the plane constraint to @a g at time @a t and modifies points to not cross the plane
  */
  void operator()( GraphType* g, double t ) {
    (void) t;
    for( auto it = g->node_begin(); it != g->node_end(); ++it ) {
      if((*it).position().z < k) {
        (*it).position().z = k;
        (*it).value().vel.z = 0;
      }
    }
  }
};

/**
  Constraint that prevents nodes from entering sphere centered at (0.5, 0.5, -0.5) and radius 0.15
*/
struct SphereConstraint: public Constraint {

  /**
    Checks the sphere constraint for all nodes in graph @a g at time @a t and modifies them so that they do not enter sphere
  */
   void operator()( GraphType* g, double t ) {
    (void) t;
    Point c = Point(0.5, 0.5, -0.5);
    float r = 0.15;
    for( auto it = g->node_begin(); it != g->node_end(); ++it ) {
      if( norm( c - (*it).position() )  < r) {
        Point q = ((*it).position() - c)/norm(c - (*it).position());
        (*it).position() = c + r*q ;

        (*it).value().vel -= inner_prod( (*it).value().vel, q )*q;
      }
    }
  }

};

/**
  Constraint that deletes any node (and its edges) that comes in contact with a sphere centered at (0.5, 0.5, -0.5) and radius 0.15
*/
struct TearConstraint : public Constraint {
  /**
  Applies the TearConstraint to graph @a g at time @a t
  */
  void operator()( GraphType* g, double t ) {
    (void) t;
    Point c = Point(0.5, 0.5, -0.5);
    float r = 0.15;
    for(auto it = g->node_begin(); it != g->node_end(); ) {
      if( norm( c - (*it).position() )  < r) {
        g->remove_node(it);
        it = g->node_begin();
      }
      else{ ++it; }
    }
  }
};


/**
  A functor that represents the combination of Constraints
*/
struct CombinedConstraint {
  std::vector<Constraint*> constraints;
  /**
  Constructs a combined constraint with a vector of constraints
  */
  CombinedConstraint(std::vector<Constraint*> c) : constraints(c) {}
  /**
  Applies all constraints this was constructed with to @a g at time @a t
  */
  void operator()(GraphType* g, double t){
    for(size_t i = 0; i < constraints.size(); i++ ){
      (*constraints[i])(g, t);
    }
  }
};

/**
@return  a Combined Constraint with @a c1 and @a c2
*/
CombinedConstraint make_combined_constraint( Constraint&& c1, Constraint&& c2 ) {
  std::vector<Constraint*> constraints;
  constraints.push_back(&c1);
  constraints.push_back(&c2);
  return CombinedConstraint( constraints );
}

/**
@return  a Combined Constraint with @a c1 and @a c2
*/
CombinedConstraint make_combined_constraint( Constraint&& c1, Constraint&& c2, Constraint&& c3 ) {
  std::vector<Constraint*> constraints;
  constraints.push_back(&c1);
  constraints.push_back(&c2);
  constraints.push_back(&c3);
  return CombinedConstraint( constraints );
}




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
double symp_euler_step(G& g, double t, double dt, F force, C constraint = C()) {
  // Compute the t+dt position

  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocitys
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  constraint( &g, t );

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
    (void) t;
    if (n.position() ==  Point(0,0,0) || n.position() ==  Point(1,0,0))
      return Point(0, 0, 0);
    Point result = Point(0,0,0);
    for( auto it = n.edge_begin(); it != n.edge_end(); ++it ) {
      double dist = norm( (*it).node2().position() - n.position() );
      Point diff = n.position() - (*it).node2().position();
      
      result += -1*(*it).value().K * diff/dist *(dist - (*it).value().L);
    }

    return result + n.value().mass * Point(0, 0, -1*grav);
  }
};



