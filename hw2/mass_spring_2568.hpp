/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <thread>
#include <cmath>

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
  Point initial;
  NodeData() : vel(0), mass(1), initial(0) {}
};


// Define the Graph type
using GraphType = Graph<NodeData,double>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using size_type = unsigned;

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

  //add the constraints to the graph
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}



/** Force object parent class.
 * @param[in] constant double indicating the constant value for the specific force
 * @param[in] n Node type indicating the specific node
 * @param[in] t double indicating time at which the force acts upon the object
 * return Point value indicating force vector
 */
class Force {
  public:
    Force(double val): constant(val) {}

    double constant;

    virtual Point operator()(Node n, double t) {
      (void)n;
      (void)t;
      return Point(0,0,0);
    }
};

/** Gravity Force child class
 * @param[in] grav default gravitational acceleration
 * @param[in] n Node type being acted on
 * @param[in] t double indicating time
 * @return gravitational force
 * @post force vector should be 0 in all components except z
 * @post z directional component should be negative
 */
class GravityForce: public Force {
  public:
    GravityForce(): Force(grav) {}

    Point operator()(Node n, double t) {
      double result = (-1)*n.value().mass*constant;
      (void)t;
      return Point(0,0,result);
    }
};


/** Mass Spring Force child class
 * @param[in] constant default spring constant, set to 100
 * @param[in] n Node type object
 * @param[in] t double indicating time
 * @return total force acted by springs along a node's edges
 * @pre n should have incident edges, i.e. n.degree() != 0
 * @pre edges should have initial length value
 * @post return value inversely related to length, i.e. as length increases,
 * force should become more negative
 */
class MassSpringForce: public Force {
  public:
    MassSpringForce(): Force(100) {}

    Point operator()(Node n, double t) {
      Point spring = Point(0,0,0);

      for(auto i = n.edge_begin(); i != n.edge_end(); ++i) {
        Edge e = *i;
        Point dist = e.node1().position() - e.node2().position();
        double length = std::sqrt(normSq(dist));
        spring += (-1)*constant*(1/length)*(length - e.value())*dist;
      }

      (void)t;
      return spring;
    }
};

/** Damping Force child class
 * @param[in] constant default damping constant, set to 1
 * @param[in] n Node type object
 * @param[in] t double indicating time
 * @return friction force created by damping
 * @pre damping constant should be greater than 0
 * @post damping force should be inversely related to velocity
 */
class DampingForce: public Force {
  public:
    DampingForce(): Force(1) {}

    Point operator()(Node n, double t) {
      (void)t;
      return (-1)*n.value().vel*constant;
    }
};


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);
    Point spring = Point(0,0,0);
    for(auto i = n.edge_begin(); i != n.edge_end(); ++i) {
      Edge e = *i;
      Point dist = e.node1().position() - e.node2().position();
      double length = std::sqrt(normSq(dist));
      spring += (-100)*(1/length)*(length - e.value())*dist;
    }
    spring.z -= n.value().mass*grav;
    (void)t;
    return spring;
  }
};


/** Combined Force functor
 * @param[in] input vector of Force pointers, pointing to Force parent or child classes
 * @param[in] n Node Type object
 * @param[in] t double indicating time
 * @return Point indicating sum of all forces acting upon the Node
 * @pre input vector should not be empty
 * @pre input vector must include pointers to Force or derived class objects
 @ @post return Point object
 */
class CombinedForce {
  public:
    CombinedForce(std::vector<Force*> input): forces(input) {}
    Point operator()(Node n, double t) {
      Point total = Point(0,0,0);
      for(unsigned i = 0; i < forces.size(); i++) {
        total += forces[i]->operator()(n,t);
      }
      return total;
    }
  private:
    std::vector<Force*> forces;
};

/** Function to create CombinedForce functor
 * @param[in] force Force or child class object, pass by reference
 * @return CombinedForce functor with input vector including all given forces
 * @pre must give three force or child objects
 * @post returns CombinedForce functor
 */
template<typename A, typename B, typename D>
CombinedForce make_combined_force(A &force1, B &force2, D &force3) {
  std::vector<Force*> inputs;

  Force *a;
  Force *b;
  Force *c;

  a = &force1;
  b = &force2;
  c = &force3;



  inputs.push_back(a);
  inputs.push_back(b);
  inputs.push_back(c);

  CombinedForce combine(inputs);

  return combine;
}



/** Constraint parent class
 * @param[in] g, Graph type object, pass by reference
 * @param[in] t, double indicating time
 * @return operations acting upon nodes of g
 * @post default constraint changes nothing
 */
class Constraint {
  public:
    Constraint() {}

    virtual void operator()(GraphType& g, double t) {
      (void)g;
      (void)t;
      return;
    }
};


/** Pin Constraint child class
 * @param[in] g, Graph type object, pass by reference
 * @param[in] t, double indicating time
 * @return fix nodes located at (0,0,0) and (1,0,0)
 * @post nodes with initial points at (0,0,0) and (1,0,0) have 0 velocity
 * @post indicated nodes should not change position at all
 * @post any other nodes are not affected
 */
class PinConstraint: public Constraint {
  public:
    PinConstraint(): Constraint() {}

    virtual void operator()(GraphType& g, double t) {
      for(auto it = g.node_begin(); it != g.node_end(); ++it) {
        Node n = *it;
        if(n.value().initial == Point(0,0,0) || n.value().initial == Point(1,0,0)) {
          n.position() = n.value().initial;
          n.value().vel = Point(0,0,0);
        }
      }
      (void)t;
      return;
    }
};

/** Plane Constraint child class
 * @param[in] g, Graph type object, pass by reference
 * @param[in] t, double indicating time
 * @return make sure nodes do not go past plane located at z = -0.75
 * @post nodes with incoming z position greater than -0.75 unchanged
 * @post nodes with incoming z position less than -0.75 have z set to -0.75
 */
class PlaneConstraint: public Constraint {
  public:
    PlaneConstraint(): Constraint() {}

    virtual void operator()(GraphType& g, double t) {
      for(auto it = g.node_begin(); it != g.node_end(); ++it) {
        Node n = *it;
        if(n.position().z < -0.75) {
          n.position().z = -0.75;
          n.value().vel.z = 0;
        }
      }
      (void)t;
      return;
    }
};


/** Sphere Constraint child class
 * @param[in] g, Graph type object, pass by reference
 * @param[in] t, double indicating time
 * @return make sure nodes do not go inside sphere located at
 * (0.5,0.5, -0.5) with radius 0.15
 * @post nodes that are greater than 0.15 units away from sphere
 * center are unchanged
 * @post nodes less than 0.15 units away from center have velocity
 * adjusted to be normal to sphere's surface
 */
class SphereConstraint: public Constraint {
  public:
    SphereConstraint(): Constraint() {}

    virtual void operator()(GraphType& g, double t) {
      for(auto it = g.node_begin(); it != g.node_end(); ++it) {
        Node n = *it;
        Point distvec = n.position() - Point(0.5,0.5,-0.5);
        double dist = std::sqrt(normSq(distvec));
        if(dist < 0.15) {
          n.position().x = (0.15/dist)*(n.position().x - 0.5) + 0.5;
          n.position().y = (0.15/dist)*(n.position().y - 0.5) + 0.5;
          n.position().z = (0.15/dist)*(n.position().z + 0.5) - 0.5;
          Point Ri = distvec/dist;
          n.value().vel = n.value().vel - (inner_prod(n.value().vel, Ri))*Ri;
        }
      }
      (void)t;
      return;
    }
};

/** Tear Constraint child class
 * @param[in] g, Graph type object, pass by reference
 * @param[in] t, double indicating time
 * @return remove nodes that fall within 0.15 of (0.5,0.5,-0.5)
 * @post nodes that are greater than 0.15 units away from center
 * unchanged
 * @post nodes within distance are removed, along with incident edges
 * @post nodes are permanently removed, nodes size and edges size reduced
 */
class TearConstraint: public Constraint {
  public:
    TearConstraint(): Constraint() {}

    virtual void operator()(GraphType& g, double t) {
      auto it = g.node_begin();
      while(it != g.node_end()) {
        Node n = *it;
        Point distvec = n.position() - Point(0.5,0.5,-0.5);
        double dist = std::sqrt(normSq(distvec));
        if(dist < 0.15) {
         it = g.remove_node(it);
        } else {
         ++it;
        }
      }
      (void)t;
      return;
    }
};


/** Combined Constraints functor
 * @param[in] inputs, vector of Constraint pointers
 * @param[in] g, Graph type object, pass by reference
 * @param[in] t, double indicating time
 * @return adjusted graph
 * @post invoke all constraints in input vector
 */
class CombinedConstraints {
  public:
    CombinedConstraints(std::vector<Constraint*> inputs): constraints(inputs) {}

    void operator()(GraphType& g, double t) {
      for(int i = 0; i < (int)constraints.size(); i++) {
        constraints[i]->operator()(g,t);
      }
      return;
    }
  private:
    std::vector<Constraint*> constraints;
};

/** Function to make combined constraints functor
 * @param[in] constraint, 3 constraint or child objects
 * @return CombinedConstraints functor with input vector
 * including all given constraints
 * @post Functor must have all given constraints
 */
template<typename H, typename J, typename K>
CombinedConstraints make_combined_constraints(H &constraint1, J &constraint2, K &constraint3) {

  std::vector<Constraint*> inputs;
  Constraint* a;
  Constraint* b;
  Constraint* c;

  a = &constraint1;
  b = &constraint2;
  c = &constraint3;

  inputs.push_back(a);
  inputs.push_back(b);
  inputs.push_back(c);

  CombinedConstraints combine(inputs);

  return combine;
}

//--functionality_0
//--Passed all tests!
//--END

//--design_0
//--Well designed!
//--END

//--style_0
//--Good coding style!
//--END

//--documentation_0
//--good
//--END


