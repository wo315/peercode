/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <algorithm>
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
static constexpr double damping_constant = 0.002;

/**@struct NodeData
 * @brief An object to hold node-specific data
 *  Allows each of our nodes to hold a number of node-specific attributes.
 *  Easily accessible via Node::value() function in the graph class.
 */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  bool pinned;
  Point pinned_position;
  /**Constructor for NodeData object*/
  NodeData() : vel(0), mass(1.0), pinned(false), pinned_position(0) {}
};

struct EdgeData {
  double L;     //< Rest length
  double K;     //< Spring constant
  // Constructor
  EdgeData() : L(1), K(100) {}
};


// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/**@brief Change a graph's nodes according to a step of the symplectic Euler
 *        method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G is the type of graph to operate with.
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


/**@brief New update function which utilizes a constraint object in addition to a force object.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G is the type of graph to operate with.
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is the type of constraint to apply at each step.
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

  //Apply constraints
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


//
// Constraints
//
/**@class Constraint
 * @brief A template for constraints
 */
class Constraint {
public:
    /** @function operator()
     * @param[in]     g   Graph
     * @param[in]     t   time
     * Virtual operator that does nothing
     */
    virtual void operator()(GraphType& g, double t){
        (void)g;
        (void)t;
    }
};

/** @class PinConstraint
 * @brief Ensures pinned points maintain the same position.
 *
 */
class PinConstraint : public Constraint {
public:
    /** @function operator()
     * @param[in,out]     g   Graph
     * @param[in]         t   time
     *
     * Iterates through all nodes in the graph  to ensure that the
     * pinned points maintain the same position. Manually resets each
     * pinned point to its pinned location, and sets velocity to zero.
     */
    virtual void operator()(GraphType& g, double t){
        (void)t;
        // Keep velocities of points at 0,0,0 and 1,0,0 equal to zero?
        auto first = g.node_begin();
        auto last = g.node_end();
        while (first != last){
            if ((*first).value().pinned == true){
                (*first).position() = (*first).value().pinned_position;
                (*first).value().vel = Point(0,0,0);
            }

            ++first;
        }
    }
};

/** @class PlaneConstraint
 * @brief Ensures that all points stay above the plane z = -0.75.
 *
 */
class PlaneConstraint : public Constraint {
public:
    PlaneConstraint(){
        plane_ = -0.75;
    }
    /** @function operator()
     * @param[in,out]     g   Graph
     * @param[in]         t   time
     * Iterates through all nodes in the graph  to ensure that
     * all nodes stay above the plane z = -0.75. All nodes that fall
     * below this plane are reset to lie on the plane, and their
     * velocities are set to zero.
     */
    virtual void operator()(GraphType& g, double t){
        (void)t;
        auto first = g.node_begin();
        auto last = g.node_end();
        while (first != last){
            if ((*first).position().z < plane_) {
                // "Fix this node"
                (*first).position().z = plane_;
                (*first).value().vel.z = 0.0;
            }
            ++first;
        }
    }

private:
    float plane_;
};


/** @class SphereConstraint
 * @brief Mimicks a blockaiding sphere in the graph.
 *
 */
class SphereConstraint : public Constraint {
public:
    SphereConstraint(){
        center_ = Point(0.5,0.5,-0.5);
        radius_ = 0.15;
    }

    /** @function operator()
     * @param[in,out]     g   Graph
     * @param[in]         t   time
     * Iterates through all nodes in the graph, and to each applies
     * behavior consistant with a blockaiding sphere. All nodes that
     * find themselves inside the sphere are reset to the closest
     * point on the exterior of the sphere, with their velocity
     * set to zero.
     */
    virtual void operator()(GraphType& g, double t){
        (void)t;
        auto first = g.node_begin();
        auto last = g.node_end();
        while (first != last){
            if (norm((*first).position() - center_) < radius_) {
                // "Fix this node"
                Point radial = (*first).position() - center_;
                float scal_dist = norm(radial);
                radial = radial * radius_/scal_dist;

                (*first).position() = radial + center_;
                Point Ri = ((*first).position() - center_) /
                    norm((*first).position() - center_);

                (*first).value().vel -=  dot((*first).value().vel,Ri)*Ri;
            }
            ++first;
        }
    }

protected:
    Point center_;
    float radius_;
};

/** @class RemovedSphereConstraint
 * @brief Removes all nodes that touch or enter a specified sphere
 *
 */
class RemovedSphereConstraint : public SphereConstraint {
public:
    RemovedSphereConstraint(){
        // Call the base constructor to initialize variables
        SphereConstraint();
    }

    /**@function operator()
     * @param[in,out]     g   Graph
     * @param[in]         t   time
     * Iterates through all nodes in the graph, and to each applies
     * behavior consistant with a "black hole" sphere that "eats"
     * a.k.a. removes all nodes that touch it, or find themselves
     * inside of it.
     */
    virtual void operator()(GraphType& g, double t){
        (void)t;
        auto first = g.node_begin();
        auto last = g.node_end();
        while (first != last){
            if (norm((*first).position() - center_) < radius_) {
                // "Fix this node" by removing it
                g.remove_node(*first);
            }
            ++first;
        }
    }
};


/** @struct CombinedConstraints
 * @brief Functor object to combine constarints
 *
 */
struct CombinedConstraints {
    /** Constructor using a vector of constraints */
    CombinedConstraints(std::vector<Constraint*> c){
        my_constraints_ = c;
    }

    /**@function operator()
     * @param[in,out]     g   Graph
     * @param[in]         t   time
     * Iterates through all constraints, applying each in turn.
     */
  void operator()(GraphType& g, double t) {
      (void)t;
      // Apply constraints
      std::vector<Constraint*>::iterator first = my_constraints_.begin();
      std::vector<Constraint*>::iterator last = my_constraints_.end();
      while (first != last) {
          // Apply individual constraint
          (*(*first))(g,t);
          ++first;
      }
  }

private:
    std::vector<Constraint*> my_constraints_;
};


/**@function make_combined_constraint
 * @param[in]     one   First Constraint
 * @param[in]     two   Second Constraint
 * @tparam[in]    con1 Type of first constraint
 * @tparam[in]    con2 Type of second constraint
 * @return returns a CombinedConstraint Functor object with
 * two forces.
 *@pre Constraints are both children of the Constraint class
 */
template <typename con1, typename con2>
CombinedConstraints make_combined_constraint(con1 one, con2 two){
    std::vector<Constraint*> Constraints; // Vector of pointers so we don't lose info.
    Constraints.push_back(&one); // Add pointer
    Constraints.push_back(&two); // Add pointer
    return CombinedConstraints(Constraints);
}

/**@function make_combined_constraint
 * @param[in]     one   First Constraint
 * @param[in]     two   Second Constraint
 * @tparam[in]    con1 Type of first constraint
 * @tparam[in]    con2 Type of second constraint
 * @tparam[in]    con3 Type of third constraint
 * @return returns a CombinedConstraint Functor object with
 * three forces.
 * @pre Constraints are all children of the Constraint class
 */
template <typename con1, typename con2, typename con3>
CombinedConstraints make_combined_constraint(con1 one, con2 two, con3 three){
    std::vector<Constraint*> Constraints; // Vector of pointers so we don't lose info.
    Constraints.push_back(&one); // Add pointer
    Constraints.push_back(&two); // Add pointer
    Constraints.push_back(&three);
    return CombinedConstraints(Constraints);
}


//
// Forces
//

/** @class Force
 * @brief A template for Forces
 * This is a legitimate force, but always returns a zero force.
 */
class Force {
public:

    /**@function operator()
     * @brief A zero force.
     * @param[in,out]     n   Node
     * @param[in]         t   time
     * @return A Point representing the zero vector.
     */
    virtual Point operator()(Node n, double t){
        (void)n;
        (void)t;  // No-op cast to silence compiler warnings
        return Point(0,0,0);

    }
};

/** @class GravityForce
 * @brief Applies a gravity force in the -z direction to a node.
 */
class GravityForce : public Force {
public:

    /**@function operator()
     * @brief Applies a gravity force in the -z direction.
     * @param[in,out]     n   Node
     * @param[in]         t   time
     * @return Point representing a gravity force vector
     */
    virtual Point operator()(Node n, double t) {
        (void)t; // No-op cast to silence compiler warnings
        return Point(0,0, n.value().mass * -grav);
    }
};

/** @class MassSpringForce
 * @brief Edges now mimick springs, with the appropriate reaction forces.
 */
class MassSpringForce : public Force {
public:

    /**@function operator()
     * @brief Applies a spring force.
     * @param[in,out]     n   Node
     * @param[in]         t   time
     * @return Point vector representing the net spring force on node n.
     * This operator iterates through all incident edges to this node,
     * and treats each as if it was a spring. It reads in each "springs"
     * rest length and spring constant, and then calculates the reaction
     * force on node n. Summing across all nodes, it returns the net
     * spring force on node n.
     */
    virtual Point operator()(Node n, double t) {
        (void)t; // No-op cast to silence compiler warnings
        Point force = Point(0,0,0);
        GraphType::IncidentIterator first = n.edge_begin();
        GraphType::IncidentIterator last = n.edge_end();
        Point disp = Point(0,0,0);
        while (first != last){
            disp = n.position() - (*first).node2().position();
            force += -(*first).value().K * disp/norm(disp) *
                            (norm(disp) - (*first).value().L);
            ++first;
        }
        return force;
    }
};

/** @class DampingForce
 * @brief Applies a damping force to a node.
 */
class DampingForce : public Force {
public:
    /**@function operator()
     * @brief Applies a damping force.
     * @param[in,out]     n   Node
     * @param[in]         t   time
     * @return Point representing the damping force on node n.
     * This operator() applies a constant damping force relative to
     * the velocity of the node. Force is always in the opposite
     * direction to the velocity of the node.
     */
    virtual Point operator()(Node n, double t) {
        (void)t; // No-op cast to silence compiler warnings
        return -damping_constant * n.value().vel;
    }
};

/** @struct CombinedForce
 * @brief Creates a combined force from a given series of individual forces.
 */
struct CombinedForce {
    /**@function CombinedForce
     * @brief Default constructor
     */
    CombinedForce(){
    }

    /**@function CombinedForce
     * @brief Constructor
     * @param[in]     a   Vector of force object pointers
     * @param[in]     t   time
     * Generates a CombinedForce object using the
     * provided vector of pointers to forces.
     */
    CombinedForce(std::vector<Force*> a){
        my_forces = a;
    }

    /**@function operator()
     * @brief Combines all forces
     * @param[in,out]     n   Node
     * @param[in]         t   time
     * @return Net force vector on node n
     * This operator() iterates through all forces provided
     * in the CombinedForce constructor, summing each force
     * and returning the net resulting force on node n.
     */
  template <typename NODE>
  Point operator()(NODE n, double t) {
      (void)t;
      // Need constraints to make this act like Problem1Force
      Point outputforce = Point(0,0,0);
      std::vector<Force*>::iterator first = my_forces.begin();
      std::vector<Force*>::iterator last = my_forces.end();
      while (first != last) {
          outputforce += (*(*first))(n,t);
          ++first;
      }
      return outputforce;
  }

private:
    std::vector<Force*> my_forces;
};

/**@function make_combined_force
 * @brief Combines two forces.
 * @param[in]     one         First force
 * @param[in]     two         Second force
 * @tparam[in]    force1  Type of first force
 * @tparam[in]    force2  Type of second force
 * @return A CombinedForce object
 * Generates a combined force object using the provided forces.
 */
//--design_1
//--argument variables are destroyed after function execution if pass by value
//--START
template <typename force1, typename force2>
CombinedForce make_combined_force(force1 one, force2 two){
    std::vector<Force*> Forces; // Vector of pointers so we don't lose info.
    Forces.push_back(&one); // Add pointer
    Forces.push_back(&two); // Add pointer
    return CombinedForce(Forces);
}
//--END

/**@function make_combined_force
 * @brief Combines three
 * @param[in]     one          First force
 * @param[in]     two          Second force
 * @param[in]     three     Three force
 * @tparam[in]    force1   Type of first force
 * @tparam[in]    force2   Type of second force
 * @tparam[in]    force3   Type of thirdforce
 * @return A CombinedForce object
 * Generates a combined force object using the provided forces.
 */
template <typename force1, typename force2, typename force3>
CombinedForce make_combined_force(force1 one, force2 two, force3 three){
    std::vector<Force*> Forces;
    Forces.push_back(&one);
    Forces.push_back(&two);
    Forces.push_back(&three);
    return CombinedForce(Forces); // Pass the vector of pointers
}

/** @struct Problem1Force
 * @brief Basic mass-spring force plus gravity force.
 */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */

  /**@function operator()
   * @brief Applies mass-spring forces and gravity forces.
   * @param[in,out]     n   Node
   * @param[in]         t   time
   * @return Net force vector on node n
   * This operator() iterates through all nodes and applies a
   * mass-spring force and a gravity force to each.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // Points at (0,0,0) and (1,0,0) never move. Return zero force
    if (n.position() == Point(0,0,0) ||
        n.position() == Point(1,0,0)){
        return Point(0,0,0);
    }

    // Each node has a number of incident edges
    // Each edge
    Point force = Point(0,0, n.value().mass * -grav);

    GraphType::IncidentIterator first = n.edge_begin();
    GraphType::IncidentIterator last = n.edge_end();
    Point disp = Point(0,0,0);
    // Changed implementation to work with EdgeData structure
    while (first != last){
        disp = n.position() - (*first).node2().position();
        force += -(*first).value().K * disp/norm(disp) *
                        (norm(disp) - (*first).value().L);
        ++first;
    }
    (void)t;
      return force;
  }
};
