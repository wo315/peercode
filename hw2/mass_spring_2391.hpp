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
static constexpr double k_const = 100;
static constexpr double l_const = 0.01;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  bool pin1;      //true if it's position == Point(0,0,0)
                  //and needs to be pined in PinConstraint
  bool pin2;      //true if it's position == Point(1,0,0)
                  // and needs to be pined in PinConstraint
    /**
     @brief: constructor of NodeData
     */
  NodeData() : vel(0), mass(1), pin1(false), pin2(false) {}
};

struct EdgeData {
    double L;
    double K;
    /**
     @brief: constructor of EdgeData
     */
    EdgeData():L(1),K(100){}
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
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
      if (n.position()!=Point(0,0,0) && n.position()!=Point(1,0,0))
          n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (n.position()!=Point(0,0,0) && n.position()!=Point(1,0,0))
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
      //if the nodes are pined, fix them at the original position and return zero force
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
        return Point(0,0,0);
    Point f_grav = n.value().mass*Point(0,0,-grav);
    Point f_spr = Point(0,0,0);
    Point xi = n.position();
      //loop over all incident edges of n, calculate spring force
    for (auto it = n.edge_begin(); it!=n.edge_end(); ++it)
    {
//            f_spr += -(*it).value().K*(n.position()-(*it).node2().position())
//            *((*it).length()-(*it).value().L)/((*it).length());
          double K = (*it).value().K;
          double L = (*it).value().L;
          Point xj = (*it).node2().position();
          double l2 = (*it).length();
          f_spr += (-K * (xi - xj) / l2) * (l2 - L);
    }
    return f_spr + f_grav;


  }
};

/** Base force class for HW2 #1. */
class Force
{
public:
    /**
        operator() of Force, returns a zero force
     @param[in] n node the force to apply to
     @param[t] time t
     @return zero force
     */
    virtual Point operator()(Node n, double t)
    {
        (void) n;
        (void) t;
        return Point(0,0,0);

    }
    /**
        destructor of force
     */
    virtual ~Force(){};
};

/**
 Gravity Force class for calculating the gravity
 */
class GravityForce: public Force
{
public:
    /**
        operator() of GravityForce, returns a gravity force
     @param[in] n node the force to apply to
     @param[t] time t
     @return Gravity Force
     */
    virtual Point operator()(Node n, double t)
    {
        (void) t;
        Point f_grav = n.value().mass*Point(0,0,-grav);
        return f_grav;
    }
    /**
     copy constructor
     @param[in] copy_c a GravityForce object
     */
    GravityForce(const GravityForce& copy_c){(void) copy_c;};
    //default constructor
    GravityForce(){};
};

/**
 Mass Spring Force class for calculating the spring force
 */
class MassSpringForce: public Force
{
public:
    /**
        operator() of MassSpringForce, returns a spring force
     @param[in] n node the force to apply to
     @param[t] time t
     @return Spring Force
     */
    virtual Point operator()(Node n, double t)
    {
        (void) t;
        Point f_spr = Point(0,0,0);
        Point xi = n.position();
        //loop over all incident edges of n
        for (auto it = n.edge_begin(); it!=n.edge_end(); ++it)
        {
//            f_spr += -(*it).value().K*(n.position()-(*it).node2().position())
//            *((*it).length()-(*it).value().L)/((*it).length());
            double K = (*it).value().K;
            double L = (*it).value().L;
            Point xj = (*it).node2().position();
            double l2 = (*it).length();
            f_spr += (-K * (xi - xj) / l2) * (l2 - L);
        }
        return f_spr;
    }
    /**
     copy constructor
     @param[in] copy_c a MassSpringForce object
     */
    MassSpringForce(const MassSpringForce& copy_c){(void) copy_c;};
    //default constructor
    MassSpringForce(){};
};
/**
 DampingForce class for calculating the damping force
 */
class DampingForce: public Force
{
    double c;
public:
    /**
     constructor
     @param c_ set the value of variable c
     */
    DampingForce(double c_):c(c_){};
    /**
     copy constructor
     @param[in] copy_c a DampingForce object
     */
    DampingForce(const DampingForce& copy_c):c(copy_c.c){};
    DampingForce():c(0){};
    /**
        operator() of DampingForce, returns a damping force
     @param[in] n node the force to apply to
     @param[t] time t
     @return Damping Force
     */
    virtual Point operator()(Node n, double t)
    {
        (void) t;
        Point f_damp = (-c)*n.value().vel;
//        std::cout<<"damp Force "<<std::endl;
        return f_damp;
    }
};

/**
 CombinedForce functor to sum/combine different forces
 */
class CombinedForce
{
    std::vector<Force*> force_vec;
    public:
    /**
     constructor
     @param[in] vec set  force_vec = vec
     */
        CombinedForce(const std::vector<Force*>& vec):force_vec(vec){}
    /**
        operator() of CombinedForce, returns a combined force
     @param[in] n node the force to apply to
     @param[t] time t
     @return combined Force
     */
    Point operator()(Node n, double t)
    {
        (void) t;
        Point sum_force = Point(0,0,0);
        //loop over the vector to add/combine forces
        for (unsigned i = 0; i<force_vec.size(); i++)
        {
            sum_force += force_vec[i]->operator()(n, t);
        }
        return sum_force;
    }
    /**
     Destructor, release all memory used on heap
     */
    ~CombinedForce()
    {
        //loop over the vector to release heap memory used
        for (unsigned i = 0; i<force_vec.size(); i++)
        {
            delete force_vec[i];
        }
    }
};
/**
 combine three forces of potentially different types in arbitrary order
 @param force1 the first force
 @param force2 the second force
 @param force3 the third force
 @return a CombinedForce object that can be used to sum all the three forces
 */
template <typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 force1, F2 force2, F3 force3)
{

    Force* f1 = new F1(force1);
    Force* f2 = new F2(force2);
    Force* f3 = new F3(force3);
    std::vector<Force*> vec {f1, f2, f3};
    return CombinedForce(vec);
}
/**
 combine two forces of potentially different types in arbitrary order
 @param force1 the first force
 @param force2 the second force
 @return a CombinedForce object that can be used to sum all the two forces
 */
template <typename F1, typename F2>
CombinedForce make_combined_force(F1 force1, F2 force2)
{

    Force* f1 = new F1(force1);
    Force* f2 = new F2(force2);
    std::vector<Force*> vec {f1, f2};
    return CombinedForce(vec);
}

/**
 Base Constraint class
 */
class Constraint
{
    public:
    /**
        operator() of Constraint,
     @param[in] gra the graph to be applied constraints to
     @param[t] time t
     */
    virtual void operator()(GraphType& gra, double t){
        (void) gra;
        (void) t;
    }
    virtual ~Constraint(){};
};

/**
 PinConstraint class: fix the position of nodes with initial position P ==Point(0,0,0) || P==Point(1,0,0)
 */
class PinConstraint: public Constraint
{
    public:
    /**
        operator() of Constraint,
     @param[in] gra the graph to be applied constraints to
     @param[t] time t
     */
        void operator()(GraphType& gra, double t)
        {
            (void) t;
            //loop over all nodes in the graph, and fix the nodes with
            //initial position P ==Point(0,0,0) || P==Point(1,0,0)
            for (auto it = gra.node_begin(); it!=gra.node_end();++it)
            {
                if((*it).value().pin1)
                    (*it).position() =  Point(0,0,0);
                if((*it).value().pin2)
                    (*it).position() =  Point(1,0,0);
            }
        }
    /**
     copy constructor
     @param[in] copy_c a PinConstraint object
     */
    PinConstraint(const PinConstraint& copy_c){(void) copy_c;};
    //default constructor
    PinConstraint(){};

};
/**
 PlaneConstraint class: move the nodes violate the limit to the closest point on plane z=-0.75
 */
class PlaneConstraint: public Constraint
{
    double limit_plane = -0.75;
    public:
    /**
     copy constructor
     @param[in] copy_c a PlaneConstraint object
     */
        PlaneConstraint(const PlaneConstraint& copy_c):limit_plane(copy_c.limit_plane){};
    //default constructor
        PlaneConstraint():limit_plane(-0.75){};
    /**
        operator() of Constraint,
     @param[in] gra the graph to be applied constraints to
     @param[t] time t
     */
        void operator()(GraphType& gra, double t)
        {
            (void) t;
            //loop over all nodes in the graph, and move the nodes violates the constraint
            //to the closest point on plane z=-0.75
            for (auto it = gra.node_begin(); it!=gra.node_end();++it)
            {
                if ((*it).position().z<limit_plane)
                {
                    (*it).value().vel.z = 0;
                    (*it).position().z = limit_plane;

                }
            }
        }
};
/**
 SphereConstraint class
 */
class SphereConstraint: public Constraint
{
    Point center_c = Point(0.5,0.5,-0.5);
    double radius = 0.15;
    public:
    /**
     copy constructor
     @param[in] copy_c a SphereConstraint object
     */
    SphereConstraint(const SphereConstraint& copy_c):center_c(copy_c.center_c),radius(copy_c.radius){};
    //default constructor
    SphereConstraint():center_c(Point(0.5,0.5,-0.5)),radius(0.15){};
    /**
        operator() of Constraint,
     @param[in] gra the graph to be applied constraints to
     @param[t] time t
     */
        void operator()(GraphType& gra, double t)
        {
            (void) t;
            //loop over all the nodes in the graph
            for (auto it = gra.node_begin(); it!=gra.node_end();++it)
            {
                Point xi = (*it).position();
                Point Ri = (xi-center_c)/norm(xi-center_c);
                //if the node violates the constraint
                if (norm(xi-center_c)<radius)
                {
                    (*it).value().vel = (*it).value().vel-inner_prod((*it).value().vel,Ri)*Ri;
                    (*it).position() = center_c + radius*(xi-center_c)/norm(xi-center_c);
                }
            }
        }
};
/**
 SphereRmConstraint class
 */
class SphereRmConstraint: public Constraint
{
    Point center_c = Point(0.5,0.5,-0.5);
    double radius = 0.15;
    public:
    /**
     copy constructor
     @param[in] copy_c a SphereRmConstraint object
     */
        SphereRmConstraint(const SphereRmConstraint& copy_c):center_c(copy_c.center_c),radius(copy_c.radius){};
    //default constructor
        SphereRmConstraint():center_c(Point(0.5,0.5,-0.5)),radius(0.15){};
    /**
        operator() of Constraint,
     @param[in] gra the graph to be applied constraints to
     @param[t] time t
     */
        void operator()(GraphType& gra, double t)
        {
            (void) t;
            //loop over all the nodes in the graph
            for (auto it = gra.node_begin(); it!=gra.node_end();++it)
            {
                Point xi = (*it).position();
                //if the node violates the constraint
                if (norm(xi-center_c)<radius)
                {
                    gra.remove_node(*it);
                }
            }

        }

};

/**
 CombinedConstraints class, combines constraints together
 */
class CombinedConstraints
{
    std::vector<Constraint*> vec;
    public:
    /**
     constructor
     @param[in] vec_ set vec=vec_
     */
        CombinedConstraints(const std::vector<Constraint*>& vec_):vec(vec_){};
    /**
        operator() of Constraint,
     @param[in] gra the graph to be applied constraints to
     @param[t] time t
     */
        void operator()(GraphType& gra, double t)
        {
            (void) t;
            // loop over the vector _vec_
            for (unsigned i = 0; i<vec.size(); i++)
            {
                vec[i]->operator()(gra, t);
            }
        }
    /**
     destructor, release all heap memories used
     */
        ~CombinedConstraints()
        {
            //loop over _vec_ to release all heap memories used
            for (unsigned i = 0; i<vec.size(); i++)
            {
                delete vec[i];
            }
        }

};

/**
 combine three constraints of potentially different types in arbitrary order
 @param cnst1 the first constraint
 @param cnst2 the second constraint
 @param cnst3 the third constraint
 @return a CombinedConstraints object that can be used to combine all constraints
 */
template <typename C1, typename C2, typename C3>
CombinedConstraints make_combined_constraints(C1 cnst1, C2 cnst2, C3 cnst3)
{

    Constraint* c1 = new C1(cnst1);
    Constraint* c2 = new C2(cnst2);
    Constraint* c3 = new C3(cnst3);
    std::vector<Constraint*> vec {c1, c2, c3};
    return CombinedConstraints(vec);
}
/**
 combine two constraints of potentially different types in arbitrary order
 @param cnst1 the first constraint
 @param cnst2 the second constraint
 @return a CombinedConstraints object that can be used to combine all constraints
 */
template <typename C1, typename C2>
CombinedConstraints make_combined_constraints(C1 cnst1, C2 cnst2)
{

    Constraint* c1 = new C1(cnst1);
    Constraint* c2 = new C2(cnst2);
    std::vector<Constraint*> vec {c1, c2};
    return CombinedConstraints(vec);
}


 /**a new symp_euler_step that can takes constraint as an argument
 Change a graph's nodes according to a step of the symplectic Euler
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
  *@tparam C is a combined contraint applied to graph g
  */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint)
{
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
    //--design_1
    //--should apply constraint only once
    //--START
    constraint(g, t);
    //--END

      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
    }

    return t + dt;

}
