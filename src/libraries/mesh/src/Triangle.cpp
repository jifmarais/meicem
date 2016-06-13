#include <assert.h>
#include "Triangle.hpp"

Triangle::Triangle()
{
    Node p {0.0, 0.0, 0.0};
    set(p , p, p);
}

Triangle::Triangle(Node n1,
                   Node n2,
                   Node n3)
{
   set(n1, n2, n3);
}

Triangle::~Triangle()
{
    //dtor
}

void Triangle::set(Node n1,
                   Node n2,
                   Node n3)
{
    m_nodes[0] = n1;
    m_nodes[1] = n2;
    m_nodes[2] = n3;
}

bool Triangle::operator==(const Triangle& rhs) const
{
    return ( m_nodes[0] == rhs.n1() &&
             m_nodes[1] == rhs.n2() &&
             m_nodes[2] == rhs.n3() );
}

bool Triangle::operator!=(const Triangle& rhs) const
{
    return not ( this->operator==(rhs) );
}

bool operator==(const Triangle& lhs, const Triangle& rhs)
{
    return lhs.operator==(rhs);
}

bool operator!=(const Triangle& lhs, const Triangle& rhs)
{
    return lhs.operator!=(rhs);
}

double Triangle::area() const
{
    return 0.5*Node::cross(m_nodes[1]-m_nodes[0], m_nodes[2]-m_nodes[0]).magnitude();
}

Node Triangle::centre() const
{
    return ( m_nodes[0] + m_nodes[1] + m_nodes[2] ) / 3.0;
}

Node Triangle::toSimplex(const Node& p) const
{
    // From "The Method of Moments in Electromagnetics" chapter 9, p. 262.
    // The equations there work, but for special cases it does not and I
    // changed it as below.
    // TODO: I'm not quite happy with this implementation, aspecially alpha.
    // I'm also not super happy with mapping of x, y, z to alpha, beta and gamma.

    Node p1mp2 = m_nodes[0] - m_nodes[1];
    Node p3mp1 = m_nodes[2] - m_nodes[0];

    double beta  = ( ((p.x() - m_nodes[0].x())*p1mp2.z()*p1mp2.y()) +
                     ((p.y() - m_nodes[0].y())*p1mp2.z()*p1mp2.x()) +
                 (2.0*(p.z() - m_nodes[0].z())*p1mp2.x()*p1mp2.y()) ) /
                    ( (p3mp1.x()*p1mp2.z()*p1mp2.y()) +
                      (p3mp1.y()*p1mp2.z()*p1mp2.x()) +
                  (2.0*p3mp1.z()*p1mp2.x()*p1mp2.y()) );
    double alpha = 0;
    if (p1mp2.x() != 0.0)
    {
        alpha = (beta*p3mp1.x() + m_nodes[0].x() - p.x()) / p1mp2.x();
    }
    else
    {
        alpha = (beta*p3mp1.y() + m_nodes[0].y() - p.y()) / p1mp2.y();
    }

    Node BarryCentricNode; // alpha <-> y, beta <-> z, gamma <-> x (origin)
    BarryCentricNode.set( 1- alpha - beta, alpha, beta);
    return BarryCentricNode;
}

//Node Triangle::toSimplex(const Node& p) const
//{
//    // From "The Method of Moments in Electromagnetics" chapter 9, p. 262.

//    double A   = m_nodes[1].x() - m_nodes[0].x();
//    double B   = m_nodes[2].x() - m_nodes[0].x();
//    double C   = m_nodes[0].x() - p.x();
//    double DpG = (m_nodes[1].y() - m_nodes[0].y()) + (m_nodes[1].z() - m_nodes[0].z());
//    double EpH = (m_nodes[2].y() - m_nodes[0].y()) + (m_nodes[2].z() - m_nodes[0].z());
//    double FpI = (m_nodes[0].y() - p.y()) + (m_nodes[0].z() - p.z());

//    Node BarryCentricNode; // alpha <-> y, beta <-> z, gamma <-> x (origin)
//    double alpha = (B*FpI - C*EpH)/(A*EpH - B*DpG);
//    double beta  = (A*FpI - C*DpG)/(B*DpG - A*EpH);
//    BarryCentricNode.set( 1- alpha - beta, alpha, beta);
//    return BarryCentricNode;
//}

Node Triangle::normal() const
{
    Node normalVector;
    normalVector = Node::cross( (m_nodes[1]-m_nodes[0]).norm(), (m_nodes[2]-m_nodes[0]).norm() );
    return normalVector.norm();
}

Node Triangle::fromSimplex(const Node& p) const
{
    assert(p.x() + p.y() + p.z() == 1);

    double x = m_nodes[0].x()*p.x() + m_nodes[1].x()*p.y() + m_nodes[2].x()*p.z();
    double y = m_nodes[0].y()*p.x() + m_nodes[1].y()*p.y() + m_nodes[2].y()*p.z();
    double z = m_nodes[0].z()*p.x() + m_nodes[1].z()*p.y() + m_nodes[2].z()*p.z();

    Node globalCoordinate;
    globalCoordinate.set(x, y, z);
    return globalCoordinate;
}

Node Triangle::operator[](unsigned index) const
{
    assert(index < 3);
    return m_nodes[index];
}

Node Triangle::n1() const
{
    return m_nodes[0];
}

Node Triangle::n2() const
{
    return m_nodes[1];
}

Node Triangle::n3() const
{
    return m_nodes[2];
}
