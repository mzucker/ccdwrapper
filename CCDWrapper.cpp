#include "CCDWrapper.h"
#include <assert.h>
#include <string.h>
#include <sstream>
#include <algorithm>

namespace ccdw {

static inline vec3& ccd2mz(ccd_vec3_t* c) { 
    return *((vec3*)c); 
}

static inline const vec3& ccd2mz(const ccd_vec3_t* c) { 
    return *((const vec3*)c); 
}

static inline ccd_vec3_t* mz2ccd(vec3& v) { 
    return (ccd_vec3_t*)(&v); 
}

static inline vec3 vabs(const vec3& v) {
    return vec3(::fabs(v[0]), ::fabs(v[1]), ::fabs(v[2]));
}

    /*
static inline vec3 vsgn(const vec3& v) {
    return vec3(ccdSign(v[0]), ccdSign(v[1]), ccdSign(v[2]));
}
    */
    
static inline vec3 vmin(const vec3& a, 
                        const vec3& b) {
    return vec3(std::min(a[0], b[0]),
                std::min(a[1], b[1]),
                std::min(a[2], b[2]));
}

static inline vec3 vmax(const vec3& a, 
                        const vec3& b) {
    return vec3(std::max(a[0], b[0]),
                std::max(a[1], b[1]),
                std::max(a[2], b[2]));
}
 
static inline int vargmin(const vec3& a) {
    if (a[0] < a[1]) { 
        // a1 not min
        if (a[0] < a[2]) {
            return 0;
        } else {
            return 2;
        }
    } else { 
        // a0 not min
        if (a[1] < a[2]) {
            return 1;
        } else {
            return 2;
        }
    }
}

    /*
static inline int vargmax(const vec3& a) {
    if (a[0] > a[1]) { 
        // a1 not max
        if (a[0] > a[2]) {
            return 0;
        } else {
            return 2;
        }
    } else { 
        // a0 not max
        if (a[1] > a[2]) {
            return 1;
        } else {
            return 2;
        }
    }
}
    */

//////////////////////////////////////////////////////////////////////

void ConvexConstVisitor::visit(const Convex* c) {}
void ConvexConstVisitor::visit(const Point* c) {}
void ConvexConstVisitor::visit(const Line* c) {}
void ConvexConstVisitor::visit(const Box* c) {}
void ConvexConstVisitor::visit(const Cylinder* c) {}
void ConvexConstVisitor::visit(const DilatedConvex* c) {}
void ConvexConstVisitor::visit(const TransformedConvex* c) {}

ConvexConstVisitor::~ConvexConstVisitor() {}

//////////////////////////////////////////////////////////////////////

void support(const void* obj,
             const ccd_vec3_t* dir,
             ccd_vec3_t* vec) {

    const Convex* c = (const Convex*)obj;
    assert( c );
    c->support( ccd2mz(dir), ccd2mz(vec) );

}
  
void center(const void* obj,
            ccd_vec3_t* center) {

    const Convex* c = (const Convex*)obj;
    assert( c );
    
    c->center( ccd2mz(center) );

}

//////////////////////////////////////////////////////////////////////

Convex::~Convex() {}

std::string Convex::description() const {
    std::ostringstream ostr;
    describe(ostr);
    return ostr.str();
}
    

void Convex::center(vec3& c) const {
    c = vec3(0,0,0);
}

bool Convex::isDilated() const {
    return false;
}

void Convex::accept(ConvexConstVisitor* v) const {
    v->visit(this);
}

//////////////////////////////////////////////////////////////////////

Box::Box(): extents(1) {}

Box::Box(const vec3& e): extents(e) {}

Box::~Box() {}

void Box::describe(std::ostream& ostr) const {
    ostr << "Box(" << extents.transpose() << ")";
}

bool Box::contains(const vec3& p, vec3* pc) const {

    // inside if abs(p[i]) <= 0.5*extents[i] for all i. 
    // outside if abs(p[i]) > 0.5*extents[i] for some i
    // if inside, closest point on boundary is based on distance

    const vec3 a = vabs(p);
    const vec3 b = ccd_real_t(0.5)*extents;

    // bounds minus point
    const vec3 c = b-a;

    // index of smallest (most negative) excess -- negative means outside box
    int i = vargmin(c);

    // if smallest is not negative, we are inside
    bool inside = (c[i] >= 0);

    if (pc) {
        if (inside) {
            // copy original point but force axis closest to box onto boundary
            *pc = p;
            (*pc)[i] = ccdSign(p[i]) * b[i];
        } else {
            // clamp to box boundary
            *pc = vmax(-b, vmin(p, b));
        }
    }

    return inside;

}
  

/*
  void Box::closest(vec3& v) const {

  bool inside = true;

  vec3 absdists(0), vout(0);
  int minaxis = 0;
    
  for (int axis=0; axis<3; ++axis) {
  if (v[axis] < -0.5*extents[axis]) {
  inside = false;
  v[axis] = -0.5*extents[axis];
  } else if (v[axis] > 0.5*extents[axis]) {
  inside = false;
  v[axis] = 0.5*extents[axis];
  } else {
  if (v[axis] < 0) {
  absdists[axis] = 0.5*extents[axis] + v[axis];
  vout[axis] = -0.5*extents[axis];
  } else {
  absdists[axis] = 0.5*extents[axis] - v[axis];
  vout[axis] = 0.5*extents[axis];
  }
  if (absdists[axis] < absdists[minaxis]) {
  minaxis = axis;
  }
  }
  }

  if (inside) { 
  v[minaxis] = vout[minaxis];
  }

  };
*/

ccd_real_t Box::maxDist() const {
    return 0.5*extents.norm();
}
  
void Box::support(const vec3& dir, vec3& s) const {

    s = vec3( ccdSign(dir[0]) * extents[0] * 0.5,
              ccdSign(dir[1]) * extents[1] * 0.5,
              ccdSign(dir[2]) * extents[2] * 0.5 );

}

void Box::accept(ConvexConstVisitor* v) const {
    v->visit(this);
}
    

/*
  static inline vec3 cwiseproduct(const vec3& a, const vec3& b) {
  return vec3(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
  }
*/

  
//////////////////////////////////////////////////////////////////////


Point::Point() {}

Point::~Point() {}

void Point::describe(std::ostream& ostr) const {
    ostr << "Point()";
}

ccd_real_t Point::maxDist() const {
    return 0;
}

void Point::support(const vec3& dir, vec3& s) const {

    s = vec3(0,0,0);

}


bool Point::contains(const vec3& p, vec3* pc) const {
    if (pc) { *pc = vec3(0,0,0); }
    return p.norm() == 0;
}

void Point::accept(ConvexConstVisitor* v) const {
    v->visit(this);
}
    

//////////////////////////////////////////////////////////////////////

Line::Line(): length(1) {}

Line::Line(ccd_real_t l): length(l) {}
  
Line::~Line() {}

void Line::describe(std::ostream& ostr) const {
    ostr << "Line(" << length << ")";
}

/*
  void Line::closest(vec3& v) const {

  if (v[2] < -0.5*length) {
  v[2] = -0.5*length;
  } else if (v[2] > 0.5*length) {
  v[2] = 0.5*length;
  }
  v[0] = v[1] = 0;

  }
*/

ccd_real_t Line::maxDist() const {
    return 0.5*length;
}

void Line::support(const vec3& dir, vec3& s) const {
    s = vec3( 0, 0, ccdSign(dir[2]) * length * 0.5 );
}

    
bool Line::contains(const vec3& p, vec3* pc) const {

    if (pc) {
        *pc = vec3(0, 0, 
                   std::max(-ccd_real_t(0.5)*length, 
                            std::min(ccd_real_t(0.5)*length, p.z())));
    }

    return ( p.x() == 0 && 
             p.y() == 0 && 
             fabs(p.z()) <= 0.5*length );

}

void Line::accept(ConvexConstVisitor* v) const {
    v->visit(this);
}
    

//////////////////////////////////////////////////////////////////////


Cylinder::Cylinder(): length(1), radius(0.5) {}

Cylinder::Cylinder(ccd_real_t l, ccd_real_t r): length(l), radius(r) {}
  
Cylinder::~Cylinder() {}

void Cylinder::describe(std::ostream& ostr) const {
    ostr << "Cylinder(" << length << ", " << radius << ")";
}

ccd_real_t Cylinder::maxDist() const {
    return sqrt(0.25*length*length + radius*radius);
}

void Cylinder::support(const vec3& dir, vec3& s) const {
    //vec2_t<ccd_real_t> p = dir.trunc();
    Eigen::Vector2d p(dir[0], dir[1]);
    p *= (radius / p.norm());
    s = vec3( p[0], p[1], ccdSign(dir[2]) * length * 0.5 );
}


bool Cylinder::contains(const vec3& p, vec3* pc) const {

    const ccd_real_t r2 = p.x()*p.x() + p.y()*p.y();
    const ccd_real_t az = fabs(p.z());
    
    const ccd_real_t hl = ccd_real_t(0.5)*length;
      
    bool inside = ( r2 <= radius*radius &&
                    az <= hl );

    if (pc) {

        ccd_real_t r = sqrt(r2);

        if (inside) {

            ccd_real_t dr = radius - r;
            ccd_real_t dz = az < hl;

            if (dr < dz) {

                // snap to radius
                *pc = vec3(p.x()*radius/r,
                           p.y()*radius/r,
                           p.z());

            } else {

                // snap to top/bottom
                *pc = vec3(p.x(), p.y(), p.z() < ccd_real_t(0) ? -hl : hl);

            }

        } else {

            *pc = vec3(p.x()*radius/r, 
                       p.y()*radius/r, 
                       std::max(-hl, std::min(p.z(), hl)));

        }
    }

    return inside;

}

void Cylinder::accept(ConvexConstVisitor* v) const {
    v->visit(this);
}

//////////////////////////////////////////////////////////////////////

DilatedConvex::DilatedConvex(): child(0), dilation(0) {}

DilatedConvex::DilatedConvex(const Convex* c, ccd_real_t r): 
    child(c), dilation(r) 
{
    assert(child);
    assert(dilation > 0);
}

DilatedConvex::~DilatedConvex() {}

void DilatedConvex::describe(std::ostream& ostr) const {
    ostr << "DilatedConvex(";
    child->describe(ostr);
    ostr << ", " << dilation << ")";
}

/*
  void DilatedConvex::closest(vec3& v) const {
  assert(child);
  assert(dilation > 0);
  vec3 vc;
  child->closest(vc);
  vec3 dir = v - vc;
  v = vc + dir * (dilation/dir.norm());
  }
*/


ccd_real_t DilatedConvex::maxDist() const {
    assert(child);
    assert(dilation > 0);
    return child->maxDist() + dilation;
}

void DilatedConvex::support(const vec3& dir, vec3& s) const {
    assert(child);
    assert(dilation > 0);
    child->support(dir, s);
    s += dir * (dilation / dir.norm());
}


bool DilatedConvex::isDilated() const {
    return true;
}

bool DilatedConvex::contains(const vec3& p, vec3* pc) const {

    vec3 cpc;

    bool inside_child = child->contains(p, &cpc);

    if (inside_child && !pc) {
        return true;
    }

    vec3 cvec = p-cpc;

    ccd_real_t l2 = cvec.dot(cvec);

    bool inside = inside_child || l2 < dilation*dilation;

    if (pc) {
        cvec *= dilation / sqrt(l2);
        if (inside_child) {
            *pc -= cvec;
        } else {
            *pc += cvec;
        }
    }

    return inside;

}

void DilatedConvex::accept(ConvexConstVisitor* v) const {
    v->visit(this);
}
    
//////////////////////////////////////////////////////////////////////

TransformedConvex::TransformedConvex(): child(0) {}
    
TransformedConvex::TransformedConvex(const Convex* c, const Transform3& x):
    child(c), xform(x)
{
    assert(child);
}

TransformedConvex::~TransformedConvex() {}

void TransformedConvex::describe(std::ostream& ostr) const {
    ostr << "TransformedConvex(";
    child->describe(ostr);
    ostr << ", " << "XFORM" << ")";
}

/*
  void TransformedConvex::closest(vec3& v) const {
  assert(child);
  v = xform.transformInv(v);
  child->closest(v);
  v = xform.transformFwd(v);
  }
*/

ccd_real_t TransformedConvex::maxDist() const {
    assert(child);
    return child->maxDist();
}

void TransformedConvex::support(const vec3& dir, vec3& s) const {
    assert(child);
    child->support( xform.rotation().inverse()*dir, s );
    s = xform * s;
}

void TransformedConvex::center(vec3& c) const {
    assert(child);
    child->center(c);
    c = xform * c;
}


bool TransformedConvex::isDilated() const {
    return child->isDilated();
}

bool TransformedConvex::contains(const vec3& p, vec3* pc) const {

    bool inside = child->contains(xform.transformInv(p), pc);

    if (pc) {
        *pc = xform.transformFwd(*pc);
    }

    return inside;

}

void TransformedConvex::accept(ConvexConstVisitor* v) const {
    v->visit(this);
}


//////////////////////////////////////////////////////////////////////

DilatedConvex* sphere(ccd_real_t radius) {
    return new DilatedConvex(new Point(), radius);
}

DilatedConvex* capsule(ccd_real_t length, ccd_real_t radius) {
    return new DilatedConvex(new Line(length), radius);
}

TransformedConvex* capsule(const vec3& p0,
                           const vec3& p1, 
                           ccd_real_t radius) {

    vec3 z = p1-p0;
    ccd_real_t length = z.norm();

    vec3 mid = 0.5*(p1+p0);

    return transform(capsule(length, radius),
                     Transform3(quatFromOneVector(z), mid));
    
}

TransformedConvex* cylinder(const vec3& p0,
                            const vec3& p1, 
                            ccd_real_t radius) {
    
    vec3 z = p1-p0;
    ccd_real_t length = z.norm();

    vec3 mid = 0.5*(p1+p0);
    
    return transform(new Cylinder(length, radius),
                     Transform3(quatFromOneVector(z), mid));
    
}

TransformedConvex* transform(const Convex* c, const Transform3& xform) {

    const TransformedConvex* t = dynamic_cast<const TransformedConvex*>(c);

    if (t) {
        return transform(t->child, xform * t->xform);
    } else {
        return new TransformedConvex(c, xform);
    }

}

Convex* dilate(const Convex* c, ccd_real_t radius) {

    const TransformedConvex* t = dynamic_cast<const TransformedConvex*>(c);
    
    if (t) {
        return transform(dilate(t->child, radius), t->xform);
    } else {
        const DilatedConvex* d = dynamic_cast<const DilatedConvex*>(c);
        if (d) {
            return dilate(d->child, radius + d->dilation);
        } else {
            return new DilatedConvex(c, radius);
        }
    }

}

//////////////////////////////////////////////////////////////////////

Report::Report(const Convex* a, const Convex* b): 
    c1(a), c2(b), flags(0), distance(0), algorithm(NUM_ALGORITHMS) {}

Checker::Checker() {
    CCD_INIT(&ccd);
    ccd.support1 = ccdw::support;
    ccd.support2 = ccdw::support;
    ccd.center1 = ccdw::center;
    ccd.center2 = ccdw::center;
    ccd.max_iterations = 1000;
    algorithm = ALGORITHM_MPR;
}

bool Checker::intersect(const Convex* c1, const Convex* c2, 
                        ccd_real_t dmin) const {

    return query(QUERY_INTERSECT, 0, c1, c2, dmin);

}
  
bool Checker::separate(Report& report,
                       const Convex* c1, const Convex* c2, 
                       ccd_real_t dmin) const {

    return query(QUERY_SEPARATION, &report, c1, c2, dmin);

}
  
bool Checker::penetration(Report& report,
                          const Convex* c1, const Convex* c2,
                          ccd_real_t dmin) const {
    
    return query(QUERY_PENETRATION, &report, c1, c2, dmin);

}
  
bool Checker::query(QueryType qtype, Report* report,
                    const Convex* orig_c1, const Convex* orig_c2,
                    ccd_real_t dmin) const {

    assert( orig_c1 );
    assert( orig_c2 );
    assert( dmin >= 0 );

    if (report) {
        *report = Report(orig_c1, orig_c2);
    }

    // do bounding sphere test
    vec3 ctr1, ctr2;
    orig_c1->center(ctr1);
    orig_c2->center(ctr2);

    ccd_real_t r1 = orig_c1->maxDist();
    ccd_real_t r2 = orig_c2->maxDist();

    ccd_real_t d = dmin + r1 + r2;

    if ((ctr2-ctr1).norm()  > d) {
        return false;
    }

    const Convex* c1 = orig_c1;
    const Convex* c2 = orig_c2;

    // do relative transform if necessary
    TransformedConvex* transformed = 0;

    const TransformedConvex* t1 = dynamic_cast<const TransformedConvex*>(orig_c1);
    const TransformedConvex* t2 = dynamic_cast<const TransformedConvex*>(orig_c2);

    if (t1 && t2) {
        transformed = transform(t2->child, t1->xform.inverse() * t2->xform);
        c1 = t1->child;
        c2 = transformed;
    }



    AlgorithmType actual_algorithm = algorithm;
    if (qtype == QUERY_SEPARATION) {
        actual_algorithm = ALGORITHM_GJK;
    }

    Convex* dilated1 = 0;
    Convex* dilated2 = 0;

    if (dmin) {
        c1 = dilated1 = dilate(c1, 0.5*dmin);
        c2 = dilated2 = dilate(c2, 0.5*dmin);
    }

    bool intersect = false;

    switch (qtype) {
    case QUERY_INTERSECT:
        if (actual_algorithm == ALGORITHM_GJK) {
            intersect = ccdGJKIntersect(c1, c2, &ccd);
        } else {
            intersect = ccdMPRIntersect(c1, c2, &ccd);
        }
        if (report && intersect) { 
            report->flags = INTERSECT;
        }
        break;
    case QUERY_SEPARATION: {
        vec3 sep;
        assert(actual_algorithm == ALGORITHM_GJK);
        intersect = (ccdGJKSeparate(c1, c2, &ccd, mz2ccd(sep)) == 0);
        if (report && intersect) {

            report->flags = INTERSECT | HAVE_SEPARATION;

            vec3 s1, s2;
            c1->support( sep, s1);
            c2->support(-sep, s2);

            report->pos1 = s1;
            report->pos2 = s2;

            report->distance = sep.norm();
            report->direction = sep / report->distance;

        }
        break;
    }
    case QUERY_PENETRATION: {
        ccd_real_t depth;
        vec3 dir, pos;
        if (actual_algorithm == ALGORITHM_GJK) {
            intersect = (ccdGJKPenetration(c1, c2, &ccd, &depth, 
                                           mz2ccd(dir), mz2ccd(pos)) == 0);
        } else {
            intersect = (ccdMPRPenetration(c1, c2, &ccd, &depth, 
                                           mz2ccd(dir), mz2ccd(pos)) == 0);
        }
        if (report && intersect) {
            report->flags = INTERSECT | HAVE_SEPARATION | HAVE_POSITION;
            report->distance = depth;
            report->direction = dir;
            report->pos1 = pos + (0.5*depth)*dir;
            report->pos2 = pos - (0.5*depth)*dir;
        }
        break;
    }
    default:
        break;
    };

    if (transformed && report) {
        report->pos1 = t1->xform * report->pos1;
        report->pos2 = t1->xform * report->pos2;
        report->direction = t1->xform.rotation() * report->direction;
    }

    if (dmin && report) {
        report->pos1 -= 0.5*dmin*report->direction;
        report->pos2 += 0.5*dmin*report->direction;
        report->distance -= dmin;
    }

    if (report) {
        report->algorithm = actual_algorithm;
    }

    delete dilated1;
    delete dilated2;

    delete transformed;
    
    return intersect;

}

//////////////////////////////////////////////////////////////////////

static inline ccd_real_t frac(size_t i, size_t n) {
    return ccd_real_t(i) / ccd_real_t(n-1);
}

void cubePoints(size_t n, std::vector<vec3>& points) {

    points.clear();

    
    for (int side=0; side<6; ++side) {

        int ax0 = (side + 0) % 3;
        int ax1 = (side + 1) % 3;
        int ax2 = (side + 2) % 3;

        int sign = (side / 3) ? -1 : 1;
      
        vec3 p;
      
        p[ax0] = sign;

        for (size_t u=0; u<n; ++u) {
            p[ax1] = sign*(frac(u, n)*2-1);
            for (size_t v=0; v<n; ++v) {
                p[ax2] = sign*(frac(v, n)*2-1);
                points.push_back(p);
            }
        }

    }

}


}
