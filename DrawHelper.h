#ifndef _CCDWRAPPERDRAW_H_
#define _CCDWRAPPERDRAW_H_

#include "glstuff.h"
#include "CCDWrapper.h"

namespace ccdw {

class DrawHelper: public ConvexConstVisitor {
public:
    
    GLUquadric* quadric;
    
    int slices;
    int sstacks;
    int cstacks;

    double dilation;
    
    DrawHelper();
    virtual ~DrawHelper();
    
    GLUquadric* getQuadric();

    virtual void visit(const Point* c);
    virtual void visit(const Line* c);
    virtual void visit(const Box* c);
    virtual void visit(const Cylinder* c);
    virtual void visit(const DilatedConvex* c);
    virtual void visit(const TransformedConvex* c);

    void render(const Convex* c);
    
};

}

#endif
