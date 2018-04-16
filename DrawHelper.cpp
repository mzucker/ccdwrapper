#include "DrawHelper.h"

namespace ccdw {
    
DrawHelper::DrawHelper(): 
    quadric(0),
    slices(32),
    sstacks(16),
    cstacks(1),
    dilation(0)
{}

DrawHelper::~DrawHelper() {
    if (quadric) { 
        gluDeleteQuadric(quadric);
    }
}

GLUquadric* DrawHelper::getQuadric() {
    if (!quadric) {
        quadric = gluNewQuadric();
    }
    return quadric;
}

void DrawHelper::visit(const Box* b) {

    Box3 box(-0.5*b->extents, 0.5*b->extents);

    if (dilation) {

        glstuff::draw_round_box(box, dilation, slices, sstacks);
      
    } else {

        glstuff::draw_box( box );

    }

}

void DrawHelper::visit(const Point* p) {

    if (dilation) {
        gluSphere(getQuadric(), dilation, slices, sstacks);
    } else {
        glPushAttrib(GL_ENABLE_BIT);
        glDisable(GL_LIGHTING);
        glBegin(GL_POINTS);
        glVertex3f(0,0,0);
        glEnd();
        glPopAttrib();
    }

}

void DrawHelper::visit(const Line* l) {

    vec3 p0(0, 0, -0.5*l->length);
    vec3 p1(0, 0,  0.5*l->length);

    if (dilation) {
        glstuff::draw_capsule(getQuadric(), p0, p1, dilation, 
                              slices, sstacks, cstacks);
    } else {
        glPushAttrib(GL_ENABLE_BIT);
        glDisable(GL_LIGHTING);
        glBegin(GL_LINES);
        glstuff::vertex(p0);
        glstuff::vertex(p1);
        glEnd();
        glPopAttrib();
    }

}

void DrawHelper::visit(const Cylinder* c) {

    vec3 p0(0, 0, -0.5*c->length);
    vec3 p1(0, 0,  0.5*c->length);

    if (dilation) {
        // TODO: deal with half-torii
    } else {
        glstuff::draw_cylinder(getQuadric(), p0, p1, c->radius, 
                               slices, cstacks);
    }

}

void DrawHelper::visit(const DilatedConvex* d) {
    assert(d->child);
    double old_dilation = dilation;
    dilation += d->dilation;
    render(d->child);
    dilation = old_dilation;
}

void DrawHelper::visit(const TransformedConvex* t) {
    assert(t->child);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glstuff::mult_transform(t->xform);
    render(t->child);
    glPopMatrix();
}

void DrawHelper::render(const Convex* c) {
    c->accept(this);
}
    
}

